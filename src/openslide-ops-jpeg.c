/*
 *  OpenSlide, a library for reading whole slide image files
 *
 *  Copyright (c) 2007-2008 Carnegie Mellon University
 *  All rights reserved.
 *
 *  OpenSlide is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, version 2.
 *
 *  OpenSlide is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with OpenSlide. If not, see <http://www.gnu.org/licenses/>.
 *
 *  Linking OpenSlide statically or dynamically with other modules is
 *  making a combined work based on OpenSlide. Thus, the terms and
 *  conditions of the GNU General Public License cover the whole
 *  combination.
 */

/*
 * Part of this file is:
 *
 * Copyright (C) 1994-1996, Thomas G. Lane.
 * This file is part of the Independent JPEG Group's software.
 * For conditions of distribution and use, see the accompanying README file.
 */

#include "config.h"

#include <glib.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <jpeglib.h>
#include <jerror.h>
#include <inttypes.h>

#include <sys/types.h>   // for off_t ?

#include "openslide-private.h"
#include "openslide-cache.h"
#include "openslide-tilehelper.h"

struct one_jpeg {
  FILE *f;
  int64_t file_size;

  int32_t mcu_starts_count;
  int64_t *mcu_starts;
  int64_t *unreliable_mcu_starts;

  int32_t tile_width;
  int32_t tile_height;

  int32_t width;
  int32_t height;

  char *comment;
};

struct layer {
  struct one_jpeg **layer_jpegs; // count given by jpeg_w * jpeg_h

  // total size (not premultiplied by scale_denom)
  int64_t pixel_w;
  int64_t pixel_h;

  int32_t jpegs_across;       // how many distinct jpeg files across?
  int32_t jpegs_down;         // how many distinct jpeg files down?

  // the size of image (0,0), which is used to find the jpeg we want
  // from a given (x,y) (not premultiplied)
  int32_t image00_w;
  int32_t image00_h;

  int32_t scale_denom;
  double no_scale_denom_downsample;  // layer0_w div non_premult_pixel_w
};

struct jpegops_data {
  int32_t jpeg_count;
  struct one_jpeg *all_jpegs;

  // layer_count is in the osr struct
  struct layer *layers;

  // cache
  struct _openslide_cache *cache;

  // thread stuff, for background search of restart markers
  GMutex *restart_marker_mutex;
  GThread *restart_marker_thread;
  boolean restart_marker_thread_should_terminate;
};



/*
 * Source manager for doing fancy things with libjpeg and restart markers,
 * initially copied from jdatasrc.c from IJG libjpeg.
 */
struct my_src_mgr {
  struct jpeg_source_mgr pub;   /* public fields */

  JOCTET *buffer;               /* start of buffer */
  int buffer_size;
};

static void init_source (j_decompress_ptr cinfo) {
  /* nothing to be done */
}

static boolean fill_input_buffer (j_decompress_ptr cinfo) {
  /* everything is done already */
  return TRUE;
}


static void skip_input_data (j_decompress_ptr cinfo, long num_bytes) {
  struct my_src_mgr *src = (struct my_src_mgr *) cinfo->src;

  src->pub.next_input_byte += (size_t) num_bytes;
  src->pub.bytes_in_buffer -= (size_t) num_bytes;
}


static void term_source (j_decompress_ptr cinfo) {
  /* nothing to do */
}

static void jpeg_random_access_src (j_decompress_ptr cinfo, FILE *infile,
				    int64_t header_stop_position,
				    int64_t start_position, int64_t stop_position) {
  struct my_src_mgr *src;

  if (cinfo->src == NULL) {     /* first time for this JPEG object? */
    cinfo->src = (struct jpeg_source_mgr *)
      (*cinfo->mem->alloc_small) ((j_common_ptr) cinfo, JPOOL_PERMANENT,
				  sizeof(struct my_src_mgr));
  }

  src = (struct my_src_mgr *) cinfo->src;
  src->pub.init_source = init_source;
  src->pub.fill_input_buffer = fill_input_buffer;
  src->pub.skip_input_data = skip_input_data;
  src->pub.resync_to_restart = jpeg_resync_to_restart; /* use default method */
  src->pub.term_source = term_source;

  // compute size of buffer and allocate
  src->buffer_size = header_stop_position + stop_position - start_position;
  src->pub.bytes_in_buffer = src->buffer_size;
  src->buffer = g_slice_alloc(src->buffer_size);

  src->pub.next_input_byte = src->buffer;

  // read in the 2 parts
  //  g_debug("reading from %" PRId64, start_position);
  rewind(infile);
  fread(src->buffer, header_stop_position, 1, infile);
  fseeko(infile, start_position, SEEK_SET);
  fread(src->buffer + header_stop_position,
	stop_position - start_position, 1, infile);

  // change the final byte to EOI
  g_return_if_fail(src->buffer[header_stop_position] != 0xFF);
  g_return_if_fail(src->buffer[src->buffer_size - 2] == 0xFF);
  src->buffer[src->buffer_size - 1] = JPEG_EOI;
}

static bool is_zxy_successor(int64_t pz, int64_t px, int64_t py,
			     int64_t z, int64_t x, int64_t y) {
  //  g_debug("p_zxy: (%" PRId64 ",%" PRId64 ",%" PRId64 "), zxy: (%"
  //	  PRId64 ",%" PRId64 ",%" PRId64 ")",
  //	  pz, px, py, z, x, y);
  if (z == pz + 1) {
    return x == 0 && y == 0;
  }
  if (z != pz) {
    return false;
  }

  // z == pz

  if (y == py + 1) {
    return x == 0;
  }
  if (y != py) {
    return false;
  }

  // y == py

  return x == px + 1;
}

static guint int64_hash(gconstpointer v) {
  int64_t i = *((const int64_t *) v);
  return i ^ (i >> 32);
}

static gboolean int64_equal(gconstpointer v1, gconstpointer v2) {
  return *((int64_t *) v1) == *((int64_t *) v2);
}

static void int64_free(gpointer data) {
  g_slice_free(int64_t, data);
}

static void layer_free(gpointer data) {
  //  g_debug("layer_free: %p", data);

  struct layer *l = data;

  //  g_debug("g_free(%p)", (void *) l->layer_jpegs);
  g_free(l->layer_jpegs);
  g_slice_free(struct layer, l);
}

static void print_wlmap_entry(gpointer key, gpointer value,
			      gpointer user_data) {
  int64_t k = *((int64_t *) key);
  struct layer *v = (struct layer *) value;

  g_debug("%" PRId64 " -> ( pw: %" PRId64 ", ph: %" PRId64
	  ", jw: %" PRId32 ", jh: %" PRId32 ", scale_denom: %" PRId32
	  ", img00_w: %" PRId32 ", img00_h: %" PRId32 ", no_scale_denom_downsample: %g )",
	  k, v->pixel_w, v->pixel_h, v->jpegs_across, v->jpegs_down, v->scale_denom, v->image00_w, v->image00_h, v->no_scale_denom_downsample);
}

static void generate_layers_into_map(GSList *jpegs,
				     int32_t jpegs_across, int32_t jpegs_down,
				     int64_t pixel_w, int64_t pixel_h,
				     int32_t image00_w, int32_t image00_h,
				     int64_t layer0_w,
				     GHashTable *width_to_layer_map) {
  // JPEG files can give us 1/1, 1/2, 1/4, 1/8 downsamples, so we
  // need to create 4 layers per set of JPEGs

  int32_t num_jpegs = jpegs_across * jpegs_down;

  int scale_denom = 1;
  while (scale_denom <= 8) {
    // create layer
    struct layer *l = g_slice_new0(struct layer);
    l->jpegs_across = jpegs_across;
    l->jpegs_down = jpegs_down;
    l->pixel_w = pixel_w;
    l->pixel_h = pixel_h;
    l->scale_denom = scale_denom;
    l->image00_w = image00_w;
    l->image00_h = image00_h;
    l->no_scale_denom_downsample = (double) layer0_w / (double) pixel_w;

    // create array and copy
    l->layer_jpegs = g_new(struct one_jpeg *, num_jpegs);
    //    g_debug("g_new(struct one_jpeg *) -> %p", (void *) l->layer_jpegs);
    GSList *jj = jpegs;
    for (int32_t i = 0; i < num_jpegs; i++) {
      g_assert(jj);
      l->layer_jpegs[i] = (struct one_jpeg *) jj->data;
      jj = jj->next;
    }

    // put into map
    int64_t *key = g_slice_new(int64_t);
    *key = l->pixel_w / l->scale_denom;

    //    g_debug("insert %" PRId64 ", scale_denom: %d", *key, scale_denom);
    g_hash_table_insert(width_to_layer_map, key, l);

    scale_denom <<= 1;
  }
}

static GHashTable *create_width_to_layer_map(int32_t count,
					     struct _openslide_jpeg_fragment **fragments,
					     struct one_jpeg *jpegs) {
  int64_t prev_z = -1;
  int64_t prev_x = -1;
  int64_t prev_y = -1;

  GSList *layer_jpegs_tmp = NULL;
  int64_t l_pw = 0;
  int64_t l_ph = 0;

  int32_t img00_w = 0;
  int32_t img00_h = 0;

  int64_t layer0_w = 0;

  // int* -> struct layer*
  GHashTable *width_to_layer_map = g_hash_table_new_full(int64_hash,
							 int64_equal,
							 int64_free,
							 layer_free);

  // go through the fragments, accumulating to layers
  for (int32_t i = 0; i < count; i++) {
    struct _openslide_jpeg_fragment *fr = fragments[i];
    struct one_jpeg *oj = jpegs + i;

    // the fragments MUST be in sorted order by z,x,y
    g_assert(is_zxy_successor(prev_z, prev_x, prev_y,
			      fr->z, fr->x, fr->y));

    // special case for first layer
    if (prev_z == -1) {
      prev_z = 0;
      prev_x = 0;
      prev_y = 0;
    }

    // save first image dimensions
    if (fr->x == 0 && fr->y == 0) {
      img00_w = oj->width;
      img00_h = oj->height;
    }

    // accumulate size
    if (fr->y == 0) {
      l_pw += oj->width;
    }
    if (fr->x == 0) {
      l_ph += oj->height;
    }

    //    g_debug(" pw: %" PRId64 ", ph: %" PRId64, l_pw, l_ph);

    // accumulate to layer
    layer_jpegs_tmp = g_slist_prepend(layer_jpegs_tmp, oj);

    // is this the end of this layer? then flush
    if (i == count - 1 || fragments[i + 1]->z != fr->z) {
      layer_jpegs_tmp = g_slist_reverse(layer_jpegs_tmp);

      // save layer0 width
      if (fr->z == 0) {
	layer0_w = l_pw;
      }

      generate_layers_into_map(layer_jpegs_tmp, fr->x + 1, fr->y + 1,
			       l_pw, l_ph,
			       img00_w, img00_h,
			       layer0_w,
			       width_to_layer_map);

      // clear for next round
      l_pw = 0;
      l_ph = 0;
      img00_w = 0;
      img00_h = 0;

      while (layer_jpegs_tmp != NULL) {
	layer_jpegs_tmp = g_slist_delete_link(layer_jpegs_tmp, layer_jpegs_tmp);
      }
    }

    // update prevs
    prev_z = fr->z;
    prev_x = fr->x;
    prev_y = fr->y;
  }

  return width_to_layer_map;
}


static void init_optimization(FILE *f,
			      int32_t *mcu_starts_count,
			      int64_t **mcu_starts) {
  struct jpeg_decompress_struct cinfo;
  struct jpeg_error_mgr jerr;

  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_decompress(&cinfo);
  rewind(f);
  jpeg_stdio_src(&cinfo, f);
  jpeg_read_header(&cinfo, TRUE);
  jpeg_start_decompress(&cinfo);

  int32_t MCUs = cinfo.MCUs_per_row * cinfo.MCU_rows_in_scan;
  *mcu_starts_count = MCUs / cinfo.restart_interval;
  *mcu_starts = g_new(int64_t, *mcu_starts_count);

  // init all to -1
  for (int32_t i = 0; i < *mcu_starts_count; i++) {
    (*mcu_starts)[i] = -1;
  }

  // set the first entry
  (*mcu_starts)[0] = ftello(f) - cinfo.src->bytes_in_buffer;

  jpeg_destroy_decompress(&cinfo);
}

static uint8_t find_next_ff_marker(FILE *f,
				   uint8_t *buf_start,
				   uint8_t **buf,
				   int buf_size,
				   int64_t file_size,
				   int64_t *after_marker_pos,
				   int *bytes_in_buf) {
  //g_debug("bytes_in_buf: %d", *bytes_in_buf);
  int64_t file_pos = ftello(f);
  boolean last_was_ff = false;
  *after_marker_pos = -1;
  while (true) {
    if (*bytes_in_buf == 0) {
      // fill buffer
      *buf = buf_start;
      int bytes_to_read = MIN(buf_size, file_size - file_pos);

      //g_debug("bytes_to_read: %d", bytes_to_read);
      size_t result = fread(*buf, bytes_to_read, 1, f);
      if (result == 0) {
	return 0;
      }

      file_pos += bytes_to_read;
      *bytes_in_buf = bytes_to_read;
    }

    // special case where the last time ended with FF
    if (last_was_ff) {
      //g_debug("last_was_ff");
      uint8_t marker = (*buf)[0];
      (*buf)++;
      (*bytes_in_buf)--;
      *after_marker_pos = file_pos - *bytes_in_buf;
      return marker;
    }

    // search for ff
    uint8_t *ff = memchr(*buf, 0xFF, *bytes_in_buf);
    if (ff == NULL) {
      // keep searching
      *bytes_in_buf = 0;
    } else {
      // ff found, advance buffer to consume everything including ff
      int offset = ff - *buf + 1;
      *bytes_in_buf -= offset;
      *buf += offset;
      g_assert(*bytes_in_buf >= 0);

      if (*bytes_in_buf == 0) {
	last_was_ff = true;
      } else {
	(*bytes_in_buf)--;
	(*buf)++;
	*after_marker_pos = file_pos - *bytes_in_buf;
	return ff[1];
      }
    }
  }
}

static void compute_mcu_start(FILE *f,
			      int64_t *mcu_starts,
			      int64_t *unreliable_mcu_starts,
			      int64_t file_size,
			      int64_t target) {
  if (mcu_starts[target] != -1) {
    // already done
    return;
  }

  g_assert(target != 0); // first item is always filled

  // check the unreliable_mcu_starts store first,
  // and use it if valid
  int64_t offset = -1;
  if (unreliable_mcu_starts != NULL) {
    offset = unreliable_mcu_starts[target];
  }

  if (offset != -1) {
    uint8_t buf[2];
    fseeko(f, offset - 2, SEEK_SET);

    size_t result = fread(buf, 2, 1, f);
    if (result == 0 ||
	buf[0] != 0xFF || buf[1] < 0xD0 || buf[1] > 0xD7) {
      g_warning("Restart marker not found in expected place");
    } else {
      mcu_starts[target] = offset;
      return;
    }
  }


  // otherwise, walk backwards, to find the first non -1 offset
  int64_t first_good = target - 1;
  while (mcu_starts[first_good] == -1) {
    first_good--;
  }
  //  g_debug("target: %d, first_good: %d", target, first_good);

  // now search for the new restart markers
  fseeko(f, mcu_starts[first_good], SEEK_SET);

  uint8_t buf_start[4096];
  uint8_t *buf = buf_start;
  int bytes_in_buf = 0;
  while (first_good < target) {
    int64_t after_marker_pos;
    uint8_t b = find_next_ff_marker(f, buf_start, &buf, 4096,
				    file_size,
				    &after_marker_pos,
				    &bytes_in_buf);
    g_assert(after_marker_pos > 0);
    //g_debug("after_marker_pos: %" PRId64, after_marker_pos);

    // EOI?
    if (b == JPEG_EOI) {
      // we're done
      break;
    } else if (b >= 0xD0 && b < 0xD8) {
      // marker
      mcu_starts[1 + first_good++] = after_marker_pos;
    }
  }
}


static void read_from_one_jpeg (struct one_jpeg *jpeg,
				uint32_t *dest,
				int32_t x, int32_t y,
				int32_t scale_denom) {
  //  g_debug("read_from_one_jpeg: %p, dest: %p, x: %d, y: %d, scale_denom: %d", (void *) jpeg, (void *) dest, x, y, scale_denom);

  // figure out where to start the data stream
  int32_t tile_y = y / jpeg->tile_height;
  int32_t tile_x = x / jpeg->tile_width;

  int32_t stride_in_tiles = jpeg->width / jpeg->tile_width;

  int64_t mcu_start = tile_y * stride_in_tiles + tile_x;

  int64_t stop_position;

  compute_mcu_start(jpeg->f,
		    jpeg->mcu_starts, jpeg->unreliable_mcu_starts,
		    jpeg->file_size,
		    mcu_start);
  if (jpeg->mcu_starts_count == mcu_start + 1) {
    // EOF
    stop_position = jpeg->file_size;
  } else {
    compute_mcu_start(jpeg->f,
		      jpeg->mcu_starts, jpeg->unreliable_mcu_starts,
		      jpeg->file_size,
		      mcu_start + 1);
    stop_position = jpeg->mcu_starts[mcu_start + 1];
  }

  // begin decompress
  struct jpeg_decompress_struct cinfo;
  struct jpeg_error_mgr jerr;

  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_decompress(&cinfo);

  jpeg_random_access_src(&cinfo, jpeg->f,
			 jpeg->mcu_starts[0],
			 jpeg->mcu_starts[mcu_start],
			 stop_position);

  jpeg_read_header(&cinfo, FALSE);
  cinfo.scale_denom = scale_denom;
  cinfo.image_width = jpeg->tile_width;  // cunning
  cinfo.image_height = jpeg->tile_height;

  jpeg_start_decompress(&cinfo);

  //g_debug("output_width: %d", cinfo.output_width);
  //g_debug("output_height: %d", cinfo.output_height);

  // allocate scanline buffers
  JSAMPARRAY buffer =
    g_slice_alloc(sizeof(JSAMPROW) * cinfo.rec_outbuf_height);
  gsize row_size =
    sizeof(JSAMPLE)
    * cinfo.output_width
    * 3;  // output components
  for (int i = 0; i < cinfo.rec_outbuf_height; i++) {
    buffer[i] = g_slice_alloc(row_size);
    //g_debug("buffer[%d]: %p", i, buffer[i]);
  }

  // decompress
  while (cinfo.output_scanline < cinfo.output_height) {
    JDIMENSION rows_read = jpeg_read_scanlines(&cinfo,
					       buffer,
					       cinfo.rec_outbuf_height);
    //g_debug("just read scanline %d", cinfo.output_scanline - rows_read);
    //g_debug(" rows read: %d", rows_read);
    int cur_buffer = 0;
    while (rows_read > 0) {
      // copy a row
      int32_t i;
      for (i = 0; i < (int32_t) cinfo.output_width; i++) {
	dest[i] = 0xFF000000 |                          // A
	  buffer[cur_buffer][i * 3 + 0] << 16 | // R
	  buffer[cur_buffer][i * 3 + 1] << 8 |  // G
	  buffer[cur_buffer][i * 3 + 2];        // B
      }

      // advance everything 1 row
      rows_read--;
      cur_buffer++;
      dest += cinfo.output_width;
    }
  }

  //  g_debug("pixels wasted: %llu", pixels_wasted);

  // free buffers
  for (int i = 0; i < cinfo.rec_outbuf_height; i++) {
    g_slice_free1(row_size, buffer[i]);
  }
  g_slice_free1(sizeof(JSAMPROW) * cinfo.rec_outbuf_height, buffer);

  // last thing, stop jpeg
  struct my_src_mgr *src = (struct my_src_mgr *) cinfo.src;   // sorry
  g_slice_free1(src->buffer_size, src->buffer);
  jpeg_destroy_decompress(&cinfo);
}


static void tilereader_read(void *tilereader_data,
			    uint32_t *dest, int64_t src_x, int64_t src_y) {
  struct layer *l = tilereader_data;

  src_y *= l->scale_denom;
  int32_t file_y = src_y / l->image00_h;
  int64_t origin_src_segment_y = file_y * (int64_t) l->image00_h;
  int32_t start_in_src_segment_y = src_y - origin_src_segment_y;

  src_x *= l->scale_denom;
  int32_t file_x = src_x / l->image00_w;
  int64_t origin_src_segment_x = file_x * (int64_t) l->image00_w;
  int32_t start_in_src_segment_x = src_x - origin_src_segment_x;

  int file_number = file_y * l->jpegs_across + file_x;
  g_assert(file_number < l->jpegs_across * l->jpegs_down);
  struct one_jpeg *jpeg = l->layer_jpegs[file_y * l->jpegs_across + file_x];

  read_from_one_jpeg(jpeg, dest,
		     start_in_src_segment_x, start_in_src_segment_y,
		     l->scale_denom);
}

static void read_region(openslide_t *osr, uint32_t *dest,
			int64_t x, int64_t y,
			int32_t layer,
			int64_t w, int64_t h) {
  //  g_debug("jpeg ops read_region: x: %" PRId64 ", y: %" PRId64 ", layer: %d, w: %" PRId64 ", h: %" PRId64 "",
  //	  x, y, layer, w, h);

  struct jpegops_data *data = osr->data;

  // get the layer
  struct layer *l = data->layers + layer;
  int32_t scale_denom = l->scale_denom;
  double rel_downsample = l->no_scale_denom_downsample;
  //  g_debug("layer: %d, rel_downsample: %g, scale_denom: %d",
  //	  layer, rel_downsample, scale_denom);

  // figure out tile dimensions
  int64_t tw = l->layer_jpegs[0]->tile_width / scale_denom;
  int64_t th = l->layer_jpegs[0]->tile_height / scale_denom;

  int64_t ds_x = x / rel_downsample / scale_denom;
  int64_t ds_y = y / rel_downsample / scale_denom;

  int64_t end_x = ds_x + w;
  int64_t end_y = ds_y + h;

  // check bounds
  int64_t pw = l->pixel_w / scale_denom;
  int64_t ph = l->pixel_h / scale_denom;
  if (end_x >= pw) {
    end_x = pw - 1;
  }
  if (end_y >= ph) {
    end_y = ph - 1;
  }

  g_mutex_lock(data->restart_marker_mutex);
  _openslide_read_tiles(ds_x, ds_y,
			end_x, end_y, 0, 0, w, h, layer, tw, th,
			tilereader_read, l,
			dest, data->cache);
  g_mutex_unlock(data->restart_marker_mutex);
}


static void destroy(openslide_t *osr) {
  struct jpegops_data *data = osr->data;

  // tell the thread to finish and wait
  g_mutex_lock(data->restart_marker_mutex);
  data->restart_marker_thread_should_terminate = true;
  g_mutex_unlock(data->restart_marker_mutex);
  g_thread_join(data->restart_marker_thread);

  // each jpeg in turn
  for (int32_t i = 0; i < data->jpeg_count; i++) {
    struct one_jpeg *jpeg = data->all_jpegs + i;

    fclose(jpeg->f);
    g_free(jpeg->mcu_starts);
    g_free(jpeg->unreliable_mcu_starts);
    g_free(jpeg->comment);
  }

  // each layer in turn
  for (int32_t i = 0; i < osr->layer_count; i++) {
    struct layer *l = data->layers + i;

    //    g_debug("g_free(%p)", (void *) l->layer_jpegs);
    g_free(l->layer_jpegs);
  }

  // the JPEG array
  g_free(data->all_jpegs);

  // the layer array
  g_free(data->layers);

  // the cache
  _openslide_cache_destroy(data->cache);

  // the mutex
  g_mutex_free(data->restart_marker_mutex);

  // the structure
  g_slice_free(struct jpegops_data, data);
}

static void get_dimensions(openslide_t *osr, int32_t layer,
			   int64_t *w, int64_t *h) {
  struct jpegops_data *data = osr->data;

  // check bounds
  if (layer >= osr->layer_count) {
    *w = 0;
    *h = 0;
    return;
  }

  struct layer *l = data->layers + layer;
  *w = l->pixel_w / l->scale_denom;
  *h = l->pixel_h / l->scale_denom;

  //  g_debug("dimensions of layer %" PRId32 ": (%" PRId32 ",%" PRId32 ")", layer, *w, *h);
}

static const char* get_comment(openslide_t *osr) {
  struct jpegops_data *data = osr->data;
  return data->all_jpegs[0].comment;
}

static struct _openslide_ops jpeg_ops = {
  .read_region = read_region,
  .destroy = destroy,
  .get_dimensions = get_dimensions,
  .get_comment = get_comment,
};


static void init_one_jpeg(struct one_jpeg *onej,
			  struct _openslide_jpeg_fragment *fragment) {
  FILE *f = onej->f = fragment->f;
  struct jpeg_decompress_struct cinfo;
  struct jpeg_error_mgr jerr;

  // file size
  fseeko(f, 0, SEEK_END);
  onej->file_size = ftello(f);

  // optimization
  g_assert((fragment->mcu_starts_count == 0 && fragment->mcu_starts == NULL) ||
	   (fragment->mcu_starts_count != 0 && fragment->mcu_starts != NULL));
  init_optimization(fragment->f,
		    &onej->mcu_starts_count, &onej->mcu_starts);
  g_assert((fragment->mcu_starts_count == 0)
	   || (fragment->mcu_starts_count == onej->mcu_starts_count));
  onej->unreliable_mcu_starts = fragment->mcu_starts;

  // init jpeg
  rewind(f);

  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_decompress(&cinfo);
  jpeg_stdio_src(&cinfo, f);

  // extract comment
  jpeg_save_markers(&cinfo, JPEG_COM, 0xFFFF);
  jpeg_read_header(&cinfo, FALSE);
  if (cinfo.marker_list) {
    // copy everything out
    char *com = g_strndup((const gchar *) cinfo.marker_list->data,
			  cinfo.marker_list->data_length);
    // but only really save everything up to the first '\0'
    onej->comment = g_strdup(com);
    g_free(com);
  }
  jpeg_save_markers(&cinfo, JPEG_COM, 0);  // stop saving

  // save dimensions
  jpeg_calc_output_dimensions(&cinfo);
  onej->width = cinfo.output_width;
  onej->height = cinfo.output_height;

  //  g_debug(" w: %d, h: %d", cinfo.output_width, cinfo.output_height);

  // save "tile" dimensions
  jpeg_start_decompress(&cinfo);
  onej->tile_width = onej->width /
    (cinfo.MCUs_per_row / cinfo.restart_interval);
  onej->tile_height = onej->height / cinfo.MCU_rows_in_scan;

  //  g_debug("jpeg \"tile\" dimensions: %dx%d", onej->tile_width, onej->tile_height);

  // destroy jpeg
  jpeg_destroy_decompress(&cinfo);
}

static gint width_compare(gconstpointer a, gconstpointer b) {
  int64_t w1 = *((const int64_t *) a);
  int64_t w2 = *((const int64_t *) b);

  g_assert(w1 >= 0 && w2 >= 0);

  return (w1 < w2) - (w1 > w2);
}

static void get_keys(gpointer key, gpointer value,
		     gpointer user_data) {
  GList *keys = *((GList **) user_data);
  keys = g_list_prepend(keys, key);
  *((GList **) user_data) = keys;
}

static void verify_mcu_starts(struct jpegops_data *data) {
  g_debug("verifying mcu starts");

  int32_t current_jpeg = 0;
  int32_t current_mcu_start = 1;

  while(current_jpeg < data->jpeg_count) {
    struct one_jpeg *oj = data->all_jpegs + current_jpeg;

    int64_t offset = oj->mcu_starts[current_mcu_start];
    g_assert(offset != -1);
    fseeko(oj->f, offset - 2, SEEK_SET);
    g_assert(getc(oj->f) == 0xFF);
    int marker = getc(oj->f);
    g_assert(marker >= 0xD0 && marker <= 0xD7);

    current_mcu_start++;
    if (current_mcu_start >= oj->mcu_starts_count) {
      current_mcu_start = 1;
      current_jpeg++;
      g_debug("done verifying jpeg %d", current_jpeg);
    }
  }
}

static gpointer restart_marker_thread_func(gpointer d) {
  struct jpegops_data *data = d;

  int32_t current_jpeg = 0;
  int32_t current_mcu_start = 0;

  while(current_jpeg < data->jpeg_count) {
    g_mutex_lock(data->restart_marker_mutex);

    // check for exit
    if (data->restart_marker_thread_should_terminate) {
      g_mutex_unlock(data->restart_marker_mutex);
      break;
    }

    //    g_debug("current_jpeg: %d, current_mcu_start: %d",
    //	    current_jpeg, current_mcu_start);

    struct one_jpeg *oj = data->all_jpegs + current_jpeg;
    compute_mcu_start(oj->f, oj->mcu_starts, oj->unreliable_mcu_starts,
		      oj->file_size,
		      current_mcu_start);

    current_mcu_start++;
    if (current_mcu_start >= oj->mcu_starts_count) {
      current_mcu_start = 0;
      current_jpeg++;
    }

    g_mutex_unlock(data->restart_marker_mutex);
  }

  //g_debug("restart_marker_thread_func done!");
  return NULL;
}

void _openslide_add_jpeg_ops(openslide_t *osr,
			     int32_t count,
			     struct _openslide_jpeg_fragment **fragments) {
  //  g_debug("count: %d", count);
  //  for (int32_t i = 0; i < count; i++) {
    //    struct _openslide_jpeg_fragment *frag = fragments[i];
    //    g_debug("%d: file: %p, x: %d, y: %d, z: %d",
    //	    i, (void *) frag->f, frag->x, frag->y, frag->z);
  //  }

  if (osr == NULL) {
    // free now and return
    for (int32_t i = 0; i < count; i++) {
      fclose(fragments[i]->f);
      g_free(fragments[i]->mcu_starts);
      g_slice_free(struct _openslide_jpeg_fragment, fragments[i]);
    }
    g_free(fragments);
    return;
  }

  g_assert(osr->data == NULL);


  // allocate private data
  struct jpegops_data *data = g_slice_new0(struct jpegops_data);
  osr->data = data;

  // load all jpegs (assume all are useful)
  data->jpeg_count = count;
  data->all_jpegs = g_new0(struct one_jpeg, count);
  for (int32_t i = 0; i < data->jpeg_count; i++) {
    g_debug("init JPEG %d", i);
    init_one_jpeg(&data->all_jpegs[i], fragments[i]);
  }

  // create map from width to layers, using the fragments
  GHashTable *width_to_layer_map = create_width_to_layer_map(count,
							     fragments,
							     data->all_jpegs);

  //  g_hash_table_foreach(width_to_layer_map, print_wlmap_entry, NULL);

  // delete all the fragments
  for (int32_t i = 0; i < count; i++) {
    g_slice_free(struct _openslide_jpeg_fragment, fragments[i]);
  }
  g_free(fragments);

  // get sorted keys
  GList *layer_keys = NULL;
  g_hash_table_foreach(width_to_layer_map, get_keys, &layer_keys);
  layer_keys = g_list_sort(layer_keys, width_compare);

  //  g_debug("number of keys: %d", g_list_length(layer_keys));


  // populate the layer_count
  osr->layer_count = g_hash_table_size(width_to_layer_map);

  // load into data array
  data->layers = g_new(struct layer, g_hash_table_size(width_to_layer_map));
  GList *tmp_list = layer_keys;

  int i = 0;

  //  g_debug("copying sorted layers");
  while(tmp_list != NULL) {
    // get a key and value
    struct layer *l = g_hash_table_lookup(width_to_layer_map, tmp_list->data);

    //    print_wlmap_entry(tmp_list->data, l, NULL);

    // copy
    struct layer *dest = data->layers + i;
    *dest = *l;    // shallow copy

    // manually free some things, because of that shallow copy
    g_hash_table_steal(width_to_layer_map, tmp_list->data);
    int64_free(tmp_list->data);  // key
    g_slice_free(struct layer, l); // shallow deletion of layer

    // consume the head and continue
    tmp_list = g_list_delete_link(tmp_list, tmp_list);
    i++;
  }

  // init cache
  data->cache = _openslide_cache_create(_OPENSLIDE_USEFUL_CACHE_SIZE);

  // unref the hash table
  g_hash_table_unref(width_to_layer_map);

  // init background thread for finding restart markers
  data->restart_marker_mutex = g_mutex_new();
  data->restart_marker_thread = g_thread_create(restart_marker_thread_func,
						data,
						TRUE,
						NULL);

  // for debugging
  /*
  g_thread_join(data->restart_marker_thread);
  verify_mcu_starts(data);
  */

  // set ops
  osr->ops = &jpeg_ops;
}
