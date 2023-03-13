/*
 *  OpenSlide, a library for reading whole slide image files
 *
 *  Copyright (c) 2007-2015 Carnegie Mellon University
 *  Copyright (c) 2011 Google, Inc.
 *  Copyright (c) 2022 Benjamin Gilbert
 *  All rights reserved.
 *
 *  OpenSlide is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, version 2.1.
 *
 *  OpenSlide is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with OpenSlide. If not, see
 *  <http://www.gnu.org/licenses/>.
 *
 */

/*
 * DICOM (.dcm) support
 */

#include "openslide-private.h"
#include "openslide-decode-jpeg.h"
#include "openslide-decode-sqlite.h"
#include "openslide-hash.h"

#include <glib.h>
#include <glib-object.h>
#include <gio/gio.h>
#include <string.h>
#include <errno.h>

#include <dicom/dicom.h>

#define FREEF(F, V) { \
    if (V) { \
        (F)((V)); \
        (V) = NULL; \
    } \
}

/* The SOP we check for.
 */
#define VLWholeSlideMicroscopyImageStorage "1.2.840.10008.5.1.4.1.1.77.1.6"

/*
#define DEBUG
 */

enum image_format {
  FORMAT_UNKNOWN,
  FORMAT_JPEG,
};

struct dicom_ops_data {
  char *dirname;
  int tile_size;
};

struct dicom_file {
  char *filename;

  GMutex lock;
  DcmFilehandle *filehandle;
  DcmDataSet *metadata;
  DcmDataSet *meta;
  DcmBOT *bot;
};

struct level {
  struct _openslide_level base;
  struct _openslide_grid *grid;

  enum image_format image_format;
  int image_width;
  int image_height;
  int tile_w;
  int tile_h;
  guint32 num_frames;
  int tiles_across;
  int tiles_down;

  struct dicom_file *file;
};

struct associated {
  struct _openslide_associated_image base;

  const char *name;
  struct dicom_file *file;
};

// a set of allowed image types for a class of image
struct allowed_types { 
  const char ***types;
  int n_types;
};

// the ImageTypes we allow for pyr levels
static const char *original_types[4] = {
  "ORIGINAL", "PRIMARY", "VOLUME", "NONE"
};
static const char *resampled_types[4] = {
  "DERIVED", "PRIMARY", "VOLUME", "RESAMPLED"
};
static const char **level_type_strings[] = {
  original_types,
  resampled_types
};

static const struct allowed_types level_types = {
  level_type_strings, 
  sizeof(level_type_strings) / sizeof(level_type_strings[0])
};

// the ImageTypes we allow for associated images
static const char *label_types[4] = {
  "ORIGINAL", "PRIMARY", "LABEL", "NONE",
};
static const char *overview_types[4] = {
  "ORIGINAL", "PRIMARY", "OVERVIEW", "NONE",
};
static const char **associated_type_strings[] = {
  label_types,
  overview_types
};
static const struct allowed_types associated_types = {
  associated_type_strings, 
  sizeof(associated_type_strings) / sizeof(associated_type_strings[0])
};

G_DEFINE_AUTOPTR_CLEANUP_FUNC(DcmFilehandle, dcm_filehandle_destroy)
G_DEFINE_AUTOPTR_CLEANUP_FUNC(DcmDataSet, dcm_dataset_destroy)

static void set_gerror_from_dcm_error(GError **err, DcmError **dcm_error)
{
  char *msg = g_strdup_printf("libdicom %s: %s - %s",
    dcm_error_code_str(dcm_error_code(*dcm_error)),
    dcm_error_summary(*dcm_error),
    dcm_error_message(*dcm_error));
  g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED, "%s", msg);
  g_free(msg);
  dcm_error_clear(dcm_error);
}

#ifdef DEBUG
static void print_file(struct dicom_file *f) {
  printf("file:\n" );
  printf("  filename = %s\n", f->filename);
  printf("  filehandle = %p\n", f->filehandle);
  printf("  metadata = %p\n", f->metadata);
  printf("  meta = %p\n", f->meta);
  printf("  bot = %p\n", f->bot);
}

static void print_level(struct level *l) {
  printf("level:\n" );
  print_file(l->file);
  printf("  base.downsample = %g\n", l->base.downsample);
  printf("  grid = %p\n", l->grid);
  printf("  format = %d\n", l->image_format);
  printf("  image_width = %d\n", l->image_width);
  printf("  image_height = %d\n", l->image_height);
  printf("  tile_w = %d\n", l->tile_w);
  printf("  tile_h = %d\n", l->tile_h);
  printf("  num_frames = %d\n", l->num_frames);
  printf("  tiles_across = %d\n", l->tiles_across);
  printf("  tiles_down = %d\n", l->tiles_down);
}

static void print_frame(DcmFrame *frame) {
  printf("value = %p\n", dcm_frame_get_value(frame));
  printf("length = %u bytes\n", dcm_frame_get_length(frame));
  printf("rows = %u\n", dcm_frame_get_rows(frame));
  printf("columns = %u\n", dcm_frame_get_columns(frame));
  printf("samples per pixel = %u\n",
         dcm_frame_get_samples_per_pixel(frame));
  printf("bits allocated = %u\n", dcm_frame_get_bits_allocated(frame));
  printf("bits stored = %u\n", dcm_frame_get_bits_stored(frame));
  printf("high bit = %u\n", dcm_frame_get_high_bit(frame));
  printf("pixel representation = %u\n",
         dcm_frame_get_pixel_representation(frame));
  printf("planar configuration = %u\n",
         dcm_frame_get_planar_configuration(frame));
  printf("photometric interpretation = %s\n",
         dcm_frame_get_photometric_interpretation(frame));
  printf("transfer syntax uid = %s\n",
         dcm_frame_get_transfer_syntax_uid(frame));
}
#endif /*DEBUG*/

static bool dicom_detect(const char *filename,
                         struct _openslide_tifflike *tl, 
                         GError **err) {
  // reject TIFFs
  if (tl) {
    g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED,
                "Is a TIFF file");
    return false;
  }

  // verify existence
  GError *tmp_err = NULL;
  if (!_openslide_fexists(filename, &tmp_err)) {
    if (tmp_err != NULL) {
      g_propagate_prefixed_error(err, tmp_err, "Testing whether file exists: ");
    } else {
      g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED,
                  "File does not exist");
    }
    return false;
  }

  // should be able to open as a DICOM
  DcmError *dcm_error = NULL;
  g_autoptr(DcmFilehandle) filehandle = 
    dcm_filehandle_create_from_file(&dcm_error, filename);
  if (filehandle == NULL) {
    set_gerror_from_dcm_error(err, &dcm_error);
    return false;
  }

  // should be able to open as a DICOM
  g_autoptr(DcmDataSet) meta = 
    dcm_filehandle_read_file_meta(&dcm_error, filehandle);
  if (!meta) {
    set_gerror_from_dcm_error(err, &dcm_error);
    return false;
  }

  return true;
}

static void dicom_file_destroy(struct dicom_file *f) {
  FREEF(dcm_filehandle_destroy, f->filehandle);
  FREEF(dcm_dataset_destroy, f->meta);
  FREEF(dcm_dataset_destroy, f->metadata);
  FREEF(dcm_bot_destroy, f->bot);
  g_mutex_clear(&f->lock);
  g_free(f->filename);
  g_slice_free(struct dicom_file, f);
}

static bool get_tag_int(DcmDataSet *dataset, const char *keyword, int *result) {
  uint32_t tag = dcm_dict_tag_from_keyword(keyword);
  DcmElement *element = dcm_dataset_get(NULL, dataset, tag);
  int64_t value;
  if (!dcm_element_get_value_integer(NULL, element, 0, &value)) {
    return false;
  }
  *result = value;
  return true;
}

static bool get_tag_str(DcmDataSet *dataset, 
                        const char *keyword, 
                        int index, 
                        const char **result) {
  uint32_t tag = dcm_dict_tag_from_keyword(keyword);
  DcmElement *element = dcm_dataset_get(NULL, dataset, tag);
  return dcm_element_get_value_string(NULL, element, index, result);
}

static struct dicom_file *dicom_file_new(char *filename, GError **err) {
  struct dicom_file *f = g_slice_new0(struct dicom_file);
  f->filename = g_strdup(filename);

  // should be able to open as a DICOM
  DcmError *dcm_error = NULL;
  f->filehandle = dcm_filehandle_create_from_file(&dcm_error, f->filename);
  if (!f->filehandle) {
    set_gerror_from_dcm_error(err, &dcm_error);
    dicom_file_destroy(f);
    return NULL;
  }
  g_mutex_init(&f->lock);

  f->meta = dcm_filehandle_read_file_meta(&dcm_error, 
                                          f->filehandle);
  if (!f->meta) {
    set_gerror_from_dcm_error(err, &dcm_error);
    dicom_file_destroy(f);
    return NULL;
  }

  const char *sop;
  get_tag_str(f->meta, "MediaStorageSOPClassUID", 0, &sop);
  if (strcmp(sop, VLWholeSlideMicroscopyImageStorage) != 0) {
    g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED,
                "Not a WSI DICOM");
    dicom_file_destroy(f);
    return NULL;
  }

  return f;
}

static void level_destroy(struct level *l) {
  _openslide_grid_destroy(l->grid);
  FREEF(dicom_file_destroy, l->file);
  g_slice_free(struct level, l);
}

static void destroy(openslide_t *osr) {
  struct dicom_ops_data *data = osr->data;
  g_free(data->dirname);
  g_slice_free(struct dicom_ops_data, data);

  for (int32_t i = 0; i < osr->level_count; i++) {
    level_destroy((struct level *) osr->levels[i]);
  }
  g_free(osr->levels);
}

static bool read_tile(openslide_t *osr,
                      cairo_t *cr,
                      struct _openslide_level *level,
                      int64_t tile_col, int64_t tile_row,
                      void *arg G_GNUC_UNUSED,
                      GError **err) {
  struct dicom_ops_data *data = osr->data;
  struct level *l = (struct level *) level;
  int32_t tile_size = data->tile_size;

#ifdef DEBUG
  printf("read_tile: tile_col = %" PRIu64 ", tile_row = %" PRIu64 "\n", 
         tile_col, 
         tile_row);
  printf("read_tile level:\n");
  print_level(l);
#endif /*DEBUG*/

  // cache
  g_autoptr(_openslide_cache_entry) cache_entry = NULL;
  guint32 *tiledata = _openslide_cache_get(osr->cache,
                                            level, tile_col, tile_row,
                                            &cache_entry);
  if (!tiledata) {
    g_auto(_openslide_slice) box =
      _openslide_slice_alloc(tile_size * tile_size * 4);
    guint32 frame_number = 1 + tile_col + l->tiles_across * tile_row;
    if (frame_number < 1 || frame_number > l->num_frames) {
      g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED,
                  "Frame number out of range");
      return false;
    }

    g_mutex_lock(&l->file->lock);

    DcmError *dcm_error = NULL;
    DcmFrame *frame = dcm_filehandle_read_frame(&dcm_error,
                                                l->file->filehandle, 
                                                l->file->metadata,
                                                l->file->bot, 
                                                frame_number);

    g_mutex_unlock(&l->file->lock);

    if (frame == NULL) {
      set_gerror_from_dcm_error(err, &dcm_error);
      return false;
    }

    const char *frame_value = dcm_frame_get_value(frame);
    uint32_t frame_length = dcm_frame_get_length(frame);

#ifdef DEBUG
    print_frame(frame);
#endif /*DEBUG*/

    if (!_openslide_jpeg_decode_buffer(frame_value,
                                       frame_length,
                                       box.p,
                                       dcm_frame_get_columns(frame),
                                       dcm_frame_get_rows(frame),
                                       err)) {
      return false;
    }

    // clip, if necessary
    if (!_openslide_clip_tile(box.p,
                              tile_size, tile_size,
                              l->base.w - tile_col * tile_size,
                              l->base.h - tile_row * tile_size,
                              err)) {
      return false;
    }

    // put it in the cache
    tiledata = _openslide_slice_steal(&box);
    _openslide_cache_put(osr->cache,
			 level, tile_col, tile_row,
			 tiledata, tile_size * tile_size * 4,
			 &cache_entry);
  }

  // draw it
  g_autoptr(cairo_surface_t) surface =
    cairo_image_surface_create_for_data((unsigned char *) tiledata,
                                        CAIRO_FORMAT_ARGB32,
                                        tile_size, tile_size, tile_size * 4);
  cairo_set_source_surface(cr, surface, 0, 0);
  cairo_paint(cr);

  return true;
}

static bool paint_region(openslide_t *osr, 
                         cairo_t *cr,
                         int64_t x, 
                         int64_t y,
                         struct _openslide_level *level,
                         int32_t w, 
                         int32_t h,
                         GError **err) {
  struct dicom_ops_data *data = osr->data;
  struct level *l = (struct level *) level;

#ifdef DEBUG
  printf("paint_region: x = %" PRId64 ", y = %" PRId64 ", w = %d, h = %d\n", 
         x, y, w, h);
  printf("paint_region level:\n");
  print_level(l);
#endif /*DEBUG*/

  return _openslide_grid_paint_region(l->grid, 
                                      cr, 
                                      data,
                                      x / l->base.downsample,
                                      y / l->base.downsample,
                                      level, 
                                      w, 
                                      h,
                                      err);
}

static const struct _openslide_ops dicom_ops = {
  .paint_region = paint_region,
  .destroy = destroy,
};

static bool is_type(struct dicom_file *f, const struct allowed_types *types) {
  // ImageType must be one of the combinations we accept
  int match = -1;
  for (int i = 0; i < types->n_types; i++) {
    bool found_difference = false;

    for (int j = 0; j < 4; j++) {
      const char *type;
      if (!get_tag_str(f->metadata, "ImageType", j, &type)) {
        return false;
      }

      if (strcmp(type, types->types[i][j]) != 0) {
        found_difference = true;
        break;
      }
    }

    if (!found_difference) {
      match = i;
      break;
    }
  }

  return match >= 0;
}

static bool associated_get_argb_data(struct _openslide_associated_image *_img,
                                uint32_t *dest,
                                GError **err) {
  struct associated *a = (struct associated *) _img;

  DcmError *dcm_error = NULL;
  DcmFrame *frame = dcm_filehandle_read_frame(&dcm_error,
                                              a->file->filehandle, 
                                              a->file->metadata,
                                              a->file->bot, 
                                              1);
  if (frame == NULL) {
    set_gerror_from_dcm_error(err, &dcm_error);
    return false;
  }

#ifdef DEBUG
  print_frame(frame);
#endif /*DEBUG*/

  return _openslide_jpeg_decode_buffer(dcm_frame_get_value(frame),
                                       dcm_frame_get_length(frame),
                                       dest,
                                       a->base.w,
                                       a->base.h,
                                       err);
}

static void associated_destroy(struct _openslide_associated_image *_img) {
  struct associated *a = (struct associated *) _img;
  FREEF(dicom_file_destroy, a->file);
  g_slice_free(struct associated, a);
}

static const struct _openslide_associated_image_ops dicom_associated_ops = {
  .get_argb_data = associated_get_argb_data,
  .destroy = associated_destroy
};

static bool read_whole_file(struct dicom_file *f) {
  // try to read the rest of the dicom file 
  if (!f->metadata) {
    f->metadata = dcm_filehandle_read_metadata(NULL, f->filehandle);
  }
  if (!f->metadata) {
    return false;
  }

  if (!f->bot) {
    f->bot = dcm_filehandle_read_bot(NULL, f->filehandle, f->metadata);
  }
  if (!f->bot) {
    /* Try to build the BOT instead.
     */
    f->bot = dcm_filehandle_build_bot(NULL, f->filehandle, f->metadata);
  }
  if (!f->bot) {
    return false;
  }

  return true;
}

static struct associated *associated_new(struct dicom_file *f) {
  // we use BOT to get the label and overview frames
  if (!read_whole_file(f)) { 
    return NULL;
  }

  // ImageType must be one of the combinations we accept
  if (!is_type(f, &associated_types)) {
    return NULL;
  }

  const char *name;
  if (!get_tag_str(f->metadata, "ImageType", 2, &name)) {
    return NULL;
  }

  int image_width;
  int image_height;
  if (!get_tag_int(f->metadata, "TotalPixelMatrixColumns", &image_width) ||
    !get_tag_int(f->metadata, "TotalPixelMatrixRows", &image_height)) {
    return NULL;
  }

  // set this once we know we will succeed ... this passes ownership to the
  // level
  struct associated *a = g_slice_new0(struct associated);
  a->file = f;
  a->name = name;
  a->base.ops = &dicom_associated_ops;
  a->base.w = image_width;
  a->base.h = image_height;

#ifdef DEBUG
  printf("associated_new: %s\n", a->name);
#endif /*DEBUG*/

  return a;
}

static struct level *level_new(struct dicom_file *f) {
  if (!read_whole_file(f)) { 
    return NULL;
  }

  struct level *l = g_slice_new0(struct level);

  if (!get_tag_int(f->metadata, "TotalPixelMatrixColumns", &l->image_width) ||
    !get_tag_int(f->metadata, "TotalPixelMatrixRows", &l->image_height) ||
    !get_tag_int(f->metadata, "Columns", &l->tile_w) ||
    !get_tag_int(f->metadata, "Rows", &l->tile_h)) {
    level_destroy(l);
    return NULL;
  }

  // ImageType must be one of the combinations we accept
  if (!is_type(f, &level_types)) {
    level_destroy(l);
    return NULL;
  }

  // we only allow square tiles
  if (l->tile_w != l->tile_h) {
    level_destroy(l);
    return NULL;
  }

  // set this once we know we will succeed ... this passes ownership to the
  // level
  l->file = f;
  l->num_frames = dcm_bot_get_num_frames(l->file->bot);
  l->base.w = l->image_width;
  l->base.h = l->image_height;
  l->base.tile_w = l->tile_w;
  l->base.tile_h = l->tile_h;
  l->tiles_across = (l->base.w / l->tile_w) + !!(l->base.w % l->tile_w);
  l->tiles_down = (l->base.h / l->tile_h) + !!(l->base.h % l->tile_h);

  return l;
}

static int find_levels(gpointer key, gpointer value, gpointer user_data) {
  const char *filename = (const char *) key;
  struct dicom_file *f = (struct dicom_file *) value;
  GHashTable *level_hash = (GHashTable *) user_data;

  struct level *l = level_new(f);
  if (l) {
    g_hash_table_insert(level_hash, (char *) filename, l);
    return true;
  }
  else {
    return false;
  }
}

static void find_largest(gpointer key G_GNUC_UNUSED,
                         gpointer value, 
                         gpointer user_data) {
  struct level *l = (struct level *) value;
  struct level **largest = (struct level **) user_data;
  struct dicom_file *f = l->file;

  if (!get_tag_int(f->metadata, "TotalPixelMatrixColumns", &l->image_width)) {
    return;
  }

  if (!*largest ||
    l->image_width > (*largest)->image_width) {
    *largest = l;
  }
}

static int remove_bad_level(gpointer key G_GNUC_UNUSED, 
                            gpointer value, 
                            gpointer user_data) {
  struct level *l = (struct level *) value;
  struct dicom_file *f = l->file;
  const char *slide_id = (const char *) user_data;

  // true to remove this file
  const char *this_slide_id;
  return !get_tag_str(f->metadata, "SeriesInstanceUID", 0, &this_slide_id) ||
         strcmp(slide_id, this_slide_id) != 0;
}

static int remove_bad_dicom(gpointer key G_GNUC_UNUSED, 
                            gpointer value, 
                            gpointer user_data) {
  struct dicom_file *f = (struct dicom_file *) value;
  const char *slide_id = (const char *) user_data;

  // we might not have read all the metadta for this dicom yet
  if (!f->metadata) {
    f->metadata = dcm_filehandle_read_metadata(NULL, f->filehandle);
  }
  if (!f->metadata) {
    // a broken file, so true to remove it
    return true;
  }

  // true to remove this file
  const char *this_slide_id;
  return !get_tag_str(f->metadata, "SeriesInstanceUID", 0, &this_slide_id) ||
         strcmp(slide_id, this_slide_id) != 0;
}

static void set_downsample(gpointer key G_GNUC_UNUSED, 
                           gpointer value, 
                           gpointer user_data) {
  struct level *l = (struct level *) value;
  const struct level *largest = (struct level *) user_data;

  // need to compute downsample now we have largest
  int downsample = largest->image_width / l->image_width;
  l->base.downsample = downsample;
}

static int add_level_array(gpointer key G_GNUC_UNUSED, 
                           gpointer value, 
                           gpointer user_data) {
  struct level *l = (struct level *) value;
  GPtrArray *level_array = (GPtrArray *) user_data;

  g_ptr_array_add(level_array, l);

  return true;
}

#ifdef DEBUG
static void print_file_from_hash(gpointer key G_GNUC_UNUSED, 
                                 gpointer value, 
                                 gpointer user_data G_GNUC_UNUSED) {
  struct dicom_file *f = (struct dicom_file *) value;

  print_file(f);
}

static void print_level_from_hash(gpointer key G_GNUC_UNUSED, 
                                  gpointer value, 
                                  gpointer user_data G_GNUC_UNUSED) {
  struct level *l = (struct level *) value;

  print_level(l);
}
#endif /*DEBUG*/

static int find_associated(gpointer key G_GNUC_UNUSED, 
                            gpointer value, 
                            gpointer user_data) {
  struct dicom_file *f = (struct dicom_file *) value;
  openslide_t *osr = (openslide_t *) user_data; 

  struct associated *a = associated_new(f);
  if (a) {
    // FIXME ... this will use LABEL and OVERVIEW ... is this OK?
    g_hash_table_insert(osr->associated_images, g_strdup(a->name), a);
    return true;
  }
  else {
    return false;
  }
}

static gint compare_level_downsamples(const void *a, const void *b) {
  const struct level *aa = *((const struct level **) a);
  const struct level *bb = *((const struct level **) b);

  return aa->base.downsample - bb->base.downsample;
}

static void make_grid(gpointer key G_GNUC_UNUSED, 
                      gpointer value, 
                      gpointer user_data) {
  struct level *l = (struct level *) value;
  openslide_t *osr = (openslide_t *) user_data; 

  l->grid = _openslide_grid_create_simple(osr,
                                          l->tiles_across, 
                                          l->tiles_down,
                                          l->tile_w, 
                                          l->tile_h,
                                          read_tile);
}

static void add_properties(openslide_t *osr, const struct level *l) {
  // why not
  double mmpp = l->image_width;
  g_hash_table_insert(osr->properties,
                      g_strdup(OPENSLIDE_PROPERTY_NAME_MPP_X),
                      _openslide_format_double(mmpp * 1000));
  g_hash_table_insert(osr->properties,
                      g_strdup(OPENSLIDE_PROPERTY_NAME_MPP_Y),
                      _openslide_format_double(mmpp * 1000));
}

static void file_hash(const struct dicom_file *f,
                      struct _openslide_hash *quickhash1) {
  uint32_t n = dcm_dataset_count(f->metadata);
  for(uint32_t i = 0; i < n; i++) {
      DcmElement *element = dcm_dataset_get(NULL, f->metadata, i);
      if (element) {
        uint32_t tag = dcm_element_get_tag(element);
        _openslide_hash_data(quickhash1, &tag, sizeof(tag));
        // need to get the vr class
        //_openslide_hash_data(quickhash1, data, datalen);
      }
  }
}

static bool dicom_open(openslide_t *osr, 
                       const char *filename,
                       struct _openslide_tifflike *tl G_GNUC_UNUSED,
                       struct _openslide_hash *quickhash1, 
                       GError **err) {
  g_autofree char *dirname = g_path_get_dirname(filename);

  g_autoptr(GDir) dir = g_dir_open(dirname, 0, err);
  if (!dir) {
    return false;
  }

  g_autoptr(GHashTable) dicom_file_hash =
    g_hash_table_new_full(g_str_hash, 
                          g_str_equal, 
                          g_free,
                          (GDestroyNotify) dicom_file_destroy);

  // open all DICOM files that look like parts of a slide image and 
  // get the file metadata
  const char *name;
  while ((name = g_dir_read_name(dir))) {
    char *filename = g_build_path("/", dirname, name, NULL);

#ifdef DEBUG
    printf("trying to open: %s ...\n", filename);
#endif /*DEBUG*/
    GError *local_error = NULL;
    struct dicom_file *f = dicom_file_new(filename, &local_error);
    if (!f) {
      g_free(filename);
#ifdef DEBUG
      printf( "open failed: %s\n", local_error->message);
      g_error_free(local_error);
#endif /*DEBUG*/
      continue;
    }

    g_hash_table_insert(dicom_file_hash, filename, f);
  }

#ifdef DEBUG
  printf("found WSI DICOM files:\n");
  g_hash_table_foreach(dicom_file_hash, print_file_from_hash, NULL);
#endif /*DEBUG*/

  // pull out the subset of DICOM files that look like pyramid levels
  g_autoptr(GHashTable) level_hash =
    g_hash_table_new_full(g_str_hash,
                          g_str_equal, 
                          g_free,
                          (GDestroyNotify) level_destroy);
  g_hash_table_foreach_steal(dicom_file_hash, find_levels, level_hash);

#ifdef DEBUG
  printf("found pyr levels DICOM files:\n");
  g_hash_table_foreach(level_hash, print_level_from_hash, NULL);
#endif /*DEBUG*/

  // we can have several slides in one directory -- find the largest pyramid
  // layer and pick that as the slide image we are opening
  const struct level *largest = NULL;
  g_hash_table_foreach(level_hash, find_largest, &largest);
  if (!largest) {
    return false;
  }

  const char *slide_id;
  if (!get_tag_str(largest->file->metadata, 
                   "SeriesInstanceUID", 
                   0, 
                   &slide_id)) {
    return false;
  }

  // throw away all files which don't have this slide_id
  g_hash_table_foreach_remove(level_hash, remove_bad_level, (char *) slide_id);
  g_hash_table_foreach_remove(dicom_file_hash, 
                              remove_bad_dicom, 
                              (char *) slide_id);

  // compute the downsample for each sublevel
  g_hash_table_foreach(level_hash, set_downsample, (struct level *) largest);

  // make the tile cache
  g_hash_table_foreach(level_hash, make_grid, osr);

#ifdef DEBUG
  printf("\nfinal pyr levels:\n");
  g_hash_table_foreach(level_hash, print_level_from_hash, NULL);
#endif /*DEBUG*/

  // now sort levels by downsample to level array
  int32_t level_count = g_hash_table_size(level_hash);
  if (level_count == 0) {
    g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED,
                "Couldn't find any tiles");
    return false;
  }

  g_autoptr(GPtrArray) level_array = 
    g_ptr_array_new_full(level_count, (GDestroyNotify) level_destroy);
  g_hash_table_foreach_steal(level_hash, add_level_array, level_array);
  g_ptr_array_sort(level_array, compare_level_downsamples); 

  // steal any associated images from the remaining set
  g_hash_table_foreach_steal(dicom_file_hash, find_associated, osr);

  // take props from the largest pyr layer
  add_properties(osr, largest);

  // perhaps add all the top-level metadata we've loaded for the largest level
  file_hash(largest->file, quickhash1);

  g_assert(osr->data == NULL);
  g_assert(osr->levels == NULL);

  // build ops data
  struct dicom_ops_data *data = g_slice_new0(struct dicom_ops_data);
  data->dirname = g_strdup(dirname);
  gsize count;

  /* g_ptr_array_steal() is only available in recent glib ... hack around 
   * this for older versions.
   */
#if GLIB_CHECK_VERSION(2, 64, 0)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
  osr->levels = (struct _openslide_level **) g_ptr_array_steal(level_array, 
                                                               &count);
#pragma GCC diagnostic pop
#else
  osr->levels = (struct _openslide_level **) level_array->pdata;
  level_array->pdata = NULL;
  level_array->len = 0;
#endif

  osr->level_count = count;
  data->tile_size = osr->levels[0]->tile_w;
  osr->data = data;
  osr->ops = &dicom_ops;

#ifdef DEBUG
  printf("sorted levels:\n");
  for (gsize i = 0; i < count; i++)
    printf("%zd: downsample = %g\n", i, osr->levels[i]->downsample);
#endif /*DEBUG*/

  return true;
}

const struct _openslide_format _openslide_format_dicom = {
  .name = "dicom",
  .vendor = "dicom",
  .detect = dicom_detect,
  .open = dicom_open,
};
