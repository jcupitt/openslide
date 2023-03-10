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
 * Sakura (svslide) support
 *
 * quickhash comes from a selection of metadata fields, the binary header
 * blob, and the lowest-resolution level
 *
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
  DcmFilehandle *filehandle;
  DcmDataSet *metadata;
  DcmDataSet *file_metadata;
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

  struct dicom_file *file;
};

struct associated_image {
  struct _openslide_associated_image base;
  char *filename;
  DcmFilehandle *filehandle;
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
  g_autoptr(DcmFilehandle) filehandle = 
    dcm_filehandle_create_from_file(NULL, filename);
  if (filehandle == NULL) {
    g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED,
                "File is not a DICOM");
    return false;
  }

  return true;
}

static void dicom_file_destroy(struct dicom_file *f) {
  FREEF(dcm_filehandle_destroy, f->filehandle);
  FREEF(dcm_dataset_destroy, f->file_metadata);
  FREEF(dcm_dataset_destroy, f->metadata);
  FREEF(dcm_bot_destroy, f->bot);
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

static struct dicom_file *dicom_file_new(char *filename) {
  struct dicom_file *f = g_slice_new0(struct dicom_file);
  f->filename = g_strdup(filename);

  f->filehandle = dcm_filehandle_create_from_file(NULL, f->filename);
  if (!f->filehandle) {
    dicom_file_destroy(f);
    return NULL;
  }

  f->file_metadata = dcm_filehandle_read_file_metadata(NULL, f->filehandle);
  if (!f->file_metadata) {
    dicom_file_destroy(f);
    return NULL;
  }

  const char *sop;
  get_tag_str(f->file_metadata, "MediaStorageSOPClassUID", 0, &sop);
  if (strcmp(sop, "VLWholeSlideMicroscopyImageStorage") != 0) {
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
                      void *arg,
                      GError **err) {
  struct dicom_ops_data *data = osr->data;
  struct level *l = (struct level *) level;
  int32_t tile_size = data->tile_size;

  // cache
  g_autoptr(_openslide_cache_entry) cache_entry = NULL;
  uint32_t *tiledata = _openslide_cache_get(osr->cache,
                                            level, tile_col, tile_row,
                                            &cache_entry);
  if (!tiledata) {
    g_auto(_openslide_slice) box =
      _openslide_slice_alloc(tile_size * tile_size * 4);

    // read tile
    /*
    GError *tmp_err = NULL;
    if (!read_image(box.p, tile_col, tile_row, l->base.downsample,
                    data->focal_plane, tile_size, stmt, &tmp_err)) {
      if (g_error_matches(tmp_err, OPENSLIDE_ERROR,
                          OPENSLIDE_ERROR_NO_VALUE)) {
        // no such tile
        g_clear_error(&tmp_err);
        return true;
      } else {
        g_propagate_error(err, tmp_err);
        return false;
      }

      dcm_get_frame
      then call _openslide_jpeg_decompress,
    }
     */

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

/*
static bool get_associated_image_data(struct _openslide_associated_image *_img,
                                      uint32_t *dest,
                                      GError **err) {
  struct associated_image *img = (struct associated_image *) _img;

  // read data
  const void *buf;
  int buflen;
  // read_image_from_dicom

  // decode it
  return _openslide_jpeg_decode_buffer(buf, buflen, dest,
                                       img->base.w, img->base.h, err);
}

static void destroy_associated_image(struct _openslide_associated_image *_img) {
  struct associated_image *img = (struct associated_image *) _img;

  g_free(img->filename);
  if(img->filehandle) {
    dcm_filehandle_destroy(img->filehandle);
    img->filehandle = NULL;
  }

  g_slice_free(struct associated_image, img);
}

static const struct _openslide_associated_image_ops dicom_associated_ops = {
  .get_argb_data = get_associated_image_data,
  .destroy = destroy_associated_image,
};

static bool add_associated_image(openslide_t *osr,
                                 const char *filename,
                                 const char *name,
                                 GError **err) {
  // read data
  const void *buf;
  int buflen;

  // read dimensions from JPEG header
  int32_t w, h;
  if (!_openslide_jpeg_decode_buffer_dimensions(buf, buflen, &w, &h, err)) {
    return false;
  }

  // create struct
  struct associated_image *img = g_slice_new0(struct associated_image);
  img->base.ops = &dicom_associated_ops;
  img->base.w = w;
  img->base.h = h;
  img->filename = g_strdup(filename);

  // add it
  g_hash_table_insert(osr->associated_images, g_strdup(name), img);

  return true;
}
 */

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

  return match > 0;
}

static struct level *level_new(struct dicom_file *f) {
  // try to read the rest of the dicom file as a pyr level
  f->metadata = dcm_filehandle_read_metadata(NULL, f->filehandle);
  if (!f->metadata) {
    return NULL;
  }

  f->bot = dcm_filehandle_read_bot(NULL, f->filehandle, f->metadata);
  if (!f->bot) {
    /* Try to build the BOT instead.
     */
    f->bot = dcm_filehandle_build_bot(NULL, f->filehandle, f->metadata);
  }
  if (!f->bot) {
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
  l->base.w = l->image_width;
  l->base.h = l->image_height;
  l->base.tile_w = l->tile_w;
  l->base.tile_h = l->tile_h;

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

static void find_largest(gpointer key, gpointer value, gpointer user_data) {
  const char *filename = (const char *) key;
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

static int remove_bad_level(gpointer key, gpointer value, gpointer user_data) {
  const char *filename = (const char *) key;
  struct level *l = (struct level *) value;
  struct dicom_file *f = l->file;
  const char *slide_id = (const char *) user_data;

  const char *this_slide_id;

  return get_tag_str(f->metadata, "SeriesInstanceUID", 0, &this_slide_id) &&
         strcmp(slide_id, this_slide_id) != 0;
}

static int remove_bad_dicom(gpointer key, gpointer value, gpointer user_data) {
  const char *filename = (const char *) key;
  struct dicom_file *f = (struct dicom_file *) value;
  const char *slide_id = (const char *) user_data;

  const char *this_slide_id;

  return get_tag_str(f->metadata, "SeriesInstanceUID", 0, &this_slide_id) &&
         strcmp(slide_id, this_slide_id) != 0;
}

static void set_downsample(gpointer key, gpointer value, gpointer user_data) {
  const char *filename = (const char *) key;
  struct level *l = (struct level *) value;
  const struct level *largest = (struct level *) user_data;

  // need to compute downsample now we have largest
  int downsample = largest->image_width / l->image_width;
  l->base.downsample = downsample;
}

static int add_level_array(gpointer key, gpointer value, gpointer user_data) {
  const char *filename = (const char *) key;
  struct level *l = (struct level *) value;
  GPtrArray *level_array = (GPtrArray *) user_data;

  g_ptr_array_add(level_array, l);

  return true;
}

/*
static bool find_associated(gpointer key, gpointer value, gpointer user_data) {
  const char *filename = (const char *) key;
  struct dicom_file *f = (struct dicom_file *) value;
  GHashTable *associated_hash = (GHashTable *) user_data;

  struct associated_image *a = associated_new(f);
  if (a) {
    g_hash_table_insert(associated_hash, filename, a);
    return true;
  }
  else {
    return false;
  }
}
 */

static gint compare_level_downsamples(const void *a, const void *b) {
  const struct level *aa = (const struct level *) a;
  const struct level *bb = (const struct level *) bb;

  return aa->base.downsample - bb->base.downsample;
}

static void make_grid(gpointer key, gpointer value, gpointer user_data) {
  const char *filename = (const char *) key;
  struct level *l = (struct level *) value;
  openslide_t *osr = (openslide_t *) user_data; 

  int64_t tiles_across = (l->base.w / l->tile_w) + !!(l->base.w % l->tile_w);
  int64_t tiles_down = (l->base.h / l->tile_h) + !!(l->base.h % l->tile_h);

  l->grid = _openslide_grid_create_simple(osr,
                                          tiles_across, 
                                          tiles_down,
                                          l->tile_w, 
                                          l->tile_h,
                                          read_tile);
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

    struct dicom_file *f = dicom_file_new(filename);
    if (!f) {
      g_free(filename);
      continue;
    }

    g_hash_table_insert(dicom_file_hash, filename, f);
  }

  // pull out the subset of DICOM files that look like pyramid levels
  g_autoptr(GHashTable) level_hash =
    g_hash_table_new_full(g_str_hash,
                          g_str_equal, 
                          g_free,
                          (GDestroyNotify) level_destroy);

  g_hash_table_foreach_steal(dicom_file_hash, find_levels, level_hash);

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

  /*
  // add properties
  add_properties(osr, db, unique_table_name);

  // add associated images
  // errors are non-fatal
  add_associated_image(osr, db, filename, "label",
                       "SELECT Image FROM SVScannedImageDataXPO JOIN "
                       "SVSlideDataXPO ON SVSlideDataXPO.m_labelScan = "
                       "SVScannedImageDataXPO.OID", NULL);
  add_associated_image(osr, db, filename, "macro",
                       "SELECT Image FROM SVScannedImageDataXPO JOIN "
                       "SVSlideDataXPO ON SVSlideDataXPO.m_overviewScan = "
                       "SVScannedImageDataXPO.OID", NULL);
  add_associated_image(osr, db, filename, "thumbnail",
                       "SELECT ThumbnailImage FROM SVHRScanDataXPO JOIN "
                       "SVSlideDataXPO ON SVHRScanDataXPO.ParentSlide = "
                       "SVSlideDataXPO.OID", NULL);
    */

  g_assert(osr->data == NULL);
  g_assert(osr->levels == NULL);

  // build ops data
  struct dicom_ops_data *data = g_slice_new0(struct dicom_ops_data);
  data->dirname = g_strdup(dirname);
  gsize count;
  osr->levels = (struct _openslide_level **) g_ptr_array_steal(level_array, 
                                                               &count);
  osr->level_count = count;
  data->tile_size = osr->levels[0]->tile_w;
  osr->data = data;
  osr->ops = &dicom_ops;

  return true;
}

const struct _openslide_format _openslide_format_dicom = {
  .name = "dicom",
  .vendor = "dicom",
  .detect = dicom_detect,
  .open = dicom_open,
};
