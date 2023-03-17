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
#include "openslide-hash.h"

#include <glib.h>
#include <glib-object.h>
#include <string.h>
#include <errno.h>

#include <dicom/dicom.h>

/*
#define DEBUG
 */

enum image_format {
  FORMAT_UNKNOWN,
  FORMAT_JPEG,
};

struct dicom_ops_data {
  char *dirname;
  int64_t tile_size;
};

typedef struct dicom_file {
  char *filename;

  GMutex lock;
  DcmFilehandle *filehandle;
  DcmDataSet *metadata;
  DcmDataSet *meta;
  DcmBOT *bot;
} DicomFile;

typedef struct level {
  struct _openslide_level base;
  struct _openslide_grid *grid;

  enum image_format image_format;
  int64_t image_width;
  int64_t image_height;
  int64_t tile_w;
  int64_t tile_h;
  guint32 num_frames;
  int64_t tiles_across;
  int64_t tiles_down;

  DicomFile *file;
} Level;

struct associated {
  struct _openslide_associated_image base;

  const char *name;
  DicomFile *file;
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

/* The DICOM UIDs and fields we check.
 */
static const char VLWholeSlideMicroscopyImageStorage[] = 
  "1.2.840.10008.5.1.4.1.1.77.1.6";
static const char SeriesInstanceUID[] = "SeriesInstanceUID";
static const char TotalPixelMatrixColumns[] = "TotalPixelMatrixColumns";
static const char TotalPixelMatrixRows[] = "TotalPixelMatrixRows";
static const char Columns[] = "Columns";
static const char Rows[] = "Rows";

G_DEFINE_AUTOPTR_CLEANUP_FUNC(DcmFilehandle, dcm_filehandle_destroy)
G_DEFINE_AUTOPTR_CLEANUP_FUNC(DcmDataSet, dcm_dataset_destroy)

static void set_gerror_from_dcm_error(GError **err, DcmError **dcm_error)
{
  g_autofree char *msg = g_strdup_printf("libdicom %s: %s - %s",
    dcm_error_code_str(dcm_error_code(*dcm_error)),
    dcm_error_summary(*dcm_error),
    dcm_error_message(*dcm_error));
  g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED, "%s", msg);
  dcm_error_clear(dcm_error);
}

static void set_dcm_error_from_gerror(DcmError **dcm_error, GError **err)
{
  dcm_error_set(dcm_error, DCM_ERROR_CODE_INVALID,
    g_quark_to_string((*err)->domain),
    "%s", (*err)->message);
  g_clear_error(err);
}

#ifdef DEBUG
static void print_file(DicomFile *f) {
  printf("file:\n" );
  printf("  filename = %s\n", f->filename);
  printf("  filehandle = %p\n", f->filehandle);
  printf("  metadata = %p\n", f->metadata);
  printf("  meta = %p\n", f->meta);
  printf("  bot = %p\n", f->bot);
}

static void print_level(Level *l) {
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

static void *dicom_openslide_vfs_open(DcmError **dcm_error, void *client) {
  const char *filename = (const char *) client;

  GError *err = NULL;
  struct _openslide_file *file = _openslide_fopen(filename, &err);
  if (!file) {
    set_dcm_error_from_gerror(dcm_error, &err);
    return NULL;
  }

  return file;
}

static int dicom_openslide_vfs_close(DcmError **dcm_error G_GNUC_UNUSED, 
                                     void *data) {
  struct _openslide_file *file = (struct _openslide_file *) data;
  _openslide_fclose(file);
  return 0;
}

static int64_t dicom_openslide_vfs_read(DcmError **dcm_error G_GNUC_UNUSED, 
                                        void *data, 
                                        char *buffer, 
                                        int64_t length) {
  struct _openslide_file *file = (struct _openslide_file *) data;
  // openslide VFS has no error return for read()
  return _openslide_fread(file, buffer, length);
}

static int64_t dicom_openslide_vfs_seek(DcmError **dcm_error, 
                                        void *data, 
                                        int64_t offset, 
                                        int whence) {
  struct _openslide_file *file = (struct _openslide_file *) data;

  GError *err = NULL;
  if (!_openslide_fseek(file, offset, whence, &err)) {
    set_dcm_error_from_gerror(dcm_error, &err);
    return -1;
  }

  // libdicom uses lseek(2) semantics, so it must always return the new file
  // pointer
  off_t new_position = _openslide_ftell(file, &err);
  if (new_position < 0) {
    set_dcm_error_from_gerror(dcm_error, &err);
  }

  return new_position;
}

static const DcmIO dicom_io_funcs = {
  .open = dicom_openslide_vfs_open,
  .close = dicom_openslide_vfs_close,
  .read = dicom_openslide_vfs_read,
  .seek = dicom_openslide_vfs_seek,
};

static DcmFilehandle *dicom_open_openslide_vfs(const char *filename, 
                                               GError **err) {
  DcmFilehandle *result;
  DcmError *dcm_error = NULL;
  result = dcm_filehandle_create(&dcm_error, 
                                 &dicom_io_funcs, 
                                 (void *) filename);
  if (!result) {
    set_gerror_from_dcm_error(err, &dcm_error);
    return NULL;
  }

  return result;
}

static bool dicom_detect(const char *filename,
                         struct _openslide_tifflike *tl, 
                         GError **err) {
  // reject TIFFs
  if (tl) {
    g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED,
                "Is a TIFF file");
    return false;
  }

  // should be able to open as a DICOM
  g_autoptr(DcmFilehandle) filehandle = dicom_open_openslide_vfs(filename, err);
  if (!filehandle) {
    return false;
  }

  DcmError *dcm_error = NULL;
  g_autoptr(DcmDataSet) meta = 
    dcm_filehandle_read_file_meta(&dcm_error, filehandle);
  if (!meta) {
    set_gerror_from_dcm_error(err, &dcm_error);
    return false;
  }

  return true;
}

static void dicom_file_destroy(DicomFile *f) {
  dcm_filehandle_destroy(f->filehandle);
  dcm_dataset_destroy(f->meta);
  dcm_dataset_destroy(f->metadata);
  dcm_bot_destroy(f->bot);
  g_mutex_clear(&f->lock);
  g_free(f->filename);
  g_slice_free(DicomFile, f);
}

G_DEFINE_AUTOPTR_CLEANUP_FUNC(DicomFile, dicom_file_destroy)

static bool get_tag_int(DcmDataSet *dataset, 
                        const char *keyword, 
                        int64_t *result) {
  uint32_t tag = dcm_dict_tag_from_keyword(keyword);
  DcmElement *element = dcm_dataset_get(NULL, dataset, tag);
  return dcm_element_get_value_integer(NULL, element, 0, result);
}

static bool get_tag_str(DcmDataSet *dataset, 
                        const char *keyword, 
                        int index, 
                        const char **result) {
  uint32_t tag = dcm_dict_tag_from_keyword(keyword);
  DcmElement *element = dcm_dataset_get(NULL, dataset, tag);
  return dcm_element_get_value_string(NULL, element, index, result);
}

static DicomFile *dicom_file_new(char *filename, GError **err) {
  g_autoptr(DicomFile) f = g_slice_new0(DicomFile);

  f->filename = g_strdup(filename);
  f->filehandle = dicom_open_openslide_vfs(filename, err);
  if (!f->filehandle) {
    return NULL;
  }

  g_mutex_init(&f->lock);

  DcmError *dcm_error = NULL;
  f->meta = dcm_filehandle_read_file_meta(&dcm_error, 
                                          f->filehandle);
  if (!f->meta) {
    set_gerror_from_dcm_error(err, &dcm_error);
    return NULL;
  }

  const char *sop;
  get_tag_str(f->meta, "MediaStorageSOPClassUID", 0, &sop);
  if (!g_str_equal(sop, VLWholeSlideMicroscopyImageStorage)) {
    g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED,
                "Not a WSI DICOM");
    return NULL;
  }

  return g_steal_pointer(&f);
}

static void level_destroy(Level *l) {
  _openslide_grid_destroy(l->grid);
  if (l->file) {
    dicom_file_destroy(l->file);
  }
  g_slice_free(Level, l);
}

G_DEFINE_AUTOPTR_CLEANUP_FUNC(Level, level_destroy)

static void destroy(openslide_t *osr) {
  struct dicom_ops_data *data = osr->data;
  g_free(data->dirname);
  g_slice_free(struct dicom_ops_data, data);

  for (int32_t i = 0; i < osr->level_count; i++) {
    level_destroy((Level *) osr->levels[i]);
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
  Level *l = (Level *) level;
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
                  "Frame number out of range %d - %d",
                  1,
                  l->num_frames);
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
  Level *l = (Level *) level;

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

static bool is_type(DicomFile *f, const struct allowed_types *types) {
  // ImageType must be one of the combinations we accept
  int match = -1;
  for (int i = 0; i < types->n_types; i++) {
    bool found_difference = false;

    for (int j = 0; j < 4; j++) {
      const char *type;
      if (!get_tag_str(f->metadata, "ImageType", j, &type)) {
        return false;
      }

      if (!g_str_equal(type, types->types[i][j])) {
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
  dicom_file_destroy(a->file);
  g_slice_free(struct associated, a);
}

static const struct _openslide_associated_image_ops dicom_associated_ops = {
  .get_argb_data = associated_get_argb_data,
  .destroy = associated_destroy
};

static bool read_whole_file(DicomFile *f) {
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

static struct associated *associated_new(DicomFile *f) {
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

  int64_t image_width;
  int64_t image_height;
  if (!get_tag_int(f->metadata, TotalPixelMatrixColumns, &image_width) ||
    !get_tag_int(f->metadata, TotalPixelMatrixRows, &image_height)) {
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

static Level *level_new(DicomFile *f) {
  if (!read_whole_file(f)) { 
    return NULL;
  }

  g_autoptr(Level) l = g_slice_new0(Level);

  if (!get_tag_int(f->metadata, TotalPixelMatrixColumns, &l->image_width) ||
    !get_tag_int(f->metadata, TotalPixelMatrixRows, &l->image_height) ||
    !get_tag_int(f->metadata, Columns, &l->tile_w) ||
    !get_tag_int(f->metadata, Rows, &l->tile_h)) {
    return NULL;
  }

  // ImageType must be one of the combinations we accept
  if (!is_type(f, &level_types)) {
    return NULL;
  }

  // we only allow square tiles
  if (l->tile_w != l->tile_h) {
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

  return g_steal_pointer(&l);
}

static int find_levels(gpointer key, gpointer value, gpointer user_data) {
  const char *filename = (const char *) key;
  DicomFile *f = (DicomFile *) value;
  GHashTable *level_hash = (GHashTable *) user_data;

  Level *l = level_new(f);
  if (l) {
    g_hash_table_insert(level_hash, (char *) filename, l);
    return true;
  }
  else {
    return false;
  }
}

static int remove_bad_slide_id(gpointer key G_GNUC_UNUSED, 
                               gpointer value, 
                               gpointer user_data) {
  DicomFile *f = (DicomFile *) value;
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
  return !get_tag_str(f->metadata, SeriesInstanceUID, 0, &this_slide_id) ||
         !g_str_equal(slide_id, this_slide_id);
}

static int add_level_array(gpointer key G_GNUC_UNUSED, 
                           gpointer value, 
                           gpointer user_data) {
  Level *l = (Level *) value;
  GPtrArray *level_array = (GPtrArray *) user_data;

  g_ptr_array_add(level_array, l);

  return true;
}

#ifdef DEBUG
static void print_file_from_hash(gpointer key G_GNUC_UNUSED, 
                                 gpointer value, 
                                 gpointer user_data G_GNUC_UNUSED) {
  DicomFile *f = (DicomFile *) value;

  print_file(f);
}

static void print_level_from_hash(gpointer key G_GNUC_UNUSED, 
                                  gpointer value, 
                                  gpointer user_data G_GNUC_UNUSED) {
  Level *l = (Level *) value;

  print_level(l);
}
#endif /*DEBUG*/

static int find_associated(gpointer key G_GNUC_UNUSED, 
                            gpointer value, 
                            gpointer user_data) {
  DicomFile *f = (DicomFile *) value;
  openslide_t *osr = (openslide_t *) user_data; 

  struct associated *a = associated_new(f);
  if (a) {
    // map dicom associated image names to openslide names
    char *openslide_name = NULL;
    if (g_str_equal(a->name, "LABEL")) {
      openslide_name = "label";
    }
    else if (g_str_equal(a->name, "OVERVIEW")) {
      openslide_name = "macro";
    }

    if (openslide_name) {
      g_hash_table_insert(osr->associated_images, 
                          g_strdup(openslide_name), 
                          a);
      return true;
    }
  }

  return false;
}

static gint compare_level_width(const void *a, const void *b) {
  const Level *aa = *((const Level **) a);
  const Level *bb = *((const Level **) b);

  return bb->image_width - aa->image_width;
}

static void make_grid(gpointer key G_GNUC_UNUSED, 
                      gpointer value, 
                      gpointer user_data) {
  Level *l = (Level *) value;
  openslide_t *osr = (openslide_t *) user_data; 

  l->grid = _openslide_grid_create_simple(osr,
                                          l->tiles_across, 
                                          l->tiles_down,
                                          l->tile_w, 
                                          l->tile_h,
                                          read_tile);
}

static void add_properties(openslide_t *osr, const Level *l) {
  // why not
  double mmpp = l->image_width;
  g_hash_table_insert(osr->properties,
                      g_strdup(OPENSLIDE_PROPERTY_NAME_MPP_X),
                      _openslide_format_double(mmpp * 1000));
  g_hash_table_insert(osr->properties,
                      g_strdup(OPENSLIDE_PROPERTY_NAME_MPP_Y),
                      _openslide_format_double(mmpp * 1000));
}

static void file_hash(const DicomFile *f,
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
    g_autofree char *filename = g_build_path("/", dirname, name, NULL);

#ifdef DEBUG
    printf("trying to open: %s ...\n", filename);
#endif /*DEBUG*/
    GError *local_error = NULL;
    DicomFile *f = dicom_file_new(filename, &local_error);
    if (!f) {
#ifdef DEBUG
      printf( "open failed: %s\n", local_error->message);
      g_error_free(local_error);
#endif /*DEBUG*/
      continue;
    }

    g_hash_table_insert(dicom_file_hash, g_steal_pointer(&filename), f);
  }

#ifdef DEBUG
  printf("found WSI DICOM files:\n");
  g_hash_table_foreach(dicom_file_hash, print_file_from_hash, NULL);
#endif /*DEBUG*/

  // the filename we were passed should be one of these WSI files ... get the
  // slide id from that
  DicomFile *start = g_hash_table_lookup(dicom_file_hash, filename); 
  if (!start) {
    g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED, 
        "%s is not a DICOM WSI", filename);
    return false;
  }

  if (!read_whole_file(start)) { 
    return NULL;
  }

  const char *slide_id;
  if (!get_tag_str(start->metadata, SeriesInstanceUID, 0, &slide_id)) {
    return false;
  }

  // throw away all files which don't have this slide_id
  g_hash_table_foreach_remove(dicom_file_hash, 
                              remove_bad_slide_id, 
                              (char *) slide_id);

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

  // make the rendering grid
  g_hash_table_foreach(level_hash, make_grid, osr);

#ifdef DEBUG
  printf("\nfinal pyr levels:\n");
  g_hash_table_foreach(level_hash, print_level_from_hash, NULL);
#endif /*DEBUG*/

  // now sort levels by image_width to make the level array
  guint level_count = g_hash_table_size(level_hash);
  if (level_count == 0) {
    g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED,
                "Couldn't find any tiles");
    return false;
  }

  g_autoptr(GPtrArray) level_array = 
    g_ptr_array_new_full(level_count, (GDestroyNotify) level_destroy);
  g_hash_table_foreach_steal(level_hash, add_level_array, level_array);
  g_ptr_array_sort(level_array, compare_level_width); 

  // steal any associated images from the remaining set
  g_hash_table_foreach_steal(dicom_file_hash, find_associated, osr);

  g_assert(osr->data == NULL);
  g_assert(osr->levels == NULL);

  // build ops data
  struct dicom_ops_data *data = g_slice_new0(struct dicom_ops_data);
  data->dirname = g_strdup(dirname);

  guint count = level_array->len;
  // we need to use g_ptr_array_free() to steal the array for openslide, but
  // level_array is an autoptr ... add an extra ref
  g_ptr_array_ref(level_array);
  osr->levels = (struct _openslide_level **) g_ptr_array_free(level_array, 
                                                              false);
  osr->level_count = count;
  data->tile_size = osr->levels[0]->tile_w;
  osr->data = data;
  osr->ops = &dicom_ops;

  // take props from the largest pyr layer
  const Level *largest = (const Level *) osr->levels[0];
  add_properties(osr, largest);

  // perhaps add all the top-level metadata we've loaded for the largest level
  file_hash(largest->file, quickhash1);

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
