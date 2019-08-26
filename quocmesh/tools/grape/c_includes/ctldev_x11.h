/*
 * Copyright (c) 2001-2014 AG Rumpf, INS, Universitaet Bonn                      *
 *                                                                               *
 * The contents of this file are subject to the terms of the Common Development  *
 * and Distribution License Version 1.0 (the "License"); you may not use       *
 * this file except in compliance with the License. You may obtain a copy of     *
 * the License at http://www.opensource.org/licenses/CDDL-1.0                    *
 *                                                                               *
 * Software distributed under the License is distributed on an "AS IS" basis,  *
 * WITHOUT WARRANTY OF ANY KIND, either expressed or implied.                    *
 */

#include "grape.h"   /* precendence necessary for Class Complex in
      differential geometry applications */
#define NEED_EVENTS
#include <X11/Xlibint.h>  /* includes <X11/Xlib.h> */
#include <X11/Xutil.h>
#include <X11/Xresource.h>
#include <X11/cursorfont.h>
#define XK_XKB_KEYS
#include <X11/keysym.h>
#ifdef G_X11_SPACEBALL
#include <X11/extensions/XI.h>
#include <X11/extensions/XInput.h>
#endif
#include "vislst_x11.h"

/* grape colormap for control device: 16 read/only and 14 r/w cells for
   lightsources and surface properties */
#define NO_GR_COLORCELLS (30)
#define NO_BASE_COLS 16

typedef struct gr_col_map {
  int read_only;    /*  = 0 color cell is read/write */
  unsigned long pixel;    /* pixel value of server's colormap */
} GR_COL_MAP;

/* struct to keep mouse in ruler, slider, ... */
typedef struct ms_grab {
  int grabbed;      /* = 1, if mouse is grabbed */
  double xmax, xmin, ymax, ymin;  /* in area */
} MS_GRABBED;

typedef struct g_control_info_x11{
/* some variables used for grape control device driver programming */
  RESOLUTION_TYPE resolution;  /* current dimensions in pixels of window */
  double aspectratio;    /* current ratio of width and height in pixel */
  int textheight;    /* textheight in pixel */
  double x_grid_resolution;   /* no. of patches in x-direction (min. 12) */
  double y_grid_resolution;   /* no. of patches in y-direction (min. 25) */
  double x_pixel_per_patch;   /* no. of pixels/patch in x-dir  (min. 28) */
  double y_pixel_per_patch;   /* no. of pixels/patch in y-dir  (min. 15) */
  double x_curr_pos, y_curr_pos;  /* current x-, y-position for draw */
  double transx, transy;    /* patch size in pixel (x-, y-dir.) */
  int current_color;  /* type of current color (e.g. BACKGROUND_COLOR) */
  int text_color;  /* type of current text color (e.g. TEXT_COLOR) */
  int free_color_index;  /* next index in Grape's Control Color map */
  int antialias;    /* antialias available */
  int doublebuffer;    /* doublebuffering available */
/* variables for special X11 window handling */
  VisListEntry *vis;    /* all information about current visual */
  Display *display;    /* connected to X Display */
  int screen;      /* server may have multiple screens */
  Window window;    /* Window id */
  GC gr_context;    /* Graphics Context */
  XFontStruct *font_info;  /* metric information for an entire font */

  GR_COL_MAP gr_color_map[NO_GR_COLORCELLS];
  MS_GRABBED mausi;

  double event_oldx, event_oldy;
  unsigned int buttons;

  Atom insert_prop;
  char *pasted_text;
  int paste_next;
  double paste_x, paste_y;

  char *copied_text;
  Time last_event_time;

#ifdef G_X11_SPACEBALL
  XDevice *spaceball_device;
  int spaceball_motion_type;
  int spaceball_button_press_type, spaceball_button_release_type;
  struct {
    int axis[6];
    int button[9];
    Time time, time_diff;
  } spaceball_data;
#endif
} G_CONTROL_INFO_X11;

extern int g_x11_ctl_x_min_size, g_x11_ctl_y_min_size, g_x11_ctl_y_std_size;

/*  general grape functions implemented in drvctl_x11.c */
void g_x11_ctl_update(void),
  g_x11_ctl_clear(void),
  g_x11_ctl_set_color(int),
  g_x11_ctl_set_text_color(int),
  g_x11_ctl_set_rgb_color(const VEC3),
  g_x11_ctl_get_color(int *),
  g_x11_ctl_change_color(int, const VEC3),
  g_x11_ctl_moveto(double, double),
  g_x11_ctl_drawto(double, double),
  g_x11_ctl_relative_drawto(double, double),
  g_x11_ctl_draw_string(double, double, const char *),
  g_x11_ctl_get_mouse(double *, double *),
  g_x11_ctl_fix_mouse(double, double, double, double, double, double),
  g_x11_ctl_set_mouse(double, double),
  g_x11_ctl_free_mouse(void),
  g_x11_ctl_rectangle(double, double, double, double),
  g_x11_ctl_circle(double, double, double),
  g_x11_ctl_get_control_size(double *, double *),
  g_x11_ctl_wait(double),
  g_x11_ctl_clear_event_queue(void),
  g_x11_ctl_reinit(void),
  g_x11_ctl_set_clip_rects(GRECT_LIST *),
  g_x11_ctl_copy_string(const char *),
  g_x11_ctl_get_spaceball_data(double *, double *, double *, double *,
             double *, double *, double *, int *);

int g_x11_ctl_attribute(int, int, int *),
  g_x11_ctl_menu_button(void),
  g_x11_ctl_middle_button(void),
  g_x11_ctl_read_string(double, double, char *, double),
  g_x11_ctl_get_free_color_index(void),
  g_x11_ctl_handle_events(EVENT *, int);

double  g_x11_ctl_stringwidth(const char *),
  g_x11_ctl_stringheight(const char *),
  g_x11_ctl_get_unit_per_pixel(void);

G_CLIP_DESCRIPTOR *g_x11_ctl_kind_of_clipping(void);
