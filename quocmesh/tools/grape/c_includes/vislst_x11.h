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

/* list of open displays, devices, visuals
 */

/* define debugging mode in all sources using this header file */
/* #define DEBUG */
/* sometimes changed to DEBUG_OLD, if boring */
#define NUM_STIPPLES 16
#define MAX_SHORT_UNSIGNED 65535.0

typedef struct VISLISTENTRY {
  struct VISLISTENTRY *next;
  Visual *visual;
  Display *display;
  char *disp_name;    /* to identify identical displays */
  int screen;
  unsigned int depth;
  int class;
  Pixmap stip[NUM_STIPPLES];
  int color_cube_max[3];
  int color_cube_mult[3];
  int size_of_col_cube;
  unsigned long *pixels;
  unsigned long white, black;
  Colormap colormap;
} VisListEntry;

/* this is just to print the visual class intelligibly */
static char *visual_class[] = {
  "StaticGray", "GrayScale", "StaticColor", "PseudoColor",
  "TrueColor", "DirectColor"
};

#ifdef DEBUG
/* to print the window events on stderr */
static char *event_names[] = {
"", "",
"KeyPress",    "KeyRelease",
"ButtonPress",    "ButtonRelease",
"MotionNotify",   "EnterNotify",    "LeaveNotify",
"FocusIn",    "FocusOut",    "KeymapNotify",
"Expose",    "GraphicsExpose",  "NoExpose",
"VisibilityNotify",  "CreateNotify",    "DestroyNotify",
"UnmapNotify",    "MapNotify",
"MapRequest",
"ReparentNotify",
"ConfigureNotify",
"ConfigureRequest",
"GravityNotify",
"ResizeRequest",
"CirculateNotify",
"CirculateRequest",
"PropertyNotify",
"SelectionClear",
"SelectionRequest",
"SelectionNotify",
"ColormapNotify",
"ClientMessage",
"MappingNotify" };
#endif
