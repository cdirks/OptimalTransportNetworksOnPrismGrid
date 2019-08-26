/*----------------------------------------------------------------------*\
 **                  **
 ** drvctl_x11.c              **
 **                  **
 ** Description: primitive control device handling routines    **
 **    for X11 display, g_ctl_x11_<functions> are sorted  **
 **    according to list in 'ctldev_x11.c'      **
 **    g_control_<function> are unknown to grape    **
 **                  **
 ** Copyright (C) 1991  by Universitaet Bonn, Freiburg      **
 **      Institut fuer Angewandte Mathematik    **
 **      Sonderforschungsbereich 256      **
 **      D-53115 Bonn, D-79104 Freiburg i.Br.    **
 **                  **
 ** Authors:  Alfred Schmidt, Thomas Kelterbach, Thomas Mackeben  **
 */
static char rcsid[]="$Id: drvctl_x11.c,v 2.40 2000/09/05 12:31:28 wu bn $";

#include "ctldev_x11.h" /* inclusive  grape.h, X11's, vislst_x11.h */
#include <X11/Xatom.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

extern G_CONTROL_INFO_X11 *current_ctl_info_x11;

#define INFO current_ctl_info_x11

void g_x11_ctl_update()
{
  XFlush(INFO->display);
}


void g_x11_ctl_clear()
{
  XClearWindow(INFO->display,INFO->window);
}


void g_control_initialize_color()
     /* this function is called only by new-instance on ControlX11
  and allocates the read only colors */
{
  XColor exact_def;
  int i;

  char *base_col[NO_BASE_COLS], *bw_col[NO_BASE_COLS];

  base_col[0] = "black"; base_col[1] = "red"; base_col[2] = "green";
  base_col[3] = "yellow"; base_col[4] = "blue"; base_col[5] = "magenta";
  base_col[6] = "cyan"; base_col[7] = "white";
  base_col[8]  = "white";  /* FOREGROUND_COLOR */
  base_col[9]  = "gray80";  /* BACKGROUND_COLOR (was gray40)*/
  base_col[10] = "white";  /* PRESSED_LOWER_COLOR and (gray80) (Umrandung links oben)*/
  base_col[11] = "gray20";  /* PRESSED_UPPER_COLOR for shading (gray20) (Umrandung rechts unten) */
  base_col[12] = "gray20";  /* TEXT_COLOR (gray90)*/
  base_col[13] = "gray20";  /* LINE_COLOR (gray10)*/
  base_col[14] = "goldenrod";  /* LIGHT_ON_COLOR for selected light */
  base_col[15] = "unused";  /* not used and for debugging */

  /* if color allocation fails, simply try black and white */
  bw_col[0] = "black"; bw_col[1] = "white"; bw_col[2] = "white";
  bw_col[3] = "white"; bw_col[4] = "white"; bw_col[5] = "black";
  bw_col[6] = "white"; bw_col[7] = "white";
  bw_col[8]  = "white";    /* FOREGROUND_COLOR */
  bw_col[9]  = "black";    /* BACKGROUND_COLOR */
  bw_col[10] = "white";    /* PRESSED_LOWER_COLOR and */
  bw_col[11] = "black";    /* PRESSED_UPPER_COLOR for shading */
  bw_col[12] = "white";    /* TEXT_COLOR */
  bw_col[13] = "black";    /* LINE_COLOR */
  bw_col[14] = "white";    /* LIGHT_ON_COLOR for selected light */
  bw_col[15] = "black";    /* not used and for debugging */

  for ( i=0; i < NO_BASE_COLS; i++ ) {
#ifndef DEBUG
    if( strcmp(base_col[i],"unused") == 0) continue;
#endif
    if( ! XParseColor (INFO->vis->display, INFO->vis->colormap,
           base_col[i], &exact_def)) {
      printf("g_control_initialize_color : %s not in data base, trying %s\n",
       base_col[i], bw_col[i]);
      if( ! XParseColor (INFO->vis->display, INFO->vis->colormap,
       bw_col[i], &exact_def)) {
  printf("g_control_initialize_color : %s not in", bw_col[i]);
  printf(" data base, giving up!\n");
  exit(-1); }
    }
    if( XAllocColor(INFO->vis->display, INFO->vis->colormap, &exact_def)
       && (INFO->vis->depth > 1) ) {
      INFO->gr_color_map[i].read_only = 1;
      INFO->gr_color_map[i].pixel = exact_def.pixel;
#ifdef DEBUG_OLD
      printf("g_control_initialize_color : %s has rgb values %f %f %f\n",
       base_col[i], INFO->gr_color_map[i].rgb[0],
       INFO->gr_color_map[i].rgb[1], INFO->gr_color_map[i].rgb[2]);
#endif
    }
    else {
#ifdef DEBUG_OLD
      printf("g_control_initialize_color : no cell for %s, using %s\n",
       base_col[i], bw_col[i]);
#endif
      if( strcmp(bw_col[i],"black") == 0 )
  INFO->gr_color_map[i].pixel = INFO->vis->black;
      else INFO->gr_color_map[i].pixel = INFO->vis->white;
    }
  }
  INFO->free_color_index = NO_BASE_COLS;
}

unsigned long g_control_get_backgroundpixel()
     /* this function is called only by new-instance on ControlX11 */
{
  return( INFO->gr_color_map[BACKGROUND_COLOR].pixel);
}

void g_x11_ctl_set_color(int col)
{        /* avoid unneccessary X calls */
  if( col < 0 || INFO->free_color_index <= col) {
    fprintf(stderr,"g_x11_ctl_set_color : index %d out of range\n",col);
    return; }
  if (INFO->gr_color_map[INFO->current_color].pixel !=
      INFO->gr_color_map[col].pixel)
    XSetForeground(INFO->display, INFO->gr_context,
       INFO->gr_color_map[col].pixel);
  INFO->current_color = col;
}

void g_x11_ctl_set_text_color(int col)
{
  INFO->text_color = col;     /* the foreground pixel is set in draw_string */
}


void g_x11_ctl_set_rgb_color(const VEC3 rgb)
     /* this function is replaced by 'set_color' and 'change_color' */
{
  printf("Please use 'get-free-color-index', "
   "'change-color' and 'set-color'.\n");
  printf("to set rgb color %f %f %f\n", rgb[0], rgb[1], rgb[2]);
  printf("g_x11_ctl_set_grb_color : function "
   "is no longer available, sorry.\n");
}


void g_x11_ctl_get_color(int *col)
{
  *col = INFO->current_color;
}


void g_x11_ctl_change_color(int col, const VEC3 rgb)
{
  XColor exact_def;

  if( col < NO_BASE_COLS || INFO->free_color_index <= col) {
    fprintf(stderr,"g_x11_ctl_set_color : index %d out of range,",col);
    fprintf(stderr," using black.\n");
    XSetForeground(INFO->display, INFO->gr_context, INFO->vis->black);
    return; }

  /* if we really need that we have to use XSetForeground (...col..)
   * or set_color, because INFO->current_color is assumed to be
   * the GC's foreground. Of course the color has to be changed before
   * that.
   */
  /* INFO->current_color = col; */

  switch(INFO->vis->class) {

  case DirectColor:
  case PseudoColor:    /* read/write cells */
  case GrayScale:
    if (INFO->gr_color_map[col].read_only == 0) {
      exact_def.flags = DoRed | DoGreen | DoBlue;
      exact_def.red = (unsigned short) (rgb[0] * MAX_SHORT_UNSIGNED);
      exact_def.green = (unsigned short) (rgb[1] * MAX_SHORT_UNSIGNED);
      exact_def.blue = (unsigned short) (rgb[2] * MAX_SHORT_UNSIGNED);
      exact_def.pixel = INFO->gr_color_map[col].pixel;
      XStoreColor(INFO->display, INFO->vis->colormap, &exact_def);
      break;
    }
  case TrueColor:    /* otherwise use read/only stategy */
  case StaticColor:
  case StaticGray:
    if ((INFO->vis->depth > 1)) {
      /* if we free white_pixel or black_pixel, we cause BadAccess,
       * so we don't. (That is because black&white were not allocateded with
       * XAllocColor sometimes).
       */
      if (INFO->gr_color_map[col].pixel != INFO->vis->white &&
    INFO->gr_color_map[col].pixel != INFO->vis->black)
  XFreeColors(INFO->vis->display, INFO->vis->colormap,
        &(INFO->gr_color_map[col].pixel), 1, 0);
      exact_def.red = (unsigned short) (rgb[0] * MAX_SHORT_UNSIGNED);
      exact_def.green = (unsigned short) (rgb[1] * MAX_SHORT_UNSIGNED);
      exact_def.blue = (unsigned short) (rgb[2] * MAX_SHORT_UNSIGNED);
      if( XAllocColor(INFO->vis->display, INFO->vis->colormap, &exact_def)) {
  INFO->gr_color_map[col].pixel = exact_def.pixel;
  break;
      }
    }
  default:      /* simply black&white */
    if (rgb[0]+rgb[1]+rgb[2] > 1.5)
      INFO->gr_color_map[col].pixel = INFO->vis->white;
    else
      INFO->gr_color_map[col].pixel = INFO->vis->black;
    break;
  }

  g_x11_ctl_set_color (col);
}


void g_x11_ctl_moveto(double x, double y)
{
  INFO->x_curr_pos = x;
  INFO->y_curr_pos = y;
}


void g_x11_ctl_drawto(double x, double y)
{
  int start[2], end[2];

  start[0] = (int) (INFO->transx * (INFO->x_curr_pos) + .5);
  start[1] = (int) (INFO->transy * (INFO->y_curr_pos) + .5);
  end[0]   = (int) (INFO->transx * x + .5);
  end[1]   = (int) (INFO->transy * y + .5);
  XDrawLine(INFO->display,INFO->window,INFO->gr_context,
      start[0],start[1],end[0],end[1]);
  INFO->x_curr_pos = x;
  INFO->y_curr_pos = y;
}


void g_x11_ctl_relative_drawto(double x, double y)
{
  int start[2], end[2];

  start[0] = (int) (INFO->transx * (INFO->x_curr_pos) + .5);
  start[1] = (int) (INFO->transy * (INFO->y_curr_pos) + .5);
  end[0] = (int) (INFO->transx * (INFO->x_curr_pos += x)+ .5);
  end[1] = (int) (INFO->transy * (INFO->y_curr_pos += y)+ .5);
  XDrawLine(INFO->display,INFO->window,INFO->gr_context,
      start[0],start[1],end[0],end[1]);
  return;
}


void g_x11_ctl_draw_string(double x, double y, const char *str)
{
  int X, Y, current_color = INFO->current_color;

  X = (int) (INFO->transx * x+ .5);
  Y = (int) (INFO->transy * y+ .5);
  g_x11_ctl_set_color(INFO->text_color);
  XDrawString(INFO->display,INFO->window,INFO->gr_context,X,Y,
        str,(int)strlen(str));
  /* INFO->current_color is assumed to be the GC's foreground. */
  g_x11_ctl_set_color(current_color);
}


void g_x11_ctl_get_mouse(double *x, double *y)
{
  Window root, child;
  int r_x, r_y, win_x, win_y;
  unsigned int buttons;

  while( ! XQueryPointer(INFO->display,INFO->window,
       &root, &child, &r_x, &r_y, &win_x, &win_y, &buttons));

  *x = (double) (((double) win_x) / INFO->transx);
  *y = (double) (((double) win_y) / INFO->transy);

  if(INFO->mausi.grabbed) {
    if( *x > INFO->mausi.xmax) *x = INFO->mausi.xmax;
    if( *x < INFO->mausi.xmin) *x = INFO->mausi.xmin;
    /* very special handling of y coordinate is inherited by sgi's gl driver */
    if( *y > INFO->mausi.ymin) *y = INFO->mausi.ymin;
    if( *y < INFO->mausi.ymax) *y = INFO->mausi.ymax;
  }
}


void g_x11_ctl_fix_mouse(double x, double xmin, double xmax,
       double y, double ymin, double ymax)
{
#ifdef DEBUG_OLD
  printf("fix-mouse: x =%9.6f, xmax = %7.4f, xmin = %7.4f\n",x,xmax,xmin);
  printf("           y =%9.6f, ymax = %7.4f, ymin = %7.4f\n",y,ymax,ymin);
#endif
  XDefineCursor(INFO->display, INFO->window,
    XCreateFontCursor(INFO->display, XC_bottom_side));
  /* store values of area in which mouse pointer has to be enclosed */
  INFO->mausi.grabbed = 1;
  INFO->mausi.xmin = xmin; INFO->mausi.xmax = xmax;
  INFO->mausi.ymin = ymin; INFO->mausi.ymax = ymax;
  /* and put the mouse pointer into the desired initial position */
  XWarpPointer(INFO->display, None, INFO->window, 0,0,0,0,
         (int) (INFO->transx * x), (int) (INFO->transy * y) );
}

void g_x11_ctl_free_mouse(void)
{
  INFO->mausi.grabbed = 0;
}

void g_x11_ctl_set_mouse( double x, double y)
{
  XEvent event;
  long event_mask = MotionNotify;

  XWarpPointer(INFO->display, None, INFO->window, 0,0,0,0,
         (int) (INFO->transx * x), (int) (INFO->transy * y) );

  while(XCheckMaskEvent(INFO->display,event_mask,&event))
    ;
}


void g_x11_ctl_rectangle(double x1, double y1, double x2, double y2)
{
  double X1,X2,Y1,Y2;
  int llx, lly;

  X1 = (x1 > x2) ? x2 : x1;
  X2 = (x1 < x2) ? x2 : x1;
  Y1 = (y1 > y2) ? y2 : y1;
  Y2 = (y1 < y2) ? y2 : y1;

  X2 *= INFO->transx;
  Y2 *= INFO->transy;
  X1 *= INFO->transx;
  Y1 *= INFO->transy;
  llx = (int)(X1 +.5);
  lly = (int)(Y1 + .5);
  X2 = X2 - llx + 1;
  Y2 = Y2 - lly + 1;
  XFillRectangle(INFO->display, INFO->window, INFO->gr_context,
     llx, lly, (int)(X2 +.5), (int)(Y2 +.5));
}

void g_x11_ctl_circle(double x, double y, double r)
{
  double X, Y, RX, RY;
  X = INFO->transx * g_fabs(x-r); /* Grape math function for */
  Y = INFO->transy * g_fabs(y-r); /* machine compatability */
  RX = INFO->transx * r * 2.0;
  RY = INFO->transy * r * 2.0;
  XFillArc(INFO->display, INFO->window, INFO->gr_context,
     (int)(X+.5), (int)(Y+.5), (unsigned int)(RX +.5),
     (unsigned int)(RY+.5), 0, 360*64);
}


void g_x11_ctl_get_control_size(double *x, double *y)
{
  *x = INFO->x_grid_resolution;
  *y = INFO->y_grid_resolution;
}


void g_x11_ctl_wait(double wait_time)
{
  fprintf (stderr, "g_x11_ctl_wait called with argument %g\n", wait_time);
}


void g_x11_ctl_clear_event_queue(void)
{
  XEvent event;
  long event_mask = (long) 0xffffffff;
  while(XCheckMaskEvent(INFO->display,event_mask,&event))
    ;
}


void g_x11_ctl_reinit(void)
     /*  supply info structure variables with values */
{
  XWindowAttributes att;

  XGetWindowAttributes(INFO->display, INFO->window, &att);
  INFO->resolution.x = att.width;
  INFO->resolution.y = att.height;
  if (INFO->resolution.y > 0)
    INFO->aspectratio = INFO->resolution.x / (double) INFO->resolution.y;
  else
    INFO->aspectratio = 1.0;

  if (INFO->resolution.x >= g_x11_ctl_x_min_size * 12) {
    INFO->x_grid_resolution = INFO->resolution.x
      / (double) g_x11_ctl_x_min_size;
    INFO->x_pixel_per_patch = g_x11_ctl_x_min_size;
  }
  if (INFO->resolution.y >= g_x11_ctl_y_std_size * 26) {
    INFO->y_grid_resolution = INFO->resolution.y
      / (double) g_x11_ctl_y_std_size;
    INFO->y_pixel_per_patch = g_x11_ctl_y_std_size;
  }
  else if (INFO->resolution.y >= g_x11_ctl_y_min_size * 26) {
    INFO->y_grid_resolution = 26.;
    INFO->y_pixel_per_patch = INFO->resolution.y / 26.;
  } else { /* sometimes allowed by server */
    INFO->y_grid_resolution = INFO->resolution.y
      / (double) g_x11_ctl_y_min_size;
    INFO->y_pixel_per_patch = g_x11_ctl_y_min_size;
  }

  INFO->transx = INFO->x_pixel_per_patch;
  INFO->transy = INFO->y_pixel_per_patch;

  INFO->textheight = INFO->font_info->ascent + INFO->font_info->descent;
#ifdef DEBUG
  printf("dbg: g_x11_ctl_reinit: Texthoehe %d\n",INFO->textheight);
#endif

  INFO->current_color = BACKGROUND_COLOR;
  XSetForeground(INFO->display, INFO->gr_context,
     INFO->gr_color_map[BACKGROUND_COLOR].pixel);
  INFO->text_color = TEXT_COLOR;

  INFO->antialias = 0;    /* no antialiasing available */
  INFO->doublebuffer = 0;  /* no doublebuffering available */
}


int g_x11_ctl_attribute(int mode, int type, int *data)
{

  /***************************************\
   *  switch for different attributes  *
   \***************************************/

  switch (type) {

    /***********************************\
     *  doublebuffing handling    *
     \***********************************/

  case G_DOUBLEBUFFER:
    switch (mode) {
    case G_MODE_GET:
      *(int *)data = 0;
      break;
    case G_MODE_SET:  break;
    case G_MODE_SWAP:  break;
    default:    printf("g_x11_ctl_attribute : Unknown mode %d for type %d\n",
           mode, G_DOUBLEBUFFER); return(0);
    }
    break;

  case G_HAS_DOUBLEBUFFER:
    switch (mode) {
    case G_MODE_GET:
      *(int *)data = 0;
      break;
    default:    printf("g_x11_ctl_attribute : Unknown mode %d for type %d\n",
           mode, G_HAS_DOUBLEBUFFER); return(0);
    }
    break;
  default:  printf("g_x11_ctl_attribute : mode %d for unknown type %d\n",
           mode, type); return(0);
  }
  return(1);      /* successful operation */
}


int g_x11_ctl_menu_button()
{
  Window root, child;
  int r_x,r_y,win_x,win_y;
  unsigned int buttons;

  if (XQueryPointer(INFO->display,INFO->window, &root, &child, &r_x,
        &r_y, &win_x, &win_y, &buttons)) {

    return ( (int) (buttons & Button1Mask) );
  }
  else return (0);
}


int g_x11_ctl_middle_button()
{
  Window root,child;
  int r_x,r_y,win_x,win_y;
  unsigned int buttons;
  if (XQueryPointer(INFO->display,INFO->window, &root, &child, &r_x,
        &r_y, &win_x, &win_y, &buttons))
    return ( (int) (buttons & Button2Mask) );
  else return (0);
}

double g_x11_ctl_stringwidth(const char *),
  g_x11_ctl_stringheight(const char *);

int g_x11_ctl_read_string(double x, double y, char *str, double max)
     /* max = display length of string in number of patches in a row */
{
  XEvent event;
  char ch;
  double llength, slength, dheight, dlength, x_unit, y_unit, curr_length;
  int ende = 0, str_length = 0;
  int curr_ptr;      /* begin of displayed string */
  short val = 0;
  char *str_disp;    /* part of string displayed */

  x_unit = 1.0 / INFO->transx;  /* height of pixel in x-dir */
  y_unit = 1.0 / INFO->transy;
  dheight = g_x11_ctl_stringheight(str);
  curr_ptr = 0;
  if(max <= 0.0) max = 12.0;  /* max number of patches in a row */
  str_disp = str;
  str[0] = '_'; str[1] = '\0';

  g_x11_ctl_draw_string(x,y,str);

  XDefineCursor(INFO->display, INFO->window,
    XCreateFontCursor(INFO->display, XC_question_arrow));
  /* input text */
  while(!ende) {
    XNextEvent(INFO->display,&event);
    switch(event.type) {
    case ButtonPress:
      /* caution!!! tk      if(val != 0) */
      ende = 1;
      break;
    case KeyPress:
      if(XLookupString((XKeyEvent *) &event,&ch,1,NULL,NULL)){
  val = ch;
  if(val >= 32 && val <= 127 && ((str_length + 1) < MAX_STRING_LENGTH)) {
    str[str_length] = val;
    str_length++;
    str[str_length] = '_'; str[str_length+1] = '\0';
    curr_length = g_x11_ctl_stringwidth(str);
    if(curr_length > max) { /* scroll string to left on screen */
      curr_ptr++;
      str_disp = &(str[curr_ptr]);
    }
  }
  else {
    switch(val) {
    case 13:    /* RETURN */
      ende = 1;
      break;
    case 27:    /* ESCAPE */
      str_length = 0;
            XNextEvent(INFO->display,&event);
            if(XLookupString((XKeyEvent *) &event,&ch,1,NULL,NULL)) val = ch ;
      ende = 2;
      break;
    case 8:    /* DELETE */
      if(str_length > 0) {
        llength = g_x11_ctl_stringwidth(str_disp);
        if(curr_ptr > 0) { /* scroll string to right on screen */
    curr_ptr--;
    str_disp = &(str[curr_ptr]);
        }
        str_length--;
        str[str_length] = '_'; str[str_length+1] = '\0';
        /* wipe out last char on screen */
        slength = g_x11_ctl_stringwidth(str_disp);
        g_x11_ctl_set_color(BACKGROUND_COLOR);
        g_x11_ctl_rectangle(x+slength, y + 2.0*y_unit,
          x + llength + x_unit, y - y_unit - dheight);
      }
    }
  }
  dlength = g_x11_ctl_stringwidth(str_disp);
  g_x11_ctl_set_color(BACKGROUND_COLOR);
  g_x11_ctl_rectangle(x, y + 2.0 * y_unit,
          x + dlength + x_unit, y - y_unit - dheight);
  g_x11_ctl_draw_string(x,y,str_disp);
  break;
      }
    }
  }
  str[str_length] = '\0';
  if(ende == 2)
    return(/*CONTROL_ESCAPE*/ 1);
  else
    return(/*CONTROL_NOTHING*/ 0);
}


int g_x11_ctl_get_free_color_index()
{
  XColor exact_def, unused;
  unsigned long plane[1], pixel[1];
  int ifree;

  ifree = INFO->free_color_index;
  if (ifree >= NO_GR_COLORCELLS) return(0);

  if (  INFO->vis->class == PseudoColor || INFO->vis->class == DirectColor ||
      INFO->vis->class == GrayScale ) {
    if( XAllocColorCells(INFO->display, INFO->vis->colormap, False,
       plane, 0, pixel, 1) ) {
      INFO->gr_color_map[ifree].pixel = pixel[0];
      INFO->gr_color_map[ifree].read_only = 0; /* r/w cell available */
      return(INFO->free_color_index++);
    }
  }
  INFO->gr_color_map[ifree].read_only = 1; /* use read/only strategy */
  if( ! XAllocNamedColor(INFO->vis->display, INFO->vis->colormap, "black",
       &exact_def, &unused) ) {
    printf("g_x11_ctl_get_free_color_index : black not in");
    printf(" data base, giving up!\n");
    exit(-1);
  }
  INFO->gr_color_map[ifree].pixel = exact_def.pixel;
  return(INFO->free_color_index++);
}

static Bool interesting_event(Display *diplay, XEvent *event, XPointer arg)
{
  if (event->xany.display == INFO->display &&
      event->xany.window == INFO->window)
    return TRUE;
  switch (event->type) {
  case Expose:
  case ConfigureNotify:
    return TRUE;
  default:
    break;
  }
  return FALSE;
}

static Bool next_window_event(Display *diplay, XEvent *event, XPointer arg)
{
  if (event->xany.display == INFO->display &&
      event->xany.window == INFO->window)
    return TRUE;
  return FALSE;
}

/* prefer updating current state of window */
static Bool window_management_event(Display *diplay, XEvent *event,
            XPointer arg)
{
  if (event->xany.display == INFO->display &&
      event->xany.window == INFO->window)
    switch (event->type) {
    case Expose:
    case ConfigureNotify:
      return TRUE;
    default:
      break;
    }
  return FALSE;
}

struct swallow {
  Display *display;
  Window window;
};

static Bool swallow_redraw_event (Display *display, XEvent *event,
          XPointer arg)
{
  struct swallow *data = (void *)arg;
  return event->xany.display == data->display &&
    /* grdev opens lots of nested windows ... ? */
    /* event->xany.window == data->window && */
    event->xany.window != INFO->window &&
    (event->xany.type == Expose ||
     event->xany.type == ConfigureNotify);
}

static void swallow_redraw (Display *display, Window window)
{
  XEvent event;
  struct swallow data;
  data.display = display;
  data.window = window;
  while (XCheckIfEvent(INFO->display, &event, swallow_redraw_event,
           (void *)&data))
    /* fprintf (stderr, "%s\n", event.xany.type == Expose ? "Expose" : "CN")*/
    ;
}

/* taken from handle.h */
GRECT *grect_alloc (double , double , double , double);
void *grect_list_free       (GRECT_LIST *);
GRECT_LIST *grect_list_append (GRECT_LIST *, GRECT *);
GRECT_LIST *grect_unpack      (GRECT_LIST *, GRECT *);
GRECT_LIST *grect_pack        (GRECT_LIST *);
GRECT_LIST *grect_cut (GRECT_LIST *, GRECT *);
GRECT_LIST *grect_list_del_intersects (GRECT_LIST *, GRECT * );
GRECT_LIST *grect_cut_list (GRECT_LIST **, GRECT_LIST *, GRECT *);
GRECT_LIST *grect_list_append (GRECT_LIST *, GRECT *);
GRECT_LIST *grect_list_free_member( GRECT_LIST *list, GRECT_LIST *member );
int grect_count_members( GRECT_LIST *list );

static int grect_list_count_entries (GRECT_LIST *list)
{
  int count = 0;
  for (; list; list = list->next)
    count++;
  return count;
}

static void add_redraw_rect (GRECT_LIST **list, XEvent *event)
{
  if (grect_list_count_entries (*list) > 32) {
    grect_list_free (*list);
    *list = grect_list_append
      (NULL, grect_alloc (0.0, 0.0,
         INFO->resolution.x / (double)INFO->transx,
         INFO->resolution.y / (double)INFO->transy));
  } else {
    GRECT rect;
    rect.x1 = event->xexpose.x / INFO->transx;
    rect.y1 = event->xexpose.y / INFO->transy;
    rect.x2 = (event->xexpose.x+event->xexpose.width) / INFO->transx;
    rect.y2 = (event->xexpose.y+event->xexpose.height) / INFO->transy;
    /* *list = grect_unpack (*list, &rect);  too buggy */
    *list = grect_list_append (*list, grect_alloc (rect.x1, rect.y1,
               rect.x2, rect.y2));
  }
}

#if 0
static void get_redraw_rect (GRECT_LIST **list, EVENT *event)
{
  GRECT rect = *(*list)->rect;
  *list = grect_list_del_intersects (*list, &rect);
  event->x      = rect.x1;
  event->y      = rect.y1;
  event->deltax = rect.x2 - rect.x1;
  event->deltay = rect.y2 - rect.y1;
}
#endif

static int check_redraw (GRECT_LIST **redraw_list, int resized,
       EVENT *g_event)
{
  XEvent event;

  memset (&event, 0, sizeof (XEvent));

  while (XCheckIfEvent(INFO->display, &event, window_management_event, NULL))
    switch (event.type) {
    case Expose:
      add_redraw_rect (redraw_list, &event);
      break;
    case ConfigureNotify:
      resized = TRUE;
      *redraw_list = grect_list_free (*redraw_list);
      g_x11_ctl_reinit ();
      g_event->x      = event.xconfigure.x / INFO->transx;
      g_event->y      = event.xconfigure.y / INFO->transy;
      g_event->deltax = event.xconfigure.width / INFO->transx;
      g_event->deltay = event.xconfigure.height / INFO->transy;
      break;
    }

  if (resized) {
    g_event->what = evResize;
    g_event->who  = evForCtl;
    /* g_event->causer = self; */
    return 2;
  }

  if (*redraw_list) {
    g_event->what = evRedraw;
    g_event->who  = evForCtl;
    /* g_event->causer = self; */
    /* get_redraw_rect (redraw_list, g_event); */
    g_event->rect_list = *redraw_list;
    *redraw_list = NULL;
    return 1;
  }

  return 0;
}

static void get_graph_events (void)
{
  XEvent event;
  memset (&event, 0, sizeof (XEvent));

  /* look for expose event if different Grape Output display (assumed next
   * visual) and send them to the Control window */
  if (INFO->vis->next)
    while (XCheckTypedEvent(INFO->vis->next->display, Expose, &event))
      if (!XSendEvent(INFO->display, INFO->window, False, ExposureMask,
          &event))
  fprintf(stderr,
    "g_x11_ctl_handle_events : Cannot "
    "send event from output window\n");
}

int g_x11_ctl_handle_events (EVENT *g_event, int wait)
{
  int ignore_event;

  do {
    XEvent event;

    get_graph_events ();

    ignore_event = FALSE;

    /* prevent false events ... */
    memset (g_event, 0, sizeof (EVENT));
    memset ( &event, 0, sizeof (XEvent));

    /* prefer updating current state of window */
    {
      GRECT_LIST *redraw_list = NULL;
      if (check_redraw (&redraw_list, FALSE, g_event))
  return 1;
    }

    if (INFO->pasted_text) {
      g_event->what = evNormalKey;
      g_event->which = INFO->pasted_text [INFO->paste_next++];
      g_event->who = evForCtl;
      g_event->x = INFO->paste_x;
      g_event->y = INFO->paste_y;
      g_event->state = 0;

      if (!INFO->pasted_text [INFO->paste_next]) {
  INFO->pasted_text = g_strfree (INFO->pasted_text);
  INFO->paste_next = 0;
      }
      return 1;
    }

    if (INFO->buttons) {
      event.xany.window = INFO->window;
      g_event->what = evMouseHold;
      if (!XCheckIfEvent(INFO->display, &event, interesting_event, NULL))
  event.type = LASTEvent;     /* prevent false events ... */
    } else {
      g_event->what = evNothing;
      if (wait)
  XNextEvent(INFO->display, &event);
      else
  if (!XCheckIfEvent(INFO->display, &event, interesting_event, NULL))
    event.type = LASTEvent;
    }

    if (event.xany.window == INFO->window) {
      g_event->who  = evForCtl;
      /* g_event->causer = self; */
    } else
      g_event->who  = evForGrph;

    switch(event.type) {
    case Expose:
      if (event.xexpose.window == INFO->window) {
  GRECT_LIST *redraw_list = NULL;
  add_redraw_rect (&redraw_list, &event);
  if (!check_redraw (&redraw_list, FALSE, g_event))
    ignore_event = TRUE;
      }
      else {
  if (event.xexpose.count != 0)
    ignore_event = TRUE;
  else {
    swallow_redraw (event.xany.display, event.xany.window);
    g_event->what = evRedraw;
  }
      }
      break;

    case ConfigureNotify:
      if (event.xconfigure.window == INFO->window) {
  GRECT_LIST *redraw_list = NULL;
  g_x11_ctl_reinit ();
  g_event->x      = event.xconfigure.x / INFO->transx;
  g_event->y      = event.xconfigure.y / INFO->transy;
  g_event->deltax = event.xconfigure.width / INFO->transx;
  g_event->deltay = event.xconfigure.height / INFO->transy;
  if (!check_redraw (&redraw_list, TRUE, g_event))
    ignore_event = TRUE;
      } else {
  swallow_redraw (event.xany.display, event.xany.window);
  g_event->what = evRedraw;
      }
      break;

    case KeyPress:
      {
  char ch;
  KeySym keysym;
  /* number of returned characters */
  int numret = XLookupString((XKeyEvent *) &event,&ch,1,&keysym,NULL);

  g_event->what = evCtrlKey;
  g_event->x     = event.xkey.x / INFO->transx;
  g_event->y     = event.xkey.y / INFO->transy;
  g_event->state = event.xkey.state;

  switch (keysym) {
#ifdef XK_KP_Enter
  case XK_KP_Enter:
#endif
#ifdef XK_ISO_Enter
  case XK_ISO_Enter:
#endif
#ifdef XK_3270_Enter
  case XK_3270_Enter:
#endif
  case XK_Return:
    g_event->which = ckReturn;
    break;
  case XK_BackSpace:
    g_event->which = ckBackspace;
    break;
#ifdef XK_KP_Tab
  case XK_KP_Tab:
#endif
  case XK_Tab:
    g_event->which = ckTab;
    break;
#if defined (XK_ISO_Left_Tab) || defined (XK_3270_BackTab)
#ifdef XK_3270_BackTab
  case XK_3270_BackTab:
#endif
#ifdef XK_ISO_Left_Tab
  case XK_ISO_Left_Tab:
#endif
    g_event->which = ckTab;
    g_event->state |= sfShiftKey;
    break;
#endif
  case XK_Escape:
    g_event->which = ckEscape;
    break;
  case XK_Delete:
    g_event->which = ckDelete;
    break;
  case XK_Home:
    g_event->which = ckHome;
    break;
  case XK_Left:
    g_event->which = ckLeftArrow;
    break;
  case XK_Up:
    g_event->which = ckUpArrow;
    break;
  case XK_Right:
    g_event->which = ckRightArrow;
    break;
  case XK_Down:
    g_event->which = ckDownArrow;
    break;
  case XK_Page_Up:
    g_event->which = ckPageUp;
    break;
  case XK_Page_Down:
    g_event->which = ckPageDown;
    break;
  case XK_End:
    g_event->which = ckEnd;
    break;
  case XK_Begin:
    g_event->which = ckBegin;
    break;
  case XK_Insert:
    XConvertSelection (INFO->display, XA_PRIMARY, XA_STRING,
           INFO->insert_prop, INFO->window,
           event.xkey.time);
    INFO->paste_x = event.xkey.x / INFO->transx;
    INFO->paste_y = event.xkey.y / INFO->transy;
    ignore_event = TRUE;
    break;
  case XK_F1:
    g_event->which = ckF1;
    break;
  case XK_F2:
    g_event->which = ckF2;
    break;
  case XK_F3:
    g_event->which = ckF3;
    break;
  case XK_F4:
    g_event->which = ckF4;
    break;
  case XK_F5:
    g_event->which = ckF5;
    break;
  case XK_F6:
    g_event->which = ckF6;
    break;
  case XK_F7:
    g_event->which = ckF7;
    break;
  case XK_F8:
    g_event->which = ckF8;
    break;
  case XK_F9:
    g_event->which = ckF9;
    break;
  case XK_F10:
    g_event->which = ckF10;
    break;
  case XK_F11:
    g_event->which = ckF11;
    break;
  case XK_F12:
    g_event->which = ckF12;
    break;

  default:
    if (numret) {
      if (isprint(ch)) {
        g_event->which = ch;
        g_event->what  = evNormalKey;
      } else if (ch && iscntrl (ch)) {
        g_event->which = ch + '@';
        g_event->what  = evNormalKey;
      } else
        ignore_event = TRUE;
    } else
      ignore_event = TRUE;
  }
  if (!ignore_event)
    INFO->last_event_time = event.xkey.time;

  break;
      }

    case SelectionNotify:
      if (event.xselection.property == INFO->insert_prop &&
    event.xselection.target == XA_STRING) {
  Atom actual_type_return;
  int actual_format_return;
  unsigned long nitems_return, bytes_after_return;
  unsigned char *prop_return = NULL;

  if (Success == XGetWindowProperty
      (INFO->display, INFO->window, INFO->insert_prop,
       0, LONG_MAX / 4, True, XA_STRING,
       &actual_type_return, &actual_format_return,
       &nitems_return, &bytes_after_return,
       &prop_return) &&
      actual_type_return == XA_STRING &&
      actual_format_return == CHAR_BIT &&
      nitems_return && prop_return)
    INFO->pasted_text = g_strcat (INFO->pasted_text,
          (char *)prop_return);
  if (prop_return)
    XFree (prop_return);
      }
      ignore_event = TRUE;
      break;

    case SelectionRequest:
      {
  XEvent sevent;
  memset (&sevent, 0, sizeof sevent);
  sevent.xselection.type = SelectionNotify;
  sevent.xselection.requestor = event.xselectionrequest.requestor;
  sevent.xselection.selection = event.xselectionrequest.selection;
  sevent.xselection.target = event.xselectionrequest.target;
  sevent.xselection.property = None;
  sevent.xselection.time = event.xselectionrequest.time;
  if (event.xselectionrequest.selection == XA_PRIMARY &&
      event.xselectionrequest.target == XA_STRING) {
    XChangeProperty (event.xselectionrequest.display,
         event.xselectionrequest.requestor,
         event.xselectionrequest.property,
         XA_STRING, CHAR_BIT, PropModeReplace,
         ( unsigned char * ) ( INFO->copied_text ),
         strlen (INFO->copied_text));
    sevent.xselection.property = event.xselectionrequest.property;
  }
  XSendEvent(event.xselectionrequest.display,
       event.xselectionrequest.requestor, False,
       0, &sevent);
      }
      ignore_event = TRUE;
      break;

    case ButtonPress:
    case ButtonRelease:
      {
  unsigned old_buttons = INFO->buttons;

  g_event->x    = event.xbutton.x / INFO->transx;
  g_event->y    = event.xbutton.y / INFO->transy;

  g_event->state  = event.xbutton.state;
  g_event->which  = 1<<(7+event.xbutton.button);

  INFO->buttons = (g_event->state ^ g_event->which) & sfSomeMouse;

  if (event.type == ButtonPress) {
    g_event->what = evMousePressed;
    if (INFO->buttons && !old_buttons)
      XDefineCursor(INFO->display, INFO->window,
        XCreateFontCursor(INFO->display,
              XC_clock)); /* busy mode */
  }
  else {
    g_event->what = evMouseReleased;
    g_event->deltax = g_event->x- INFO->event_oldx;
    g_event->deltay = g_event->y- INFO->event_oldy;
    if (old_buttons && !INFO->buttons)
      XDefineCursor(INFO->display, INFO->window,
        XCreateFontCursor(INFO->display,
              XC_hand2)); /* ready mode */
  }
  INFO->event_oldx = g_event->x;
  INFO->event_oldy = g_event->y;
      }

      if (!ignore_event)
  INFO->last_event_time = event.xbutton.time;

      break;

    case MotionNotify:
      if (INFO->buttons) {
  g_event->what = evMouseDragged;
  g_event->x = event.xmotion.x / INFO->transx;
  g_event->y = event.xmotion.y / INFO->transy;
  g_event->state = g_event->which = INFO->buttons;

  while (XCheckIfEvent(INFO->display, &event, next_window_event, NULL))
    if (event.type == MotionNotify) {
      g_event->x = event.xmotion.x / INFO->transx;
      g_event->y = event.xmotion.y / INFO->transy;
    } else {
      XPutBackEvent (INFO->display, &event);
      break;
    }

  g_event->deltax = g_event->x - INFO->event_oldx;
  g_event->deltay = g_event->y - INFO->event_oldy;
  INFO->event_oldx = g_event->x;
  INFO->event_oldy = g_event->y;

  if (!g_event->deltax && !g_event->deltay)
    g_event->what = evMouseHold;
      }
      else
  ignore_event = TRUE;

      if (!ignore_event)
  INFO->last_event_time = event.xmotion.time;

      break;

    default:
#ifdef G_X11_SPACEBALL
      {
  int which =
    ((INFO->spaceball_motion_type &&
      event.type == INFO->spaceball_motion_type) ? 1 : 0) +
    ((INFO->spaceball_button_press_type &&
      event.type == INFO->spaceball_button_press_type) ? 2 : 0) +
    ((INFO->spaceball_button_release_type &&
      event.type == INFO->spaceball_button_release_type) ? 3 : 0);

  if (which) {
    if (which == 1) {
      int event_type = event.type;
      XEvent event1;
      while (XCheckIfEvent(INFO->display, &event1, next_window_event,
         NULL))
        if (event.type == event_type) {
    event = event1;
        } else {
    XPutBackEvent (INFO->display, &event1);
    break;
        }
    }
    g_event->what = evSpaceball;
    if (which == 1) {
      /* see X11/extensions/XInput.h */
      XDeviceMotionEvent *extevent = (XDeviceMotionEvent *)&event;
      int axis;
      for (axis = 0; axis < 6; axis++)
        INFO->spaceball_data.axis[axis] =
    (axis >= extevent->first_axis &&
     axis < extevent->first_axis + extevent->axes_count) ?
    extevent->axis_data[axis] : 0;
      INFO->spaceball_data.time_diff =
        INFO->spaceball_data.time ?
        extevent->time - INFO->spaceball_data.time : 0;
      INFO->spaceball_data.time = extevent->time;
      break;
    } else {
      XDeviceButtonEvent *extevent = (XDeviceButtonEvent *)&event;
      int axis;
      for (axis = 0; axis < 6; axis++)
        INFO->spaceball_data.axis[axis] =
    (axis >= extevent->first_axis &&
     axis < extevent->first_axis + extevent->axes_count) ?
    extevent->axis_data[axis] : 0;
      INFO->spaceball_data.button[extevent->button] = 3 - which;
      INFO->spaceball_data.time_diff =
        INFO->spaceball_data.time ?
        extevent->time - INFO->spaceball_data.time : 0;
      INFO->spaceball_data.time = extevent->time;
      if (extevent->button && which == 2) {
        g_event->what = evNormalKey;
        g_event->which = "w012lstm"[extevent->button-1];
        g_event->state = extevent->state;
        g_event->x = extevent->x / INFO->transx;
        g_event->y = extevent->y / INFO->transy;
      } else
        ignore_event = TRUE;
      break;
    }
  }
      }
#endif

      if (INFO->buttons) {
  g_sleep_usec (1000000/30);
  g_event->x = INFO->event_oldx;
  g_event->y = INFO->event_oldy;
  g_event->deltax = 0.0;
  g_event->deltay = 0.0;
  g_event->state = g_event->which = INFO->buttons;
      }
      else
  g_event->what = evNothing;
      break;
    }
  } while (ignore_event);

  return 1;
}

#ifdef G_X11_SPACEBALL
void g_x11_ctl_get_spaceball_data(double *period,
          double *tx, double *ty, double *tz,
          double *rx, double *ry, double *rz,
          int *button)
{
  const double scale = 1e2, rot_scale = 1e1;
  *period = INFO->spaceball_data.time_diff * 1e0;
  *tx = INFO->spaceball_data.axis[0]*scale;
  *ty = INFO->spaceball_data.axis[1]*scale;
  *tz = INFO->spaceball_data.axis[2]*scale;
  *rx = INFO->spaceball_data.axis[3]*rot_scale;
  *ry = INFO->spaceball_data.axis[4]*rot_scale;
  *rz = INFO->spaceball_data.axis[5]*rot_scale;
  *button = INFO->spaceball_data.button[0];
}
#endif

void g_x11_ctl_copy_string (const char *string)
{
  INFO->copied_text = g_strchange (INFO->copied_text, string);
  if (TRUE) { /* newline fuer Mehrzeiler?? */
    char *pos;
    if ((pos = strrchr (INFO->copied_text, '\n')) && pos[1] != '\0')
      INFO->copied_text = g_strcat (INFO->copied_text, "\n");
  }
  if (!g_strempty (INFO->copied_text)) {
    XSetSelectionOwner(INFO->display, XA_PRIMARY, INFO->window,
           INFO->last_event_time);
#if 0
    /* wenn's nicht geklappt hat, merkt man das eh bald */
    if (XGetSelectionOwner(INFO->display, XA_PRIMARY) != INFO->window) {
      fprintf (stderr, "/* We didn't get the selection */\n");
    }
#endif

    /* zusaetzlich (fuer alte clients?) einfuegen in den cut-buffer ring */
    XRotateBuffers (INFO->display, 1);
    XStoreBytes (INFO->display, INFO->copied_text, strlen (INFO->copied_text));
  }
}

double g_x11_ctl_stringwidth(const char *name)
{
  /* returns width of string in multiples of patch width */
  return XTextWidth (INFO->font_info, name, (int)strlen(name)) /
    INFO->x_pixel_per_patch;
}


double g_x11_ctl_stringheight(const char *name)
{
  /* returns height of string in multiples of patch height */
  return INFO->textheight / INFO->y_pixel_per_patch;
}


double g_x11_ctl_get_unit_per_pixel(void)
{
  return 1.0 / INFO->x_pixel_per_patch;
}

#define G_11_MAXRECTS 1000

G_CLIP_DESCRIPTOR *g_x11_ctl_kind_of_clipping( void )
{
  static G_CLIP_DESCRIPTOR g_descript;
  g_descript.which           = g_ctEnableRects;
  g_descript.max_no_of_rects = G_11_MAXRECTS;

  return(&g_descript);
}

void g_x11_ctl_set_clip_rects (GRECT_LIST *rect)
{
  if (!rect)
    XSetClipMask( INFO->display, INFO->gr_context, None);
  else {
    GRECT_LIST *runner;
    int counter = 0;
    XRectangle g_x11_rect[G_11_MAXRECTS];
    for (runner = rect; runner; runner = runner->next) {
      if (counter >= G_11_MAXRECTS) {
  fprintf(stderr, "Too many (>%d) Clip-Rectangles. Ignoring last.\n",
    G_11_MAXRECTS);
  break;
      }
      g_x11_rect[ counter].x     = (short) (runner->rect->x1*INFO->transx+.5);
      g_x11_rect[ counter].y     = (short) (runner->rect->y1*INFO->transy+.5);
      g_x11_rect[ counter].width = (unsigned short)
  ((runner->rect->x2*INFO->transx+ .5))- g_x11_rect[ counter].x +1;
      g_x11_rect[ counter].height= (unsigned short)
  ((runner->rect->y2)*INFO->transy+ .5)- g_x11_rect[ counter].y +1;
      counter++;
    }
    XSetClipRectangles (INFO->display, INFO->gr_context, 0, 0,
      g_x11_rect, counter, Unsorted);
  }
}
