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

#include "grapeInterface1d.h"

#ifdef USE_EXTERNAL_GRAPE

TRIANG1D* quocmesh_convert_to_triang1d(qc::ScalarArray<double, qc::QC_1D>* triangData, const char* dataName)
{
  int n = triangData->getNumX();
  TRIANG1D* triang = reinterpret_cast<TRIANG1D*> (GRAPE (Triang1d, "new-instance") (dataName));
  GRAPE (triang, "list-alloc") (n);
  for (int i = 0; i < n; ++i) {
    triang->x [i] = static_cast<double> (i) / (n-1);
    triang->y [i] = triangData->get (i);
    triang->z [i] = 0;
  }
  triang->number_of_points = n;

  return triang;
}

TRIANG1D* triang1d_rescaled_disp ()
{
  TRIANG1D* self = reinterpret_cast<TRIANG1D*> (START_METHOD (G_INSTANCE));
  ASSURE (self, "", END_METHOD (NULL));


  static GROUP* group = NULL;
  static double scale = 1;

  if (!group)
    group = reinterpret_cast<GROUP*>
      (new_item (Group, I_Name, "Triang1d::scaled-disp",
     I_RSizeX, 1.0, I_SizeY, 1.0 + 1.25,
     I_FillMode, mfNorthWest,
     I_Border, bfBorder | bfTitle,
     I_Item, new_item (Function_Ruler, I_Name, "scale",
           I_FillMode, mfNorthWest,
           I_Var, &scale, dfDouble,
           I_DefaultValue, 1.0,
           I_MinMax, -10.0, 10.0,
           I_StepSize, 0.1,
           I_Offset, 0.0,
           I_Scale, 1.0,
           I_RSizeX, 1.0,
           I_End),
     I_End));

  MANAGER* mgr = reinterpret_cast<MANAGER*> (GRAPE (Manager, "get-stdmgr") ());
  if (GRAPE (mgr, "new-handle") (triang1d_rescaled_disp, 1)) GRAPE (mgr, "add-inter") (group);

  TRIANG1D* temp = reinterpret_cast<TRIANG1D*> (GRAPE (Triang1d, "new-instance") ("temp"));
  GRAPE (temp, "list-growto") (self->number_of_points);
  temp->number_of_points = self->number_of_points;

  for (int i = 0; i < temp->number_of_points; ++i) {
    temp->x [i] = self->x [i];
    temp->y [i] = scale * self->y [i];
    temp->z [i] = 0;
  }

  GRAPE (temp, "display") ();
  GRAPE (temp, "delete") ();

  END_METHOD (self);
}

TRIANG1D* triang1d_height_disp ()
{
  TRIANG1D* self = reinterpret_cast<TRIANG1D*> (START_METHOD (G_INSTANCE));
  ASSURE (self, "", END_METHOD (NULL));

  GRAPHICDEVICE* dev = reinterpret_cast<GRAPHICDEVICE*> (GRAPE (GraphicDevice, "get-stddev") ());

  static GROUP* group = NULL;
  static double scale = 1;

  if (!group)
    group = reinterpret_cast<GROUP*>
      (new_item (Group, I_Name, "Triang1d::height-disp",
     I_RSizeX, 1.0, I_SizeY, 1.0 + 1.25,
     I_FillMode, mfNorthWest,
     I_Border, bfBorder | bfTitle,
     I_Item, new_item (Function_Ruler, I_Name, "scale",
           I_FillMode, mfNorthWest,
           I_Var, &scale, dfDouble,
           I_DefaultValue, 1.0,
           I_MinMax, -10.0, 10.0,
           I_StepSize, 0.1,
           I_Offset, 0.0,
           I_Scale, 1.0,
           I_RSizeX, 1.0,
           I_End),
     I_End));

  MANAGER* mgr = reinterpret_cast<MANAGER*> (GRAPE (Manager, "get-stdmgr") ());
  if (GRAPE (mgr, "new-handle") (triang1d_height_disp, 1)) GRAPE (mgr, "add-inter") (group);

  if (dev->grid_patch == G_GRID) {
    VEC3 p;

    g_vec3_set (p, self->x [0], scale * self->y [0], 0);
    dev->move (p);

    for (int i = 1; i < self->number_of_points; ++i) {

      g_vec3_set (p, self->x [i], scale * self->y [i], 0);
      dev->draw (p);
    }

    for (int i = self->number_of_points - 1; i >= 0; --i) {

      g_vec3_set (p, self->x [i], 0, 0);
      dev->draw (p);
    }

    g_vec3_set (p, self->x [0], scale * self->y [0], 0);
    dev->draw (p);
  }
  else {

    VEC3 p [4]; VEC3 n;
    for (int i = 1; i < self->number_of_points; ++i) {

      if (self->y [i-1] * self->y [i] >= 0.0) {

  g_vec3_set (p [0], self->x [i-1], 0, 0);
  g_vec3_set (p [1], self->x [i-1], scale * self->y [i-1], 0);
  g_vec3_set (p [2], self->x [i], scale * self->y [i], 0);
  g_vec3_set (p [3], self->x [i], 0, 0);

  if (p [1][2] != 0.0) g_vec3_get_normal_to_plane_quietly (n, p [0], p [1], p [3]);
  else if (p [2][2] != 0.0) g_vec3_get_normal_to_plane_quietly (n, p [0], p [2], p [3]);
  else g_vec3_set (n, 0.0, 0.0, 0.0);

  dev->begin_patch ();
  dev->patch_normal (n);
  dev->patch_vertex (p [0]);
  dev->patch_vertex (p [1]);
  dev->patch_vertex (p [2]);
  dev->patch_vertex (p [3]);
  dev->end_patch ();
      }
      else {
  double alpha = abs (self->y [i-1]) / (abs (self->y [i-1]) + abs (self->y [i])), beta = 1.0 - alpha;

  g_vec3_set (p [0], self->x [i-1], 0, 0);
  g_vec3_set (p [1], self->x [i-1], scale * self->y [i-1], 0);
  g_vec3_set (p [2], beta * self->x [i-1] + alpha * self-> x [i], 0, 0);
  g_vec3_get_normal_to_plane_quietly (n, p [0], p [1], p [2]);

  dev->begin_patch ();
  dev->patch_normal (n);
  dev->patch_vertex (p [0]);
  dev->patch_vertex (p [1]);
  dev->patch_vertex (p [2]);
  dev->end_patch ();

  g_vec3_set (p [0], self->x [i], 0, 0);
  g_vec3_set (p [1], self->x [i], scale * self->y [i], 0);
  g_vec3_set (p [2], beta * self->x [i-1] + alpha * self-> x [i], 0, 0);
  g_vec3_get_normal_to_plane_quietly (n, p [0], p [1], p [2]);

  dev->begin_patch ();
  dev->patch_normal (n);
  dev->patch_vertex (p [0]);
  dev->patch_vertex (p [1]);
  dev->patch_vertex (p [2]);
  dev->end_patch ();
      }
    }
  }

  END_METHOD (self);
}

TRIANG1D* triang1d_color_disp ()
{
  TRIANG1D* self = reinterpret_cast<TRIANG1D*> (START_METHOD (G_INSTANCE));
  ASSURE (self, "", END_METHOD (NULL));

  GRAPHICDEVICE* dev = reinterpret_cast<GRAPHICDEVICE*> (GRAPE (GraphicDevice, "get-stddev") ());

  static GROUP* group = NULL;
  static COLORBAR* colorbar = NULL;
  static double scale = 1;

  if (!group) {
    colorbar = reinterpret_cast<COLORBAR*> (GRAPE (Colorbar, "get-stdcolorbar") (triang1d_color_disp, "Triang1d::color-disp"));
    group = reinterpret_cast<GROUP*>
      (new_item (Group, I_Name, "Triang1d::color-disp",
     I_RSizeX, 1.0, I_SizeY, 2.0 + 1.25,
     I_FillMode, mfNorthWest,
     I_Border, bfBorder | bfTitle,
     I_Item, new_item (Function_Ruler, I_Name, "scale",
           I_FillMode, mfNorthWest,
           I_Var, &scale, dfDouble,
           I_DefaultValue, 1.0,
           I_MinMax, -10.0, 10.0,
           I_StepSize, 0.1,
           I_Offset, 0.0,
           I_Scale, 1.0,
           I_RSizeX, 1.0,
           I_End),
     I_Item, GRAPE (colorbar, "get-button") (),
     I_End));
  }

  MANAGER* mgr = reinterpret_cast<MANAGER*> (GRAPE (Manager, "get-stdmgr") ());
  if (GRAPE (mgr, "new-handle") (triang1d_color_disp, 1)) GRAPE (mgr, "add-inter") (group);

  if (dev->grid_patch == G_GRID) {
    VEC3 p; VEC3 c; VEC3 co;
    dev->attribute (G_MODE_GET, G_LINE_COLOR, co);

    g_get_color_for_value (colorbar, self->y [0], c);
    g_vec3_set (p, self->x [0], scale * self->y [0], 0);
    dev->attribute (G_MODE_SET, G_LINE_COLOR, c);
    dev->move (p);

    for (int i = 1; i < self->number_of_points; ++i) {

      g_get_color_for_value (colorbar, 0.5 * (self->y [i-1] + self->y [i]), c);
      dev->attribute (G_MODE_SET, G_LINE_COLOR, c);
      dev->move (p); // !
      g_vec3_set (p, self->x [i], scale * self->y [i], 0);
      dev->draw (p);
    }

    dev->attribute (G_MODE_SET, G_LINE_COLOR, co);

    for (int i = self->number_of_points - 1; i >= 0; --i) {

      g_vec3_set (p, self->x [i], 0, 0);
      dev->draw (p);
    }

    g_vec3_set (p, self->x [0], scale * self->y [0], 0);
    dev->draw (p);
  }
  else if (dev->grid_patch == G_PATCH) {

    int cm, on = ON;
    dev->attribute (G_MODE_GET, G_COLOR_MATERIAL, &cm);
    dev->attribute (G_MODE_SET, G_COLOR_MATERIAL, &on);

    VEC3 p [4]; VEC3 n; VEC3 c [2];
    for (int i = 1; i < self->number_of_points; ++i) {

      if (self->y [i-1] * self->y [i] >= 0.0) {

  g_vec3_set (p [0], self->x [i-1], 0, 0);
  g_vec3_set (p [1], self->x [i-1], scale * self->y [i-1], 0);
  g_vec3_set (p [2], self->x [i], scale * self->y [i], 0);
  g_vec3_set (p [3], self->x [i], 0, 0);
  g_get_color_for_value (colorbar, self->y [i-1], c [0]);
  g_get_color_for_value (colorbar, self->y [i], c [1]);

  if (p [1][2] != 0.0) g_vec3_get_normal_to_plane_quietly (n, p [0], p [1], p [3]);
  else if (p [2][2] != 0.0) g_vec3_get_normal_to_plane_quietly (n, p [0], p [2], p [3]);
  else g_vec3_set (n, 0.0, 0.0, 0.0);

  dev->begin_patch ();
  dev->patch_normal (n);
  dev->patch_color (c [0]);
  dev->patch_vertex (p [0]);
  dev->patch_vertex (p [1]);
  dev->patch_color (c [1]);
  dev->patch_vertex (p [2]);
  dev->patch_vertex (p [3]);
  dev->end_patch ();
      }
      else {
  double alpha = abs (self->y [i-1]) / (abs (self->y [i-1]) + abs (self->y [i])), beta = 1.0 - alpha;

  g_vec3_set (p [0], self->x [i-1], 0, 0);
  g_vec3_set (p [1], self->x [i-1], scale * self->y [i-1], 0);
  g_vec3_set (p [2], beta * self->x [i-1] + alpha * self-> x [i], 0, 0);
  g_vec3_get_normal_to_plane_quietly (n, p [0], p [1], p [2]);
  g_get_color_for_value (colorbar, self->y [i-1], c [0]);
  g_get_color_for_value (colorbar, 0, c [1]);

  dev->begin_patch ();
  dev->patch_normal (n);
  dev->patch_color (c [0]);
  dev->patch_vertex (p [0]);
  dev->patch_vertex (p [1]);
  dev->patch_color (c [1]);
  dev->patch_vertex (p [2]);
  dev->end_patch ();

  g_vec3_set (p [0], self->x [i], 0, 0);
  g_vec3_set (p [1], self->x [i], scale * self->y [i], 0);
  g_vec3_set (p [2], beta * self->x [i-1] + alpha * self-> x [i], 0, 0);
  g_vec3_get_normal_to_plane_quietly (n, p [0], p [1], p [2]);
  g_get_color_for_value (colorbar, self->y [i], c [0]);
  g_get_color_for_value (colorbar, 0, c [1]);

  dev->begin_patch ();
  dev->patch_normal (n);
  dev->patch_color (c [0]);
  dev->patch_vertex (p [0]);
  dev->patch_vertex (p [1]);
  dev->patch_color (c [1]);
  dev->patch_vertex (p [2]);
  dev->end_patch ();
      }
    }

    dev->attribute (G_MODE_SET, G_COLOR_MATERIAL, &cm);
  }
  else {

    int handle = g_colorbar_get_texture (colorbar, dev);
    g_texture_set_attribute (dev, G_TEXTURE_MODE, handle, G_TEXTURE_MODULATE_MODE);
    dev->begin_texture (handle);

    VEC3 p [4]; VEC3 n; VEC3 c [2];
    for (int i = 1; i < self->number_of_points; ++i) {

      if (self->y [i-1] * self->y [i] >= 0.0) {

  g_vec3_set (p [0], self->x [i-1], 0, 0);
  g_vec3_set (p [1], self->x [i-1], scale * self->y [i-1], 0);
  g_vec3_set (p [2], self->x [i], scale * self->y [i], 0);
  g_vec3_set (p [3], self->x [i], 0, 0);
  g_colorbar_get_texture_coord (colorbar, dev, self->y [i-1], c [0]);
  g_colorbar_get_texture_coord (colorbar, dev, self->y [i], c [1]);

  if (p [1][2] != 0.0) g_vec3_get_normal_to_plane_quietly (n, p [0], p [1], p [3]);
  else if (p [2][2] != 0.0) g_vec3_get_normal_to_plane_quietly (n, p [0], p [2], p [3]);
  else g_vec3_set (n, 0.0, 0.0, 0.0);

  dev->begin_patch ();
  dev->patch_normal (n);
  dev->texture_coord (c [0]);
  dev->patch_vertex (p [0]);
  dev->patch_vertex (p [1]);
  dev->texture_coord (c [1]);
  dev->patch_vertex (p [2]);
  dev->patch_vertex (p [3]);
  dev->end_patch ();
      }
      else {
  double alpha = abs (self->y [i-1]) / (abs (self->y [i-1]) + abs (self->y [i])), beta = 1.0 - alpha;

  g_vec3_set (p [0], self->x [i-1], 0, 0);
  g_vec3_set (p [1], self->x [i-1], scale * self->y [i-1], 0);
  g_vec3_set (p [2], beta * self->x [i-1] + alpha * self-> x [i], 0, 0);
  g_vec3_get_normal_to_plane_quietly (n, p [0], p [1], p [2]);
  g_colorbar_get_texture_coord (colorbar, dev, self->y [i-1], c [0]);
  g_colorbar_get_texture_coord (colorbar, dev, 0, c [1]);

  dev->begin_patch ();
  dev->patch_normal (n);
  dev->texture_coord (c [0]);
  dev->patch_vertex (p [0]);
  dev->patch_vertex (p [1]);
  dev->texture_coord (c [1]);
  dev->patch_vertex (p [2]);
  dev->end_patch ();

  g_vec3_set (p [0], self->x [i], 0, 0);
  g_vec3_set (p [1], self->x [i], scale * self->y [i], 0);
  g_vec3_set (p [2], beta * self->x [i-1] + alpha * self-> x [i], 0, 0);
  g_vec3_get_normal_to_plane_quietly (n, p [0], p [1], p [2]);
  g_colorbar_get_texture_coord (colorbar, dev, self->y [i], c [0]);
  g_colorbar_get_texture_coord (colorbar, dev, 0, c [1]);

  dev->begin_patch ();
  dev->patch_normal (n);
  dev->texture_coord (c [0]);
  dev->patch_vertex (p [0]);
  dev->patch_vertex (p [1]);
  dev->texture_coord (c [1]);
  dev->patch_vertex (p [2]);
  dev->end_patch ();
      }
    }

    dev->end_texture ();
  }

  GRAPE (colorbar, "display") ();

  END_METHOD (self);
}

void addMethodsAndProjects1d()
{
  GRAPE (Triang1d, "add-method") ("rescaled-disp", triang1d_rescaled_disp);
  GRAPE (Triang1d, "add-method") ("height-disp", triang1d_height_disp);
  GRAPE (Triang1d, "add-method") ("color-disp", triang1d_color_disp);
}

void initStartGrape(TRIANG1D* triang)
{
  addMethodsAndProjects1d();

  MANAGER*  mgr;
  mgr = reinterpret_cast<MANAGER*>(GRAPE(Manager,"get-stdmgr")());

  GRAPE(mgr,"handle")(triang);
}

#endif
