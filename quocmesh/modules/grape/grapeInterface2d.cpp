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

#include "grapeInterface2d.h"

#ifdef USE_EXTERNAL_GRAPE

/* structure for the function: it contains
  the original function data and an estimator */
struct functionComponents
{
  qc::ScalarArray<double, qc::QC_2D>* data;
  qc::Estimator2d<double>* estimator;
  double (*get_element_estimate) (HELEMENT2D *, void *function_data);
  void (*get_vertex_estimate) (HELEMENT2D *, double *, void *function_data);
  std::vector<qc::ScalarArray<double, qc::QC_2D>*>* data_seq;
  std::vector<qc::Estimator2d<double>*>* estimator_seq;
};

struct functionComponentsMV
{
  aol::MultiVector<double>* data;
  std::vector<aol::MultiVector<double>*>* data_seq;
};

struct MESH_INFO {
  qc::ScalarArray<double, qc::QC_2D>* data;
  char *filename;

  double t;
};

/* two vectors for managing the elements2d */
vector<HELEMENT2D*> elements_used;
vector<HELEMENT2D*> elements_free;


/* vertex inheritance structures for the cubes - first: the weights */
/* ATTENTION: For using the vinherit structures the vertex-estimator
  must be initialized!!! */
static double weightOne[1] = {1.0};
static double weightTwo[2] = {0.5,0.5};
static double weightFour[4] = {0.25,0.25,0.25,0.25};

/* second: the structures itself */
int indices0[1] = {0}; VINHERIT vinh0 = {1,indices0,weightOne};
int indices1[2] = {0,1}; VINHERIT vinh1 = {2,indices1,weightTwo};
int indices2[1] = {1}; VINHERIT vinh2 = {1,indices2,weightOne};
int indices3[2] = {3,0}; VINHERIT vinh3 = {2,indices3,weightTwo};
int indices4[4] = {0,1,2,3}; VINHERIT vinh4 = {4,indices4,weightFour};
int indices5[2] = {1,2}; VINHERIT vinh5 = {2,indices5,weightTwo};
int indices6[1] = {3}; VINHERIT vinh6 = {1,indices6,weightOne};
int indices7[2] = {2,3}; VINHERIT vinh7 = {2,indices7,weightTwo};
int indices8[1] = {2}; VINHERIT vinh8 = {1,indices8,weightOne};

/* third: the arrays of vinherit structures for the child-elements */
G_CONST VINHERIT vinhChild0[4] = {vinh0,vinh1,vinh4,vinh3};
G_CONST VINHERIT vinhChild1[4] = {vinh1,vinh2,vinh5,vinh4};
G_CONST VINHERIT vinhChild2[4] = {vinh3,vinh4,vinh7,vinh6};
G_CONST VINHERIT vinhChild3[4] = {vinh4,vinh5,vinh8,vinh7};



/* local coordinates for the vertices of the square */
static G_CONST double s_c0[2] = {0.,0.};
static G_CONST double s_c1[2] = {1.,0.};
static G_CONST double s_c2[2] = {1.,1.};
static G_CONST double s_c3[2] = {0.,1.};

/* the coord-system */
static G_CONST double *square_coord_system[4] = {s_c0,s_c1,s_c2,s_c3};

/* The functions for the description of the elements */
static void        square_free_element   (HELEMENT2D*);
static int         square_boundary       (HELEMENT2D*,int);
static int         square_world_to_coord (HELEMENT2D*,double G_CONST*,double*);
static int         square_check_inside   (HELEMENT2D*,double G_CONST*);
static void        square_coord_to_world (HELEMENT2D*,double G_CONST*,double*);
static void        square_coord_of_parent(HELEMENT2D*,double*,double*);
static HELEMENT2D* square_neighbour      (HELEMENT2D*,int,int,double*,double*,MESH_ELEMENT_FLAGS);
static HELEMENT2D* get_element(GENMESH2D* gmesh);

static void scalar_function(ELEMENT2D*,int,double[],double[],void*);
static void vector_function(ELEMENT2D*,int,double[],double[],void*);
static void square_f_el_info (ELEMENT2D*,F_EL_INFO2D *,void*);


/* now the DESCRIPTION of the elements, here simple squares */
static HELEMENT2D_DESCRIPTION square_description;

/* Extend copy methods to also copy user data and functions */
GENMESH2D *genmesh2d_my_softcopy (GENMESH2D *copy)
{
  GENMESH2D *self;

  self = reinterpret_cast<GENMESH2D*> (START_METHOD (G_INSTANCE));
  ASSURE (self, "", END_METHOD (NULL));

  if ((copy = reinterpret_cast<GENMESH2D*> (GRAPE(self, "^softcopy")(copy)))) {
    copy->dimension_of_world = self->dimension_of_world;

    /* function pointers should be set in new-instance */
    copy->first_child = self->first_child;
    copy->next_child = self->next_child;
    copy->first_macro = self->first_macro;
    copy->next_macro = self->next_macro;
    copy->select_child = self->select_child;

    copy->copy_element = self->copy_element;
    copy->free_element = self->free_element;
    copy->complete_element = self->complete_element;
    copy->set_time = self->set_time;
    copy->get_time = self->get_time;

    copy->max_level = self->max_level;
    copy->level_of_interest = self->level_of_interest;

    /* function pointers should be set in new-instance */
    copy->get_geometry_vertex_estimate = self->get_geometry_vertex_estimate;
    copy->get_geometry_element_estimate = self->get_geometry_element_estimate;

    copy->get_lens_element_estimate = self->get_lens_element_estimate;

    copy->threshold = self->threshold;

    /* extension */
    copy->user_data = self->user_data;
    copy->f_data = self->f_data;
  }

  END_METHOD (copy);
}

/* fill the description with some useful things */
void fillDescription(void)
{
  square_description.dindex = 0;
  square_description.number_of_vertices = 4;
  square_description.dimension_of_coord = 2;
  square_description.coord = square_coord_system;
  square_description.parametric_degree = 1;

  square_description.world_to_coord = square_world_to_coord;
  square_description.coord_to_world = square_coord_to_world;
  square_description.check_inside = square_check_inside;
  square_description.neighbour = square_neighbour;
  square_description.boundary = square_boundary;
  square_description.coord_of_parent = square_coord_of_parent;
}

/* filling an element with some useful informations */
static HELEMENT2D* fillHElement2D(HELEMENT2D* helement, GENMESH2D* gmesh)
{
  helement->mesh = gmesh;
  helement->eindex = 0;
  helement->descr = reinterpret_cast<ELEMENT2D_DESCRIPTION*>(&square_description);
  helement->present = static_cast<MESH_ELEMENT_FLAGS> (hefVertex | hefVindex | hefVinh);
  helement->user_data = NULL;
  helement->parent = NULL;
  helement->vinh = NULL;
  helement->ref_rule = -1;
  helement->level = 0;
  helement->has_children = 1;
  return(helement);
}

/* filling the function data with some useful informations with the NEW fC-Structure */
void fillScalarDataFunction(F_HDATA2D* scal_function, functionComponents* fC)
{
  scal_function->name                 = const_cast<char*>("Skalare Daten");
  scal_function->dimension_of_value   = 1;
  scal_function->continuous_data      = 1;
  scal_function->f                    = scalar_function;
  scal_function->last                 = NULL;
  scal_function->next                 = NULL;
  scal_function->get_bounds           = NULL;
  scal_function->threshold            = 0.0;
  scal_function->hp_threshold         = 0.0;
  scal_function->hp_maxlevel          = 0;
  scal_function->f_el_info            = square_f_el_info;
  scal_function->get_vertex_estimate  = fC->get_vertex_estimate;
  scal_function->get_element_estimate = fC->get_element_estimate;
  scal_function->function_data        = fC;
}


/* filling the function data with some useful informations */
void fillVectorDataFunction(F_HDATA2D* vec_function, functionComponentsMV* fC)
{
  vec_function->name                 = const_cast<char*>("Vektordaten");
  vec_function->dimension_of_value   = fC->data->numComponents();
  vec_function->continuous_data      = 1;
  vec_function->f                    = vector_function;
  vec_function->last                 = NULL;
  vec_function->next                 = NULL;
  vec_function->get_bounds           = NULL;
  vec_function->threshold            = 0.0;
  vec_function->hp_threshold         = 0.0;
  vec_function->hp_maxlevel          = 0;
  vec_function->f_el_info            = square_f_el_info;
  vec_function->get_vertex_estimate  = NULL;
  vec_function->get_element_estimate = NULL;
  vec_function->function_data        = fC;
}

/* **********************************************************
*
*  the functions of the description
* ********************************************************** */

static void square_coord_of_parent(HELEMENT2D* /*helem*/,double* /*p*/,double* /*q*/)
{
  fprintf(stderr,"coord_of_parent: not implemented yet!\n");

  return;
}

static HELEMENT2D *square_neighbour(HELEMENT2D * /*helem*/, int /*np*/,int /*flag*/,
            double */*coord*/, double */*xyz*/,
            MESH_ELEMENT_FLAGS /*fl*/)
{
  fprintf(stderr,"neighbour: not implemented yet!\n");

  return NULL;
}

// this is an auxiliary function for testing whether a
// point with index ind is an inner or a boundary point.
// N is the number of elements in one line
int is_boundary(int ind, int N)
{
  // the indices of the elements of the undermost line are < N,
  // on the left side they are multiple of N,
  // on the right side they are multiple - 1,
  // on the ceiling line the index is bigger or equal N*(N-1)
  if ( (ind<N) || ((ind%N)==0) || ((ind%N)==4) || (ind>=N*(N-1)) )
  {
    return(1);
  }
  else return(0);
}

static int square_boundary(HELEMENT2D *helement, int np)
{
  // first get the number of vertices in one line
  int N = (1<<helement->mesh->max_level)+1;
  // the given index is interpreted as the starting point
  // of the edge
  int ind0 = helement->vindex[np];
  int ind1 = helement->vindex[(np+1) % 4];

  if (is_boundary(ind0,N) && is_boundary(ind1,N))
    return(-1);
  else
    return(0);
}

static int square_world_to_coord(HELEMENT2D* helement,double G_CONST* xyz,double* coord)
{
  // simply subtract the coordinates of the first vertex of the element
  // but first get the absolute width/height of the element
  double h = helement->vertex[1][0] - helement->vertex[0][0];
  coord[0] = (xyz[0] - helement->vertex[0][0])/h;
  coord[1] = (xyz[1] - helement->vertex[0][1])/h;
  return(square_check_inside(helement,coord));
}


static void square_coord_to_world(HELEMENT2D* helement,double G_CONST* coord,double* xyz)
{
  // simply add the coordinates of the first vertex of the element
  // but first get the absolute width/height of the element
  double h = helement->vertex[1][0] - helement->vertex[0][0];
  xyz[0] = coord[0]*h + helement->vertex[0][0];
  xyz[1] = coord[1]*h + helement->vertex[0][1];
  return;
}


static int square_check_inside(HELEMENT2D */*helement*/, double G_CONST*coord)
{
  if ((coord[0]>=0) && (coord[0]<=1) && (coord[1]>=0) && (coord[1]<=1))
    return(-1);
  else      // point is not inside => finde the index of the separating face
  {
    if (coord[0]<0) return(3);
    if (coord[0]>1) return(1);
    if (coord[1]<0) return(0);
    if (coord[1]>1) return(2);
  }
  cout<<"Check_inside: Paradoxically the point is neither inside nor outside...???"<<endl;
  return(-1);
}


/* **********************************************************
*
*    the functions of the HMESH
*
* ********************************************************** */


static HELEMENT2D* square_first_child(HELEMENT2D* parent,MESH_ELEMENT_FLAGS /*fl*/)
{

  int level = parent->level;

  /* still not the finest level? */
  if (level < parent->mesh->max_level)
  {
    // increase the level
    ++level;

    /* for the first child the initialized coordinates and vindex
      has to be updated */

    // calc the vertical and horizontal number-distance
    int dx = (parent->vindex[1]-parent->vindex[0])/2;
    int dy = (parent->vindex[3]-parent->vindex[0])/2;

    // get an element to fill
    HELEMENT2D* elementret = get_element(parent->mesh);
    elementret->vindex[0] = parent->vindex[0];
    elementret->vindex[1] = parent->vindex[0]+dx;
    elementret->vindex[2] = parent->vindex[0]+dx+dy;
    elementret->vindex[3] = parent->vindex[0]+dy;

    /* now the coordinates */
    double h = 1. / (static_cast<double>(1<<level));
    double x = parent->vertex[0][0];
    double y = parent->vertex[0][1];

    elementret->vertex[0][0] = x; elementret->vertex[0][1] = y; elementret->vertex[0][2] = 0.;
    elementret->vertex[1][0] = x+h; elementret->vertex[1][1] = y; elementret->vertex[1][2] = 0.;
    elementret->vertex[2][0] = x+h; elementret->vertex[2][1] = y+h; elementret->vertex[2][2] = 0.;
    elementret->vertex[3][0] = x; elementret->vertex[3][1] = y+h; elementret->vertex[3][2] = 0.;

    /* the number of the children is saved in the eindex-entry,
      running from 0 to 3, 0 means: first child */
    elementret->eindex = 0;
    elementret->parent = parent;
    elementret->level = level;
    elementret->vinh = vinhChild0;
    if (level == parent->mesh->max_level) elementret->has_children=0;

    return(reinterpret_cast<HELEMENT2D*>(elementret));
  }
  else
  {
    return NULL;
  }
}



static HELEMENT2D* square_next_child(HELEMENT2D* helement,MESH_ELEMENT_FLAGS /*fl*/)
{
  int level = helement->level;

  // level >0, there are four children, if its already the fourth one => return NULL
  if ((level>0) && (helement->eindex < 3))
  {
    /* now adjust the necessary entries of the element-structure */
    ++(helement->eindex);          // increase number of the element

    // calc the vertical and horizontal number-distance
    int dx = helement->vindex[1]-helement->vindex[0];
    int dy = helement->vindex[3]-helement->vindex[0];

    // calc the coord-distance
    double h = helement->vertex[1][0] - helement->vertex[0][0];

    // z-coord = 0
    helement->vertex[0][2] = helement->vertex[1][2] = 0.;
    helement->vertex[2][2] = helement->vertex[3][2] = 0.;

    switch(helement->eindex)
    {
      int i;
      // eine Position weiter rechts
      case 1:
          for (i=0; i<4; ++i) helement->vindex[i] += dx;
          for (i=0; i<4; ++i) helement->vertex[i][0] += h;
          helement->vinh = vinhChild1;
        break;
      // links oben
      case 2:
          for (i=0; i<4; ++i) helement->vindex[i] += (dy-dx);
          for (i=0; i<4; ++i)
          {
            helement->vertex[i][0] -= h;
            helement->vertex[i][1] += h;
          }
          helement->vinh = vinhChild2;
        break;
      // rechts oben
      case 3:
          for (i=0; i<4; ++i) helement->vindex[i] += dx;
          for (i=0; i<4; ++i) helement->vertex[i][0] += h;
          helement->vinh = vinhChild3;
        break;
      default:
        throw aol::UnimplementedCodeException ( "Unhandled switch case!", __FILE__, __LINE__ );
    };

    return(reinterpret_cast<HELEMENT2D*>(helement));

  }
  else
  {
    square_free_element(helement);
    return(NULL);
  }

}


HELEMENT2D* square_copy_element(HELEMENT2D* helement,MESH_ELEMENT_FLAGS /*fl*/)
{
  // get an element to fill
  HELEMENT2D* elementret = get_element(helement->mesh);

  // copy the entries of the element
  elementret->mesh = helement->mesh;
  elementret->eindex = helement->eindex;
  elementret->descr = helement->descr;
  elementret->present = helement->present;
  elementret->user_data = helement->user_data;
  elementret->parent = helement->parent;
  elementret->vinh = helement->vinh;
  elementret->ref_rule = helement->ref_rule;
  elementret->level = helement->level;
  elementret->has_children = helement->has_children;

  int i,j;
  for (i=0; i<4; ++i)
    elementret->vindex[i] = helement->vindex[i];

  for (i=0; i<4; ++i)
    for (j=0; j<2; ++j)
      elementret->vertex[i][j] = helement->vertex[i][j];

  return(helement);
}

static HELEMENT2D* square_select_child(ELEMENT2D* parent,
  double* parent_coord,double* child_coord,MESH_ELEMENT_FLAGS fl)
{
  if (   (parent_coord[0] >= 0.) && (parent_coord[0] <= 1.) &&
        (parent_coord[1] >= 0.) && (parent_coord[1] <= 1.) )
  {
    // get a copy of the parent element
    HELEMENT2D* child_el = square_copy_element(parent, fl);
    //return child_el;

    // adjust the entries of the parent structure
    ++child_el->level;

    // calc the vertical and horizontal number-distance
    int dx = parent->vindex[1]-parent->vindex[0];
    int dy = parent->vindex[3]-parent->vindex[0];

    // calc the coord-distance
    //double h = 1. / ((double) (1<<(child_el->level)) );
    double h = child_el->vertex[1][0] - child_el->vertex[0][0];


    // get the number of the child element (0..3)
    int child_number = 0;
    if (parent_coord[0] > 0.5) ++child_number;
    if (parent_coord[1] > 0.5) child_number += 2;

    // now, depending on the child_number, adjust the coords and
    // vindex-numbers

    // links vorne
    if (child_number==0)
    {
      child_el->vinh = vinhChild0;
      child_el->vindex[1] -= dx; child_el->vindex[2] -= dx+dy; child_el->vindex[3] -= dy;

      child_el->vertex[1][0] -= h; child_el->vertex[2][0] -= h;
      child_el->vertex[2][1] -= h; child_el->vertex[3][1] -= h;
      child_coord [0] = parent_coord [0] * 2.0;
      child_coord [1] = parent_coord [1] * 2.0;
    }
    // rechts vorne
    if (child_number==1)
    {
      child_el->vinh = vinhChild1;
      child_el->vindex[0] += dx; child_el->vindex[2] -= dy; child_el->vindex[3] += dx-dy;

      child_el->vertex[0][0] += h; child_el->vertex[2][1] -= h;
      child_el->vertex[3][0] += h; child_el->vertex[3][1] -= h;
      child_coord [0] = (parent_coord [0] - 0.5) * 2.0;
      child_coord [1] = parent_coord [1] * 2.0;
    }
    // links hinten
    if (child_number==2)
    {
      child_el->vinh = vinhChild2;
      child_el->vindex[0] += dy; child_el->vindex[1] += dy-dx; child_el->vindex[2] -= dx;

      child_el->vertex[0][1] += h; child_el->vertex[1][0] -= h;
      child_el->vertex[1][1] += h; child_el->vertex[2][0] -= h;
      child_coord [0] = parent_coord [0] * 2.0;
      child_coord [1] = (parent_coord [1] - 0.5) * 2.0;
    }
    // rechts hinten
    if (child_number==3)
    {
      child_el->vinh = vinhChild3;
      child_el->vindex[0] += dx+dy; child_el->vindex[1] += dy; child_el->vindex[3] += dx;

      child_el->vertex[0][0] += h; child_el->vertex[0][1] += h;
      child_el->vertex[1][1] += h; child_el->vertex[3][0] += h;
      child_coord [0] = (parent_coord [0] - 0.5) * 2.0;
      child_coord [1] = (parent_coord [1] - 0.5) * 2.0;
    }

    return child_el;

  }
  else return NULL;    // not a child element
}

static HELEMENT2D* square_first_macro(GENMESH2D* hmesh,MESH_ELEMENT_FLAGS /*fl*/)
{
  // get an element to fill
  HELEMENT2D* elementret = get_element(hmesh);

  /* for the first macro the initialized coordinates and vindex
    has to be updated, eindex is already correct */
  /* first calc the indices of the edges: */
  int d = hmesh->max_level;
  int pd = 1<<d;          // 2^d
  int p2d = 1<<(2*d);      // 2^(2d)

  elementret->vindex[0] = 0; elementret->vindex[1] = pd;
  elementret->vindex[3] = pd+p2d; elementret->vindex[2] = p2d+2*pd;

  /* now the coordinates */
  elementret->vertex[0][0] = 0.; elementret->vertex[0][1] = 0.; elementret->vertex[0][2] = 0.;
  elementret->vertex[1][0] = 1.; elementret->vertex[1][1] = 0.; elementret->vertex[1][2] = 0.;
  elementret->vertex[2][0] = 1.; elementret->vertex[2][1] = 1.; elementret->vertex[2][2] = 0.;
  elementret->vertex[3][0] = 0.; elementret->vertex[3][1] = 1.; elementret->vertex[3][2] = 0.;
  elementret->level = 0;
  elementret->vinh = NULL;


  return(reinterpret_cast<HELEMENT2D*>(elementret));
}


static HELEMENT2D* square_next_macro(HELEMENT2D* helement,MESH_ELEMENT_FLAGS /*fl*/)
{
  /* There is only one macro-element => return NULL */
  square_free_element(helement);
  return NULL;
}


static void square_free_element(HELEMENT2D *el)
{
  // put the element into the free-element-list
  elements_free.push_back(el);
  elements_used.pop_back();
  //cout<<"----------- Element freigeworden, Size of free-List: "<<elements_free.size()<<endl;

  return;
}


static int square_set_time(GENMESH2D* mesh,double t)
{
  (reinterpret_cast<MESH_INFO*>(mesh->user_data))->t = t;

  F_HDATA2D* fun;
  for (fun = reinterpret_cast<F_HDATA2D*> (mesh->f_data); fun->last; fun = reinterpret_cast<F_HDATA2D*>(fun->last)) ; // Goto first function

  // Interpolate all functions
  for (; fun; fun = reinterpret_cast<F_HDATA2D*> (fun->next)) {

    if (fun->dimension_of_value == 1) {

      int i = static_cast<int> (floor (t));
      int j = i+1;

      double jfac = t - i;
      double ifac = 1 - jfac;

      functionComponents* fC = reinterpret_cast<functionComponents*>(fun->function_data);
      int n = fC->data_seq->size ();

      if (i < 0) { i = j = 0; ifac = 1; jfac = 0; }
      if (j >= n) { i = j = n-1; ifac = 1; jfac = 0; }

      qc::ScalarArray<double, qc::QC_2D>& idata = *((*fC->data_seq) [i]);
      qc::ScalarArray<double, qc::QC_2D>& jdata = *((*fC->data_seq) [j]);
      qc::ScalarArray<double, qc::QC_2D>& data = *(fC->data);

      for (int a = 0; a < data.size (); ++a)
  data [a] = ifac * idata [a] + jfac * jdata [a];

      if (!fC->estimator) continue;

      qc::Estimator2d<double>& iestimator = *((*fC->estimator_seq) [i]);
      qc::Estimator2d<double>& jestimator = *((*fC->estimator_seq) [j]);
      qc::Estimator2d<double>& estimator = *(fC->estimator);

      for (int a = 0; a < estimator.size (); ++a)
  estimator [a] = ifac * iestimator [a] + jfac * jestimator [a];
    }
    else {
      // Currently no time-dependency for vector valued functions
    }
  }

  return TRUE;
}

/*

why unused??

static int square_get_time(GENMESH2D* mesh, double* t, double* t1, double* t2)
{
  *t = ((MESH_INFO*)mesh->user_data)->t;
  *t1 = floor (*t); *t2 = ceil (*t);

  functionComponents* fC = (functionComponents*)(((F_HDATA2D*)mesh->f_data)->function_data);
  int n = fC->data_seq->size ();
  if (*t1 < 0) *t1 = *t2 = 0;
  if (*t2 >= n) *t1 = *t2 = n-1;

  return TRUE;
}
*/

GENMESH2D* genmesh2d_get_times (double* t1, double* t2)
{
  GENMESH2D* self = reinterpret_cast<GENMESH2D*> (START_METHOD( G_INSTANCE ));
  ASSURE(self,"",END_METHOD(NULL));

  *t1 = 0;

  functionComponents* fC = reinterpret_cast<functionComponents*>(reinterpret_cast<F_HDATA2D*>(self->f_data)->function_data);
  *t2 = static_cast<int>(fC->data_seq->size () - 1);

  END_METHOD (self);
}

/* ***************************************************************
*
*  Now the functions for managing the used elements in a vector
*
* *************************************************************** */

HELEMENT2D* get_element(GENMESH2D* gmesh)
{
  // the necessary pointers
  HELEMENT2D* elementret;

  // no free element available => allocate a new one
  // and put it on the stack
  if (elements_free.size() == 0)
  {
    elementret = reinterpret_cast<HELEMENT2D*>(mem_alloc(sizeof(HELEMENT2D)));
    elementret->vindex = reinterpret_cast<int*>(mem_alloc(4*sizeof(int)));
    elementret->vertex = reinterpret_cast<double**>(mem_alloc(4*sizeof(double*)));
    for (int i=0; i<4; ++i) {
      elementret->vertex[i] = reinterpret_cast<double*>(mem_alloc(3*sizeof(double)));
        for (int j=0; j<3; ++j) elementret->vertex[i][j] = 0;
    }
    //cout<<"------------ Neues Element deklariert, ";
  }
  else
  {
    // ok, there is a free element available, so take it
    // out of the vector
    elementret = elements_free.back();
    elements_free.pop_back();
    //cout<<"------------ Freies Element verfuegbar, ";
  }
  elements_used.push_back(elementret);
  // initialize and return it
  //  cout<<" S.o. used: "<<elements_used.size()<<
  //      ", S.o. free: "<<elements_free.size()<<endl;
  fillHElement2D(elementret, gmesh);
  return(reinterpret_cast<HELEMENT2D*>(elementret));
}


/* ***************************************************************
*
*  Now the functions itself
*
* *************************************************************** */

// first: the easy case of a scalar function
static void scalar_function(HELEMENT2D* helement,
          int        index,
          double     coord[],
          double     val[],
          void*      function_data)
{

  //qcScalarArray<double, qc::QC_2D>* meshData = (qcScalarArray<double, qc::QC_2D>*)function_data;
  qc::ScalarArray<double, qc::QC_2D>* meshData = reinterpret_cast<functionComponents*>(function_data)->data;

  // zi,si are the indices of the array2d
  int i,ind;
  double f[4];


  if(coord)
  {
    double komplx = 1.-coord[0];
    double t1,t2;

    for (i=0; i<4; ++i)
    {
      ind = helement->vindex[i];
      f[i] = (*meshData)[ind];
    }

    t1 = coord[0]*f[2] + komplx*f[3];
    t2 = coord[0]*f[1] + komplx*f[0];

    val[0] = coord[1]*t1 + (1.-coord[1])*t2;
  }
  else
  {
    ind = helement->vindex[index];
    val[0] = (*meshData)[ind];
  }

  return;
}


// second: the case of a n-dimensional vector-valued function
static void vector_function(HELEMENT2D* helement,
          int        index,
          double     coord[],
          double     val[],
          void*      function_data)
{

  aol::MultiVector<double>* meshDataVec = reinterpret_cast<functionComponentsMV*>(function_data)->data;

  int dim = meshDataVec->numComponents ();    // correct (!) dimension of the multivector
  int i,j;
  int ind[4];
  double f[4], t1,t2,komplx,komply;

  // get the indices of the four vertices
  for (i=0; i<4; ++i)
    ind[i] = helement->vindex[i];

  if(coord)
  {
    komplx = 1.-coord[0];
    komply = 1.-coord[1];

    // bilinear interpolating of each coordinate
    for (i=0; i<dim; ++i)
    {
      for (j=0; j<4; ++j)
        f[j] = (*meshDataVec)[i][ind[j]];

      t1 = coord[0]*f[2] + komplx*f[3];
      t2 = coord[0]*f[1] + komplx*f[0];

      val[i] = coord[1]*t1 + komply*t2;
    }
  }
  else
  {
    ind[0] = helement->vindex[index];
    for (i=0; i<dim; ++i)
      val[i] = (*meshDataVec) [i][ind[0]];
  }

  return;
}



static void square_f_el_info(ELEMENT2D* /*el*/,F_EL_INFO2D* hel_info,void* /*data*/)
{
  hel_info->polynomial_degree = 1;

  return;
}





/* ***************************************************************
*
*  Now the ESTIMATORS of the function-data
*
* *************************************************************** */

// the element-estimator calcs the maximum of the error values
// of the five inner vertices (one level higher)
double calcElementError (HELEMENT2D *el, void *function_data)
{
  int level = el->level;
  int max_level = el->mesh->max_level;

  // finest level => no estimator necessary
  if (level == max_level)
  {
    return 0;
  }
  else
  {
    int mx,my,ind,halfstep;
    int zwei_d = (1<<max_level) + 1;

    // We now have to get the scalar2d-indices of the left upper vertex
    // of the element. Therefore calc the x- and y-indices-distance.
    ind = el->vindex[0];
    halfstep = (el->vindex[1] - ind) >> 1;
    mx = ind / zwei_d;
    my = ind % zwei_d;

    qc::Estimator2d<double>* Est2d = reinterpret_cast<qc::Estimator2d<double>*>(reinterpret_cast<functionComponents*>(function_data)->estimator);

    // get the maximum of the five inner points (one level finer) of the element
    double e = Est2d->maxChildNodeValue(mx,my,halfstep);

    return(e);
  }
  return 0;

}


// the vertex-estimator reads the error out of the saturated error structure
void calcVertexError (HELEMENT2D *el, double *error, void *function_data)
{
  qc::Estimator2d<double>* Est2d = reinterpret_cast<qc::Estimator2d<double>*>(reinterpret_cast<functionComponents*>(function_data)->estimator);

  for (int i=0; i<4; ++i)
    error[i] = (*Est2d)[ el->vindex[i] ];

}






/* *********************************************************************
*
*    Taetaetae: the convert-functions
*    expects a QuocMesh-Grid and provides a GenMesh
*    the other argument is either a ScalarArray<QC_2D> or a MultiVector
*
* ********************************************************************* */

// First data-independent entries are filled, args are the mesh and its depth
void fillGMesh(GENMESH2D* gmesh, int d)
{
  gmesh->first_child            = square_first_child;
  gmesh->next_child             = square_next_child;
  gmesh->select_child             = square_select_child;
  gmesh->first_macro            = square_first_macro;
  gmesh->next_macro             = square_next_macro;

  gmesh->free_element           = square_free_element;
  gmesh->max_level              = d;
  gmesh->level_of_interest      = d;
  gmesh->max_vindex             = ((1<<d)+1)*((1<<d)+1);
  gmesh->max_eindex             = 1<<(2*d);

  gmesh->dimension_of_world     = 2;
  gmesh->max_dimension_of_coord = 2;
  gmesh->max_number_of_vertices = 4;
  gmesh->set_time               = square_set_time;
  gmesh->get_geometry_vertex_estimate  = NULL;
  gmesh->get_geometry_element_estimate = NULL;

  MESH_INFO* mesh_info = reinterpret_cast<MESH_INFO*> (mem_alloc(sizeof(MESH_INFO)));
  mesh_info->t = 0;
  gmesh->user_data = mesh_info;
}



/* ********************************************************************
    function tests, whether a scalar-array is quadratic or not,
    and provides (if its quadratic) the depth of it!
  ********************************************************************* */
int depthOfScalarArray(qc::ScalarArray<double, qc::QC_2D>* data)
{
  int N = data->getNumX ();
  int d = qc::logBaseTwo (N);
  if (data->getNumY () != N)
    throw aol::Exception ("Image not square", __FILE__, __LINE__);
  if ((1 << d) + 1 != N)
    throw aol::Exception ("dimension no power of two", __FILE__, __LINE__);

  // everything is ok (be lucky)
  return(d);
}

int depthOfMultiVector(aol::MultiVector<double>* data)
{
  int N = (data->getTotalSize()) / (data->numComponents());
  int d = qc::logBaseTwo(N) / 2;
  int erg = (1 << d) + 1;
  if ( erg*erg != N)
    throw aol::Exception ("dimension no power of two", __FILE__, __LINE__);

  // everything is ok (be lucky)
  return(d);
}



/* *********************************************************************
* the convert-routine for scalar data, without own estimator
* the estimator is generated automatically here
* TODO: The estimator must be generated
* ********************************************************************* */

GENMESH2D* quocmesh_convert_to_gmesh2d(qc::ScalarArray<double, qc::QC_2D>* meshData,
                                       const char* dataName)
{
  GENMESH2D* gmesh;

  // initialize the gmesh-structure
  gmesh   = reinterpret_cast<GENMESH2D*>(GRAPE(GenMesh2d,"new-instance")("Adaptive QuocMesh"));
  // fill the description of one element
  fillDescription();
  // Fill the gmesh-structure
  int d = depthOfScalarArray(meshData);
  fillGMesh(gmesh,d);

  // allocate memory for the data structure
  F_HDATA2D* data_function = reinterpret_cast<F_HDATA2D*>(mem_alloc(sizeof(F_HDATA2D)));
  gmesh->f_data = reinterpret_cast<GENMESH_FDATA*>(data_function);

  // -------------------------------------------------------------------------------
  // --------------------------- TODO: Define an error estimator and saturate it
  // -------------------------------------------------------------------------------


  // fill the data structure
  functionComponents* fC = reinterpret_cast<functionComponents*> (mem_alloc(sizeof(functionComponents)));
  fC->data = new qc::ScalarArray<double, qc::QC_2D> (*meshData);
  fC->estimator = NULL;
  fC->get_element_estimate = NULL;
  fC->get_vertex_estimate = NULL;

  fC->data_seq = new std::vector<qc::ScalarArray<double, qc::QC_2D>*>;
  fC->data_seq->push_back (meshData);
  fC->estimator_seq = new std::vector<qc::Estimator2d<double>*>;
  fC->estimator_seq->push_back (NULL);

  fillScalarDataFunction(data_function, fC);
  data_function->name = const_cast<char*>(dataName);

  return(gmesh);
}

/* *********************************************************************
* the convert-routine for scalar data, with own estimator, here you
* can use the estimator2d structure with a makeSaturatedErrorArray-call
* before, or you deliver NULL and no estimator is used
* ********************************************************************* */

GENMESH2D* quocmesh_convert_to_gmesh2d(qc::ScalarArray<double, qc::QC_2D>* meshData,
                                       qc::Estimator2d<double>* errorArray, const char* dataName)
{

  cerr<<"In convert ..."<<endl;

  GENMESH2D* gmesh;

  // initialize the gmesh-structure
  gmesh   = reinterpret_cast<GENMESH2D*>(GRAPE(GenMesh2d,"new-instance")("Adaptive QuocMesh"));
  // fill the description of one element
  fillDescription();
  // Fill the gmesh-structure
  int d = depthOfScalarArray(meshData);
  int dE = depthOfScalarArray(errorArray);
  if (d != dE)
    throw aol::Exception ("ScalarArray and Estimator-Array have NOT the same depth!", __FILE__, __LINE__);
  fillGMesh(gmesh,d);

  // allocate memory for the data structure
  F_HDATA2D* data_function = reinterpret_cast<F_HDATA2D*>(mem_alloc(sizeof(F_HDATA2D)));
  gmesh->f_data = reinterpret_cast<GENMESH_FDATA*>(data_function);

  // fill the data structure
  /** REMARK: memory leak: memory of functionComponents is not freed
              (but only a few bytes)
      */
  functionComponents* fC = reinterpret_cast<functionComponents*>(mem_alloc(sizeof(functionComponents)));
  fC->data = new qc::ScalarArray<double, qc::QC_2D> (*meshData);
  fC->estimator = new qc::Estimator2d<double> (*errorArray);
  fC->get_element_estimate = &calcElementError;
  fC->get_vertex_estimate = &calcVertexError;

  fC->data_seq = new std::vector<qc::ScalarArray<double, qc::QC_2D>*>;
  fC->data_seq->push_back (meshData);
  fC->estimator_seq = new std::vector<qc::Estimator2d<double>*>;
  fC->estimator_seq->push_back (errorArray);

  fillScalarDataFunction(data_function, fC);
  data_function->name = const_cast<char*>(dataName);

  return(gmesh);
}

/* *********************************************************************
* the convert-routine for vector valued data
* ********************************************************************* */

GENMESH2D* quocmesh_convert_to_gmesh2d(aol::MultiVector<double>* meshData,
                                       const char* dataName)
{
  GENMESH2D* gmesh;

  // initialize the gmesh-structure
  gmesh   = reinterpret_cast<GENMESH2D*>(GRAPE(GenMesh2d,"new-instance")("Adaptive QuocMesh"));
  // fill the description of one element
  fillDescription();

  // Fill the gmesh-structure
  int d = depthOfMultiVector(meshData);
  fillGMesh(gmesh,d);

  // allocate memory for the data structure
  F_HDATA2D* data_function = reinterpret_cast<F_HDATA2D*>(mem_alloc(sizeof(F_HDATA2D)));
  gmesh->f_data = reinterpret_cast<GENMESH_FDATA*>(data_function);

  // fill the data structure
  /** REMARK: memory leak: memory of functionComponents is not freed
              (but only a few bytes)
      */
  functionComponentsMV* fC = reinterpret_cast<functionComponentsMV*>(mem_alloc(sizeof(functionComponentsMV)));
  fC->data = new aol::MultiVector<double> (*meshData);

  fC->data_seq = new std::vector<aol::MultiVector<double>*>;
  fC->data_seq->push_back (meshData);

  // fill the data structure
  fillVectorDataFunction(data_function, fC);
  data_function->name = const_cast<char*>(dataName);

  return(gmesh);
}

GENMESH2D *reload_data( ) {
  GENMESH2D *self;
  self = reinterpret_cast<GENMESH2D*> (START_METHOD( G_INSTANCE ));
  ASSURE(self,"",END_METHOD(NULL));

  cerr << "reloading data... " << reinterpret_cast<MESH_INFO*>(self->user_data)->filename;

  reinterpret_cast<MESH_INFO*>(self->user_data)->data->load(reinterpret_cast<MESH_INFO*>(self->user_data)->filename );

  END_METHOD( self );
}

void add_reload_button( GENMESH2D *mesh, qc::ScalarArray<double, qc::QC_2D> &data, const char *filename ) {
  reinterpret_cast<MESH_INFO*>(mesh->user_data)->filename = const_cast<char*>(filename);
  reinterpret_cast<MESH_INFO*>(mesh->user_data)->data = &data;

  GRAPE(GenMesh2d,"add-method")("reload-data-send", reload_data );

  MANAGER*  mgr;
  mgr = reinterpret_cast<MANAGER*>(GRAPE(Manager,"get-stdmgr")());
  GRAPE(mgr,"add-inter")( reinterpret_cast<BUTTON*>(new_item( Button,
                    I_Name, "reload data",
                    I_Size, 3.5, 1.0,
                    I_Instance, mesh,
                    I_Method, "reload-data-send",
                    I_End )));

}

/* ---------------------------------------------------------------------------------------
* ***** ADDING NEW FUNCTIONS *******************************
* --------------------------------------------------------------------------------------- */

/* *********************************************************************
* add scalar data
* ********************************************************************* */

void addScalarData(GENMESH2D* gmesh, qc::ScalarArray<double, qc::QC_2D>* meshData, const char* dataName)
{
    // allocate memory for the data structures
    F_HDATA2D* data_function = reinterpret_cast<F_HDATA2D*>(mem_alloc(sizeof(F_HDATA2D)));

    // fill the data structure
    functionComponents* fC = reinterpret_cast<functionComponents*>(mem_alloc(sizeof(functionComponents)));
    fC->estimator = NULL;
    fC->get_element_estimate = NULL;
    fC->get_vertex_estimate = NULL;

    fC->data_seq = new std::vector<qc::ScalarArray<double, qc::QC_2D>*>;
    fC->data_seq->push_back (meshData);
    fC->estimator_seq = new std::vector<qc::Estimator2d<double>*>;
    fC->estimator_seq->push_back (NULL);

    fC->data = new qc::ScalarArray<double, qc::QC_2D> (*meshData);
    fillScalarDataFunction(data_function, fC);
    data_function->name = const_cast<char*>(dataName);

    // now add the function via GRAPE
    GRAPE(gmesh,"add-function")(data_function);
}


/* *********************************************************************
* add timestep with scalar data
* ********************************************************************* */


void addTimestep (GENMESH2D* gmesh, qc::ScalarArray<double, qc::QC_2D>* meshData, const char* dataName, qc::Estimator2d<double>* errorArray)
{
  F_HDATA2D* fun;
  for (fun = reinterpret_cast<F_HDATA2D*>(gmesh->f_data); fun->last; fun = reinterpret_cast<F_HDATA2D*>(fun->last)) ; // Goto first function

  // Find correct function
  for (; fun && strcmp (fun->name, dataName); fun = reinterpret_cast<F_HDATA2D*>(fun->next)) ;

  if (strcmp (fun->name, dataName))
    throw aol::Exception ("Function not found", __FILE__, __LINE__);

  functionComponents* fC = reinterpret_cast<functionComponents*>(fun->function_data);

  fC->data_seq->push_back (meshData);
  fC->estimator_seq->push_back (errorArray);
}

/* *********************************************************************
* add timestep with vector data
* ********************************************************************* */

void addTimestep (GENMESH2D* gmesh, aol::MultiVector<double>* meshData, const char* dataName)
{
  F_HDATA2D* fun;
  for (fun = reinterpret_cast<F_HDATA2D*>(gmesh->f_data); fun->last; fun = reinterpret_cast<F_HDATA2D*>(fun->last)) ; // Goto first function

  // Find correct function
  for (; fun && strcmp (fun->name, dataName); fun = reinterpret_cast<F_HDATA2D*>(fun->next)) ;

  if (strcmp (fun->name, dataName))
    throw aol::Exception ("Function not found", __FILE__, __LINE__);

  functionComponentsMV* fC = reinterpret_cast<functionComponentsMV*>(fun->function_data);

  fC->data_seq->push_back (meshData);
}


/* *********************************************************************
* add scalar data with an qcEstimator2d
* ********************************************************************* */

void addScalarDataWithEstimator(GENMESH2D* gmesh, qc::ScalarArray<double, qc::QC_2D>* meshData, const char* dataName,
                                qc::Estimator2d<double>* errorArray)
{
    // allocate memory for the data structures
    F_HDATA2D* data_function = reinterpret_cast<F_HDATA2D*>(mem_alloc(sizeof(F_HDATA2D)));

    // fill the data structure
    functionComponents* fC = reinterpret_cast<functionComponents*>(mem_alloc(sizeof(functionComponents)));
    fC->estimator = new qc::Estimator2d<double> (*errorArray);
    fC->get_element_estimate = &calcElementError;
    fC->get_vertex_estimate = &calcVertexError;

    fC->data_seq = new std::vector<qc::ScalarArray<double, qc::QC_2D>*>;
    fC->data_seq->push_back (meshData);
    fC->estimator_seq = new std::vector<qc::Estimator2d<double>*>;
    fC->estimator_seq->push_back (errorArray);

    fC->data = new qc::ScalarArray<double, qc::QC_2D> (*meshData);
    fillScalarDataFunction(data_function, fC);
    data_function->name = const_cast<char*>(dataName);

    // now add the function via GRAPE
    GRAPE(gmesh,"add-function")(data_function);
}




/* *********************************************************************
* add vector data
* ********************************************************************* */

void addVectorData(GENMESH2D* gmesh, aol::MultiVector<double>* meshData, const char* dataName)
{
    // allocate memory for the data structures
    F_HDATA2D* data_function = reinterpret_cast<F_HDATA2D*>(mem_alloc(sizeof(F_HDATA2D)));

    // fill the data structure
    /** REMARK: memory leak: memory of functionComponents is not freed
  (but only a few bytes)
    */
    functionComponentsMV* fC = reinterpret_cast<functionComponentsMV*>(mem_alloc(sizeof(functionComponentsMV)));
    fC->data = new aol::MultiVector<double> (*meshData);

    fC->data_seq = new std::vector<aol::MultiVector<double>*>;
    fC->data_seq->push_back (meshData);

    // fill the data structure
    fillVectorDataFunction(data_function, fC);
    data_function->name = const_cast<char*>(dataName);

    // now add the function via GRAPE
    GRAPE(gmesh,"add-function")(data_function);
}

void addMethodsAndProjects2d ()
{
  g_project_add(const_cast<char*>("uif-m2"));

  GRAPE (GenMesh2d, "delete-method") ("softcopy");
  GRAPE (GenMesh2d, "add-method") ("softcopy", genmesh2d_my_softcopy);
  GRAPE (GenMesh2d, "add-method") ("get-times", genmesh2d_get_times);
}

void initStartGrape(GENMESH2D* gmesh, const char* dataName)
{
  addMethodsAndProjects2d();

  MANAGER*  mgr;
  mgr = reinterpret_cast<MANAGER*>(GRAPE(Manager,"get-stdmgr")());
  GRAPE(gmesh,"select-function")("default",dataName);

  GRAPE(mgr,"handle")(gmesh);
}

void initStartGrape(GENMESH2D* gmesh)
{
  addMethodsAndProjects2d();

  MANAGER*  mgr;
  mgr = reinterpret_cast<MANAGER*>(GRAPE(Manager,"get-stdmgr")());

  GRAPE(mgr,"handle")(gmesh);
}

void initStartGrapeTime(GENMESH2D* gmesh)
{
  addMethodsAndProjects2d();

  MANAGER*  mgr;
  mgr = reinterpret_cast<MANAGER*>(GRAPE(Manager,"get-stdmgr")());

  TIMESCENE* tsc = reinterpret_cast<TIMESCENE*> (GRAPE (TimeScene, "new-instance") ("sequence"));
  tsc->dynamic = reinterpret_cast<TREEOBJECT*> (gmesh);
  tsc->object = reinterpret_cast<TREEOBJECT*> (GRAPE (gmesh, "softcopy") (NULL));

  GRAPE(mgr,"handle")(tsc);
}

#endif
