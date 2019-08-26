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

#include "grapeInterface3d.h"

#ifdef USE_EXTERNAL_GRAPE

#include "float.h"

/* the next two lines represent a really wonderful trick:
  these lines already exists in Quoc.cpp, but if they exist
  only there, this program produces an error:
  `undefined reference to `aol::ZOTrait<short>::zero'
  but with the lines here, everything works, wonderful! */
/*const short aol::ZOTrait<short>::one = 1;
const short aol::ZOTrait<short>::zero = 0;*/

// Avoid name collisions between 2d and 3d interface
namespace {

  // Here we create a forward instantiation of all the functions involved
  void _n_e_v_e_r_C_a_l_l_M_e_() {
    quocmesh_convert_to_gmesh3d(static_cast<qc::ScalarArray<float, qc::QC_3D>*>(NULL),
        static_cast<char*>(NULL));
    quocmesh_convert_to_gmesh3d(static_cast<qc::ScalarArray<double, qc::QC_3D>*>(NULL),
        static_cast<char*>(NULL));

    quocmesh_convert_to_gmesh3d(static_cast<qc::ScalarArray<float, qc::QC_3D>*>(NULL),
                static_cast<qc::Estimator3d<float>*>(NULL),
                static_cast<char*>(NULL));
    quocmesh_convert_to_gmesh3d(static_cast<qc::ScalarArray<double, qc::QC_3D>*>(NULL),
                static_cast<qc::Estimator3d<double>*>(NULL),
                static_cast<char*>(NULL));

    quocmesh_convert_to_gmesh3d(static_cast<aol::MultiVector<float>*>(NULL),
                static_cast<char*>(NULL));
    quocmesh_convert_to_gmesh3d(static_cast<aol::MultiVector<double>*>(NULL),
                static_cast<char*>(NULL));

    addScalarData(static_cast<GENMESH3D*>(NULL),
      static_cast<qc::ScalarArray<float, qc::QC_3D>*>(NULL),
      static_cast<char*>(NULL));
    addScalarData(static_cast<GENMESH3D*>(NULL),
      static_cast<qc::ScalarArray<double, qc::QC_3D>*>(NULL),
      static_cast<char*>(NULL));

    addScalarDataWithEstimator(static_cast<GENMESH3D*>(NULL),
             static_cast<qc::ScalarArray<float, qc::QC_3D>*>(NULL),
             static_cast<char*>(NULL),
             static_cast<qc::Estimator3d<float>*>(NULL));
    addScalarDataWithEstimator(static_cast<GENMESH3D*>(NULL),
             static_cast<qc::ScalarArray<double, qc::QC_3D>*>(NULL),
             static_cast<char*>(NULL),
             static_cast<qc::Estimator3d<double>*>(NULL));

    addVectorData(static_cast<GENMESH3D*>(NULL),
                  static_cast<aol::MultiVector<float>*>(NULL),
                  static_cast<char*>(NULL));
    addVectorData(static_cast<GENMESH3D*>(NULL),
                  static_cast<aol::MultiVector<double>*>(NULL),
                  static_cast<char*>(NULL));
  }


#define class GRAPE_class
#define CLASS GRAPE_CLASS
#define GRAPE_METHOD_RETURN ___GRAPE_METHOD_RETURN
#define GRAPE_METHOD_CODE ___GRAPE_METHOD_CODE
extern "C" {
  typedef struct {
    GENMESH3D_STRUCT;
    MATRIX44 sort_matrix;
    int sort_order[8];
  } QUOCMESH3D;
}
#undef GRAPE_METHOD_RETURN
#undef GRAPE_METHOD_CODE
#undef class
#undef CLASS

GRAPE_CLASS *QuocMesh3d = NULL;

// offsets for the vertices of the child elements in one cube
// needed for the sorted first_ and next_child-methods
const int childOffsetsX[8] = { 0,1,0,1,0,1,0,1 };
const int childOffsetsY[8] = { 0,0,1,1,0,0,1,1 };
const int childOffsetsZ[8] = { 0,0,0,0,1,1,1,1 };


/* structure for the function: it contains
  the original function data and an estimator */
template <typename REAL>
struct functionComponents
{
  qc::ScalarArray<REAL, qc::QC_3D>* data;
  qc::Estimator3d<REAL>* estimator;
  double (*get_element_estimate) (HELEMENT3D *, void *function_data);
  void (*get_vertex_estimate) (HELEMENT3D *, double *, void *function_data);
  std::vector<qc::ScalarArray<REAL, qc::QC_3D>*>* data_seq;
  std::vector<qc::Estimator3d<REAL>*>* estimator_seq;
};

struct MESH_INFO {

  double t;
};

/* Extend copy methods to also copy user data and functions */
GENMESH3D *genmesh3d_my_softcopy (GENMESH3D *copy)
{
  GENMESH3D *self;

  self = reinterpret_cast<GENMESH3D *>(START_METHOD (G_INSTANCE));
  ASSURE (self, "", END_METHOD (NULL));

  if ((copy = reinterpret_cast<GENMESH3D *>(GRAPE(self, "^softcopy")(copy)))) {

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
    copy->threshold = self->threshold;

    /* function pointers should be set in new-instance */
    copy->get_geometry_vertex_estimate = self->get_geometry_vertex_estimate;
    copy->get_geometry_element_estimate = self->get_geometry_element_estimate;

    copy->get_lens_element_estimate = self->get_lens_element_estimate;

    /* extension */
    copy->user_data = self->user_data;
    copy->f_data = self->f_data;
  }

  END_METHOD (copy);
}

/* two vectors for managing the elements3d */
vector<HELEMENT3D*> elements_used;
vector<HELEMENT3D*> elements_free;

/* vertex inheritance structures for the cubes */
static double weightOne[1] = {1.0};
static double weightTwo[2] = {0.5,0.5};
static double weightFour[4] = {0.25,0.25,0.25,0.25};
static double weightEight[8] = {0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125};

/* the ground */
int indices0[1] = {0}; VINHERIT vinh0 = {1,indices0,weightOne};
int indices1[2] = {1,0}; VINHERIT vinh1 = {2,indices1,weightTwo};
int indices2[1] = {1}; VINHERIT vinh2 = {1,indices2,weightOne};
int indices3[2] = {0,3}; VINHERIT vinh3 = {2,indices3,weightTwo};
int indices4[4] = {0,1,2,3}; VINHERIT vinh4 = {4,indices4,weightFour};
int indices5[2] = {2,1}; VINHERIT vinh5 = {2,indices5,weightTwo};
int indices6[1] = {3}; VINHERIT vinh6 = {1,indices6,weightOne};
int indices7[2] = {3,2}; VINHERIT vinh7 = {2,indices7,weightTwo};
int indices8[1] = {2}; VINHERIT vinh8 = {1,indices8,weightOne};
/* the plane in the middle */
int indices9[2] = {0,4}; VINHERIT vinh9 = {2,indices9,weightTwo};
int indices10[4] = {0,1,5,4}; VINHERIT vinh10 = {4,indices10,weightFour};
int indices11[2] = {1,5}; VINHERIT vinh11 = {2,indices11,weightTwo};
int indices12[4] = {0,3,7,4}; VINHERIT vinh12 = {4,indices12,weightFour};
int indices13[8] = {0,1,2,3,4,5,6,7}; VINHERIT vinh13 = {8,indices13,weightEight};
int indices14[4] = {1,2,6,5}; VINHERIT vinh14 = {4,indices14,weightFour};
int indices15[2] = {3,7}; VINHERIT vinh15 = {2,indices15,weightTwo};
int indices16[4] = {2,3,7,6}; VINHERIT vinh16 = {4,indices16,weightFour};
int indices17[2] = {2,6}; VINHERIT vinh17 = {2,indices17,weightTwo};
/* the ceiling plane */
int indices18[1] = {4}; VINHERIT vinh18 = {1,indices18,weightOne};
int indices19[2] = {4,5}; VINHERIT vinh19 = {2,indices19,weightTwo};
int indices20[1] = {5}; VINHERIT vinh20 = {1,indices20,weightOne};
int indices21[2] = {4,7}; VINHERIT vinh21 = {2,indices21,weightTwo};
int indices22[4] = {4,5,6,7}; VINHERIT vinh22 = {4,indices22,weightFour};
int indices23[2] = {5,6}; VINHERIT vinh23 = {2,indices23,weightTwo};
int indices24[1] = {7}; VINHERIT vinh24 = {1,indices24,weightOne};
int indices25[2] = {7,6}; VINHERIT vinh25 = {2,indices25,weightTwo};
int indices26[1] = {6}; VINHERIT vinh26 = {1,indices26,weightOne};

G_CONST VINHERIT vinhChild0[8] = {vinh0,vinh1,vinh4,vinh3,vinh9,vinh10,vinh13,vinh12};
G_CONST VINHERIT vinhChild1[8] = {vinh1,vinh2,vinh5,vinh4,vinh10,vinh11,vinh14,vinh13};
G_CONST VINHERIT vinhChild2[8] = {vinh3,vinh4,vinh7,vinh6,vinh12,vinh13,vinh16,vinh15};
G_CONST VINHERIT vinhChild3[8] = {vinh4,vinh5,vinh8,vinh7,vinh13,vinh14,vinh17,vinh16};
G_CONST VINHERIT vinhChild4[8] = {vinh9,vinh10,vinh13,vinh12,vinh18,vinh19,vinh22,vinh21};
G_CONST VINHERIT vinhChild5[8] = {vinh10,vinh11,vinh14,vinh13,vinh19,vinh20,vinh23,vinh22};
G_CONST VINHERIT vinhChild6[8] = {vinh12,vinh13,vinh16,vinh15,vinh21,vinh22,vinh25,vinh24};
G_CONST VINHERIT vinhChild7[8] = {vinh13,vinh14,vinh17,vinh16,vinh22,vinh23,vinh26,vinh25};

/* local coordinates for the coordsystem of the cube */
static G_CONST double s_c0[3] = {0.,0.,0.};
static G_CONST double s_c1[3] = {1.,0.,0.};
static G_CONST double s_c2[3] = {1.,0.,1.};
static G_CONST double s_c3[3] = {0.,0.,1.};
static G_CONST double s_c4[3] = {0.,1.,0.};
static G_CONST double s_c5[3] = {1.,1.,0.};
static G_CONST double s_c6[3] = {1.,1.,1.};
static G_CONST double s_c7[3] = {0.,1.,1.};

/* the coord-system */
static G_CONST double *cube_coords[8] = {s_c0,s_c1,s_c2,s_c3,s_c4,s_c5,s_c6,s_c7};
/* each polygon has 4 vertices */
static G_CONST int poly_length[6] = {4,4,4,4,4,4};
/* defining the vertices of each polygon */
static G_CONST int poly_vert1[4] = {0,3,2,1};
static G_CONST int poly_vert2[4] = {4,5,6,7};
static G_CONST int poly_vert3[4] = {0,1,5,4};
static G_CONST int poly_vert4[4] = {1,2,6,5};
static G_CONST int poly_vert5[4] = {2,3,7,6};
static G_CONST int poly_vert6[4] = {0,4,7,3};
static G_CONST int *poly_vertex[6] =
  {poly_vert1,poly_vert2,poly_vert3,poly_vert4,poly_vert5,poly_vert6};
/* defining the neighbour-hoods */
static G_CONST int pn1[4] = {5,4,3,2};
static G_CONST int pn2[4] = {2,3,4,5};
static G_CONST int pn3[4] = {0,3,1,5};
static G_CONST int pn4[4] = {0,4,1,2};
static G_CONST int pn5[4] = {0,5,1,3};
static G_CONST int pn6[4] = {2,1,4,0};
static G_CONST int *poly_neighbour[6] = {pn1,pn2,pn3,pn4,pn5,pn6};

} // nameless namespace

/* The functions for the description of the elements */
static void        cube_free_element   (HELEMENT3D*);
static int         cube_boundary       (HELEMENT3D*,int);
static int         cube_world_to_coord (HELEMENT3D*,double G_CONST*,double*);
static int         cube_check_inside   (HELEMENT3D*,double G_CONST*);
static void        cube_coord_to_world (HELEMENT3D*,double G_CONST*,double*);
static void        cube_coord_of_parent(HELEMENT3D*,double*,double*);
static HELEMENT3D* cube_neighbour      (HELEMENT3D*,int,int,double*,double*,MESH_ELEMENT_FLAGS);
static HELEMENT3D* get_element(GENMESH3D* gmesh);
static HELEMENT3D* cube_copy_element(HELEMENT3D*,MESH_ELEMENT_FLAGS);

template <typename REAL>
static void scalar_function3d(ELEMENT3D*,int,double[],double[],void*);
template <typename REAL>
static void vector_function3d(ELEMENT3D*,int,double[],double[],void*);
static void cube_f_el_info (ELEMENT3D*,F_EL_INFO3D *,void*);
template <typename REAL>
static void cube_get_bounds_scalar (ELEMENT3D*, double *, double *, void *function_data);
static void cube_get_bounds_vector (ELEMENT3D*, double *, double *, void *function_data);

/* now the DESCRIPTION of the elements, here simple cubes */
static HELEMENT3D_DESCRIPTION cube_description;

namespace {
/* fill the description with some useful things */
void fillDescription(void)
{
  cube_description.dindex = 0;
  cube_description.number_of_vertices = 8;
  cube_description.number_of_polygons = 6;
  cube_description.polygon_length =poly_length;
  cube_description.polygon_vertex = poly_vertex;
  cube_description.polygon_neighbour = poly_neighbour;
  cube_description.dimension_of_coord = 3;
  cube_description.coord = cube_coords;
  cube_description.parametric_degree = 1;

  cube_description.world_to_coord = cube_world_to_coord;
  cube_description.coord_to_world = cube_coord_to_world;
  cube_description.check_inside = cube_check_inside;
  cube_description.neighbour = cube_neighbour;
  cube_description.boundary = cube_boundary;
  cube_description.coord_of_parent = cube_coord_of_parent;
}
}

/* filling an element with some useful informations */
static HELEMENT3D* fillHELEMENT3D(HELEMENT3D* helement, GENMESH3D* gmesh)
{
  helement->mesh = gmesh;
  helement->eindex = 0;
  helement->descr = reinterpret_cast<ELEMENT3D_DESCRIPTION*> (&cube_description);
  helement->present = static_cast<MESH_ELEMENT_FLAGS> (hefVertex | hefVindex | hefVinh);
  helement->user_data = NULL;
  helement->parent = NULL;
  helement->vinh = NULL;
  helement->ref_rule = 0;
  helement->level = 0;
  helement->has_children = 1;

  return(helement);
}

/* filling the function data with some useful informations */
template <typename REAL>
void fillScalarDataFunction3d(F_HDATA3D* scal_function, functionComponents<REAL> * fC)
{
  scal_function->name                 = const_cast<char*>("Skalare 3D-Daten");
  scal_function->dimension_of_value   = 1;
  scal_function->continuous_data      = 1;
  scal_function->f                    = scalar_function3d<REAL>;
  scal_function->last                 = NULL;
  scal_function->next                 = NULL;
  scal_function->get_bounds           = cube_get_bounds_scalar<REAL>;
  scal_function->threshold            = 0.0;
  scal_function->geometry_threshold   = 0.0;
  scal_function->hp_threshold         = 0.0;
  scal_function->hp_maxlevel          = 0;
  scal_function->f_el_info            = cube_f_el_info;
  scal_function->get_vertex_estimate  = fC->get_vertex_estimate;
  scal_function->get_element_estimate = fC->get_element_estimate;
  scal_function->function_data        = fC;
}

/* filling the function data with some useful informations */
template <typename REAL>
void fillVectorDataFunction3d(F_HDATA3D* vec_function, aol::MultiVector<REAL>* data)
{
  if (! data->allDimsEqual())
    throw aol::Exception("\n\nERROR: The components of the MultiVector are NOT of the same size!!!!!!!\n\n"
      ,__FILE__, __LINE__ );

  vec_function->name                 = const_cast<char*>("Vektordaten");
  vec_function->dimension_of_value   = (*data).numComponents();
  vec_function->continuous_data      = 1;
  vec_function->f                    = vector_function3d<REAL>;
  vec_function->last                 = NULL;
  vec_function->next                 = NULL;
  vec_function->get_bounds           = cube_get_bounds_vector;
  vec_function->threshold            = 0.0;
  vec_function->geometry_threshold   = 0.0;
  vec_function->hp_threshold         = 0.0;
  vec_function->hp_maxlevel          = 0;
  vec_function->f_el_info            = cube_f_el_info;
  vec_function->get_vertex_estimate  = NULL;
  vec_function->get_element_estimate = NULL;
  vec_function->function_data        = data;
}





/* ******************************************************************************************
*
*  the functions of the description
* ******************************************************************************************** */

static void cube_coord_of_parent(HELEMENT3D* /*helem*/,double* /*p*/,double* /*q*/)
{
  fprintf(stderr,"coord_of_parent: not implemented yet!\n");

  return;
}

static HELEMENT3D *cube_neighbour(HELEMENT3D */*helem*/, int /*np*/, int /*flag*/,
          double */*coord*/, double */*xyz*/,
          MESH_ELEMENT_FLAGS /*fl*/)
{
  fprintf(stderr,"neighbour: not implemented yet!\n");

  return NULL;
}

// this is an auxiliary function for testing whether a
// point with index ind is an inner or a boundary point.
// N is the number of elements in one line
int is_bnd3d(int ind, int N)
{
  // the elements on the undermost and on the ceiling plane
  // are boundary-elements
  if ( (ind >= N*N*(N-1)) || (ind < N*N) ) return(1);

  // not there, then calc the index mod N*N => the same case as 2d
  ind %= (N*N);

  // the indices of the elements of the undermost line are < N,
  // on the left side they are multiple of N,
  // on the right side they are multiple of N-1,
  // on the ceiling line the index is bigger or equal N*(N-1)
  if ( (ind<N) || ((ind%N)==0) || ((ind%N)==(N-1)) || (ind>=N*(N-1)) )
  {
    return(1);
  }
  else return(0);
}

static int cube_boundary(HELEMENT3D *helement, int np)
{
  // first get the number of vertices in one line
  int N = (1<<helement->mesh->max_level)+1;

  // now test whether all vertices are at the boundary or not
  if (is_bnd3d(helement->vindex[cube_description.polygon_vertex[np][0]],N) &&
      is_bnd3d(helement->vindex[cube_description.polygon_vertex[np][1]],N) &&
      is_bnd3d(helement->vindex[cube_description.polygon_vertex[np][2]],N) &&
      is_bnd3d(helement->vindex[cube_description.polygon_vertex[np][3]],N)   )
    return(-1);
  else
    return(0);
}

static int cube_world_to_coord(HELEMENT3D* helement,double G_CONST* xyz,double* coord)
{
  // simply subtract the coordinates of the first vertex of the element
  // but first get the absolute width/height of the element
  double h = helement->vertex[1][0] - helement->vertex[0][0];
  coord[0] = (xyz[0] - helement->vertex[0][0])/h;
  coord[1] = (xyz[1] - helement->vertex[0][1])/h;
  coord[2] = (xyz[2] - helement->vertex[0][2])/h;
  return(cube_check_inside(helement,coord));
}


static void cube_coord_to_world(HELEMENT3D* helement,double G_CONST* coord,double* xyz)
{
  // simply add the coordinates of the first vertex of the element
  // but first get the absolute width/height of the element
  double h = helement->vertex[1][0] - helement->vertex[0][0];
  xyz[0] = coord[0]*h + helement->vertex[0][0];
  xyz[1] = coord[1]*h + helement->vertex[0][1];
  xyz[2] = coord[2]*h + helement->vertex[0][2];
  return;
}


static int cube_check_inside(HELEMENT3D */*helement*/, double G_CONST*coord)
{
  if ((coord[0]>=0) && (coord[0]<=1) &&
      (coord[1]>=0) && (coord[1]<=1) &&
      (coord[2]>=0) && (coord[2]<=1) )
    return(-1);
  else      // point is not inside => finde the index of the separating face
  {
    if (coord[0]<0) return(5);
    if (coord[0]>1) return(3);
    if (coord[1]<0) return(2);
    if (coord[1]>1) return(4);
    if (coord[2]<0) return(0);
    if (coord[2]>1) return(1);
  }
  cout<<"cube_heck_inside: Paradoxically the point is neither inside nor outside...???"<<endl;
  return(-1);
}


/* **********************************************************
*
*    the functions of the HMESH
*
* ********************************************************** */


static HELEMENT3D* cube_first_child(HELEMENT3D* parent,MESH_ELEMENT_FLAGS /*fl*/)
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
    int dz = (parent->vindex[4]-parent->vindex[0])/2;

    // get an element to fill
    HELEMENT3D* elementret = get_element(parent->mesh);
    elementret->vindex[0] = parent->vindex[0];
    elementret->vindex[1] = parent->vindex[0]+dx;
    elementret->vindex[2] = parent->vindex[0]+dx+dy;
    elementret->vindex[3] = parent->vindex[0]+dy;
    elementret->vindex[4] = parent->vindex[0]+dz;
    elementret->vindex[5] = parent->vindex[0]+dx+dz;
    elementret->vindex[6] = parent->vindex[0]+dx+dy+dz;
    elementret->vindex[7] = parent->vindex[0]+dy+dz;

    /* now the coordinates */
    double h = 1. / (static_cast<double> (1<<level));

    double x = parent->vertex[0][0];
    double y = parent->vertex[0][1];
    double z = parent->vertex[0][2];

    elementret->vertex[0][0] = x; elementret->vertex[0][1] = y; elementret->vertex[0][2] = z;
    elementret->vertex[1][0] = x+h; elementret->vertex[1][1] = y; elementret->vertex[1][2] = z;
    elementret->vertex[2][0] = x+h; elementret->vertex[2][1] = y+h; elementret->vertex[2][2] = z;
    elementret->vertex[3][0] = x; elementret->vertex[3][1] = y+h; elementret->vertex[3][2] = z;
    elementret->vertex[4][0] = x; elementret->vertex[4][1] = y; elementret->vertex[4][2] = z+h;
    elementret->vertex[5][0] = x+h; elementret->vertex[5][1] = y; elementret->vertex[5][2] = z+h;
    elementret->vertex[6][0] = x+h; elementret->vertex[6][1] = y+h; elementret->vertex[6][2] = z+h;
    elementret->vertex[7][0] = x; elementret->vertex[7][1] = y+h; elementret->vertex[7][2] = z+h;

    /* the number of the children is saved in the eindex-entry,
      running from 0 to 7, 0 means: first child */
    elementret->eindex = 0;
    elementret->parent = parent;
    elementret->level = level;
    elementret->vinh = vinhChild0;
    if (level == parent->mesh->max_level) elementret->has_children=0;

    return reinterpret_cast<HELEMENT3D*>(elementret);
  }
  else
  {
    return NULL;
  }
}



static HELEMENT3D* cube_next_child(HELEMENT3D* helement,MESH_ELEMENT_FLAGS /*fl*/)
{
  int level = helement->level;

  // level >0, there are eight children, if its already the eighth one => return NULL
  if ((level>0) && (helement->eindex < 7))
  {
    /* now adjust the necessary entries of the element-structure */
    ++(helement->eindex);          // increase number of the element

    // calc the vertical and horizontal number-distance
    int dx = helement->vindex[1]-helement->vindex[0];
    int dy = helement->vindex[3]-helement->vindex[0];
    int dz = helement->vindex[4]-helement->vindex[0];

    // calc the coord-distance
    double h = 1. / (static_cast<double> (1<<level) );

    // now, dependent on the number of the child-element, calc
    // its next position and indices
    switch(helement->eindex)
    {
      int i;
      // eine Position weiter rechts (alles unten und vorne)
      case 1:
          // y- and z-coords remain, only translate in x-direction
          for (i=0; i<8; ++i) helement->vindex[i] += dx;
          for (i=0; i<8; ++i) helement->vertex[i][0] += h;
          helement->vinh = vinhChild1;
        break;

      // links hinten
      case 2:
          // z-coords remain, translate in x- and y-direction
          for (i=0; i<8; ++i) helement->vindex[i] += (dy-dx);
          for (i=0; i<8; ++i)
          {
            helement->vertex[i][0] -= h;
            helement->vertex[i][1] += h;
          }
          helement->vinh = vinhChild2;
        break;
      // rechts hinten
      case 3:
          // y- and z-coords remain, only translate in x-direction
          for (i=0; i<8; ++i) helement->vindex[i] += dx;
          for (i=0; i<8; ++i) helement->vertex[i][0] += h;
          helement->vinh = vinhChild3;
        break;
      // jetzt Wechsel von hinten unten rechts nach oben vorne links
      case 4:
          // translate in  every direction
          for (i=0; i<8; ++i) helement->vindex[i] += (dz-dx-dy);
          for (i=0; i<8; ++i)
          {
            helement->vertex[i][0] -= h;
            helement->vertex[i][1] -= h;
            helement->vertex[i][2] += h;
          }
          helement->vinh = vinhChild4;
        break;
        // eine Position weiter rechts (alles oben und vorne)
      case 5:
          // y- and z-coords remain, only translate in x-direction
          for (i=0; i<8; ++i) helement->vindex[i] += dx;
          for (i=0; i<8; ++i) helement->vertex[i][0] += h;
          helement->vinh = vinhChild5;
        break;
      // links hinten
      case 6:
          // z-coords remain, translate in x- andy y-direction
          for (i=0; i<8; ++i) helement->vindex[i] += (dy-dx);
          for (i=0; i<8; ++i)
          {
            helement->vertex[i][0] -= h;
            helement->vertex[i][1] += h;
          }
          helement->vinh = vinhChild6;
        break;
      // rechts hinten
      case 7:
          // y- and z-coords remain, only translate in x-direction
          for (i=0; i<8; ++i) helement->vindex[i] += dx;
          for (i=0; i<8; ++i) helement->vertex[i][0] += h;
          helement->vinh = vinhChild7;
        break;
      default:
        throw aol::UnimplementedCodeException ( "Unhandled switch case!", __FILE__, __LINE__ );
    };

    return reinterpret_cast<HELEMENT3D*> (helement);

  }
  else
  {
    cube_free_element(helement);
    return(NULL);
  }
}
static HELEMENT3D* cube_nth_child(HELEMENT3D *elementret, HELEMENT3D* parent, int n, MESH_ELEMENT_FLAGS /*fl*/)
{

//   if (el)
//     cube_free_element(el);
//   el = cube_first_child(parent, fl);
//   for (int i = 0; i < n; i++)
//     el = cube_next_child(el, fl);
//   return el;
//
//
//
  // increase the level
  int level = parent->level;
  ++level;

  // calc the vertical and horizontal number-distance
  int dx = (parent->vindex[1]-parent->vindex[0])/2;
  int dy = (parent->vindex[3]-parent->vindex[0])/2;
  int dz = (parent->vindex[4]-parent->vindex[0])/2;


  // if it is the first element get an element to fill
  if (!elementret) elementret = get_element(parent->mesh);
//     HELEMENT3D* elementret = get_element(parent->mesh);

  int vindexOffset = childOffsetsX[n]*dx + childOffsetsY[n]*dy + childOffsetsZ[n]*dz;
  elementret->vindex[0] = parent->vindex[0] + vindexOffset;
  elementret->vindex[1] = parent->vindex[0] + dx + vindexOffset;
  elementret->vindex[2] = parent->vindex[0] + dx + dy + vindexOffset;
  elementret->vindex[3] = parent->vindex[0] + dy + vindexOffset;
  elementret->vindex[4] = parent->vindex[0] + dz + vindexOffset;
  elementret->vindex[5] = parent->vindex[0] + dx + dz + vindexOffset;
  elementret->vindex[6] = parent->vindex[0] + dx + dy + dz + vindexOffset;
  elementret->vindex[7] = parent->vindex[0] + dy + dz + vindexOffset;

  /* now the coordinates */
  double h = 1. / (static_cast<double> (1<<level));

  double x = parent->vertex[0][0] + childOffsetsX[n]*h;
  double y = parent->vertex[0][1] + childOffsetsY[n]*h;
  double z = parent->vertex[0][2] + childOffsetsZ[n]*h;

  elementret->vertex[0][0] = x; elementret->vertex[0][1] = y; elementret->vertex[0][2] = z;
  elementret->vertex[1][0] = x+h; elementret->vertex[1][1] = y; elementret->vertex[1][2] = z;
  elementret->vertex[2][0] = x+h; elementret->vertex[2][1] = y+h; elementret->vertex[2][2] = z;
  elementret->vertex[3][0] = x; elementret->vertex[3][1] = y+h; elementret->vertex[3][2] = z;
  elementret->vertex[4][0] = x; elementret->vertex[4][1] = y; elementret->vertex[4][2] = z+h;
  elementret->vertex[5][0] = x+h; elementret->vertex[5][1] = y; elementret->vertex[5][2] = z+h;
  elementret->vertex[6][0] = x+h; elementret->vertex[6][1] = y+h; elementret->vertex[6][2] = z+h;
  elementret->vertex[7][0] = x; elementret->vertex[7][1] = y+h; elementret->vertex[7][2] = z+h;

  /* the number of the children is saved in the eindex-entry,
  running from 0 to 7, 0 means: first child */
  elementret->parent = parent;
  elementret->level = level;
  switch (n) { // could introduce an array...
    case 0: elementret->vinh = vinhChild0; break;
    case 1: elementret->vinh = vinhChild1; break;
    case 2: elementret->vinh = vinhChild2; break;
    case 3: elementret->vinh = vinhChild3; break;
    case 4: elementret->vinh = vinhChild4; break;
    case 5: elementret->vinh = vinhChild5; break;
    case 6: elementret->vinh = vinhChild6; break;
    case 7: elementret->vinh = vinhChild7; break;
    default: throw aol::UnimplementedCodeException ( "Unhandled switch case!", __FILE__, __LINE__ );
  }
  if (level == parent->mesh->max_level) elementret->has_children=0;

  return reinterpret_cast<HELEMENT3D*> (elementret);


}

static HELEMENT3D* cube_sorted_first_child(HELEMENT3D* parent,MESH_ELEMENT_FLAGS fl)
{
  int level = parent->level;

  /* still not the finest level? */
  if (level < parent->mesh->max_level)
  {
    HELEMENT3D *elementret = cube_nth_child(NULL, parent, reinterpret_cast<QUOCMESH3D *>(parent->mesh)->sort_order[0], fl);
    elementret->eindex = 0;
    return elementret;
  }  else  {
    return NULL;
  }
}

static HELEMENT3D* cube_sorted_next_child(HELEMENT3D* helement,MESH_ELEMENT_FLAGS fl)
{
  helement->eindex++;
//   cerr << ", eindex Anfang: "<<helement->eindex;
  if (helement->eindex < 8) {
    int tmp = helement->eindex;
    helement = cube_nth_child(helement, helement->parent, reinterpret_cast<QUOCMESH3D *>(helement->mesh)->sort_order[helement->eindex], fl);
    helement->eindex = tmp;
//     cerr << ", eindex Ende: "<<helement->eindex<<endl;
    return helement;
  }
  else
  {
    cube_free_element(helement);
    return(NULL);
  }
}

static HELEMENT3D* cube_select_child(ELEMENT3D* parent,
             double* parent_coord,double* child_coord,MESH_ELEMENT_FLAGS fl)
{
  if (   (parent_coord[0] > 0.) && (parent_coord[0] < 1.) &&
        (parent_coord[1] > 0.) && (parent_coord[1] < 1.) &&
        (parent_coord[2] > 0.) && (parent_coord[2] < 1.) )
  {
    // get a copy of the parent element
    HELEMENT3D* child_el = cube_copy_element(parent, fl);
    //return child_el;

    // adjust the entries of the parent structure
    ++child_el->level;

    // calc the vertical and horizontal number-distance
    int dx = parent->vindex[1]-parent->vindex[0];
    int dy = parent->vindex[3]-parent->vindex[0];
    int dz = parent->vindex[4]-parent->vindex[0];
    int i;

    // calc the coord-distance
    //double h = 1. / ((double) (1<<(child_el->level)) );
    double h = child_el->vertex[1][0] - child_el->vertex[0][0];


    // get the number of the child element (0..7)
    int child_number = 0;
    if (parent_coord[0] > 0.5) ++child_number;
    if (parent_coord[1] > 0.5) child_number += 2;
    if (parent_coord[2] > 0.5) child_number += 4;

    // now, depending on the child_number, adjust the coords and
    // vindex-numbers

    // the lower four child-elements have at their celing decreased z-coords
    if (child_number<4)
    {
      for (i=4; i<8; ++i)      // z-coord
      {
          child_el->vertex[i][2] -= h;
          child_el->vindex[i] -= dz;
      }
      child_coord [2] = parent_coord [2] * 2.0;
    }
    // the ceilng four child-elements have at their ground increased z-coords
    if (child_number>3)
    {
      for (i=0; i<4; ++i)      // z-coord
      {
          child_el->vertex[i][2] += h;
          child_el->vindex[i] += dz;
      }
      child_coord [2] = (parent_coord [2] - 0.5) * 2.0;
    }

    // links vorne
    if (child_number==0 || child_number==4)
    {
      if (child_number==0 ) child_el->vinh = vinhChild0; else child_el->vinh = vinhChild4;
      child_el->vindex[1] -= dx; child_el->vindex[2] -= dx+dy; child_el->vindex[3] -= dy;
      child_el->vindex[5] -= dx; child_el->vindex[6] -= dx+dy; child_el->vindex[7] -= dy;

      child_el->vertex[1][0] -= h; child_el->vertex[2][0] -= h;
      child_el->vertex[2][1] -= h; child_el->vertex[3][1] -= h;
      child_el->vertex[5][0] -= h; child_el->vertex[6][0] -= h;
      child_el->vertex[6][1] -= h; child_el->vertex[7][1] -= h;

      child_coord [0] = parent_coord [0] * 2.0;
      child_coord [1] = parent_coord [1] * 2.0;
    }
    // rechts vorne
    if (child_number==1 || child_number==5)
    {
      if (child_number==1 ) child_el->vinh = vinhChild1; else child_el->vinh = vinhChild5;
      child_el->vindex[0] += dx; child_el->vindex[2] -= dy; child_el->vindex[3] += dx-dy;
      child_el->vindex[4] += dx; child_el->vindex[6] -= dy; child_el->vindex[7] += dx-dy;

      child_el->vertex[0][0] += h; child_el->vertex[2][1] -= h;
      child_el->vertex[3][0] += h; child_el->vertex[3][1] -= h;
      child_el->vertex[4][0] += h; child_el->vertex[6][0] -= h;
      child_el->vertex[7][0] += h; child_el->vertex[7][1] -= h;

      child_coord [0] = (parent_coord [0] - 0.5) * 2.0;
      child_coord [1] = parent_coord [1] * 2.0;
    }
    // links hinten
    if (child_number==2 || child_number==6)
    {
      if (child_number==2 ) child_el->vinh = vinhChild2; else child_el->vinh = vinhChild6;
      child_el->vindex[0] += dy; child_el->vindex[1] += dy-dx; child_el->vindex[2] -= dx;
      child_el->vindex[4] += dy; child_el->vindex[5] += dy-dx; child_el->vindex[6] -= dx;

      child_el->vertex[0][1] += h; child_el->vertex[1][0] -= h;
      child_el->vertex[1][1] += h; child_el->vertex[2][0] -= h;
      child_el->vertex[4][1] += h; child_el->vertex[5][0] -= h;
      child_el->vertex[5][1] += h; child_el->vertex[6][0] -= h;

      child_coord [0] = parent_coord [0] * 2.0;
      child_coord [1] = (parent_coord [1] - 0.5) * 2.0;
    }
    // rechts hinten
    if (child_number==3 || child_number==7)
    {
      if (child_number==3 ) child_el->vinh = vinhChild3; else child_el->vinh = vinhChild7;
      child_el->vindex[0] += dx+dy; child_el->vindex[1] += dy; child_el->vindex[3] += dx;
      child_el->vindex[4] += dx+dy; child_el->vindex[5] += dy; child_el->vindex[7] += dx;

      child_el->vertex[0][0] += h; child_el->vertex[0][1] += h;
      child_el->vertex[1][1] += h; child_el->vertex[3][0] += h;
      child_el->vertex[4][0] += h; child_el->vertex[4][1] += h;
      child_el->vertex[5][1] += h; child_el->vertex[7][0] += h;

      child_coord [0] = (parent_coord [0] - 0.5) * 2.0;
      child_coord [1] = (parent_coord [1] - 0.5) * 2.0;
    }

    return child_el;

  }
  else return NULL;    // not a child element
}

int g_projecting_matrix (const MATRIX44 mat)
{
  return mat[3][0] || mat[3][1] || mat[3][2] || 1 != mat[3][3];
}

double *g_matrix44_lift_dir (VEC3 dir, G_CONST MATRIX44 mat, const VEC3 view)
{
  VEC4 view4, dir4;
  if (g_projecting_matrix (mat))
    return NULL;
  g_vec3_assign (view4, view);
  view4[3] = 0;
  g_solve4 (mat, view4, dir4);
  g_vec3_assign (dir, dir4);
  return dir;
}

static HELEMENT3D* cube_first_macro(GENMESH3D* hmesh,MESH_ELEMENT_FLAGS /*fl*/)
{

  // get an element to fill
  HELEMENT3D* elementret = get_element(hmesh);

  // set the traverse mode (sorted or fast)
  if (hmesh->access_mode & mafSorted) {
    hmesh->first_child = cube_sorted_first_child;
    hmesh->next_child = cube_sorted_next_child;

    QUOCMESH3D *self = reinterpret_cast<QUOCMESH3D *> (hmesh);
    const VEC3 view_dir = {0,0,1};
    VEC3 sort_direction;
    g_matrix44_lift_dir (sort_direction, self->sort_matrix, view_dir);
    const VEC3 gc[8] = {
      {0.0, 0.0, 0.0},
      {0.5, 0.0, 0.0},
      {0.0, 0.5, 0.0},
      {0.5, 0.5, 0.0},
      {0.0, 0.0, 0.5},
      {0.5, 0.0, 0.5},
      {0.0, 0.5, 0.5},
      {0.5, 0.5, 0.5}
    };
    double depth[8];
    int i;
    for (i = 0; i < 8; i++)
      depth[i] = g_vec3_skp(sort_direction, gc[i]);
    for (int j = 0; j < 8; j++)
      for (i = 0; i < 7-j; i++)
        if (depth[self->sort_order[i]] > depth[self->sort_order[i+1]]) {
          int tmp = self->sort_order[i];
          self->sort_order[i] = self->sort_order[i+1];
          self->sort_order[i+1] = tmp;
        }

//     for (i=0; i<8; i++)
//       cerr << self->sort_order[i]<<", ";
//     cerr << endl;

  } else {
    hmesh->first_child = cube_first_child;
    hmesh->next_child = cube_next_child;
  }


  /* for the first macro the initialized coordinates and vindex
    has to be updated, eindex is already correct */
  /* first calc the indices of the edges: */
  int d = hmesh->max_level;
  int pd = 1<<d;          // 2^d
  int p2d = 1<<(2*d);      // 2^(2d)
  int p3d = 1<<(3*d);      // 2^(3d)

  elementret->vindex[0] = 0; elementret->vindex[1] = pd;
  elementret->vindex[3] = pd+p2d; elementret->vindex[2] = p2d+2*pd;
  elementret->vindex[4] = p3d+2*p2d+1*pd; elementret->vindex[5] = p3d+2*p2d+2*pd;
  elementret->vindex[6] = p3d+3*p2d+3*pd; elementret->vindex[7] = p3d+3*p2d+2*pd;

  /* now the coordinates */
  elementret->vertex[0][0] = 0.; elementret->vertex[0][1] = 0.; elementret->vertex[0][2] = 0.;
  elementret->vertex[1][0] = 1.; elementret->vertex[1][1] = 0.; elementret->vertex[1][2] = 0.;
  elementret->vertex[2][0] = 1.; elementret->vertex[2][1] = 1.; elementret->vertex[2][2] = 0.;
  elementret->vertex[3][0] = 0.; elementret->vertex[3][1] = 1.; elementret->vertex[3][2] = 0.;
  elementret->vertex[4][0] = 0.; elementret->vertex[4][1] = 0.; elementret->vertex[4][2] = 1.;
  elementret->vertex[5][0] = 1.; elementret->vertex[5][1] = 0.; elementret->vertex[5][2] = 1.;
  elementret->vertex[6][0] = 1.; elementret->vertex[6][1] = 1.; elementret->vertex[6][2] = 1.;
  elementret->vertex[7][0] = 0.; elementret->vertex[7][1] = 1.; elementret->vertex[7][2] = 1.;
  elementret->level = 0;

  return reinterpret_cast<HELEMENT3D*> (elementret);
}

static HELEMENT3D* cube_next_macro(HELEMENT3D* helement,MESH_ELEMENT_FLAGS /*fl*/)
{
  /* There is only one macro-element => return NULL */
  cube_free_element(helement);
  return NULL;
}


static void cube_free_element(HELEMENT3D *el)
{
  // put the element into the free-element-list
  elements_free.push_back(el);
  elements_used.pop_back();
  //cout<<"----------- qc::Element freigeworden, Size of free-List: "<<elements_free.size()<<endl;

  return;
}



HELEMENT3D* cube_copy_element(HELEMENT3D* helement,MESH_ELEMENT_FLAGS /*fl*/)
{
  // get an element to fill
  HELEMENT3D* elementret = get_element(helement->mesh);

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
  for (i=0; i<8; ++i)
    elementret->vindex[i] = helement->vindex[i];

  for (i=0; i<8; ++i)
    for (j=0; j<3; ++j)
      elementret->vertex[i][j] = helement->vertex[i][j];

  return(helement);
}




/* ***************************************************************
*
*  Now the functions for managing the used elements in a vector
*
* *************************************************************** */

HELEMENT3D* get_element(GENMESH3D* gmesh)
{
  // the necessary pointers
  HELEMENT3D* elementret;

  // no free element available => allocate a new one
  // and put it on the stack
  if (elements_free.size() == 0)
  {
    elementret = reinterpret_cast<HELEMENT3D*>(mem_alloc(sizeof(HELEMENT3D)));
    elementret->vindex = reinterpret_cast<int*>(mem_alloc(8*sizeof(int)));
    elementret->vertex = reinterpret_cast<double**>(mem_alloc(8*sizeof(double*)));
    for (int i=0; i<8; ++i)
    {
      elementret->vertex[i] = reinterpret_cast<double*>(mem_alloc(3*sizeof(double)));
      for (int j=0; j<3; ++j) elementret->vertex[i][j] = 0;
    }
    //cout<<"------------ Neues qc::Element deklariert, ";
  }
  else
  {
    // ok, there is a free element available, so take it
    // out of the vector
    elementret = elements_free.back();
    elements_free.pop_back();
    //cout<<"------------ Freies qc::Element verfuegbar, ";
  }
  elements_used.push_back(elementret);
  // initialize and return it
    //cout<<" S.o. used: "<<elements_used.size()<<
    //    ", S.o. free: "<<elements_free.size()<<endl;

  fillHELEMENT3D(elementret, gmesh);
  return reinterpret_cast<HELEMENT3D*>(elementret);
}

template <typename REAL>
static int cube_set_time(GENMESH3D* mesh,REAL t)
{
  reinterpret_cast<MESH_INFO*>(mesh->user_data)->t = t;

  F_HDATA3D* fun;
  for (fun = reinterpret_cast<F_HDATA3D*>(mesh->f_data); fun->last; fun = reinterpret_cast<F_HDATA3D*>(fun->last)) ; // Goto first function

  // Interpolate all functions
  for (; fun; fun = reinterpret_cast<F_HDATA3D*>(fun->next)) {

    if (fun->dimension_of_value == 1) {

      int i = static_cast<int> (floor (t));
      int j = i+1;

      double jfac = t - i;
      double ifac = 1 - jfac;

      functionComponents<REAL>* fC = reinterpret_cast<functionComponents<REAL>*>(fun->function_data);
      int n = fC->data_seq->size ();

      if (i < 0) { i = j = 0; ifac = 1; jfac = 0; }
      if (j >= n) { i = j = n-1; ifac = 1; jfac = 0; }

      qc::ScalarArray<REAL, qc::QC_3D>& idata = *((*fC->data_seq) [i]);
      qc::ScalarArray<REAL, qc::QC_3D>& jdata = *((*fC->data_seq) [j]);
      qc::ScalarArray<REAL, qc::QC_3D>& data = *(fC->data);

      for (int a = 0; a < data.size (); ++a) {
  data [a] = ifac * idata [a] + jfac * jdata [a];
      }

      if (!fC->estimator) continue;

      qc::Estimator3d<REAL>& iestimator = *((*fC->estimator_seq) [i]);
      qc::Estimator3d<REAL>& jestimator = *((*fC->estimator_seq) [j]);
      qc::Estimator3d<REAL>& estimator = *(fC->estimator);

      for (int a = 0; a < estimator.size (); ++a)
  estimator [a] = ifac * iestimator [a] + jfac * jestimator [a];
    }
    else {
      // Currently no time-dependency for vector valued functions
    }
  }

  return TRUE;
}

#ifdef _SQUARE_GET_TIME_USED_
// produces warning: defined but not used
static int square_get_time(GENMESH3D* mesh, double* t, double* t1, double* t2)
{
  *t = ((MESH_INFO*)mesh->user_data)->t;
  *t1 = floor (*t); *t2 = ceil (*t);

  functionComponents<double>* fC = (functionComponents<double>*)(((F_HDATA3D*)mesh->f_data)->function_data);
  int n = fC->data_seq->size ();
  if (*t1 < 0) *t1 = *t2 = 0;
  if (*t2 >= n) *t1 = *t2 = n-1;

  return TRUE;
}
#endif

GENMESH3D* genmesh3d_get_times (double* t1, double* t2)
{
  GENMESH3D* self = reinterpret_cast<GENMESH3D*> (START_METHOD( G_INSTANCE ));
  ASSURE(self,"",END_METHOD(NULL));

  *t1 = 0;

  functionComponents<double>* fC = reinterpret_cast<functionComponents<double>*>(reinterpret_cast<F_HDATA3D*>(self->f_data)->function_data);
  *t2 = static_cast<int>(fC->data_seq->size () - 1);

  END_METHOD (self);
}



/* ***************************************************************
*
*  Now the functions itself
*
* *************************************************************** */


// first: the easy case of a scalar function
template <typename REAL>
static void scalar_function3d(HELEMENT3D* helement,
          int        index,
          double     coord[],
          double     val[],
          void*      function_data)
{
  //qc::ScalarArray<double, qc::QC_3D>* meshData = (qc::ScalarArray<double, qc::QC_3D>*)function_data;
  qc::ScalarArray<REAL, qc::QC_3D>* meshData = reinterpret_cast<functionComponents<REAL>*>(function_data)->data;

  // 8 function-values for the trilinear interpolation
  double f[8];

  if(coord)
  {
    double komplx = 1.-coord[0];
    double komply = 1.-coord[1];
    double t1,t2,t3,t4;

    for (int i=0; i<8; ++i)
      f[i] = (*meshData)[helement->vindex[i]];

    t1 = coord[0]*f[6] + komplx*f[7];
    t2 = coord[0]*f[5] + komplx*f[4];
    t3 = coord[0]*f[2] + komplx*f[3];
    t4 = coord[0]*f[1] + komplx*f[0];

    t1 = coord[1]*t1 + komply*t2;
    t2 = coord[1]*t3 + komply*t4;

    val[0] = coord[2]*t1 + (1.-coord[2])*t2;

  }
  else val[0] = (*meshData)[helement->vindex[index]];

  return;
}




template <typename REAL>
// second: the case of a n-dimensional vector-valued function
static void vector_function3d(HELEMENT3D* helement,
          int        index,
          double     coord[],
          double     val[],
          void*      function_data)
{

  aol::MultiVector<REAL>* meshDataVec = reinterpret_cast<aol::MultiVector<REAL>*> (function_data);
  int dim = helement->mesh->f_data->dimension_of_value;    // dimension of the multivector
  int i,j;
  double f[8], t1,t2,t3,t4,komplx,komply;
  int ind[8];

  // get the indices of the eight vertices
  for (i=0; i<8; ++i)
    ind[i] = helement->vindex[i];

  if(coord)
  {
    komplx = 1.-coord[0];
    komply = 1.-coord[1];

    // linear interpolating of each coordinate
    for (i=0; i<dim; ++i)
    {
      for (j=0; j<8; ++j)
        f[j] = (*meshDataVec)[i][ind[j]];

      t1 = coord[0]*f[6] + komplx*f[7];
      t2 = coord[0]*f[5] + komplx*f[4];
      t3 = coord[0]*f[2] + komplx*f[3];
      t4 = coord[0]*f[1] + komplx*f[0];

      t1 = coord[1]*t1 + komply*t2;
      t2 = coord[1]*t3 + komply*t4;

      val[i] = coord[2]*t1 + (1.-coord[2])*t2;
    }
  }
  else
  {
    ind[0] = helement->vindex[index];
    for (i=0; i<dim; ++i)
      val[i] = (*meshDataVec)[i][ind[0]];
  }

  return;
}






/* ***************************************************************
*
*  Now the ESTIMATORS of the function-data
*
* *************************************************************** */

// the element-estimator calcs the maximum of the error values
// of the inner vertices (one level higher)
template <typename REAL>
double calcElementError (HELEMENT3D *el, void *function_data)
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
    int mx,my,mz,ind,halfstep;
    int zwei_d = (1<<max_level) + 1;
    int zwei_d_quad = zwei_d * zwei_d;

    // We now have to get the scalar2d-indices of the left upper vertex
    // of the element. Therefore calc the x- and y-indices-distance.
    ind = el->vindex[0];
    halfstep = (el->vindex[1] - ind) >> 1;

    mx = ind % zwei_d;
    my = (ind % zwei_d_quad) / zwei_d;
    mz = ind / zwei_d_quad;

    qc::Estimator3d<double>* Est3d = reinterpret_cast<qc::Estimator3d<double>*>(reinterpret_cast<functionComponents<REAL>*>(function_data)->estimator);

    // get the maximum of the five inner points (one level finer) of the element
    double e = Est3d->maxChildNodeValue(mx,my,mz,halfstep);

    return(e);
  }
  return 0;

}


// the vertex-estimator reads the error out of the saturated error structure
template <typename REAL>
void calcVertexError (HELEMENT3D *el, double *error, void *function_data)
{
  qc::Estimator3d<double>* Est3d = reinterpret_cast<qc::Estimator3d<double>*>(reinterpret_cast<functionComponents<REAL>*>(function_data)->estimator);

  for (int i=0; i<8; ++i)
    error[i] = (*Est3d)[ el->vindex[i] ];

}


static void cube_f_el_info(ELEMENT3D* /*el*/,F_EL_INFO3D* hel_info,void* /*data*/)
{
  hel_info->polynomial_degree = 1;

  return;
}
template <typename REAL>
static void cube_get_bounds_scalar (ELEMENT3D* helement, double *min, double *max, void *function_data)
{

  // if there`s an error estimator use it
  if (reinterpret_cast<functionComponents<REAL>*> (function_data)->get_element_estimate)
  {
    // get data-array
    qc::ScalarArray<REAL, qc::QC_3D>* meshData =
      reinterpret_cast<qc::ScalarArray<REAL, qc::QC_3D>*>(reinterpret_cast<functionComponents<REAL>*> (function_data)->data);

    double a;
    *max = *min = (*meshData)[helement->vindex[0]];

    // get the min and max function-value of this element
    for (int i=1; i<8; ++i)
    {
      a = (*meshData)[helement->vindex[i]];
      if (a > *max) *max=a;
      if (a < *min) *min=a;
    }

    // get the error of the element
    a = static_cast<double> (reinterpret_cast<functionComponents<REAL>*> (function_data)->get_element_estimate(helement,function_data));

    // and add respectively subtract the error
    *max += a;
    *min -= a;
  }
  else
  {
    *max = std::numeric_limits<double>::max();
    *min = -std::numeric_limits<double>::max();
  }

  return;
}

// function returns the max and min length of the vector-valued data defined
// on one element

static void cube_get_bounds_vector (ELEMENT3D* /*helement*/, double *min, double *max, void */*function_data*/)
{
  /*int ind;    // index of the vertex
  double x;

  // get the data
  aol::MultiVector<double>* meshDataMV = (aol::MultiVector<double>*)function_data;

  // initialize min and max with first vertex
  ind = helement->vindex[0];

  *min = (*meshDataMV)[ind].norm();
  *max = (*meshDataMV)[ind].norm();

  // now get the function values of every vertex and save min and max
  for (int i=1; i<8; ++i)
  {
    ind = helement->vindex[i];
    // read the function value and adjust min and max
    x = (*meshDataMV)[ind].norm();
    if (x>*max) *max=x;
    if (x<*min) *min=x;

  }*/

  *max=std::numeric_limits<double>::max();
  *min=-std::numeric_limits<double>::max();

  return;
}


/* *********************************************************************
*
*    Taetaetae: the convert-functions
*    expects a QuocMesh-Grid and provides a GenMesh
*    the other argument is either a ScalarArray<QC_3D> or a MultiVector
*
* ********************************************************************* */

// First data-independent entries are filled, args are the mesh and the depth of it
void fillGMesh(QUOCMESH3D* gmesh, int d)
{
  gmesh->first_child            = cube_first_child;
  gmesh->next_child             = cube_next_child;
  gmesh->select_child            = cube_select_child;
  gmesh->first_macro            = cube_first_macro;
  gmesh->next_macro             = cube_next_macro;

  gmesh->free_element           = cube_free_element;
  gmesh->max_level              = d;
  gmesh->level_of_interest      = d;
  gmesh->max_vindex             = ((1<<d)+1)*((1<<d)+1)*((1<<d)+1);
  gmesh->max_dindex             = 1;
  gmesh->max_eindex             = 1<<(3*d);

  gmesh->access_mode            = static_cast<MESH_ACCESS_FLAGS> (mafNone);
  gmesh->access_capability      = static_cast<MESH_ACCESS_FLAGS> (mafNone);

  gmesh->max_dimension_of_coord = 3;
  gmesh->max_number_of_vertices = 8;
  gmesh->set_time               = cube_set_time;

  gmesh->get_geometry_vertex_estimate  = NULL;
  gmesh->get_geometry_element_estimate = NULL;

  gmesh->access_capability = mafSorted;
  g_matrix44_set_identity(gmesh->sort_matrix);
  for (int i = 0; i < 8; i++)
    gmesh->sort_order[i] = i;

  MESH_INFO* mesh_info = reinterpret_cast<MESH_INFO*> (mem_alloc(sizeof(MESH_INFO)));
  mesh_info->t = 0;
  gmesh->user_data = mesh_info;
}


/* ********************************************************************
    function tests, whether a scalar-array is quadratic or not,
    and provides (if its quadratic) the depth of it!
  ********************************************************************* */
template <typename REAL>
int depthOfScalarArray3d(qc::ScalarArray<REAL, qc::QC_3D>* data)
{
  int N = data->getNumX ();
  int d = qc::logBaseTwo (N);
  if(N == 2) d = 0;
  if ( (data->getNumY () != N) || (data->getNumZ () != N))
    throw aol::Exception ("Image not square", __FILE__, __LINE__);
  if ((1 << d) + 1 != N)
    throw aol::Exception ("ScalarArray<QC_3D>: dimension no power of two", __FILE__, __LINE__);

  // everything is ok (be lucky)
  return(d);
}

template <typename REAL>
int depthOfMultiVector3d(aol::MultiVector<REAL>* data)
{
  int N = (data->getTotalSize()) / (data->numComponents()) ;
  int d = qc::logBaseTwo(N) / 3;
  int erg = (1 << d) + 1;
  if ( erg*erg*erg != N)
    throw aol::Exception ("Multivector: dimension no power of two", __FILE__, __LINE__);

  // everything is ok (be lucky)
  return(d);
}


QUOCMESH3D *quocmesh3d_set_sort_matrix (MATRIX44 matrix)
{
  QUOCMESH3D *self;

  self = reinterpret_cast<QUOCMESH3D *> (START_METHOD(G_INSTANCE));
  ASSURE(self, "", END_METHOD(NULL));

  g_matrix44_assign (matrix, self->sort_matrix);

  END_METHOD(self);
}

static void make_class(void)
{
  if (!QuocMesh3d) {
    QuocMesh3d = reinterpret_cast<GRAPE_CLASS *> (GRAPE (GenMesh3d, "new-class")("QuocMesh3d", sizeof (QUOCMESH3D)));
    GRAPE (QuocMesh3d, "add-method") ("softcopy", genmesh3d_my_softcopy);
    GRAPE (QuocMesh3d, "add-method") ("set-sort-matrix", quocmesh3d_set_sort_matrix);
  }
}

/* *********************************************************************
* the convert-routine for scalar data
* ********************************************************************* */
template <typename REAL>
    GENMESH3D* quocmesh_convert_to_gmesh3d(qc::ScalarArray<REAL, qc::QC_3D>* meshData, const char* dataName)
{
  QUOCMESH3D* gmesh;

  make_class();
  // initialize the gmesh-structure
  gmesh = reinterpret_cast<QUOCMESH3D*> (GRAPE(QuocMesh3d,"new-instance")("Adaptive 3d-QuocMesh"));
  // fill the description of one element
  fillDescription();
  // Fill the gmesh-structure
  int d = depthOfScalarArray3d(meshData);
  fillGMesh(gmesh,d);//,grid);

  // allocate memory for the data structure
  F_HDATA3D* data_function = reinterpret_cast<F_HDATA3D*> (mem_alloc(sizeof(F_HDATA3D)));
  gmesh->f_data     = reinterpret_cast<GENMESH_FDATA*> (data_function);

  // fill the data structure
  functionComponents<REAL>* fC = reinterpret_cast<functionComponents<REAL>*> (mem_alloc(sizeof(functionComponents<REAL>)));
  fC->data = new qc::ScalarArray<REAL, qc::QC_3D> (*meshData);
  fC->estimator = NULL;
  fC->get_element_estimate = NULL;
  fC->get_vertex_estimate = NULL;

  fC->data_seq = new std::vector<qc::ScalarArray<REAL, qc::QC_3D>*>;
  fC->data_seq->push_back (meshData);
  fC->estimator_seq = new std::vector<qc::Estimator3d<REAL>*>;
  fC->estimator_seq->push_back (NULL);

  fillScalarDataFunction3d(data_function, fC);
  data_function->name = const_cast<char*>(dataName);

  return reinterpret_cast<GENMESH3D *>(gmesh);
}



/* *********************************************************************
* the convert-routine for vector valued data
* ********************************************************************* */
template <typename REAL>
GENMESH3D* quocmesh_convert_to_gmesh3d(aol::MultiVector<REAL>* meshData,
                                       const char* dataName) // = "Vektordaten" )
{
  QUOCMESH3D* gmesh;

  make_class();
  // initialize the gmesh-structure
  gmesh   = reinterpret_cast<QUOCMESH3D*> (GRAPE(QuocMesh3d,"new-instance")("Adaptive QuocMesh"));
  // fill the description of one element
  fillDescription();

  // Fill the gmesh-structure
  int d = depthOfMultiVector3d(meshData);
  fillGMesh(gmesh,d);

  // allocate memory for the data structure
  F_HDATA3D* data_function = reinterpret_cast<F_HDATA3D*> (mem_alloc(sizeof(F_HDATA3D)));
  gmesh->f_data = reinterpret_cast<GENMESH_FDATA*> (data_function);

  // fill the data structure
  fillVectorDataFunction3d(data_function, meshData);
  data_function->name = const_cast<char*>(dataName);

  return reinterpret_cast<GENMESH3D *>(gmesh);
}



/* ---------------------------------------------------------------------------------------
* NOW FOLLOWING THE CONVERT-ROUTINES WITH AN ESTIMATOR
* --------------------------------------------------------------------------------------- */

/* *********************************************************************
* the convert-routine for scalar data
* ********************************************************************* */

template <typename REAL>
GENMESH3D* quocmesh_convert_to_gmesh3d(qc::ScalarArray<REAL, qc::QC_3D>* meshData,
                                       qc::Estimator3d<REAL>* errorArray, const char* dataName) // = "Daten mit Est." )
{
  QUOCMESH3D* gmesh;

  make_class();
  // initialize the gmesh-structure
  gmesh = reinterpret_cast<QUOCMESH3D*> (GRAPE(QuocMesh3d,"new-instance")("Adaptive 3d-QuocMesh"));
  // fill the description of one element
  fillDescription();
  // Fill the gmesh-structure
  int d = depthOfScalarArray3d(meshData);
  int dE = depthOfScalarArray3d(errorArray);
  if (d != dE)
    throw aol::Exception ("ScalarArray and Estimator-Array have NOT the same depth!", __FILE__, __LINE__);

  fillGMesh(gmesh,d);

  // allocate memory for the data structure
  F_HDATA3D* data_function = reinterpret_cast<F_HDATA3D*> (mem_alloc(sizeof(F_HDATA3D)));
  gmesh->f_data     = reinterpret_cast<GENMESH_FDATA*> (data_function);

  // fill the data structure
  functionComponents<REAL>* fC = reinterpret_cast<functionComponents<REAL>*> (mem_alloc(sizeof(functionComponents<REAL>)));
  fC->data = new qc::ScalarArray<REAL, qc::QC_3D> (*meshData);

  fC->estimator = new qc::Estimator3d<REAL> (*errorArray);
  fC->get_element_estimate = &calcElementError<REAL>;
  fC->get_vertex_estimate = &calcVertexError<REAL>;

  fC->data_seq = new std::vector<qc::ScalarArray<REAL, qc::QC_3D>*>;
  fC->data_seq->push_back (meshData);
  fC->estimator_seq = new std::vector<qc::Estimator3d<REAL>*>;
  fC->estimator_seq->push_back (errorArray);

  fillScalarDataFunction3d(data_function, fC);
  data_function->name = const_cast<char*>(dataName);

  return reinterpret_cast<GENMESH3D *>(gmesh);
}



/* ---------------------------------------------------------------------------------------
* ***** ADDING NEW FUNCTIONS *******************************
* --------------------------------------------------------------------------------------- */

/* *********************************************************************
* add scalar data
* ********************************************************************* */
template <typename REAL>
    void addScalarData(GENMESH3D* gmesh, qc::ScalarArray<REAL, qc::QC_3D>* meshData, const char* dataName)
{
    // allocate memory for the data structures
    F_HDATA3D* data_function = reinterpret_cast<F_HDATA3D*> (mem_alloc(sizeof(F_HDATA3D)));

    // fill the data structure
    functionComponents<REAL>* fC = reinterpret_cast<functionComponents<REAL>*> (mem_alloc(sizeof(functionComponents<REAL>)));
    fC->estimator = NULL;
    fC->get_element_estimate = NULL;
    fC->get_vertex_estimate = NULL;

    fC->data_seq = new std::vector<qc::ScalarArray<REAL, qc::QC_3D>*>;
    fC->data_seq->push_back (meshData);
    fC->estimator_seq = new std::vector<qc::Estimator3d<REAL>*>;
    fC->estimator_seq->push_back (NULL);

    fC->data = new qc::ScalarArray<REAL, qc::QC_3D> (*meshData);
    fillScalarDataFunction3d(data_function, fC);
    data_function->name = const_cast<char*>(dataName);

    // now add the function via GRAPE
    GRAPE(gmesh,"add-function")(data_function);
}

/* *********************************************************************
* add timestep with scalar data
* ********************************************************************* */

void addTimestep (GENMESH3D* gmesh, qc::ScalarArray<double, qc::QC_3D>* meshData, const char* dataName, qc::Estimator3d<double>* errorArray)
{
  F_HDATA3D* fun;
  for (fun = reinterpret_cast<F_HDATA3D*> (gmesh->f_data); fun->last; fun = reinterpret_cast<F_HDATA3D*> (fun->last)) ; // Goto first function

  // Find correct function
  for (; fun && strcmp (fun->name, dataName); fun = reinterpret_cast<F_HDATA3D*> (fun->next)) ;

  if (strcmp (fun->name, dataName))
    throw aol::Exception ("Function not found", __FILE__, __LINE__);

  functionComponents<double>* fC = reinterpret_cast<functionComponents<double>*> (fun->function_data);

  fC->data_seq->push_back (meshData);
  fC->estimator_seq->push_back (errorArray);
}


/* *********************************************************************
* add scalar data with an qc::Estimator3d
* ********************************************************************* */

template <typename REAL>
    void addScalarDataWithEstimator(GENMESH3D* gmesh, qc::ScalarArray<REAL, qc::QC_3D>* meshData, const char* dataName,
                                qc::Estimator3d<REAL>* errorArray)
{
    // allocate memory for the data structures
    F_HDATA3D* data_function = reinterpret_cast<F_HDATA3D*> (mem_alloc(sizeof(F_HDATA3D)));

    // fill the data structure
    functionComponents<REAL>* fC = reinterpret_cast<functionComponents<REAL>*> (mem_alloc(sizeof(functionComponents<REAL>)));
    fC->estimator = new qc::Estimator3d<REAL> (*errorArray);
    fC->get_element_estimate = &calcElementError<REAL>;
    fC->get_vertex_estimate = &calcVertexError<REAL>;

    fC->data_seq = new std::vector<qc::ScalarArray<REAL, qc::QC_3D>*>;
    fC->data_seq->push_back (meshData);
    fC->estimator_seq = new std::vector<qc::Estimator3d<REAL>*>;
    fC->estimator_seq->push_back (errorArray);

    fC->data = new qc::ScalarArray<REAL, qc::QC_3D> (*meshData);
    fillScalarDataFunction3d(data_function, fC);
    data_function->name = const_cast<char*>(dataName);

    // now add the function via GRAPE
    GRAPE(gmesh,"add-function")(data_function);
}



/* *********************************************************************
* add vector data
* ********************************************************************* */
template <typename REAL>
    void addVectorData(GENMESH3D* gmesh, aol::MultiVector<REAL>* meshData, const char* dataName)
{
    // allocate memory for the data structures
    F_HDATA3D* data_function = reinterpret_cast<F_HDATA3D*> (mem_alloc(sizeof(F_HDATA3D)));

    // fill the data structure
    fillVectorDataFunction3d(data_function, meshData );
    data_function->name = const_cast<char*>(dataName);

    // now add the function via GRAPE
    GRAPE(gmesh,"add-function")(data_function);
}


void addMethodsAndProjects3d ()
{
  g_project_add(const_cast<char*>("uif-m3"));

  make_class ();
  GRAPE (GenMesh3d, "add-method") ("get-times", genmesh3d_get_times);
}

void initStartGrape(GENMESH3D* gmesh, const char* dataName)
{
  addMethodsAndProjects3d ();

  MANAGER*  mgr;
  mgr = reinterpret_cast<MANAGER*> (GRAPE(Manager,"get-stdmgr")());
  GRAPE(gmesh,"select-function")("default",dataName);

  GRAPE(mgr,"handle")(gmesh);
}

void initStartGrape(GENMESH3D* gmesh)
{
  addMethodsAndProjects3d ();

  MANAGER*  mgr;
  mgr = reinterpret_cast<MANAGER*> (GRAPE(Manager,"get-stdmgr")());

  GRAPE(mgr,"handle")(gmesh);
}

void initStartGrapeTime(GENMESH3D* gmesh)
{
  addMethodsAndProjects3d ();

  MANAGER*  mgr;
  mgr = reinterpret_cast<MANAGER*>(GRAPE(Manager,"get-stdmgr")());

  TIMESCENE* tsc = reinterpret_cast<TIMESCENE*> (GRAPE (TimeScene, "new-instance") ("sequence"));
  tsc->dynamic = reinterpret_cast<TREEOBJECT*> (gmesh);
  tsc->object = reinterpret_cast<TREEOBJECT*> (GRAPE (gmesh, "softcopy") (NULL));

  GRAPE(mgr,"handle")(tsc);
}

// Instantiation for ICC
template GENMESH3D* quocmesh_convert_to_gmesh3d (aol::MultiVector<float>*, const char*);
template GENMESH3D* quocmesh_convert_to_gmesh3d (aol::MultiVector<double>*, const char*);
template GENMESH3D* quocmesh_convert_to_gmesh3d (qc::ScalarArray<float, qc::QC_3D>*, const char*);
template GENMESH3D* quocmesh_convert_to_gmesh3d (qc::ScalarArray<double, qc::QC_3D>*, const char*);
template GENMESH3D* quocmesh_convert_to_gmesh3d (qc::ScalarArray<float, qc::QC_3D>*, qc::Estimator3d<float>*, const char*);
template GENMESH3D* quocmesh_convert_to_gmesh3d (qc::ScalarArray<double, qc::QC_3D>*, qc::Estimator3d<double>*, const char*);
template void addScalarData (GENMESH3D*, qc::ScalarArray<float, qc::QC_3D>*, const char*);
template void addScalarData (GENMESH3D*, qc::ScalarArray<double, qc::QC_3D>*, const char*);
template void addVectorData (GENMESH3D*, aol::MultiVector<float>*, const char*);
template void addVectorData (GENMESH3D*, aol::MultiVector<double>*, const char*);
template void addScalarDataWithEstimator (GENMESH3D*, qc::ScalarArray<float, qc::QC_3D>*, const char*, qc::Estimator3d<float>*);
template void addScalarDataWithEstimator (GENMESH3D*, qc::ScalarArray<double, qc::QC_3D>*, const char*, qc::Estimator3d<double>*);

#endif
