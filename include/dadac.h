#pragma once
#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <string.h>
 

#define T double
#define DATA_DIMS 0 

#ifdef USE_FLOAT32
	#define FLOAT_TYPE float
#else
	#define FLOAT_TYPE double 
#endif

#ifdef USE_INT32
	#define MY_SIZE_MAX UINT32_MAX
	#define idx_t uint32_t
#else
	#define MY_SIZE_MAX UINT64_MAX
	#define idx_t uint64_t
#endif


#define DTHR 23.92812698
#define PI_F 3.1415926f
#define ARRAY_INCREMENT 500
#define DA_DTYPE idx_t
#define NOBORDER MY_SIZE_MAX

/**********************************
 * DATA STRUCTURES FOR CLUSTERING *
 **********************************/
struct Node {
  idx_t data;
  struct Node *next;
};

struct LinkedList {
  idx_t count;
  struct Node *head;
};

struct lu_dynamicArray {
  idx_t *data;
  idx_t size;
  idx_t count;
};

struct Datapoint_info {
  FLOAT_TYPE g;
  idx_t array_idx;
  FLOAT_TYPE log_rho;
  FLOAT_TYPE log_rho_c;
  FLOAT_TYPE log_rho_err;
  idx_t kstar;
  int is_center;
  int cluster_idx;
};

struct border_t {
  FLOAT_TYPE density;
  FLOAT_TYPE error;
  idx_t idx;
};

struct SparseBorder_t {
  idx_t i;
  idx_t j;
  idx_t idx;
  FLOAT_TYPE density;
  FLOAT_TYPE error;
};

struct AdjList_t {
  idx_t count;
  idx_t size;
  struct SparseBorder_t* data;
};

struct aClusters {
  struct lu_dynamicArray centers;
  FLOAT_TYPE **border_density;
  FLOAT_TYPE **border_err;
  idx_t **border_idx;
  FLOAT_TYPE *__border_density_data;
  FLOAT_TYPE *__border_err_data;
  idx_t *__border_idx_data;
  idx_t n;
};

struct Clusters {
  int UseSparseBorders;
  struct AdjList_t *SparseBorders;
  struct lu_dynamicArray centers;
  struct border_t **borders;
  struct border_t *__borders_data;
  idx_t n;
};

struct merge_t {
  idx_t source;
  idx_t target;
  FLOAT_TYPE density;
};

typedef struct Datapoint_info Datapoint_info;
typedef struct lu_dynamicArray lu_dynamicArray;
typedef struct Clusters Clusters;
typedef struct Node Node;
typedef struct LinkedList LinkedList;
typedef struct merge_t merge_t;
typedef struct border_t border_t;
typedef struct SparseBorder_t SparseBorder_t;
typedef struct AdjList_t AdjList_t; 

void DynamicArray_allocate(lu_dynamicArray *a);
void DynamicArray_pushBack(lu_dynamicArray *a, idx_t p);
//void Clusters_allocate(Clusters *c);
void Clusters_allocate(Clusters *c, int s);
void Clusters_free(Clusters *c);

int cmp(const void *a, const void *b);
void computeRho(Datapoint_info *particles, const FLOAT_TYPE d,
                const idx_t points);
int cmpPP(const void *p1, const void *p2);
void computeCorrection(Datapoint_info *particles, int* mask, idx_t n, FLOAT_TYPE Z);

Clusters Heuristic1(Datapoint_info* dpInfo, int* mask, size_t nrows, size_t ncols);
//Clusters Heuristic1(Datapoint_info* dpInfo, int* mask, int nrows, int ncols);
void Heuristic2(Clusters* cluster, Datapoint_info* dpInfo, int* mask, size_t nrows, size_t ncols);
//void Heuristic2(Clusters* cluster, Datapoint_info* dpInfo, int* mask, size_t nrows, size_t ncols);
void Heuristic3(Clusters *cluster, Datapoint_info *particles, FLOAT_TYPE Z,int halo);
void freeDatapointArray(Datapoint_info* d, size_t n);

Datapoint_info* computeDensityFromImg(FLOAT_TYPE* vals, int* mask, int nrows, int ncols, int r);
void Delete_adjlist_element(Clusters * c, const idx_t list_idx, const idx_t el);
int merging_roles( FLOAT_TYPE dens1, FLOAT_TYPE dens1_err,
			  FLOAT_TYPE dens2, FLOAT_TYPE dens2_err,
			  FLOAT_TYPE dens_border, FLOAT_TYPE dens_border_err );
int is_a_merging( FLOAT_TYPE dens1, FLOAT_TYPE dens1_err,
			 FLOAT_TYPE dens2, FLOAT_TYPE dens2_err,
			 FLOAT_TYPE dens_border, FLOAT_TYPE dens_border_err,
			 FLOAT_TYPE Z);
