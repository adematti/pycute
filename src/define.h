#ifndef _CUTE_DEFINE_
#define _CUTE_DEFINE_

#define MAX(a, b) (((a) > (b)) ? (a) : (b)) //Maximum of two numbers
#define MIN(a, b) (((a) < (b)) ? (a) : (b)) //Minimum of two numbers
#define CLAMP(x, low, high)  (((x) > (high)) ? (high) : (((x) < (low)) ? (low) : (x))) //min(max(a,low),high)
//#define ABS(a)   (((a) < 0) ? -(a) : (a)) //Absolute value
#define MAX_ELLS 9
#define MAX_CATS 4
#define MAX_POLES 2
#define MAX_LOS 2

#ifdef _FLOAT32
typedef float histo_t;
#else
typedef double histo_t;
#endif //_FLOAT32

typedef enum {QUIET, INFO, DEBUG} VERBOSITY;
typedef enum {BIN_LIN, BIN_LOG, BIN_CUSTOM, BIN_BIN} BIN_TYPE;
typedef enum {CORR_SMU, CORR_ANGULAR, CORR_SCOS} CORR_TYPE;
typedef enum {LOS_MIDPOINT, LOS_ENDPOINT, LOS_FIRSTPOINT, LOS_CUSTOM} LOS_TYPE;
typedef enum {WEIGHT_PROD, WEIGHT_PRODSUM} WEIGHT_TYPE;
typedef enum {MULTI_ALL, MULTI_EVEN, MULTI_ODD} MULTI_TYPE;

//Box for 2PCFs
typedef struct {
  size_t n_obj;
  histo_t *pos;
  histo_t *weight;
  size_t *bin;
  size_t *index;
  size_t *visit_min;
  size_t *visit_max;
} Box; //cell

typedef struct {
  size_t n_boxes;
  Box* boxes;
} Mesh;

//Catalog
typedef struct {
  size_t n_obj;
  histo_t *pos;
  histo_t *weight;
  size_t *bin;
  histo_t *box;
} Catalog; //Catalog

//Bin
typedef struct {
  size_t n_bin;
  histo_t min,max,step,log10min,log10max;
  histo_t* edges;
  BIN_TYPE type;
} Bin;

typedef struct {
  size_t n_ells;
  MULTI_TYPE type;
} Pole;

typedef struct {
  histo_t *los;
  LOS_TYPE type;
  size_t n;
} LOS;

VERBOSITY verbose;
CORR_TYPE corr_type;
WEIGHT_TYPE weight_type;
Bin bin_main,bin_aux,bin_bin;
size_t dim_box,dim_weight,dim_pos;

#endif //_CUTE_DEFINE_
