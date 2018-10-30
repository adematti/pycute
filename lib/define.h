///////////////////////////////////////////////////////////////////////
//                                                                   //
//   Copyright 2012 David Alonso                                     //
//                                                                   //
//                                                                   //
// This file is part of CUTE.                                        //
//                                                                   //
// CUTE is free software: you can redistribute it and/or modify it   //
// under the terms of the GNU General Public License as published by //
// the Free Software Foundation, either version 3 of the License, or //
// (at your option) any later version.                               //
//                                                                   //
// CUTE is distributed in the hope that it will be useful, but       //
// WITHOUT ANY WARRANTY; without even the implied warranty of        //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU  //
// General Public License for more details.                          //
//                                                                   //
// You should have received a copy of the GNU General Public License //
// along with CUTE.  If not, see <http://www.gnu.org/licenses/>.     //
//                                                                   //
///////////////////////////////////////////////////////////////////////

#ifndef _CUTE_DEFINE_
#define _CUTE_DEFINE_

#define MAX(a, b) (((a) > (b)) ? (a) : (b)) //Maximum of two numbers
#define MIN(a, b) (((a) < (b)) ? (a) : (b)) //Minimum of two numbers
#define CLAMP(x, low, high)  (((x) > (high)) ? (high) : (((x) < (low)) ? (low) : (x))) //min(max(a,low),high)
//#define ABS(a)   (((a) < 0) ? -(a) : (a)) //Absolute value
#define MAX_ELLS 9
#define MAX_CATS 4
#define MAX_POLES 2

#ifdef _FLOAT32
typedef float histo_t;
#else
typedef double histo_t;
#endif //_FLOAT32

typedef enum {BIN_LIN, BIN_LOG, BIN_CUSTOM} BIN_TYPE;
typedef enum {CORR_SMU, CORR_ANGULAR, CORR_SCOS} CORR_TYPE;
typedef enum {LOS_MIDPOINT, LOS_ENDPOINT} LOS_TYPE;
typedef enum {MULTI_ALL, MULTI_EVEN, MULTI_ODD} MULTI_TYPE;

//Box for 2PCFs
typedef struct {
	size_t n_obj;
	histo_t *pos;
	histo_t *weight;
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
	histo_t *box;
} Catalog; //Catalog

//Bin
typedef struct {
	size_t n_bin;
	histo_t min,max,step,log10min,log10max;
	histo_t* edges;
	BIN_TYPE type;  
} Bin; //Bin

typedef struct {
	size_t n_ells;
	MULTI_TYPE type;  
} Pole; //Bin

CORR_TYPE corr_type;
LOS_TYPE los_type;
Bin bin_main,bin_aux,bin_radial;
size_t dim_box,dim_weight,dim_pos;

#endif //_CUTE_DEFINE_
