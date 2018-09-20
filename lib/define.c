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
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU //
// General Public License for more details.                          //
//                                                                   //
// You should have received a copy of the GNU General Public License //
// along with CUTE.  If not, see <http://www.gnu.org/licenses/>.     //
//                                                                   //
///////////////////////////////////////////////////////////////////////

/*********************************************************************/
//               Common global variables and macros                  //
/*********************************************************************/
#include <stdio.h>
#include <math.h>
#include "define.h"

//////////////////////////////////////

////////// Internal variables //////////
///

enum CORR_TYPE corr_type;
enum MULTI_TYPE multi_type;
Bin bin_main;
Bin bin_aux;

//Boxing variables
size_t *n_side,n_boxes,dim_box,dim_weight,dim_pos,n_ells;
histo_t *l_box;
///
////////////////////////////////////////
