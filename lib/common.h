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

#ifndef _CUTE_COMMON_
#define _CUTE_COMMON_

//General-purpose functions

histo_t wrap_phi(histo_t phi);

size_t get_bin_index(histo_t x,Bin bin);

size_t ravel_index(size_t ind[],size_t dim[],size_t n);

size_t* unravel_index(size_t ind,size_t dim[],size_t n,size_t index[]);

size_t get_dichotomy_index(histo_t x,histo_t *tab,size_t min,size_t max);

void legendre_all(histo_t mu,histo_t mu2,histo_t leg[]);

void legendre_even(histo_t mu2,histo_t leg[]);

void legendre_odd(histo_t mu,histo_t mu2,histo_t leg[]);

void timer(size_t i);

void error_mem_out(void);

void error_open_file(char *fname);

void free_catalog(Catalog cat);


//Boxes

void ind2pos(size_t index[],histo_t pos[]);

void set_meshs(Catalog *cats,Mesh *meshs,size_t n_cats);

void free_meshs(Mesh *meshs,size_t n_meshs);

//Correlators

void cross_2pcf_main(Mesh mesh1,Mesh mesh2,histo_t meanmain[],histo_t count[]);

void auto_2pcf_main(Mesh mesh1,histo_t meanmain[],histo_t count[]);

void cross_2pcf_main_aux(Mesh mesh1,Mesh mesh2,histo_t meanmain[],histo_t meanaux[],histo_t count[]);

void auto_2pcf_main_aux(Mesh mesh1,histo_t meanmain[],histo_t meanaux[],histo_t count[]);

void cross_2pcf_multi(Mesh mesh1,Mesh mesh2,histo_t meanmain[],histo_t count[],Pole pole);

void auto_2pcf_multi(Mesh mesh1,histo_t meanmain[],histo_t count[],Pole pole);

void cross_3pcf_multi(Mesh* meshs,size_t n_meshs,histo_t count[],Pole *poles);

void cross_3pcf_multi_double_los(Mesh* meshs,size_t n_meshs,histo_t count[],Pole *poles);

void cross_2pcf_multi_radial(Mesh mesh1, Mesh mesh2, histo_t count[],Pole pole,_Bool normalize);

void cross_4pcf_multi_radial(Mesh *meshs,histo_t count[],Pole *poles,_Bool normalize);

#ifdef _DEBUG
//Debug files output
void write_mesh(Mehs mesh,char *fn);

void write_catalog(Catalog cat,char *fn);

void write_histo(size_t num_histo,histo_t* hist,char *fn);
#endif //_DEBUG

#endif //_CUTE_COMMON_
