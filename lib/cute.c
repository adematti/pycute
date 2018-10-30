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

/*********************************************************************/
//                               Main                                //
/*********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#include "define.h"
#include "common.h"

static Catalog cats[MAX_CATS] = {{.n_obj=0},{.n_obj=0},{.n_obj=0},{.n_obj=0}};
static Pole poles[MAX_POLES] = {{.n_ells=0},{.n_ells=0}};
static Mesh meshs[MAX_CATS];
static size_t n_cats = MAX_CATS;

void print_num_threads()
{
	//Calculate number of threads
	size_t n_threads=0;
#pragma omp parallel
	{
#pragma omp atomic
		n_threads++;
	}
	printf(" - using %zu threads\n",n_threads);
}

void set_num_threads(size_t n_threads)
{
	omp_set_num_threads(n_threads);
#ifdef _VERBOSE
	print_num_threads();
#endif //_VERBOSE
}

void print_bin(char* mode)
{
	printf("*** Binning\n");
	Bin bin=bin_main;
	if (!strcmp(mode,"aux")) bin=bin_aux;
	else if (!strcmp(mode,"radial")) bin=bin_radial;
	printf(" - range: %.3f < %s < %.3f\n",bin.min,mode,bin.max);
	if (bin.type==BIN_LIN) printf(" - #bins(lin): %zu\n",bin.n_bin);
	else if (bin.type==BIN_LOG) printf(" - #bins(log): %zu\n",bin.n_bin);
	else if (bin.type==BIN_CUSTOM) printf(" - #bins(custom): %zu\n",bin.n_bin);
	printf(" - resolution: %.3f\n",bin.step);
}

void set_bin(char* mode,histo_t* edges,size_t n_bin,char* type)
{
	Bin bin;

	bin.min=edges[0];
	bin.max=edges[n_bin];
	bin.n_bin=n_bin;
	bin.log10min=0.;
	bin.log10max=0.;
	bin.edges=NULL;

	if (!strcmp(type,"lin")) {
		bin.type=BIN_LIN;
		bin.step=(bin.max-bin.min)/n_bin;
	}
	else if (!strcmp(type,"log")) {
		bin.type=BIN_LOG;
		bin.log10min=log10(bin.min);
		bin.log10max=log10(bin.max);
		bin.step=(bin.log10max-bin.log10min)/n_bin;
	}
	else if (!strcmp(type,"custom")) {
		bin.type=BIN_CUSTOM;
		bin.edges=edges;
		bin.step=(bin.max-bin.min)/n_bin;
	}
	else {
		bin.type=BIN_LIN;
		bin.step=(bin.max-bin.min)/n_bin;
		fprintf(stderr," - invalid main-binning type. Choices:lin, log or custom.\n");
		fprintf(stderr," - I choose linear binning.\n");
	}
	if (!strcmp(mode,"main")) bin_main=bin;
	else if (!strcmp(mode,"aux")) bin_aux=bin;
	else if (!strcmp(mode,"radial")) bin_radial=bin;
	else {
		bin_main=bin;
		fprintf(stderr," - invalid binning. Choices: main, aux, radial.\n");
		fprintf(stderr," - I choose main.\n");
	}
#ifdef _VERBOSE
	print_bin(mode);
#endif //_VERBOSE
}

void print_corr_type()
{
	if (corr_type==CORR_SMU) printf(" - corr-type: s-mu\n");
	else if (corr_type==CORR_ANGULAR) printf(" - corr-type: angular\n");
	else if (corr_type==CORR_SCOS) printf(" - corr-type: s-cos\n");
}

void set_corr_type(char* type)
{
	//printf("*** Correlation type:\n");
	if (!strcmp(type,"s-mu")) corr_type=CORR_SMU;
	else if (!strcmp(type,"angular")) corr_type=CORR_ANGULAR;
	else if (!strcmp(type,"s-cos")) corr_type=CORR_SCOS;
	else {
		corr_type=CORR_SMU;
		fprintf(stderr," - invalid correlation type. Choices: s-mu, angular or s-cos.\n");
		fprintf(stderr," - I choose s-mu.\n");
	}
#ifdef _VERBOSE
	print_corr_type();
#endif //_VERBOSE
}

void print_pole(Pole pole)
{
	printf("*** Multipoles\n");
	size_t ill;
	if (pole.type==MULTI_ALL) {
		printf(" - multi-type: all\n");
		printf(" - ells:");
		for (ill=0;ill<pole.n_ells;ill++) printf(" %zu",ill);
		printf("\n");
	}
	if (pole.type==MULTI_EVEN) {
		printf(" - multi-type: even\n");
		printf(" - ells:");
		for (ill=0;ill<pole.n_ells;ill++) printf(" %zu",2*ill);
		printf("\n");
	}
	if (pole.type==MULTI_ODD) {
		printf(" - multi-type: odd\n");
		printf(" - ells:");
		for (ill=0;ill<pole.n_ells;ill++) printf(" %zu",2*ill+1);
		printf("\n");
	}
}

void set_pole(size_t num,MULTI_TYPE type,size_t n_ells)
{
	Pole pole;
	pole.type = type;
	pole.n_ells = n_ells;
	poles[num-1] = pole;
#ifdef _VERBOSE
	print_pole(pole);
#endif //_VERBOSE	
}

void clear_poles()
{
	size_t ipole;
	for (ipole=0;ipole<MAX_POLES;ipole++) poles[ipole].n_ells = 0;
}


void print_los_type()
{
	if (los_type==LOS_MIDPOINT) printf(" - los-type: midpoint\n");
	else if (los_type==LOS_ENDPOINT) printf(" - los-type: endpoint\n");
}

void set_los_type(char* type)
{
	if (!strcmp(type,"midpoint")) los_type=LOS_MIDPOINT;
	else if (!strcmp(type,"endpoint")) los_type=LOS_ENDPOINT;
	else {
		los_type=LOS_MIDPOINT;
		fprintf(stderr," - invalid los type. Choices: middle or endpoint.\n");
		fprintf(stderr," - I choose middle.\n");
	}
#ifdef _VERBOSE
	print_los_type();
#endif //_VERBOSE
}

void print_catalogs()
{
	printf(" - #spatial dim: %zu\n",dim_box);
	printf(" - #weight dim: %zu\n",dim_weight);
	size_t icat;
	for (icat=0;icat<MAX_CATS;icat++) {
		if (cats[icat].n_obj>0) {
			printf(" - #objects in cat%zu: %zu\n",icat+1,cats[icat].n_obj);
		}
	}
}

void print_normalize(_Bool normalize)
{
	if (normalize) printf(" - normalize: yes\n");
	else printf(" - normalize: no\n");
}

void set_catalog(size_t num,histo_t *p,histo_t *w,size_t n,size_t dim_b,size_t dim_w)
{	
	Catalog cat;
	cat.pos = p;
	cat.weight = w;
	cat.n_obj = n;
	cats[num-1] = cat;
	
	dim_box=dim_b;
	dim_pos=dim_box;
	dim_weight=dim_w;

}

void set_num_catalogs()
{
	size_t icat;
	for (icat=0;icat<MAX_CATS;icat++) {
		if (cats[icat].n_obj==0) break;
	}
	n_cats = icat;
}

void clear_catalogs()
{
	size_t icat;
	for (icat=0;icat<MAX_CATS;icat++) cats[icat].n_obj = 0;
}

void run_2pcf_main(histo_t* meanmain,histo_t* count,char* corr_type,size_t n_threads)
{
	timer(0);
	set_num_catalogs();
#ifdef _VERBOSE
	printf("*** 2-point correlation function (main)\n");
	print_catalogs();
#endif //_VERBOSE
	set_corr_type(corr_type);
	set_num_threads(n_threads);
	set_meshs(cats,meshs,n_cats);
#ifdef _VERBOSE
	timer(1);
#endif //_VERBOSE
	if (n_cats==1) auto_2pcf_main(meshs[0],meanmain,count);
	else cross_2pcf_main(meshs[0],meshs[1],meanmain,count);
#ifdef _DEBUG
	if (n_cats==1) write_histo(bin_main.n_bin,count,"debug_auto.dat");
	else write_histo(bin_main.n_bin,count,"debug_12.dat");
#endif //_DEBUG
#ifdef _VERBOSE
	printf("*** Cleaning up\n");
	timer(1);
	printf("\n");
#endif //_VERBOSE
	free_meshs(meshs,n_cats);
}

void run_2pcf_main_aux(histo_t* meanmain,histo_t* meanaux,histo_t* count,char* corr_type,char *los_type,size_t n_threads)
{
	timer(0);
	set_num_catalogs();
#ifdef _VERBOSE
	printf("*** 2-point correlation function (main,aux)\n");
	print_catalogs();
#endif //_VERBOSE
	set_corr_type(corr_type);
	set_los_type(los_type);
	set_num_threads(n_threads);
	set_meshs(cats,meshs,n_cats);
#ifdef _VERBOSE
	timer(1);
#endif //_VERBOSE
	if (n_cats==1) auto_2pcf_main_aux(meshs[0],meanmain,meanaux,count);
	else cross_2pcf_main_aux(meshs[0],meshs[1],meanmain,meanaux,count);
#ifdef _DEBUG
	if (n_cats==1) write_histo(bin_main.n_bin*bin_aux.n_bin,count,"debug_auto.dat");
	else write_histo(bin_main.n_bin*bin_aux.n_bin,count,"debug_12.dat");
#endif //_DEBUG
#ifdef _VERBOSE
	printf("*** Cleaning up\n");
	timer(1);
	printf("\n");
#endif //_VERBOSE
	free_meshs(meshs,n_cats);
}

void run_2pcf_multi(histo_t *meanmain,histo_t *count,char *los_type,size_t n_threads)
{
	timer(0);
	set_num_catalogs();
#ifdef _VERBOSE
	printf("*** 2-point correlation function (multi)\n");
	print_catalogs();
#endif //_VERBOSE
	set_los_type(los_type);
	set_num_threads(n_threads);
	set_meshs(cats,meshs,n_cats);
#ifdef _VERBOSE
	timer(1);
#endif //_VERBOSE
	if (n_cats==1) auto_2pcf_multi(meshs[0],meanmain,count,poles[0]);
	else cross_2pcf_multi(meshs[0],meshs[1],meanmain,count,poles[0]);	
#ifdef _DEBUG
	if (n_cats==1) write_histo(bin_main.n_bin*n_ells,count,"debug_auto.dat");
	else write_histo(bin_main.n_bin*n_ells,count,"debug_12.dat");
#endif //_DEBUG
#ifdef _VERBOSE
	printf("*** Cleaning up\n");
	timer(1);
	printf("\n");
#endif //_VERBOSE
	free_meshs(meshs,n_cats);
}

void run_3pcf_multi(histo_t *count,char *los_type,size_t n_threads)
{
	timer(0);
	set_num_catalogs();
#ifdef _VERBOSE
	print_catalogs();
	printf("*** 3-point correlation function multipoles\n");
#endif //_VERBOSE
	set_los_type(los_type);
	set_num_threads(n_threads);
	set_meshs(cats,meshs,n_cats);
#ifdef _VERBOSE
	timer(1);
#endif //_VERBOSE
	cross_3pcf_multi(meshs,n_cats,count,poles);
#ifdef _VERBOSE
	printf("*** Cleaning up\n");
	timer(1);
	printf("\n");
#endif //_VERBOSE
	free_meshs(meshs,n_cats);
}

void run_3pcf_multi_double_los(histo_t *count,char *los_type,size_t n_threads)
{
	timer(0);
	set_num_catalogs();
#ifdef _VERBOSE
	print_catalogs();
	printf("*** 3-point correlation function multipoles with double los for cat2\n");
#endif //_VERBOSE
	set_los_type(los_type);
	set_num_threads(n_threads);
	set_meshs(cats,meshs,n_cats);
#ifdef _VERBOSE
	timer(1);
#endif //_VERBOSE
	cross_3pcf_multi_double_los(meshs,n_cats,count,poles);
#ifdef _VERBOSE
	printf("*** Cleaning up\n");
	timer(1);
	printf("\n");
#endif //_VERBOSE
	free_meshs(meshs,n_cats);
}

void run_2pcf_multi_radial(histo_t *count,_Bool normalize,char *los_type,size_t n_threads)
{
	timer(0);
	set_num_catalogs();
#ifdef _VERBOSE
	print_catalogs();
	printf("*** 2-point radial correlation function multipoles\n");
#endif //_VERBOSE
	set_los_type(los_type);
	set_num_threads(n_threads);
	set_meshs(cats,meshs,n_cats);
#ifdef _VERBOSE
	timer(1);
#endif //_VERBOSE
	cross_2pcf_multi_radial(meshs[0],meshs[1],count,poles[0],normalize);
#ifdef _VERBOSE
	printf("*** Cleaning up\n");
	timer(1);
	printf("\n");
#endif //_VERBOSE
	free_meshs(meshs,n_cats);
}

void run_4pcf_multi_radial(histo_t *count,_Bool normalize,char *los_type,size_t n_threads)
{
	timer(0);
	set_num_catalogs();
#ifdef _VERBOSE
	print_catalogs();
	printf("*** 2-point radial correlation function multipoles\n");
#endif //_VERBOSE
	set_los_type(los_type);
	set_num_threads(n_threads);
	set_meshs(cats,meshs,n_cats);
#ifdef _VERBOSE
	timer(1);
#endif //_VERBOSE
	cross_4pcf_multi_radial(meshs,count,poles,normalize);
#ifdef _VERBOSE
	printf("*** Cleaning up\n");
	timer(1);
	printf("\n");
#endif //_VERBOSE
	free_meshs(meshs,n_cats);
}


