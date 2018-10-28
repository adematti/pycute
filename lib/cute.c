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

static Catalog catalog[MAX_CATS] = {{.n_obj=0},{.n_obj=0},{.n_obj=0},{.n_obj=0}};

void print_num_threads()
{
	//Calculate number of threads
	size_t num_threads=0;
#pragma omp parallel
	{
#pragma omp atomic
		num_threads++;
	}
	printf(" - Using %zu threads\n",num_threads);
}


void set_num_threads(size_t num_threads)
{
	omp_set_num_threads(num_threads);
#ifdef _VERBOSE
	print_num_threads();
#endif //_VERBOSE
}

void print_bin(char* mode)
{
	Bin bin=bin_main;
	if (!strcmp(mode,"aux")) bin=bin_aux;
	else if (!strcmp(mode,"radial")) bin=bin_radial;
	printf(" - Range: %.3f < %s < %.3f\n",bin.min,mode,bin.max);
	if (bin.type==BIN_LIN) printf(" - #bins(lin): %zu\n",bin.n_bin);
	else if (bin.type==BIN_LOG) printf(" - #bins(log): %zu\n",bin.n_bin);
	else if (bin.type==BIN_CUSTOM) printf(" - #bins(custom): %zu\n",bin.n_bin);
	printf(" - Resolution: %.3f\n",bin.step);
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

	printf("*** Binning\n");
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
		fprintf(stderr," - Invalid main-binning type. Choices:lin, log or custom.\n");
		fprintf(stderr," - I choose linear binning.\n");
	}
	if (!strcmp(mode,"main")) bin_main=bin;
	else if (!strcmp(mode,"aux")) bin_aux=bin;
	else if (!strcmp(mode,"radial")) bin_radial=bin;
	else {
		bin_main=bin;
		fprintf(stderr," - Invalid binning. Choices: main, aux, radial.\n");
		fprintf(stderr," - I choose main.\n");
	}
#ifdef _VERBOSE
	print_bin(mode);
	printf("\n");
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
		fprintf(stderr," - Invalid correlation type. Choices: s-mu, angular or s-cos.\n");
		fprintf(stderr," - I choose s-mu.\n");
	}
#ifdef _VERBOSE
	print_corr_type();
#endif //_VERBOSE
}

void print_multi_type()
{
	size_t l;
	if (multi_type==MULTI_ALL) {
		printf(" - multi-type: all\n");
		printf(" - ells:");
		for (l=0;l<n_ells;l++) printf(" %zu",l);
		printf("\n");
	}
	if (multi_type==MULTI_EVEN) {
		printf(" - multi-type: even\n");
		printf(" - ells:");
		for (l=0;l<n_ells;l++) printf(" %zu",2*l);
		printf("\n");
	}
	if (multi_type==MULTI_ODD) {
		printf(" - multi-type: odd\n");
		printf(" - ells:");
		for (l=0;l<n_ells;l++) printf(" %zu",2*l+1);
		printf("\n");
	}
}

void set_multi_type(char* type,size_t num_ells)
{
	if (!strcmp(type,"all")) multi_type=MULTI_ALL;
	else if (!strcmp(type,"even")) multi_type=MULTI_EVEN;
	else if (!strcmp(type,"odd")) multi_type=MULTI_ODD;
	else {
		multi_type=MULTI_ALL;
		fprintf(stderr," - Invalid multipole type. Choices: all, even or odd.\n");
		fprintf(stderr," - I choose all.\n");
	}
	corr_type=CORR_SMU;
	n_ells = MIN(num_ells,MAX_ELLS);
#ifdef _VERBOSE
	print_multi_type();
#endif //_VERBOSE
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
		fprintf(stderr," - Invalid los type. Choices: middle or endpoint.\n");
		fprintf(stderr," - I choose middle.\n");
	}
#ifdef _VERBOSE
	print_los_type();
#endif //_VERBOSE
}

void print_cross(size_t ind1,size_t ind2)
{
	printf(" - Correlation: catalog %zu with catalog %zu\n",ind1,ind2);
}

void set_cross(size_t *ind1,size_t *ind2)
{
	if ((*ind1>=MAX_CATS) || (*ind2>=MAX_CATS)) {
		fprintf(stderr," - Invalid number of catalogs. Choices: 1 or 2.\n");
		fprintf(stderr," - I choose 1 - 1.\n");
		(*ind1) = 1;
		(*ind2) = 1;
	}
	if (catalog[*ind1-1].n_obj==0) fprintf(stderr," - Catalog %zu is empty.\n",*ind1);
	if (catalog[*ind2-1].n_obj==0) fprintf(stderr," - Catalog %zu is empty.\n",*ind2);
#ifdef _VERBOSE
	print_cross(*ind1,*ind2);
#endif //_VERBOSE
}

void print_catalogs()
{
	printf("*** Catalogs\n");
	printf(" - #spatial dim: %zu\n",dim_box);
	printf(" - #weight dim: %zu\n",dim_weight);
	size_t icat;
	for (icat=0;icat<MAX_CATS;icat++) {
		if (catalog[icat].n_obj>0) {
			printf(" - #objects in cat%zu: %zu\n",icat+1,catalog[icat].n_obj);
		}
	}
	printf("\n");
}

void print_normalize(_Bool normalize)
{
	if (normalize) printf(" - Normalize: yes\n");
	else printf(" - Normalize: no\n");
}

void set_catalog(size_t num,histo_t *p,histo_t *w,size_t n,size_t dim_b,size_t dim_w)
{	
	catalog[num-1].pos=p;
	catalog[num-1].weight=w;
	catalog[num-1].n_obj=n;
	
	dim_box=dim_b;
	dim_pos=dim_box;
	dim_weight=dim_w;

}

void clear_catalogs()
{
	size_t icat;
	for (icat=0;icat<MAX_CATS;icat++) catalog[icat].n_obj = 0;
}

void print_correlation(size_t ind1,size_t ind2)
{
	printf("*** Correlating\n");
	if (ind1!=ind2) printf(" - Cross-correlating\n");
	else printf(" - Auto-correlating\n");
}


void run_2pcf_main_aux(size_t ind1,size_t ind2,histo_t* meanmain,histo_t* meanaux,histo_t* count,char* corr_type,char *los_type,size_t num_threads)
{

	Box *boxes1=NULL;
	Box *boxes2=NULL;
	size_t *indices1,*indices2;
	size_t nfull1,nfull2;

#ifdef _VERBOSE
	print_catalogs();
#endif //_VERBOSE
	printf("*** 2-point correlation function (main,aux)\n");
	set_corr_type(corr_type);
	set_los_type(los_type);
	set_cross(&ind1,&ind2);
	set_num_threads(num_threads);
#ifdef _DEBUG
	write_catalog(catalog[0],"debug_cat1.dat");
	write_catalog(catalog[1],"debug_cat2.dat");
#endif //_DEBUG
	timer(0);
	if ((ind1==1)||(ind2==1)) set_catalog_box(&catalog[0]);
	if ((ind1==2)||(ind2==2)) set_catalog_box(&catalog[1]);
	if ((ind1==1)&&(ind2==1)) init_params(catalog[0],catalog[0]);
	else if ((ind1==2)&&(ind2==2)) init_params(catalog[1],catalog[1]);
	else init_params(catalog[0],catalog[1]);
	
	if ((ind1==1)||(ind2==1)) boxes1=make_boxes_from_catalog(catalog[0],&indices1,&nfull1);
	if ((ind1==2)||(ind2==2)) boxes2=make_boxes_from_catalog(catalog[1],&indices2,&nfull2);
#ifdef _VERBOSE
	timer(2);
	printf("\n");
#endif //_VERBOSE
#ifdef _DEBUG
	if ((ind1==1)||(ind2==1)) write_boxes(boxes1,"debug_box1.dat");
	if ((ind1==2)||(ind2==2)) write_boxes(boxes2,"debug_box2.dat");
#endif //_DEBUG
#ifdef _VERBOSE
	print_correlation(ind1,ind2);
#endif //_VERBOSE
	if ((ind1==1)&&(ind2==1)) auto_2pcf_main_aux(nfull1,indices1,boxes1,meanmain,meanaux,count);
	else if ((ind1==2)&&(ind2==2)) auto_2pcf_main_aux(nfull2,indices2,boxes2,meanmain,meanaux,count);
	else cross_2pcf_main_aux(nfull1,nfull2,indices1,indices2,boxes1,boxes2,meanmain,meanaux,count);
#ifdef _VERBOSE
	timer(1);
	printf("\n");
#endif //_VERBOSE	
#ifdef _DEBUG
	if ((ind1==1)&&(ind2==1)) write_histo(bin_main.n_bin*bin_aux.n_bin,count,"debug_11.dat");
	else if ((ind1==2)&&(ind2==2)) write_histo(bin_main.n_bin*bin_aux.n_bin,count,"debug_22.dat");
	else write_histo(bin_main.n_bin*bin_aux.n_bin,count,"debug_12.dat");
#endif //_DEBUG

	printf("*** Cleaning up\n");
	if ((ind1==1)||(ind2==1)) {
		free_boxes(boxes1);
		free(indices1);
	}
	if ((ind1==2)||(ind2==2)) {
		free_boxes(boxes2);
		free(indices2);
	}
	free_params();
}




void run_2pcf_main(size_t ind1,size_t ind2,histo_t* meanmain,histo_t* count,char* corr_type,size_t num_threads)
{

	Box *boxes1=NULL;
	Box *boxes2=NULL;
	size_t *indices1,*indices2;
	size_t nfull1,nfull2;

#ifdef _VERBOSE
	print_catalogs();
#endif //_VERBOSE
	printf("*** 2-point correlation function (main)\n");
	set_corr_type(corr_type);
	set_cross(&ind1,&ind2);
	set_num_threads(num_threads);
#ifdef _DEBUG
	write_catalog(catalog[0],"debug_cat1.dat");
	write_catalog(catalog[1],"debug_cat2.dat");
#endif //_DEBUG
	timer(0);
	if ((ind1==1)||(ind2==1)) set_catalog_box(&catalog[0]);
	if ((ind1==2)||(ind2==2)) set_catalog_box(&catalog[1]);
	if ((ind1==1)&&(ind2==1)) init_params(catalog[0],catalog[0]);
	else if ((ind1==2)&&(ind2==2)) init_params(catalog[1],catalog[1]);
	else init_params(catalog[0],catalog[1]);
	
	if ((ind1==1)||(ind2==1)) boxes1=make_boxes_from_catalog(catalog[0],&indices1,&nfull1);
	if ((ind1==2)||(ind2==2)) boxes2=make_boxes_from_catalog(catalog[1],&indices2,&nfull2);
#ifdef _VERBOSE
	timer(2);
	printf("\n");
#endif //_VERBOSE
#ifdef _DEBUG
	if ((ind1==1)||(ind2==1)) write_boxes(boxes1,"debug_box1.dat");
	if ((ind1==2)||(ind2==2)) write_boxes(boxes2,"debug_box2.dat");
#endif //_DEBUG
#ifdef _VERBOSE
	print_correlation(ind1,ind2);
#endif //_VERBOSE
	if ((ind1==1)&&(ind2==1)) auto_2pcf_main(nfull1,indices1,boxes1,meanmain,count);
	else if ((ind1==2)&&(ind2==2)) auto_2pcf_main(nfull2,indices2,boxes2,meanmain,count);
	else cross_2pcf_main(nfull1,nfull2,indices1,indices2,boxes1,boxes2,meanmain,count);
#ifdef _VERBOSE
	timer(1);
	printf("\n");
#endif //_VERBOSE
#ifdef _DEBUG
	if ((ind1==1)&&(ind2==1)) write_histo(bin_main.n_bin,count,"debug_11.dat");
	else if ((ind1==2)&&(ind2==2)) write_histo(bin_main.n_bin,count,"debug_22.dat");
	else write_histo(bin_main.n_bin,count,"debug_12.dat");
#endif //_DEBUG

	printf("*** Cleaning up\n");
	if ((ind1==1)||(ind2==1)) {
		free_boxes(boxes1);
		free(indices1);
	}
	if ((ind1==2)||(ind2==2)) {
		free_boxes(boxes2);
		free(indices2);
	}
	free_params();
}

void run_2pcf_multi(size_t ind1,size_t ind2,histo_t *meanmain,histo_t *count,size_t num_ells,char* multi_type,char *los_type,size_t num_threads)
{

	Box *boxes1=NULL;
	Box *boxes2=NULL;
	size_t *indices1,*indices2;
	size_t nfull1,nfull2;

#ifdef _VERBOSE
	print_catalogs();
#endif //_VERBOSE
	printf("*** 2-point correlation function multipoles\n");
	set_multi_type(multi_type,num_ells);
	set_los_type(los_type);
	set_cross(&ind1,&ind2);
	set_num_threads(num_threads);
#ifdef _DEBUG
	write_catalog(catalog[0],"debug_cat1.dat");
	write_catalog(catalog[1],"debug_cat2.dat");
#endif //_DEBUG
	timer(0);
	if ((ind1==1)||(ind2==1)) set_catalog_box(&catalog[0]);
	if ((ind1==2)||(ind2==2)) set_catalog_box(&catalog[1]);
	if ((ind1==1)&&(ind2==1)) init_params(catalog[0],catalog[0]);
	else if ((ind1==2)&&(ind2==2)) init_params(catalog[1],catalog[1]);
	else init_params(catalog[0],catalog[1]);
	
	if ((ind1==1)||(ind2==1)) boxes1=make_boxes_from_catalog(catalog[0],&indices1,&nfull1);
	if ((ind1==2)||(ind2==2)) boxes2=make_boxes_from_catalog(catalog[1],&indices2,&nfull2);
#ifdef _VERBOSE
	timer(2);
	printf("\n");
#endif //_VERBOSE	
#ifdef _DEBUG
	if ((ind1==1)||(ind2==1)) write_boxes(boxes1,"debug_box1.dat");
	if ((ind1==2)||(ind2==2)) write_boxes(boxes2,"debug_box2.dat");
#endif //_DEBUG
#ifdef _VERBOSE
	print_correlation(ind1,ind2);
#endif //_VERBOSE
	if ((ind1==1)&&(ind2==1)) auto_2pcf_multi(nfull1,indices1,boxes1,meanmain,count);
	else if ((ind1==2)&&(ind2==2)) auto_2pcf_multi(nfull2,indices2,boxes2,meanmain,count);
	else cross_2pcf_multi(nfull1,nfull2,indices1,indices2,boxes1,boxes2,meanmain,count);
#ifdef _VERBOSE
	timer(1);
	printf("\n");
#endif //_VERBOSE
#ifdef _DEBUG
	if ((ind1==1)&&(ind2==1)) write_histo(bin_main.n_bin*n_ells,count,"debug_11.dat");
	else if ((ind1==2)&&(ind2==2)) write_histo(bin_main.n_bin*n_ells,count,"debug_22.dat");
	else write_histo(bin_main.n_bin*n_ells,count,"debug_12.dat");
#endif //_DEBUG

	printf("*** Cleaning up\n");
	if ((ind1==1)||(ind2==1)) {
		free_boxes(boxes1);
		free(indices1);
	}
	if ((ind1==2)||(ind2==2)) {
		free_boxes(boxes2);
		free(indices2);
	}
	free_params();
}

void run_3pcf_multi(histo_t *count,size_t num_ells,char* multi_type,char *los_type,size_t num_threads)
{
#ifdef _VERBOSE
	print_catalogs();
#endif //_VERBOSE
	printf("*** 3-point correlation function multipoles\n");
	set_multi_type(multi_type,num_ells);
	set_los_type(los_type);
	set_num_threads(num_threads);
	cross_3pcf_multi(catalog,count);
}

void run_3pcf_multi_double_los(histo_t *count,size_t num_ells,char* multi_type,char *los_type,size_t num_threads)
{
#ifdef _VERBOSE
	print_catalogs();
#endif //_VERBOSE
	printf("*** 3-point correlation function multipoles with double los for cat2\n");
	set_multi_type(multi_type,num_ells);
	set_los_type(los_type);
	set_num_threads(num_threads);
	cross_3pcf_multi_double_los(catalog,count);
}

void run_3pcf_multi_radial(histo_t *count,size_t num_ells,char* multi_type,_Bool normalize,char *los_type,size_t num_threads)
{
#ifdef _VERBOSE
	print_catalogs();
#endif //_VERBOSE
	printf("*** 3-point radial correlation function multipoles\n");
	set_multi_type(multi_type,num_ells);
	set_los_type(los_type);
	set_num_threads(num_threads);
#ifdef _VERBOSE
	print_normalize(normalize);
#endif //_VERBOSE
	cross_3pcf_multi_radial(catalog,count,normalize);
}

void run_2pcf_multi_radial(histo_t *count,size_t num_ells,char* multi_type,_Bool normalize,char *los_type,size_t num_threads)
{
#ifdef _VERBOSE
	print_catalogs();
#endif //_VERBOSE
	printf("*** 2-point radial correlation function multipoles\n");
	set_multi_type(multi_type,num_ells);
	set_los_type(los_type);
	set_num_threads(num_threads);
#ifdef _VERBOSE
	print_normalize(normalize);
#endif //_VERBOSE
	cross_2pcf_multi_radial(catalog,count,normalize);
}

void run_4pcf_multi_radial(histo_t *count,size_t num_ells,char* multi_type,_Bool normalize,char *los_type,size_t num_threads)
{
#ifdef _VERBOSE
	print_catalogs();
#endif //_VERBOSE
	printf("*** 4-point radial correlation function multipoles\n");
	set_multi_type(multi_type,num_ells);
	set_los_type(los_type);
	set_num_threads(num_threads);
#ifdef _VERBOSE
	print_normalize(normalize);
#endif //_VERBOSE
	cross_4pcf_multi_radial(catalog,count,normalize);
}

void run_4pcf_multi(histo_t *count,size_t num_ells,char* multi_type,char *los_type,size_t num_threads)
{
#ifdef _VERBOSE
	print_catalogs();
#endif //_VERBOSE
	printf("*** 4-point correlation function multipoles\n");
	set_multi_type(multi_type,num_ells);
	set_los_type(los_type);
	set_num_threads(num_threads);
	cross_4pcf_multi(catalog,count);
}
