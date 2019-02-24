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
static LOS los[MAX_LOS] = {{.type=LOS_MIDPOINT},{.type=LOS_MIDPOINT}};
static Mesh meshs[MAX_CATS];
static size_t n_cats = MAX_CATS;

void print_num_threads()
{
	//Calculate number of threads
	size_t num_threads=0;
#pragma omp parallel
	{
#pragma omp atomic
		num_threads++;
	}
	printf(" - using %zu threads\n",num_threads);
}

void set_num_threads(size_t num_threads)
{
	if (num_threads>0) omp_set_num_threads(num_threads);
#ifdef _VERBOSE
	print_num_threads();
#endif //_VERBOSE
}

void print_bin(char* mode)
{
	printf("*** Binning\n");
	Bin bin=bin_main;
	if (!strcmp(mode,"aux")) bin=bin_aux;
	else if (!strcmp(mode,"bin")) {
		bin=bin_bin;
		printf(" - range: 0 <= %s < %zu\n",mode,bin.n_bin);
		return;
	}
	printf(" - range: %.3f < %s < %.3f\n",bin.min,mode,bin.max);
	if (bin.type==BIN_LIN) printf(" - #bins(lin): %zu\n",bin.n_bin);
	else if (bin.type==BIN_LOG) printf(" - #bins(log): %zu\n",bin.n_bin);
	else if (bin.type==BIN_CUSTOM) printf(" - #bins(custom): %zu\n",bin.n_bin);
	printf(" - resolution: %.3f\n",bin.step);
}

void set_bin(char* mode,histo_t* edges,size_t n_bin,char* type)
{
	Bin bin;

	bin.n_bin=n_bin;
	bin.log10min=0.;
	bin.log10max=0.;
	bin.edges=NULL;

	if (!strcmp(type,"log")) {
		bin.min=edges[0];
		bin.max=edges[n_bin];
		bin.type=BIN_LOG;
		bin.log10min=log10(bin.min);
		bin.log10max=log10(bin.max);
		bin.step=(bin.log10max-bin.log10min)/bin.n_bin;
		bin.edges=(histo_t *) malloc((bin.n_bin+1)*sizeof(histo_t));
		size_t ibin;
		for (ibin=0;ibin<=bin.n_bin;ibin++) bin.edges[ibin]=pow(10.,bin.step*ibin+bin.log10min);
	}
	else if (!strcmp(type,"custom")) {
		bin.min=edges[0];
		bin.max=edges[n_bin];
		bin.type=BIN_CUSTOM;
		bin.edges=edges;
		bin.step=(bin.max-bin.min)/bin.n_bin;
	}
	else if (!strcmp(type,"bin")) {
		bin.type=BIN_BIN;
	}
	else {
		bin.min=edges[0];
		bin.max=edges[n_bin];
		bin.type=BIN_LIN;
		bin.step=(bin.max-bin.min)/bin.n_bin;
		if (strcmp(type,"lin")) {
			fprintf(stderr," - invalid binning type. Choices: lin, log or custom.\n");
			fprintf(stderr," - I choose linear binning.\n");
		}
	}
	if (!strcmp(mode,"main")) bin_main=bin;
	else if (!strcmp(mode,"aux")) bin_aux=bin;
	else if (!strcmp(mode,"bin")) bin_bin=bin;
	else {
		bin_main=bin;
		fprintf(stderr," - invalid binning mode. Choices: main, aux, bin.\n");
		fprintf(stderr," - I choose main.\n");
	}
#ifdef _VERBOSE
	print_bin(mode);
#endif //_VERBOSE
}

void free_bin(Bin *bin) {
	if ((bin->type==BIN_LOG)&&(bin->n_bin>0)) free(bin->edges);
	bin->n_bin=0; 
}

void clear_bins()
{
	free_bin(&bin_main);
	free_bin(&bin_aux);
	free_bin(&bin_bin);
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

void set_pole(size_t num,char *type,size_t n_ells)
{
	Pole pole;
	if (!strcmp(type,"all")) pole.type=MULTI_ALL;
	else if (!strcmp(type,"even")) pole.type=MULTI_EVEN;
	else if (!strcmp(type,"odd")) pole.type=MULTI_ODD;
	else {
		pole.type=MULTI_ALL;
		fprintf(stderr," - invalid multipole type. Choices: all, even or odd.\n");
		fprintf(stderr," - I choose all.\n");
	}
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


void print_los(LOS l)
{
	printf("*** Line-of-sight\n");
	if (l.type==LOS_MIDPOINT) printf(" - los-type: midpoint\n");
	else if (l.type==LOS_ENDPOINT) printf(" - los-type: endpoint\n");
	else if (l.type==LOS_FIRSTPOINT) printf(" - los-type: firstpoint\n");
	else if (l.type==LOS_CUSTOM) printf(" - los-type: custom\n");
	printf(" - los-n: %zu\n",l.n);
}

void set_los(size_t num,char* type,histo_t* vec,size_t n)
{
	LOS l;
	if (!strcmp(type,"midpoint")) l.type=LOS_MIDPOINT;
	else if (!strcmp(type,"endpoint")) l.type=LOS_ENDPOINT;
	else if (!strcmp(type,"firstpoint")) l.type=LOS_FIRSTPOINT;
	else if (!strcmp(type,"custom")) {
		l.type=LOS_CUSTOM;
		l.los=vec;
	}
	else {
		l.type=LOS_MIDPOINT;
		fprintf(stderr," - invalid los type. Choices: midpoint, endpoint, firstpoint or custom.\n");
		fprintf(stderr," - I choose midpoint.\n");
	}
	l.n = n;
	los[num-1] = l;
#ifdef _VERBOSE
	print_los(l);
#endif //_VERBOSE
}

void clear_los()
{
	size_t ilos;
	for (ilos=0;ilos<MAX_LOS;ilos++) {
		los[ilos].los = NULL;
		los[ilos].n = 0;
	}
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
	if (corr_type==CORR_SCOS) set_los(0,"midpoint",NULL,0);
#ifdef _VERBOSE
	print_corr_type();
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

void set_catalog(size_t num,histo_t *p,histo_t *w,size_t *bin,size_t n,size_t dim_b,size_t dim_w)
{	
	Catalog cat;
	cat.pos = p;
	cat.weight = w;
	cat.bin = bin;
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

void run_2pcf_main(histo_t* meanmain,histo_t* count,char* corr_type,size_t num_threads)
{
	timer(0);
	set_num_catalogs();
#ifdef _VERBOSE
	printf("*** 2-point correlation function (main)\n");
	print_catalogs();
#endif //_VERBOSE
	set_corr_type(corr_type);
	set_num_threads(num_threads);
	set_meshs(cats,meshs,n_cats);
#ifdef _VERBOSE
	timer(1);
#endif //_VERBOSE
	if (n_cats==1) auto_2pcf_main(meshs[0],meanmain,count);
	else cross_2pcf_main(meshs[0],meshs[1],meanmain,count);
#ifdef _DEBUG
	if (n_cats==1) write_histo(bin_main.n_bin,count,"debug_auto.dat");
	else write_histo(bin_main.n_bin,count,"debug_cross.dat");
#endif //_DEBUG
#ifdef _VERBOSE
	printf("*** Cleaning up\n");
	timer(1);
	printf("\n");
#endif //_VERBOSE
	free_meshs(meshs,n_cats);
}

void run_2pcf_main_aux(histo_t* meanmain,histo_t* meanaux,histo_t* count,char* corr_type,size_t num_threads)
{
	timer(0);
	set_num_catalogs();
#ifdef _VERBOSE
	printf("*** 2-point correlation function (main,aux)\n");
	print_catalogs();
#endif //_VERBOSE
	set_corr_type(corr_type);
	set_num_threads(num_threads);
	set_meshs(cats,meshs,n_cats);
#ifdef _VERBOSE
	timer(1);
#endif //_VERBOSE
	if (n_cats==1) auto_2pcf_main_aux(meshs[0],meanmain,meanaux,count,los[0]);
	else cross_2pcf_main_aux(meshs[0],meshs[1],meanmain,meanaux,count,los[0]);
#ifdef _DEBUG
	if (n_cats==1) write_histo(bin_main.n_bin*bin_aux.n_bin,count,"debug_auto.dat");
	else write_histo(bin_main.n_bin*bin_aux.n_bin,count,"debug_cross.dat");
#endif //_DEBUG
#ifdef _VERBOSE
	printf("*** Cleaning up\n");
	timer(1);
	printf("\n");
#endif //_VERBOSE
	free_meshs(meshs,n_cats);
}

void run_2pcf_multi(histo_t *meanmain,histo_t *count,size_t num_threads)
{
	timer(0);
	set_num_catalogs();
#ifdef _VERBOSE
	printf("*** 2-point correlation function multipoles\n");
	print_catalogs();
#endif //_VERBOSE
	corr_type = CORR_SMU;
	set_num_threads(num_threads);
	set_meshs(cats,meshs,n_cats);
#ifdef _VERBOSE
	timer(1);
#endif //_VERBOSE
	if (n_cats==1) auto_2pcf_multi(meshs[0],meanmain,count,poles[0],los[0]);
	else cross_2pcf_multi(meshs[0],meshs[1],meanmain,count,poles[0],los[0]);	
#ifdef _DEBUG
	if (n_cats==1) write_histo(bin_main.n_bin*n_ells,count,"debug_auto.dat");
	else write_histo(bin_main.n_bin*n_ells,count,"debug_cross.dat");
#endif //_DEBUG
#ifdef _VERBOSE
	printf("*** Cleaning up\n");
	timer(1);
	printf("\n");
#endif //_VERBOSE
	free_meshs(meshs,n_cats);
}

void run_2pcf_multi_radial_legendre(histo_t *count,size_t num_threads)
{
	timer(0);
	set_num_catalogs();
#ifdef _VERBOSE
	printf("*** 2-point radial-legendre correlation function multipoles\n");
	print_catalogs();
#endif //_VERBOSE
	corr_type = CORR_SMU;
	set_num_threads(num_threads);
	set_meshs(cats,meshs,n_cats);
#ifdef _VERBOSE
	timer(1);
#endif //_VERBOSE
	cross_2pcf_multi_radial_legendre(meshs[0],meshs[1],count,poles,los[0]);
#ifdef _VERBOSE
	printf("*** Cleaning up\n");
	timer(1);
	printf("\n");
#endif //_VERBOSE
	free_meshs(meshs,n_cats);
}

void run_2pcf_multi_angular_legendre(histo_t *count,size_t num_threads)
{
	timer(0);
	set_num_catalogs();
#ifdef _VERBOSE
	printf("*** 2-point angular-legendre correlation function multipoles\n");
	print_catalogs();
#endif //_VERBOSE
	corr_type = CORR_SMU;
	set_num_threads(num_threads);
	set_meshs(cats,meshs,n_cats);
#ifdef _VERBOSE
	timer(1);
#endif //_VERBOSE
	cross_2pcf_multi_angular_legendre(meshs[0],meshs[1],count,poles,los[0]);
#ifdef _VERBOSE
	printf("*** Cleaning up\n");
	timer(1);
	printf("\n");
#endif //_VERBOSE
	free_meshs(meshs,n_cats);
}

void run_3pcf_multi(histo_t *count,size_t num_threads)
{
	timer(0);
	set_num_catalogs();
#ifdef _VERBOSE
	printf("*** 3-point correlation function multipoles\n");
	print_catalogs();
#endif //_VERBOSE
	corr_type = CORR_SMU;
	set_num_threads(num_threads);
	set_meshs(cats,meshs,n_cats);
#ifdef _VERBOSE
	timer(1);
#endif //_VERBOSE
	cross_3pcf_multi(meshs,n_cats,count,poles,los);
#ifdef _VERBOSE
	printf("*** Cleaning up\n");
	timer(1);
	printf("\n");
#endif //_VERBOSE
	free_meshs(meshs,n_cats);
}

void run_3pcf_multi_double_los(histo_t *count,size_t num_threads)
{
	timer(0);
	set_num_catalogs();
#ifdef _VERBOSE
	printf("*** 3-point correlation function multipoles with double los for cat2\n");
	print_catalogs();
	print_bin("main");
	print_bin("aux");
#endif //_VERBOSE
	corr_type = CORR_SMU;
	set_num_threads(num_threads);
	set_meshs(cats,meshs,n_cats);
#ifdef _VERBOSE
	timer(1);
#endif //_VERBOSE
	cross_3pcf_multi_double_los(meshs,n_cats,count,poles,los);
#ifdef _VERBOSE
	printf("*** Cleaning up\n");
	timer(1);
	printf("\n");
#endif //_VERBOSE
	free_meshs(meshs,n_cats);
}

void run_2pcf_multi_binned(histo_t *count,_Bool normalize,size_t num_threads)
{
	timer(0);
	set_num_catalogs();
#ifdef _VERBOSE
	printf("*** 2-point binned correlation function multipoles\n");
	print_catalogs();
	print_normalize(normalize);
#endif //_VERBOSE
	corr_type = CORR_SMU;
	set_num_threads(num_threads);
	set_meshs(cats,meshs,n_cats);
#ifdef _VERBOSE
	timer(1);
#endif //_VERBOSE
	cross_2pcf_multi_binned(meshs[0],meshs[1],count,poles[0],los[0],normalize);
#ifdef _VERBOSE
	printf("*** Cleaning up\n");
	timer(1);
	printf("\n");
#endif //_VERBOSE
	free_meshs(meshs,n_cats);
}

void run_4pcf_multi_binned(histo_t *count,_Bool normalize,size_t num_threads)
{
	timer(0);
	set_num_catalogs();
#ifdef _VERBOSE
	printf("*** 4-point binned correlation function multipoles\n");
	print_catalogs();
	print_normalize(normalize);
#endif //_VERBOSE
	corr_type = CORR_SMU;
	set_num_threads(num_threads);
	set_meshs(cats,meshs,n_cats);
#ifdef _VERBOSE
	timer(1);
#endif //_VERBOSE
	cross_4pcf_multi_binned(meshs,count,poles,los,normalize);
#ifdef _VERBOSE
	printf("*** Cleaning up\n");
	timer(1);
	printf("\n");
#endif //_VERBOSE
	free_meshs(meshs,n_cats);
}

void integrate_legendre(histo_t *count,histo_t *integral,size_t num_threads)
{
	printf("*** Integrating legendre\n");
	set_num_threads(num_threads);
	
	Pole pole = poles[0];
	size_t n_bin_tot = bin_main.n_bin*pole.n_ells;
	size_t ibin;
	for (ibin=0;ibin<n_bin_tot;ibin++) integral[ibin] = 0.;

#pragma omp for nowait schedule(dynamic)
	for (ibin=0;ibin<bin_main.n_bin;ibin++) {
		size_t ibin_aux;
		for (ibin_aux=0;ibin_aux<bin_aux.n_bin;ibin_aux++) {
			histo_t dist_aux = get_bin_mid(ibin_aux,bin_aux);
			histo_t leg[MAX_ELLS];
			legendre(dist_aux,leg,pole.type);
			histo_t weight = count[ibin_aux+bin_aux.n_bin*ibin];
			size_t ill;
			for (ill=0;ill<pole.n_ells;ill++) integral[ill+pole.n_ells*ibin] += weight*leg[ill];
		}
	}
}

void integrate_radial_legendre(histo_t *count,histo_t *integral,size_t num_threads)
{
	timer(0);
	printf("*** Integrating radial legendre\n");
	set_num_threads(num_threads);
	
	size_t n_bin_tot = bin_main.n_bin*poles[0].n_ells*bin_main.n_bin*poles[1].n_ells;
	size_t ibin;
	for (ibin=0;ibin<n_bin_tot;ibin++) integral[ibin] = 0.;

	MULTI_TYPE multi_type1 = poles[0].type;
	MULTI_TYPE multi_type2 = poles[1].type;
	size_t n_ells1 = poles[0].n_ells;
	size_t n_ells2 = poles[1].n_ells;

#pragma omp for nowait schedule(dynamic)
	for (ibin=0;ibin<bin_main.n_bin;ibin++) {
		histo_t dist_main = get_bin_mid(ibin,bin_main);
		size_t ibin_aux;
		for (ibin_aux=0;ibin_aux<bin_aux.n_bin;ibin_aux++) {
			histo_t dist_aux = get_bin_mid(ibin_aux,bin_aux);
			histo_t leg1[MAX_ELLS],leg2[MAX_ELLS];
			legendre(dist_aux,leg1,multi_type1);
			histo_t weight = count[ibin_aux+bin_aux.n_bin*ibin];
			long ibin2; //can be <0
			for (ibin2=bin_main.n_bin-1;ibin2>=0;ibin2--) {
				histo_t dist_main2 = get_bin_mid(ibin2,bin_main);
				histo_t dist_aux2 = dist_main*dist_aux/dist_main2;
				if (fabs(dist_aux2)>1.) break;
				legendre(dist_aux2,leg2,multi_type2);
				size_t ill1,ill2;
				for (ill1=0;ill1<n_ells1;ill1++) {
					histo_t tmp = weight*leg1[ill1]/dist_main2;
					for(ill2=0;ill2<n_ells2;ill2++) integral[ill2+n_ells2*(ill1+n_ells1*(ibin2+bin_main.n_bin*ibin))] += tmp*leg2[ill2];
				}
			}
		}
	} //end omp parallel
	timer(1);
}

void integrate_angular_legendre(histo_t *count,histo_t *integral,size_t num_threads)
{
	timer(0);
	printf("*** Integrating angular legendre\n");
	set_num_threads(num_threads);
	
	size_t n_bin_tot = bin_main.n_bin*poles[0].n_ells*bin_main.n_bin*poles[1].n_ells;
	size_t ibin;
	for (ibin=0;ibin<n_bin_tot;ibin++) integral[ibin] = 0.;

	MULTI_TYPE multi_type1 = poles[0].type;
	MULTI_TYPE multi_type2 = poles[1].type;
	size_t n_ells1 = poles[0].n_ells;
	size_t n_ells2 = poles[1].n_ells;

#pragma omp for nowait schedule(dynamic)
	for (ibin=0;ibin<bin_main.n_bin;ibin++) {
		histo_t fast_dist_main = get_bin_mid(ibin,bin_main);
		fast_dist_main *= fast_dist_main;
		size_t ibin_aux;
		for (ibin_aux=0;ibin_aux<bin_aux.n_bin;ibin_aux++) {
			histo_t dist_aux = get_bin_mid(ibin_aux,bin_aux);
			histo_t leg1[MAX_ELLS],leg2[MAX_ELLS],leg2bis[MAX_ELLS];
			legendre(dist_aux,leg1,multi_type1);
			histo_t weight = count[ibin_aux+bin_aux.n_bin*ibin];
			long ibin2; //can be <0
			for (ibin2=bin_main.n_bin-1;ibin2>=0;ibin2--) {
				histo_t fast_dist_main2 = get_bin_mid(ibin2,bin_main);
				fast_dist_main2 *= fast_dist_main2;
				histo_t fast_dist_aux2 = 1.-(fast_dist_main/fast_dist_main2)*(1.-dist_aux*dist_aux);
				if (fast_dist_aux2<0.) break;
				histo_t dist_aux2 = sqrt(fast_dist_aux2);
				legendre(dist_aux2,leg2,multi_type2);
				legendre(-dist_aux2,leg2bis,multi_type2);
				size_t ill1,ill2;
				for (ill1=0;ill1<n_ells1;ill1++) {
					histo_t tmp = weight*leg1[ill1]/(fast_dist_main2*dist_aux2);
					for(ill2=0;ill2<n_ells2;ill2++) integral[ill2+n_ells2*(ill1+n_ells1*(ibin2+bin_main.n_bin*ibin))] += tmp*(leg2[ill2]+leg2bis[ill2]);
				}
			}
		}
	} //end omp parallel
	timer(1);
}
