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
//                      Correlators with OpenMP                      //
/*********************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "define.h"
#include "common.h"

static histo_t fast_dist_main_min,fast_dist_main_max;
static histo_t fast_dist_aux_min,fast_dist_aux_max;
static histo_t fast_dist_radial_min,fast_dist_radial_max;

static histo_t my_sqrt(histo_t x)
{
#ifdef _FLOAT32
	return sqrtf(x);
#else
	return sqrt(x);
#endif	
}

static histo_t my_abs(histo_t x)
{
#ifdef _FLOAT32
	return fabsf(x);
#else
	return fabs(x);
#endif //_FLOAT32
}

/*
static float FastInvSqrt(float x) {
  float xhalf=0.5f * x;
  int i=*(int*)&x;         // evil floating point bit level hacking
  i=0x5f3759df - (i >> 1);  // what the fuck?
  x=*(float*)&i;
  x=x*(1.5f-(xhalf*x*x));
  return x;
}*/

static histo_t get_weight(histo_t *weight1,histo_t *weight2)
{
	size_t idim;
	histo_t weight=0;
	for(idim=0;idim<dim_weight;idim++){
		weight+=weight1[idim]*weight2[idim];
	}
	return weight;
}

static histo_t get_fast_distance_main(histo_t *pos1,histo_t *pos2)
{
	if ((corr_type==CORR_SMU)||(corr_type==CORR_SCOS)) {
		size_t idim;
		histo_t dist=0.;
		histo_t diff;
		for(idim=0;idim<dim_pos;idim++){
			diff=pos2[idim]-pos1[idim];
			dist+=diff*diff;
		}
		return dist;
	}
	else if (corr_type==CORR_ANGULAR) {
		size_t idim;
		histo_t dist=0;
		for(idim=0;idim<dim_pos;idim++){
			dist+=pos1[idim]*pos2[idim];
		}
		return 1.-my_abs(dist);
	}
	return 0.;
}

static histo_t get_distance_main(histo_t fast_dist)
{
	if ((corr_type==CORR_SMU)||(corr_type==CORR_SCOS)) return my_sqrt(fast_dist);
	else if (corr_type==CORR_ANGULAR) return acos(1.-fast_dist);
	return 0.;
}

static histo_t get_inv_distance_main(histo_t dist)
{
	if ((corr_type==CORR_SMU)||(corr_type==CORR_SCOS)) return dist*dist;
	else if (corr_type==CORR_ANGULAR) return 1.-cos(dist);
	return 0.;
}

static _Bool visit_main(histo_t fast_dist)
{
	return ((fast_dist<fast_dist_main_max)&&(fast_dist>fast_dist_main_min));
}

static void set_fast_distance_main_limit()
{
	fast_dist_main_min=get_inv_distance_main(bin_main.min);
	fast_dist_main_max=get_inv_distance_main(bin_main.max);
}

static histo_t get_fast_distance_aux(histo_t *pos1,histo_t *pos2)
{	
	size_t idim;
	histo_t dist=0.;
	histo_t norm1=0.;
	histo_t norm2=0.;
	if (corr_type==CORR_SMU) {
		histo_t diff,los;
		for(idim=0;idim<dim_pos;idim++){
			diff=pos2[idim]-pos1[idim];
			if (los_type==LOS_ENDPOINT) los=pos1[idim];
			else los=pos2[idim]+pos1[idim];
			los=pos2[idim]+pos1[idim];
			norm1+=diff*diff;
			norm2+=los*los;
			dist+=diff*los;
		}
	}
	else if (corr_type==CORR_SCOS) {
		for(idim=0;idim<dim_pos;idim++){
			norm1+=pos1[idim]*pos1[idim];
			norm2+=pos2[idim]*pos2[idim];
			dist+=pos1[idim]*pos2[idim];
		}
	}
	//return dist*dist/(norm1*norm2);
	return dist*my_abs(dist)/(norm1*norm2); //to get the sign
}


static histo_t get_distance_aux(histo_t fast_dist)
{
	if (fast_dist<0) return -my_sqrt(my_abs(fast_dist));
	return my_sqrt(fast_dist);
}

static histo_t get_inv_distance_aux(histo_t dist)
{
	//return dist*dist;
	return dist*my_abs(dist);
}

static _Bool visit_aux(histo_t fast_dist)
{
	//printf("%.3f  ",fast_dist);
	return ((fast_dist<=fast_dist_aux_max)&&(fast_dist>=fast_dist_aux_min));
}

static void set_fast_distance_aux_limit()
{
	fast_dist_aux_min=get_inv_distance_aux(bin_aux.min);
	fast_dist_aux_max=get_inv_distance_aux(bin_aux.max);
}


static histo_t get_fast_distance_radial(histo_t *pos)
{	
	size_t idim;
	histo_t dist=0.;
	for(idim=0;idim<dim_pos;idim++){
		dist+=pos[idim]*pos[idim];
	}
	//return dist*dist/(norm1*norm2);
	return dist; //to get the sign
}

static histo_t get_distance_radial(histo_t fast_dist)
{	
	return my_sqrt(fast_dist);
}

static histo_t get_inv_distance_radial(histo_t dist)
{
	return dist*dist;
}

static _Bool visit_radial(histo_t fast_dist)
{
	//printf("%.3f  ",fast_dist);
	return ((fast_dist<=fast_dist_radial_max)&&(fast_dist>=fast_dist_radial_min));
}


static void set_fast_distance_radial_limit()
{
	fast_dist_radial_min=get_inv_distance_radial(bin_radial.min);
	fast_dist_radial_max=get_inv_distance_radial(bin_radial.max);
}


static _Bool visit_box(Box box1,Box box2)
{
	if ((corr_type==CORR_SMU)||(corr_type==CORR_SCOS)) {
		size_t idim;
		for(idim=0;idim<dim_box;idim++) {
			//printf("%zu %zu %zu %zu   ",idim,box2.index[idim],box1.visit_min[idim],box1.visit_max[idim]);
			if (box2.index[idim]<box1.visit_min[idim]) return 0;
			if (box2.index[idim]>box1.visit_max[idim]) return 0;
		}
	}
	else if (corr_type==CORR_ANGULAR) {
		size_t idim=0;
		if (box2.index[idim]<box1.visit_min[idim]) return 0;
		if (box2.index[idim]>box1.visit_max[idim]) return 0;
		idim=1;
		if ((box1.visit_min[idim]<box1.visit_max[idim])&&((box2.index[idim]<box1.visit_min[idim])||(box2.index[idim]>box1.visit_max[idim]))) return 0;
		if ((box1.visit_min[idim]>box1.visit_max[idim])&&(box2.index[idim]<box1.visit_min[idim])&&(box2.index[idim]>box1.visit_max[idim])) return 0;
	}
	return 1;
}


void cross_2pcf_main(size_t nbox_full1,size_t nbox_full2,size_t *indices1,size_t *indices2,
		    Box *boxes1,Box *boxes2,
			histo_t meanmain[],histo_t count[])
{
  //////
  // 2pcf cross

	histo_t *htreadmain,*threadcount;
	set_fast_distance_main_limit();
	size_t n_bin_tot=bin_main.n_bin;
	size_t ibin,ibox1;

	for(ibin=0;ibin<n_bin_tot;ibin++) {
		meanmain[ibin]=0.;
		count[ibin]=0.;
	}

#pragma omp parallel default(none)				\
  shared(nbox_full1,nbox_full2,indices1,indices2,boxes1,boxes2,meanmain,count,n_bin_tot,bin_main,dim_pos,dim_weight) private(htreadmain,threadcount)
	{
		size_t ibin;
		htreadmain=(histo_t *) malloc(n_bin_tot*sizeof(histo_t));
		threadcount=(histo_t *) malloc(n_bin_tot*sizeof(histo_t));

		for(ibin=0;ibin<n_bin_tot;ibin++) {
			htreadmain[ibin]=0.;
			threadcount[ibin]=0.;
		}

#pragma omp for nowait schedule(dynamic)
		for(ibox1=0;ibox1<nbox_full1;ibox1++) {
			Box box1=boxes1[indices1[ibox1]];
			size_t ibox2;
			for(ibox2=0;ibox2<nbox_full2;ibox2++) {
				Box box2=boxes2[indices2[ibox2]];
				if (visit_box(box1,box2)) {
					size_t nobj1=box1.n_obj;
					size_t nobj2=box2.n_obj;
					size_t iobj1;
					for(iobj1=0;iobj1<nobj1;iobj1++) {
						histo_t *pos1=&(box1.pos[dim_pos*iobj1]);
						histo_t *weight1=&(box1.weight[dim_weight*iobj1]);
						size_t iobj2;
						for(iobj2=0;iobj2<nobj2;iobj2++) {
							histo_t *pos2=&(box2.pos[dim_pos*iobj2]);
							histo_t *weight2=&(box2.weight[dim_weight*iobj2]);
							histo_t fast_dist=get_fast_distance_main(pos1,pos2);
							if (visit_main(fast_dist)) {
								histo_t dist_main=get_distance_main(fast_dist);
								histo_t weight=get_weight(weight1,weight2);
								size_t ibin=get_bin_index(dist_main,bin_main);
								htreadmain[ibin]+=weight*dist_main;
								threadcount[ibin]+=weight;						
							}
						}
					}
				}	
			}
		} // end omp for
	
#pragma omp critical
		{
			for(ibin=0;ibin<n_bin_tot;ibin++) {
				meanmain[ibin]+=htreadmain[ibin];
				count[ibin]+=threadcount[ibin];
			}
			free(htreadmain);
			free(threadcount);
		}
	} //end omp parallel
	for(ibin=0;ibin<n_bin_tot;ibin++) {
		meanmain[ibin]/=count[ibin];
	}
}


void auto_2pcf_main(size_t nbox_full,size_t *indices,Box *boxes,
			histo_t meanmain[],histo_t count[])
{
  //////
  // 2pcf auto

	histo_t *htreadmain,*threadcount;
	set_fast_distance_main_limit();
	size_t n_bin_tot=bin_main.n_bin;
	size_t ibin;

	for(ibin=0;ibin<n_bin_tot;ibin++) {
		meanmain[ibin]=0.;
		count[ibin]=0.;
	}

#pragma omp parallel default(none)				\
  shared(nbox_full,indices,boxes,meanmain,count,n_bin_tot,bin_main,dim_pos,dim_weight) private(htreadmain,threadcount)
	{
		size_t ibin,ibox1;
		htreadmain=(histo_t *) malloc(n_bin_tot*sizeof(histo_t));
		threadcount=(histo_t *) malloc(n_bin_tot*sizeof(histo_t));

		for(ibin=0;ibin<n_bin_tot;ibin++) {
			htreadmain[ibin]=0.;
			threadcount[ibin]=0.;
		}

#pragma omp for nowait schedule(dynamic)
		for(ibox1=0;ibox1<nbox_full;ibox1++) {
			Box box1=boxes[indices[ibox1]];
			size_t nobj1=box1.n_obj;
			size_t iobj1;
			for(iobj1=0;iobj1<nobj1;iobj1++) {
				histo_t *pos1=&(box1.pos[dim_pos*iobj1]);
				histo_t *weight1=&(box1.weight[dim_weight*iobj1]);
				size_t iobj2;
				for(iobj2=iobj1+1;iobj2<nobj1;iobj2++) {
					histo_t *pos2=&(box1.pos[dim_pos*iobj2]);
					histo_t *weight2=&(box1.weight[dim_weight*iobj2]);
					histo_t fast_dist=get_fast_distance_main(pos1,pos2);
					if (visit_main(fast_dist)) {
						histo_t dist_main=get_distance_main(fast_dist);
						histo_t weight=get_weight(weight1,weight2);
						size_t ibin=get_bin_index(dist_main,bin_main);
						htreadmain[ibin]+=weight*dist_main;
						threadcount[ibin]+=weight;						
					}						
				}
			}
			size_t ibox2;
			for(ibox2=ibox1+1;ibox2<nbox_full;ibox2++) {
				Box box2=boxes[indices[ibox2]];
				if (visit_box(box1,box2)) {
					size_t nobj2=box2.n_obj;
					for(iobj1=0;iobj1<nobj1;iobj1++) {
						histo_t *pos1=&(box1.pos[dim_pos*iobj1]);
						histo_t *weight1=&(box1.weight[dim_weight*iobj1]);
						size_t iobj2;
						for(iobj2=0;iobj2<nobj2;iobj2++) {
							histo_t *pos2=&(box2.pos[dim_pos*iobj2]);
							histo_t *weight2=&(box2.weight[dim_weight*iobj2]);
							histo_t fast_dist=get_fast_distance_main(pos1,pos2);
							if (visit_main(fast_dist)) {
								histo_t dist_main=get_distance_main(fast_dist);
								histo_t weight=get_weight(weight1,weight2);
								size_t ibin=get_bin_index(dist_main,bin_main);
								htreadmain[ibin]+=weight*dist_main;
								threadcount[ibin]+=weight;						
							}
						}
					}
				}	
			}
		} // end omp for
	
#pragma omp critical
		{
			for(ibin=0;ibin<n_bin_tot;ibin++) {
				meanmain[ibin]+=htreadmain[ibin];
				count[ibin]+=threadcount[ibin];
			}
			free(htreadmain);
			free(threadcount);
		}
	} //end omp parallel
	for(ibin=0;ibin<n_bin_tot;ibin++) {
		meanmain[ibin]/=count[ibin];
	}
}



void cross_2pcf_main_aux(size_t nbox_full1,size_t nbox_full2,size_t *indices1,size_t *indices2,
		    Box *boxes1,Box *boxes2,
			histo_t meanmain[],histo_t meanaux[],histo_t count[])
{
  //////
  // 2pcf cross

	histo_t *htreadmain,*htreadaux,*threadcount;
	set_fast_distance_main_limit();
	set_fast_distance_aux_limit();
	size_t n_bin_tot=bin_main.n_bin*bin_aux.n_bin;
	size_t ibin;

	for(ibin=0;ibin<n_bin_tot;ibin++) {
		meanmain[ibin]=0.;
		meanaux[ibin]=0.;
		count[ibin]=0.;
	}

#pragma omp parallel default(none)				\
  shared(nbox_full1,nbox_full2,indices1,indices2,boxes1,boxes2,meanmain,meanaux,count,n_bin_tot,bin_main,bin_aux,dim_pos,dim_weight) private(htreadmain,htreadaux,threadcount)
	{
		size_t ibin,ibox1;
		htreadmain=(histo_t *) malloc(n_bin_tot*sizeof(histo_t));
		htreadaux=(histo_t *) malloc(n_bin_tot*sizeof(histo_t));
		threadcount=(histo_t *) malloc(n_bin_tot*sizeof(histo_t));

		for(ibin=0;ibin<n_bin_tot;ibin++) {
			htreadmain[ibin]=0.;
			htreadaux[ibin]=0.;
			threadcount[ibin]=0.;
		}

#pragma omp for nowait schedule(dynamic)
		for(ibox1=0;ibox1<nbox_full1;ibox1++) {
			Box box1=boxes1[indices1[ibox1]];
			size_t ibox2;
			for(ibox2=0;ibox2<nbox_full2;ibox2++) {
				Box box2=boxes2[indices2[ibox2]];
				if (visit_box(box1,box2)) {
					size_t nobj1=box1.n_obj;
					size_t nobj2=box2.n_obj;
					size_t iobj1;
					for(iobj1=0;iobj1<nobj1;iobj1++) {
						histo_t *pos1=&(box1.pos[dim_pos*iobj1]);
						histo_t *weight1=&(box1.weight[dim_weight*iobj1]);
						size_t iobj2;
						for(iobj2=0;iobj2<nobj2;iobj2++) {
							histo_t *pos2=&(box2.pos[dim_pos*iobj2]);
							histo_t *weight2=&(box2.weight[dim_weight*iobj2]);
							histo_t fast_dist_main=get_fast_distance_main(pos1,pos2);
							if (visit_main(fast_dist_main)) {
								histo_t fast_dist_aux=get_fast_distance_aux(pos1,pos2);
								if (visit_aux(fast_dist_aux)) {
									histo_t dist_main=get_distance_main(fast_dist_main);
									histo_t dist_aux=get_distance_aux(fast_dist_aux);
									size_t ibin=get_bin_index(dist_aux,bin_aux)+bin_aux.n_bin*get_bin_index(dist_main,bin_main);
									histo_t weight=get_weight(weight1,weight2);
									htreadmain[ibin]+=weight*dist_main;
									htreadaux[ibin]+=weight*dist_aux;
									threadcount[ibin]+=weight;
								}							
							}
						}
					}
				}	
			}
		} // end omp for
	
#pragma omp critical
		{
			for(ibin=0;ibin<n_bin_tot;ibin++) {
				meanmain[ibin]+=htreadmain[ibin];
				meanaux[ibin]+=htreadaux[ibin];
				count[ibin]+=threadcount[ibin];
			}
			free(htreadmain);
			free(htreadaux);
			free(threadcount);
		}
	} //end omp parallel
	for(ibin=0;ibin<n_bin_tot;ibin++) {
		meanmain[ibin]/=count[ibin];
		meanaux[ibin]/=count[ibin];
	}
}

void auto_2pcf_main_aux(size_t nbox_full,size_t *indices,Box *boxes,
			histo_t meanmain[],histo_t meanaux[],histo_t count[])
{
  //////
  // 2pcf auto

	histo_t *htreadmain,*htreadaux,*threadcount;
	set_fast_distance_main_limit();
	set_fast_distance_aux_limit();
	size_t n_bin_tot=bin_main.n_bin*bin_aux.n_bin;
	size_t ibin;

	for(ibin=0;ibin<n_bin_tot;ibin++) {
		meanmain[ibin]=0.;
		meanaux[ibin]=0.;
		count[ibin]=0.;
	}

#pragma omp parallel default(none)				\
  shared(nbox_full,indices,boxes,meanmain,meanaux,count,n_bin_tot,bin_main,bin_aux,dim_pos,dim_weight) private(htreadmain,htreadaux,threadcount)
	{
		size_t ibin,ibox1;
		htreadmain=(histo_t *) malloc(n_bin_tot*sizeof(histo_t));
		htreadaux=(histo_t *) malloc(n_bin_tot*sizeof(histo_t));
		threadcount=(histo_t *) malloc(n_bin_tot*sizeof(histo_t));

		for(ibin=0;ibin<n_bin_tot;ibin++) {
			htreadmain[ibin]=0.;
			htreadaux[ibin]=0.;
			threadcount[ibin]=0.;
		}

#pragma omp for nowait schedule(dynamic)
		for(ibox1=0;ibox1<nbox_full;ibox1++) {
			Box box1=boxes[indices[ibox1]];
			size_t nobj1=box1.n_obj;
			size_t iobj1;
			for(iobj1=0;iobj1<nobj1;iobj1++) {
				histo_t *pos1=&(box1.pos[dim_pos*iobj1]);
				histo_t *weight1=&(box1.weight[dim_weight*iobj1]);
				size_t iobj2;
				for(iobj2=iobj1+1;iobj2<nobj1;iobj2++) {
					histo_t *pos2=&(box1.pos[dim_pos*iobj2]);
					histo_t *weight2=&(box1.weight[dim_weight*iobj2]);
					histo_t fast_dist_main=get_fast_distance_main(pos1,pos2);
					if (visit_main(fast_dist_main)) {
						histo_t fast_dist_aux=get_fast_distance_aux(pos1,pos2);
						if (visit_aux(fast_dist_aux)) {
							histo_t dist_main=get_distance_main(fast_dist_main);
							histo_t dist_aux=get_distance_aux(fast_dist_aux);
							size_t ibin=get_bin_index(dist_aux,bin_aux)+bin_aux.n_bin*get_bin_index(dist_main,bin_main);
							histo_t weight=get_weight(weight1,weight2);
							htreadmain[ibin]+=weight*dist_main;
							htreadaux[ibin]+=weight*dist_aux;
							threadcount[ibin]+=weight;
						}							
					}
				}
			}
			size_t ibox2;
			for(ibox2=ibox1+1;ibox2<nbox_full;ibox2++) {
				Box box2=boxes[indices[ibox2]];
				if (visit_box(box1,box2)) {
					size_t nobj2=box2.n_obj;
					for(iobj1=0;iobj1<nobj1;iobj1++) {
						histo_t *pos1=&(box1.pos[dim_pos*iobj1]);
						histo_t *weight1=&(box1.weight[dim_weight*iobj1]);
						size_t iobj2;
						for(iobj2=0;iobj2<nobj2;iobj2++) {
							histo_t *pos2=&(box2.pos[dim_pos*iobj2]);
							histo_t *weight2=&(box2.weight[dim_weight*iobj2]);
							histo_t fast_dist_main=get_fast_distance_main(pos1,pos2);
							if (visit_main(fast_dist_main)) {
								histo_t fast_dist_aux=get_fast_distance_aux(pos1,pos2);
								if (visit_aux(fast_dist_aux)) {
									histo_t dist_main=get_distance_main(fast_dist_main);
									histo_t dist_aux=get_distance_aux(fast_dist_aux);
									size_t ibin=get_bin_index(dist_aux,bin_aux)+bin_aux.n_bin*get_bin_index(dist_main,bin_main);
									histo_t weight=get_weight(weight1,weight2);
									htreadmain[ibin]+=weight*dist_main;
									htreadaux[ibin]+=weight*dist_aux;
									threadcount[ibin]+=weight;
								}							
							}
						}
					}
				}	
			}
		} // end omp for
	
#pragma omp critical
		{
			for(ibin=0;ibin<n_bin_tot;ibin++) {
				meanmain[ibin]+=htreadmain[ibin];
				meanaux[ibin]+=htreadaux[ibin];
				count[ibin]+=threadcount[ibin];
			}
			free(htreadmain);
			free(htreadaux);
			free(threadcount);
		}
	} //end omp parallel
	for(ibin=0;ibin<n_bin_tot;ibin++) {
		meanmain[ibin]/=count[ibin];
		meanaux[ibin]/=count[ibin];
	}
}

void legendre(histo_t fast_dist,histo_t leg[]) {
	
	if (multi_type==MULTI_ALL) {
		histo_t dist_aux=get_distance_aux(fast_dist);
		//printf("%.3f %.3f  ",dist_aux,fast_dist);
		//if (dist_aux>0) printf("%.3f  ",dist_aux);
		legendre_all(dist_aux,my_abs(fast_dist),leg);
	}
	else if (multi_type==MULTI_EVEN) {
		//printf("%.3f ",fast_dist);
		legendre_even(my_abs(fast_dist),leg);
	}
	else if (multi_type==MULTI_ODD) {
		histo_t dist_aux=get_distance_aux(fast_dist);
		legendre_odd(dist_aux,my_abs(fast_dist),leg);
	}
}

void cross_2pcf_multi(size_t nbox_full1,size_t nbox_full2,size_t *indices1,size_t *indices2,Box *boxes1,Box *boxes2,histo_t meanmain[], histo_t count[])
{
  //////
  // 2pcf cross

	histo_t *threadmeanmain,*threadcount;
	set_fast_distance_main_limit();
	set_fast_distance_aux_limit();
	size_t n_bin_tot=bin_main.n_bin*n_ells;
	size_t ibin;
	histo_t leg[MAX_ELLS];

	for(ibin=0;ibin<n_bin_tot;ibin++) {
		meanmain[ibin]=0;
		count[ibin]=0;
	}

#pragma omp parallel default(none)				\
  shared(nbox_full1,indices1,boxes1,nbox_full2,indices2,boxes2,meanmain,count,n_bin_tot,bin_main,bin_aux,dim_pos,dim_weight,n_ells) private(threadmeanmain,threadcount,leg)
	{
		size_t ibin,ibox1;
		threadmeanmain=(histo_t *) malloc(n_bin_tot*sizeof(histo_t));
		threadcount=(histo_t *) malloc(n_bin_tot*sizeof(histo_t));

		for(ibin=0;ibin<n_bin_tot;ibin++) {
			threadmeanmain[ibin]=0;
			threadcount[ibin]=0;
		}
#pragma omp for nowait schedule(dynamic)
		for(ibox1=0;ibox1<nbox_full1;ibox1++) {
			Box box1=boxes1[indices1[ibox1]];
			size_t ibox2;
			for(ibox2=0;ibox2<nbox_full2;ibox2++) {
				Box box2=boxes2[indices2[ibox2]];
				if (visit_box(box1,box2)) {
					size_t nobj1=box1.n_obj;
					size_t nobj2=box2.n_obj;
					size_t iobj1;
					for(iobj1=0;iobj1<nobj1;iobj1++) {
						histo_t *pos1=&(box1.pos[dim_pos*iobj1]);
						histo_t *weight1=&(box1.weight[dim_weight*iobj1]);
						size_t iobj2;
						for(iobj2=0;iobj2<nobj2;iobj2++) {
							histo_t *pos2=&(box2.pos[dim_pos*iobj2]);
							histo_t fast_dist_main=get_fast_distance_main(pos1,pos2);
							if (visit_main(fast_dist_main)) {
								histo_t fast_dist_aux=get_fast_distance_aux(pos1,pos2);
								if (visit_aux(fast_dist_aux)) {
									histo_t dist_main=get_distance_main(fast_dist_main);
									histo_t weight=get_weight(weight1,&(box2.weight[dim_weight*iobj2]));
									legendre(fast_dist_aux,leg);
									size_t ill;
									for(ill=0;ill<n_ells;ill++) {
										histo_t w=weight*leg[ill];
										size_t ibin=ill+n_ells*get_bin_index(dist_main,bin_main);
										threadmeanmain[ibin]+=w*dist_main;
										threadcount[ibin]+=w;
									}
								}							
							}
						}							
					}
				}	
			}
		} // end omp for
	
#pragma omp critical
		{
			for (ibin=0;ibin<n_bin_tot;ibin++) {
				meanmain[ibin]+=threadmeanmain[ibin];
				count[ibin]+=threadcount[ibin];
			}
			free(threadmeanmain);
			free(threadcount);
		}
	} //end omp parallel
	for (ibin=0;ibin<n_bin_tot;ibin++) meanmain[ibin]/=count[ibin];
}


void auto_2pcf_multi(size_t nbox_full,size_t *indices,Box *boxes,histo_t meanmain[], histo_t count[])
{
  //////
  // 2pcf cross

	histo_t *threadmeanmain,*threadcount;
	set_fast_distance_main_limit();
	set_fast_distance_aux_limit();
	size_t n_bin_tot=bin_main.n_bin*n_ells;
	size_t ibin;
	histo_t leg[MAX_ELLS];

	for(ibin=0;ibin<n_bin_tot;ibin++) {
		meanmain[ibin]=0;
		count[ibin]=0;
	}

#pragma omp parallel default(none)				\
  shared(nbox_full,indices,boxes,meanmain,count,n_bin_tot,bin_main,bin_aux,dim_pos,dim_weight,n_ells) private(threadmeanmain,threadcount,leg)
	{
		size_t ibin,ibox1;
		threadmeanmain=(histo_t *) malloc(bin_main.n_bin*sizeof(histo_t));
		threadcount=(histo_t *) malloc(n_bin_tot*sizeof(histo_t));

		for(ibin=0;ibin<n_bin_tot;ibin++) {
			threadmeanmain[ibin]=0;
			threadcount[ibin]=0;
		}

#pragma omp for nowait schedule(dynamic)
		for(ibox1=0;ibox1<nbox_full;ibox1++) {
			Box box1=boxes[indices[ibox1]];
			size_t nobj1=box1.n_obj;
			size_t iobj1;
			for(iobj1=0;iobj1<nobj1;iobj1++) {
				histo_t *pos1=&(box1.pos[dim_pos*iobj1]);
				histo_t *weight1=&(box1.weight[dim_weight*iobj1]);
				size_t iobj2;
				for(iobj2=iobj1+1;iobj2<nobj1;iobj2++) {
					histo_t *pos2=&(box1.pos[dim_pos*iobj2]);
					histo_t fast_dist_main=get_fast_distance_main(pos1,pos2);
					if (visit_main(fast_dist_main)) {
						histo_t fast_dist_aux=get_fast_distance_aux(pos1,pos2);
						if (visit_aux(fast_dist_aux)) {
							histo_t dist_main=get_distance_main(fast_dist_main);
							histo_t weight=get_weight(weight1,&(box1.weight[dim_weight*iobj2]));
							legendre(fast_dist_aux,leg);
							size_t ill;
							for(ill=0;ill<n_ells;ill++) {
								histo_t w=weight*leg[ill];
								size_t ibin=ill+n_ells*get_bin_index(dist_main,bin_main);
								threadmeanmain[ibin]+=w*dist_main;
								threadcount[ibin]+=w;
							}
						}							
					}
				}
			}
			size_t ibox2;
			for(ibox2=ibox1+1;ibox2<nbox_full;ibox2++) {
				Box box2=boxes[indices[ibox2]];
				if (visit_box(box1,box2)) {
					size_t nobj2=box2.n_obj;
					for(iobj1=0;iobj1<nobj1;iobj1++) {
						histo_t *pos1=&(box1.pos[dim_pos*iobj1]);
						histo_t *weight1=&(box1.weight[dim_weight*iobj1]);
						size_t iobj2;
						for(iobj2=0;iobj2<nobj2;iobj2++) {
							histo_t *pos2=&(box2.pos[dim_pos*iobj2]);
							histo_t fast_dist_main=get_fast_distance_main(pos1,pos2);
							if (visit_main(fast_dist_main)) {
								histo_t fast_dist_aux=get_fast_distance_aux(pos1,pos2);
								if (visit_aux(fast_dist_aux)) {
									histo_t dist_main=get_distance_main(fast_dist_main);
									histo_t weight=get_weight(weight1,&(box2.weight[dim_weight*iobj2]));
									legendre(fast_dist_aux,leg);
									size_t ill;
									for(ill=0;ill<n_ells;ill++) {
										histo_t w=weight*leg[ill];
										size_t ibin=ill+n_ells*get_bin_index(dist_main,bin_main);
										threadmeanmain[ibin]+=w*dist_main;
										threadcount[ibin]+=w;
									}
								}							
							}
						}							
					}
				}	
			}
		} // end omp for
#pragma omp critical
		{
			for (ibin=0;ibin<n_bin_tot;ibin++) {
				meanmain[ibin]+=threadmeanmain[ibin];
				count[ibin]+=threadcount[ibin];
			}

			free(threadmeanmain);
			free(threadcount);
		}
	} //end omp parallel
	for (ibin=0;ibin<n_bin_tot;ibin++) meanmain[ibin]/=count[ibin];
}

void cross_3pcf_multi(Catalog cat[], histo_t count[])
{
  //////
  // 3pcf cross
	
	set_fast_distance_main_limit();
	set_fast_distance_aux_limit();
	size_t n_bin_tot = bin_main.n_bin*bin_main.n_bin*n_ells*n_ells;
	size_t n_bin_sec = bin_main.n_bin*n_ells;
	size_t ibin,ibin2,ibin3,ill2,ill3;
	size_t max_secs = 2;
	histo_t leg[MAX_ELLS];
	
	Catalog cat1 = cat[0];
	if (cat[2].n_obj==0) max_secs = 1; 
	
	histo_t* countsec[2];
	//printf("nells %zu\n",n_ells);
	for (ibin=0;ibin<n_bin_tot;ibin++) count[ibin]=0.;
	histo_t* threadcount;
	
#pragma omp parallel default(none)				\
  shared(count,cat,cat1,max_secs,n_bin_tot,n_bin_sec,bin_main,dim_pos,dim_weight,n_ells) private(countsec,threadcount,leg)
	{
		size_t isec,ibin,iobj1;
		for (isec=0;isec<max_secs;isec++) countsec[isec] = (histo_t*) malloc(n_bin_sec*sizeof(histo_t));
		threadcount = (histo_t*) malloc(n_bin_tot*sizeof(histo_t));
		for(ibin=0;ibin<n_bin_tot;ibin++) threadcount[ibin]=0.;
		
#pragma omp for nowait schedule(dynamic)
		for (iobj1=0;iobj1<cat1.n_obj;iobj1++) {
			histo_t* pos1=&(cat1.pos[dim_pos*iobj1]);
			histo_t weight1=cat1.weight[dim_weight*iobj1];
			if (weight1==0.) continue;
			size_t isec;
			for (isec=0;isec<max_secs;isec++) {
				Catalog cat2 = cat[isec+1];
				histo_t *count2 = countsec[isec];
				size_t ibin2,iobj2;
				for (ibin2=0;ibin2<n_bin_sec;ibin2++) count2[ibin2]=0.;
				for (iobj2=0;iobj2<cat2.n_obj;iobj2++) {
					histo_t* pos2=&(cat2.pos[dim_pos*iobj2]);
					histo_t fast_dist_main = get_fast_distance_main(pos1,pos2);
					if (visit_main(fast_dist_main)) {
						histo_t fast_dist_aux = get_fast_distance_aux(pos1,pos2);
						if (visit_aux(fast_dist_aux)) {
							size_t ibin = get_bin_index(get_distance_main(fast_dist_main),bin_main);
							legendre(fast_dist_aux,leg);
							histo_t weight=cat2.weight[dim_weight*iobj2];
							size_t ill2;
							for (ill2=0;ill2<n_ells;ill2++) count2[ill2+ibin*n_ells] += leg[ill2]*weight;
						}
					}
				}
			}
			if (max_secs==1) {
				histo_t *count2 = countsec[0];
				histo_t *count3 = countsec[0];
				size_t ibin2;
				for (ibin2=0;ibin2<bin_main.n_bin;ibin2++) {
					if (count2[ibin2*n_ells]!=0.) {
						size_t ibin3;
						for (ibin3=0;ibin3<bin_main.n_bin;ibin3++) {
							if (count3[ibin3*n_ells]!=0.) {
								size_t ill2,ill3;
								for (ill2=0;ill2<n_ells;ill2++) {
									for (ill3=ill2;ill3<n_ells;ill3++) {
										threadcount[ill3+n_ells*(ill2+n_ells*(ibin3+bin_main.n_bin*ibin2))] += weight1*count2[ill2+ibin2*n_ells]*count3[ill3+ibin3*n_ells];
									}
								}
							}
						}
					}
				}
			}
			else {
				histo_t *count2 = countsec[0];
				histo_t *count3 = countsec[1];
				size_t ibin2;
				for (ibin2=0;ibin2<bin_main.n_bin;ibin2++) {
					if (count2[ibin2*n_ells]!=0.) {
						size_t ibin3;
						for (ibin3=0;ibin3<bin_main.n_bin;ibin3++) {
							if (count3[ibin3*n_ells]!=0.) {
								size_t ill2,ill3;
								for (ill2=0;ill2<n_ells;ill2++) {
									for (ill3=0;ill3<n_ells;ill3++) {
										threadcount[ill3+n_ells*(ill2+n_ells*(ibin3+bin_main.n_bin*ibin2))] += weight1*count2[ill2+ibin2*n_ells]*count3[ill3+ibin3*n_ells];
									}
								}
							}
						}
					}
				}
			}
		}
#pragma omp critical
		{
			for(ibin=0;ibin<n_bin_tot;ibin++) count[ibin]+=threadcount[ibin];
			free(threadcount);
			for (isec=0;isec<max_secs;isec++) free(countsec[isec]);
		}
	} //end omp parallel
	if (max_secs==1) {
		for (ibin2=0;ibin2<bin_main.n_bin;ibin2++) {
			for (ibin3=0;ibin3<bin_main.n_bin;ibin3++) {
				for (ill2=0;ill2<n_ells;ill2++) {
					for (ill3=0;ill3<ill2;ill3++) {
						count[ill3+n_ells*(ill2+n_ells*(ibin3+bin_main.n_bin*ibin2))] = count[ill2+n_ells*(ill3+n_ells*(ibin2+bin_main.n_bin*ibin3))];
					}
				}
			}
		}
	}
}


void cross_3pcf_multi_double_los(Catalog cat[], histo_t count[])
{
  //////
  // 3pcf cross
	
	set_fast_distance_main_limit();
	set_fast_distance_aux_limit();
	size_t n_bin_tot = bin_main.n_bin*bin_main.n_bin*n_ells*n_ells;
	size_t n_bin_sec = bin_main.n_bin*n_ells;
	size_t ibin;
	size_t max_secs = 2;
	histo_t leg[MAX_ELLS];
	
	Catalog cat1 = cat[0];
	if (cat[2].n_obj==0) max_secs = 1; 
	
	histo_t* countsec[2];
	//printf("nells %zu\n",n_ells);
	for (ibin=0;ibin<n_bin_tot;ibin++) count[ibin]=0.;
	histo_t* threadcount;
	
#pragma omp parallel default(none)				\
  shared(count,cat,cat1,max_secs,n_bin_tot,n_bin_sec,bin_main,dim_pos,dim_weight,n_ells) private(countsec,threadcount,leg)
	{
		size_t isec,ibin,iobj1;
		for (isec=0;isec<2;isec++) countsec[isec] = (histo_t*) malloc(n_bin_sec*sizeof(histo_t));
		threadcount = (histo_t*) malloc(n_bin_tot*sizeof(histo_t));
		for(ibin=0;ibin<n_bin_tot;ibin++) threadcount[ibin]=0.;
		
#pragma omp for nowait schedule(dynamic)
		for (iobj1=0;iobj1<cat1.n_obj;iobj1++) {
			histo_t* pos1=&(cat1.pos[dim_pos*iobj1]);
			histo_t weight1=cat1.weight[dim_weight*iobj1];
			if (weight1==0.) continue;
			size_t isec;
			for (isec=0;isec<2;isec++) {
				histo_t *count2 = countsec[isec];
				size_t ibin2;
				for (ibin2=0;ibin2<n_bin_sec;ibin2++) count2[ibin2]=0.;
			}
			for (isec=0;isec<max_secs;isec++) {
				Catalog cat2 = cat[isec+1];
				histo_t *count2 = countsec[isec];
				histo_t *count3 = countsec[1];
				size_t iobj2;
				for (iobj2=0;iobj2<cat2.n_obj;iobj2++) {
					histo_t* pos2=&(cat2.pos[dim_pos*iobj2]);
					histo_t fast_dist_main = get_fast_distance_main(pos1,pos2);
					if (visit_main(fast_dist_main)) {
						histo_t fast_dist_aux_1 = get_fast_distance_aux(pos1,pos2);
						histo_t fast_dist_aux_2 = (isec==0) ? get_fast_distance_aux(pos2,pos1) : fast_dist_aux_1;
						if (visit_aux(fast_dist_aux_1) && visit_aux(fast_dist_aux_2)) {
							size_t ibin = get_bin_index(get_distance_main(fast_dist_main),bin_main);
							histo_t weight=cat2.weight[dim_weight*iobj2];
							legendre(fast_dist_aux_1,leg);
							size_t ill2;
							for (ill2=0;ill2<n_ells;ill2++) count2[ill2+ibin*n_ells] += leg[ill2]*weight;
							if (max_secs==1) {
								for (ill2=0;ill2<n_ells;ill2++) count3[ill2+ibin*n_ells] += leg[ill2]*weight;
							}
							if (isec==0) {
								legendre(fast_dist_aux_2,leg);
								for (ill2=0;ill2<n_ells;ill2++) count2[ill2+ibin*n_ells] += leg[ill2]*weight;
							}
						}
					}
				}
			}
			histo_t *count2 = countsec[0];
			histo_t *count3 = countsec[1];
			size_t ibin2;
			for (ibin2=0;ibin2<bin_main.n_bin;ibin2++) {
				if (count2[ibin2*n_ells]!=0.) {
					size_t ibin3;
					for (ibin3=0;ibin3<bin_main.n_bin;ibin3++) {
						if (count3[ibin3*n_ells]!=0.) {
							size_t ill2,ill3;
							for (ill2=0;ill2<n_ells;ill2++) {
								for (ill3=0;ill3<n_ells;ill3++) {
									threadcount[ill3+n_ells*(ill2+n_ells*(ibin3+bin_main.n_bin*ibin2))] += weight1*count2[ill2+ibin2*n_ells]*count3[ill3+ibin3*n_ells];
								}
							}
						}
					}
				}
			}
		}
#pragma omp critical
		{
			for(ibin=0;ibin<n_bin_tot;ibin++) count[ibin]+=threadcount[ibin];
			free(threadcount);
			for (isec=0;isec<2;isec++) free(countsec[isec]);
		}
	} //end omp parallel
}



void cross_3pcf_multi_radial(Catalog cat[], histo_t count[],_Bool normalize)
{
  //////
  // 3pcf cross
	
	set_fast_distance_main_limit();
	set_fast_distance_aux_limit();
	set_fast_distance_radial_limit();
	size_t n_bin_tot = bin_main.n_bin*bin_main.n_bin*n_ells*n_ells;
	size_t n_bin_sec = bin_main.n_bin*bin_radial.n_bin*n_ells;
	size_t ibin,ibin2,ibin3,ill2,ill3;
	size_t max_secs = 2;
	histo_t leg[MAX_ELLS];
	
	Catalog cat1 = cat[0];
	if (cat[2].n_obj==0) max_secs = 1; 
	
	histo_t* countsec[2];
	//printf("nells %zu\n",n_ells);
	for (ibin=0;ibin<n_bin_tot;ibin++) count[ibin]=0.;
	histo_t* threadcount;
	
#pragma omp parallel default(none)				\
  shared(count,cat,cat1,max_secs,n_bin_tot,n_bin_sec,bin_main,bin_radial,dim_pos,dim_weight,n_ells,normalize) private(countsec,threadcount,leg)
	{
		size_t isec,ibin,iobj1;
		for (isec=0;isec<max_secs;isec++) countsec[isec] = (histo_t*) malloc(n_bin_sec*sizeof(histo_t));
		threadcount = (histo_t*) malloc(n_bin_tot*sizeof(histo_t));
		for(ibin=0;ibin<n_bin_tot;ibin++) threadcount[ibin]=0.;
		
#pragma omp for nowait schedule(dynamic)
		for (iobj1=0;iobj1<cat1.n_obj;iobj1++) {
			histo_t* pos1=&(cat1.pos[dim_pos*iobj1]);
			histo_t weight1=cat1.weight[dim_weight*iobj1];
			size_t isec;
			for (isec=0;isec<max_secs;isec++) {
				Catalog cat2 = cat[isec+1];
				histo_t *count2 = countsec[isec];
				size_t ibin2;
				for (ibin2=0;ibin2<n_bin_sec;ibin2++) count2[ibin2]=0.;
				if ((isec==2) && (normalize)) {
					histo_t *norm2 = (histo_t*) malloc(bin_radial.n_bin*sizeof(histo_t));
					size_t ibin;
					for (ibin=0;ibin<bin_radial.n_bin;ibin++) norm2[ibin] = 0.;
					histo_t *fast_aux = (histo_t*) malloc(cat2.n_obj*sizeof(histo_t));
					size_t *ibinmain = (size_t*) malloc(cat2.n_obj*sizeof(size_t));
					size_t *ibinradial = (size_t*) malloc(cat2.n_obj*sizeof(size_t));
					size_t iobj2;
					for (iobj2=0;iobj2<cat2.n_obj;iobj2++) {
						fast_aux[iobj2] = -2.;
						ibinmain[iobj2] = 0;
						ibinradial[iobj2] = 0;
						histo_t* pos2=&(cat2.pos[dim_pos*iobj2]);
						histo_t fast_dist_main = get_fast_distance_main(pos1,pos2);
						histo_t fast_dist_aux = get_fast_distance_aux(pos1,pos2);
						//histo_t fast_dist_radial = fast_dist_main*fast_dist_aux;
						histo_t fast_dist_radial = get_fast_distance_radial(pos2);
						if (visit_radial(fast_dist_radial) && visit_main(fast_dist_main) && visit_aux(fast_dist_aux)) {
							fast_aux[iobj2] = fast_dist_aux;
							ibinmain[iobj2] = get_bin_index(get_distance_main(fast_dist_main),bin_main);
							ibinradial[iobj2] = get_bin_index(get_distance_radial(fast_dist_radial),bin_radial);
							norm2[ibinradial[iobj2]] += cat2.weight[dim_weight*iobj2];
						}
					}
					for (iobj2=0;iobj2<cat2.n_obj;iobj2++) {
						if (fast_aux[iobj2]>-2.) {
							ibin = ibinradial[iobj2]+bin_radial.n_bin*ibinmain[iobj2];
							legendre(fast_aux[iobj2],leg);
							histo_t weight=cat2.weight[dim_weight*iobj2]/norm2[ibinradial[iobj2]];
							size_t ill2;
							for (ill2=0;ill2<n_ells;ill2++) count2[ill2+ibin*n_ells] += leg[ill2]*weight;
						}
					}
					free(fast_aux);
					free(ibinmain);
					free(ibinradial);
				}
				else {
					size_t iobj2;		
					for (iobj2=0;iobj2<cat2.n_obj;iobj2++) {
						histo_t* pos2=&(cat2.pos[dim_pos*iobj2]);
						histo_t fast_dist_main = get_fast_distance_main(pos1,pos2);
						histo_t fast_dist_aux = get_fast_distance_aux(pos1,pos2);
						//histo_t fast_dist_radial = fast_dist_main*fast_dist_aux;
						histo_t fast_dist_radial = get_fast_distance_radial(pos2);
						if (visit_radial(fast_dist_radial) && visit_main(fast_dist_main) && visit_aux(fast_dist_aux)) {
							//printf("%zu ",get_bin_index(get_distance_aux(fast_dist_radial),bin_radial));
							ibin = get_bin_index(get_distance_radial(fast_dist_radial),bin_radial)+bin_radial.n_bin*get_bin_index(get_distance_main(fast_dist_main),bin_main);
							legendre(fast_dist_aux,leg);
							histo_t weight=cat2.weight[dim_weight*iobj2];
							size_t ill2;
							for (ill2=0;ill2<n_ells;ill2++) count2[ill2+ibin*n_ells] += leg[ill2]*weight;
						}
					}
				}
			}
			if (max_secs==1) {
				histo_t *count2 = countsec[0];
				histo_t *count3 = countsec[0];
				size_t ibin;
				for (ibin=0;ibin<bin_radial.n_bin;ibin++) {
					size_t ibin2;
					for (ibin2=0;ibin2<bin_main.n_bin;ibin2++) {
						size_t ill2;
						for (ill2=0;ill2<n_ells;ill2++) {
							histo_t tmp2 = weight1*count2[ill2+n_ells*(ibin+bin_radial.n_bin*ibin2)];
							if (tmp2!=0.) {
								size_t ibin3;
								for (ibin3=0;ibin3<bin_main.n_bin;ibin3++) {
									size_t ill3;
									for (ill3=ill2;ill3<n_ells;ill3++) {
										threadcount[ill3+n_ells*(ill2+n_ells*(ibin3+bin_main.n_bin*ibin2))] += tmp2*count3[ill3+n_ells*(ibin+bin_radial.n_bin*ibin3)];
									}
								}
							}
						}
					}
				}
			}
			else {
				histo_t *count2 = countsec[0];
				histo_t *count3 = countsec[1];
				size_t ibin;
				for (ibin=0;ibin<bin_radial.n_bin;ibin++) {
					size_t ibin2;
					for (ibin2=0;ibin2<bin_main.n_bin;ibin2++) {
						size_t ill2;
						for (ill2=0;ill2<n_ells;ill2++) {
							histo_t tmp2 = weight1*count2[ill2+n_ells*(ibin+bin_radial.n_bin*ibin2)];
							if (tmp2!=0.) {
								size_t ibin3;
								for (ibin3=0;ibin3<bin_main.n_bin;ibin3++) {
									size_t ill3;
									for (ill3=0;ill3<n_ells;ill3++) {
										threadcount[ill3+n_ells*(ill2+n_ells*(ibin3+bin_main.n_bin*ibin2))] += tmp2*count3[ill3+n_ells*(ibin+bin_radial.n_bin*ibin3)];
									}
								}
							}
						}
					}
				}
			}
		}
#pragma omp critical
		{
			for(ibin=0;ibin<n_bin_tot;ibin++) count[ibin]+=threadcount[ibin];
			free(threadcount);
			for (isec=0;isec<max_secs;isec++) free(countsec[isec]);
		}
	} //end omp parallel
	if (max_secs==1) {
		for (ibin2=0;ibin2<bin_main.n_bin;ibin2++) {
			for (ibin3=0;ibin3<bin_main.n_bin;ibin3++) {
				for (ill2=0;ill2<n_ells;ill2++) {
					for (ill3=0;ill3<ill2;ill3++) {
						count[ill3+n_ells*(ill2+n_ells*(ibin3+bin_main.n_bin*ibin2))] = count[ill2+n_ells*(ill3+n_ells*(ibin2+bin_main.n_bin*ibin3))];
					}
				}
			}
		}
	}
}

void cross_2pcf_multi_radial(Catalog cat[], histo_t count[],_Bool normalize)
{
  //////
  // 2pcf cross, aux is the s_\parallel
	
	set_fast_distance_main_limit();
	set_fast_distance_aux_limit();
	set_fast_distance_radial_limit();
	size_t n_bin_tot = bin_main.n_bin*bin_radial.n_bin*n_ells;
	size_t ibin;
	histo_t leg[MAX_ELLS];
	histo_t norm1=1.;
	histo_t *norm2 = (histo_t*) malloc(bin_radial.n_bin*sizeof(histo_t));
	
	Catalog cat1 = cat[0];
	Catalog cat2 = cat[1];
	//printf("%zu \n",cat1.n_obj);
	
	if (normalize) {
		size_t iobj;
		norm1=0;
		for (iobj=0;iobj<cat1.n_obj;iobj++) norm1 += cat1.weight[dim_weight*iobj];
		for (ibin=0;ibin<bin_radial.n_bin;ibin++) norm2[ibin] = 0.;
		for (iobj=0;iobj<cat2.n_obj;iobj++) {
			histo_t fast_dist_radial = get_fast_distance_radial(&(cat2.pos[dim_pos*iobj]));
			if (visit_radial(fast_dist_radial)) {
				size_t ibinradial = get_bin_index(get_distance_radial(fast_dist_radial),bin_radial);
				norm2[ibinradial] += cat2.weight[dim_weight*iobj];
			}
		}
	}

	for (ibin=0;ibin<n_bin_tot;ibin++) count[ibin]=0.;
	histo_t* threadcount;
	
#pragma omp parallel default(none)				\
  shared(count,cat1,cat2,n_bin_tot,bin_main,bin_radial,dim_pos,dim_weight,n_ells,normalize,norm1,norm2) private(threadcount,leg)
	{
		size_t ibin,iobj2;
		threadcount = (histo_t*) malloc(n_bin_tot*sizeof(histo_t));
		for (ibin=0;ibin<n_bin_tot;ibin++) threadcount[ibin]=0.;
		
#pragma omp for nowait schedule(dynamic)
		for (iobj2=0;iobj2<cat2.n_obj;iobj2++) {
			histo_t* pos2=&(cat2.pos[dim_pos*iobj2]);
			histo_t* weight2=&(cat2.weight[dim_weight*iobj2]);
			histo_t fast_dist_radial = get_fast_distance_radial(pos2);
			if (visit_radial(fast_dist_radial)) {
				size_t ibinradial = get_bin_index(get_distance_radial(fast_dist_radial),bin_radial);
				size_t iobj1;
				for (iobj1=0;iobj1<cat1.n_obj;iobj1++) {
					histo_t* pos1=&(cat1.pos[dim_pos*iobj1]);
					histo_t fast_dist_main = get_fast_distance_main(pos1,pos2);
					histo_t fast_dist_aux = get_fast_distance_aux(pos1,pos2);
					if (visit_main(fast_dist_main) && visit_aux(fast_dist_aux)) {
						ibin = ibinradial+bin_radial.n_bin*get_bin_index(get_distance_main(fast_dist_main),bin_main);
						legendre(fast_dist_aux,leg);
						histo_t weight=get_weight(&(cat1.weight[dim_weight*iobj1]),weight2);
						if (normalize) weight /= norm1*norm2[ibinradial];
						size_t ill;
						for (ill=0;ill<n_ells;ill++) threadcount[ill+ibin*n_ells] += leg[ill]*weight;
					}
				}
			}
		}
#pragma omp critical
		{
			for(ibin=0;ibin<n_bin_tot;ibin++) count[ibin]+=threadcount[ibin];
			free(threadcount);
		}
	} //end omp parallel
	free(norm2);
}

void cross_4pcf_multi_radial(Catalog cat[], histo_t count[], _Bool normalize)
{
  //////
  // 4pcf cross, aux is the s_\parallel
    size_t n_bin_main = bin_main.n_bin;
    size_t n_bin_radial = bin_radial.n_bin;
    size_t n_bin_2pcf = n_bin_main*n_bin_radial*n_ells;
    size_t n_bin_tot = n_bin_main*n_bin_main*n_ells*n_ells;
  	histo_t *count1 = (histo_t*) malloc(n_bin_2pcf*sizeof(histo_t));
  	cross_2pcf_multi_radial(cat,count1,0);
  	histo_t *count2 = (histo_t*) malloc(n_bin_2pcf*sizeof(histo_t));
  	cross_2pcf_multi_radial(&(cat[2]),count2,normalize);
	
	size_t ibin;
	for (ibin=0;ibin<n_bin_tot;ibin++) count[ibin]=0.;
	histo_t* threadcount;
	
#pragma omp parallel default(none)				\
  shared(count,count1,count2,n_bin_tot,n_bin_main,n_bin_radial,n_ells) private(threadcount)
	{
		size_t ibin;
		threadcount = (histo_t*) malloc(n_bin_tot*sizeof(histo_t));
		for(ibin=0;ibin<n_bin_tot;ibin++) threadcount[ibin]=0.;
		
#pragma omp for nowait schedule(dynamic)
		for (ibin=0;ibin<n_bin_radial;ibin++) {
			size_t ibin1;
			for (ibin1=0;ibin1<n_bin_main;ibin1++) {
				size_t ill1;
				for (ill1=0;ill1<n_ells;ill1++) {
					histo_t tmp1 = count1[ill1+(ibin+ibin1*n_bin_radial)*n_ells];
					if (tmp1!=0.) {
						size_t ibin2;
						for (ibin2=0;ibin2<n_bin_main;ibin2++) {
							size_t ill2;
							for (ill2=0;ill2<n_ells;ill2++) {
								threadcount[ill2+n_ells*(ill1+n_ells*(ibin2+n_bin_main*ibin1))] += tmp1*count2[ill2+(ibin+ibin2*n_bin_radial)*n_ells];
							}
						}
					}
				}
			}
		}
#pragma omp critical
		{
			for(ibin=0;ibin<n_bin_tot;ibin++) count[ibin]+=threadcount[ibin];
			free(threadcount);
		}
	} //end omp parallel
	free(count1);
	free(count2);
}

void cross_2pcf_multi_cat(Catalog cat[], histo_t count[])
{
  //////
  // 2pcf cross, aux is the s_\parallel
	
	set_fast_distance_main_limit();
	set_fast_distance_aux_limit();
	size_t n_bin_tot = bin_main.n_bin*n_ells;
	size_t ibin;
	histo_t leg[MAX_ELLS];
	
	Catalog cat1 = cat[0];
	Catalog cat2 = cat[1];
	//printf("%zu \n",cat1.n_obj);

	for (ibin=0;ibin<n_bin_tot;ibin++) count[ibin]=0.;
	histo_t* threadcount;
	
#pragma omp parallel default(none)				\
  shared(count,cat1,cat2,n_bin_tot,bin_main,bin_aux,dim_pos,dim_weight,n_ells) private(threadcount,leg)
	{
		size_t ibin,iobj1;
		threadcount = (histo_t*) malloc(n_bin_tot*sizeof(histo_t));
		for (ibin=0;ibin<n_bin_tot;ibin++) threadcount[ibin]=0.;
		
#pragma omp for nowait schedule(dynamic)
		for (iobj1=0;iobj1<cat1.n_obj;iobj1++) {
			histo_t* pos1=&(cat1.pos[dim_pos*iobj1]);
			histo_t* weight1=&(cat1.weight[dim_weight*iobj1]);
			size_t iobj2;
			for (iobj2=0;iobj2<cat2.n_obj;iobj2++) {
				histo_t* pos2=&(cat2.pos[dim_pos*iobj2]);
				histo_t fast_dist_main = get_fast_distance_main(pos1,pos2);
				histo_t fast_dist_aux = get_fast_distance_aux(pos1,pos2);
				if (visit_main(fast_dist_main) && visit_aux(fast_dist_aux)) {
					ibin = get_bin_index(get_distance_main(fast_dist_main),bin_main);
					legendre(fast_dist_aux,leg);
					histo_t weight=get_weight(weight1,&(cat2.weight[dim_weight*iobj2]));
					size_t ill;
					for (ill=0;ill<n_ells;ill++) threadcount[ill+ibin*n_ells] += leg[ill]*weight;
				}
			}
		}
#pragma omp critical
		{
			for(ibin=0;ibin<n_bin_tot;ibin++) count[ibin]+=threadcount[ibin];
			free(threadcount);
		}
	} //end omp parallel
}



void cross_4pcf_multi(Catalog cat[], histo_t count[])
{
  //////
  // 4pcf cross, aux is the s_\parallel
    
    size_t n_bin_main = bin_main.n_bin;
    size_t n_bin_2pcf = n_bin_main*n_ells;
    size_t n_bin_tot = n_bin_2pcf*n_bin_2pcf;
    
  	histo_t *count1 = (histo_t*) malloc(n_bin_2pcf*sizeof(histo_t));
  	cross_2pcf_multi_cat(cat,count1);
  	histo_t *count2 = (histo_t*) malloc(n_bin_2pcf*sizeof(histo_t));
  	cross_2pcf_multi_cat(&(cat[2]),count2);
	
	size_t ibin;
	for (ibin=0;ibin<n_bin_tot;ibin++) count[ibin]=0.;
	histo_t* threadcount;
	
#pragma omp parallel default(none)				\
  shared(count,count1,count2,n_bin_tot,n_bin_main,n_ells) private(threadcount)
	{
		size_t ibin,ibin1;
		threadcount = (histo_t*) malloc(n_bin_tot*sizeof(histo_t));
		for(ibin=0;ibin<n_bin_tot;ibin++) threadcount[ibin]=0.;
		
#pragma omp for nowait schedule(dynamic)
		for (ibin1=0;ibin1<n_bin_main;ibin1++) {
			size_t ill1;
			for (ill1=0;ill1<n_ells;ill1++) {
				histo_t tmp1 = count1[ill1+ibin1*n_ells];
				if (tmp1!=0.) {
					size_t ibin2;
					for (ibin2=0;ibin2<n_bin_main;ibin2++) {
						size_t ill2;
						for (ill2=0;ill2<n_ells;ill2++) {
							threadcount[ill2+n_ells*(ill1+n_ells*(ibin2+n_bin_main*ibin1))] += tmp1*count2[ill2+ibin2*n_ells];
						}
					}
				}
			}
		}
#pragma omp critical
		{
			for(ibin=0;ibin<n_bin_tot;ibin++) {
				count[ibin]+=threadcount[ibin];
			}
			free(threadcount);
		}
	} //end omp parallel
	free(count1);
	free(count2);
}

/*
void cross_2pcf_multi_radial(Catalog cat[], histo_t count[],_Bool normalize)
{
  //////
  // 2pcf cross, aux is the s_\parallel
	
	set_fast_distance_main_limit();
	set_fast_distance_aux_limit();
	set_fast_distance_radial_limit();
	size_t n_bin_tot = bin_main.n_bin*bin_radial.n_bin*n_ells;
	size_t ibin;
	histo_t leg[MAX_ELLS];
	histo_t norm1=1.;
	
	Catalog cat1 = cat[0];
	Catalog cat2 = cat[1];
	//printf("%zu \n",cat1.n_obj);
	
	if (normalize) {
		size_t iobj;
		norm1=0;
		for (iobj=0;iobj<cat1.n_obj;iobj++) {
			norm1 += cat1.weight[dim_weight*iobj];
		}
	}

	for (ibin=0;ibin<n_bin_tot;ibin++) count[ibin]=0.;
	histo_t* threadcount;
	
#pragma omp parallel default(none)				\
  shared(count,cat1,cat2,n_bin_tot,bin_main,bin_radial,dim_pos,dim_weight,n_ells,normalize,norm1) private(threadcount,leg)
	{
		size_t ibin,iobj1,iobj2,ill;
		threadcount = (histo_t*) malloc(n_bin_tot*sizeof(histo_t));
		for (ibin=0;ibin<n_bin_tot;ibin++) {
			threadcount[ibin]=0.;
		}
		
#pragma omp for nowait schedule(dynamic)
		for (iobj1=0;iobj1<cat1.n_obj;iobj1++) {
			histo_t* pos1=&(cat1.pos[dim_pos*iobj1]);
			histo_t* weight1=&(cat1.weight[dim_weight*iobj1]);
			if (normalize) {
				histo_t *norm2 = (histo_t*) malloc(bin_radial.n_bin*sizeof(histo_t));
				for (ibin=0;ibin<bin_radial.n_bin;ibin++) {
					norm2[ibin] = 0.;
				}
				histo_t *fast_aux = (histo_t*) malloc(cat2.n_obj*sizeof(histo_t));
				size_t *ibinmain = (size_t*) malloc(cat2.n_obj*sizeof(size_t));
				size_t *ibinradial = (size_t*) malloc(cat2.n_obj*sizeof(size_t));
				for (iobj2=0;iobj2<cat2.n_obj;iobj2++) {
					fast_aux[iobj2] = -2.;
					ibinmain[iobj2] = 0;
					ibinradial[iobj2] = 0;
					histo_t* pos2=&(cat2.pos[dim_pos*iobj2]);
					histo_t fast_dist_main = get_fast_distance_main(pos1,pos2);
					histo_t fast_dist_aux = get_fast_distance_aux(pos1,pos2);
					//histo_t fast_dist_radial = fast_dist_main*fast_dist_aux;
					histo_t fast_dist_radial = get_fast_distance_radial(pos2);
					if (visit_radial(fast_dist_radial) && visit_main(fast_dist_main) && visit_aux(fast_dist_aux)) {
						fast_aux[iobj2] = fast_dist_aux;
						ibinmain[iobj2] = get_bin_index(get_distance_main(fast_dist_main),bin_main);
						ibinradial[iobj2] = get_bin_index(get_distance_radial(fast_dist_radial),bin_radial);
						norm2[ibinradial[iobj2]] += cat2.weight[dim_weight*iobj2];
					}
				}
				for (iobj2=0;iobj2<cat2.n_obj;iobj2++) {
					if (fast_aux[iobj2]>-2.) {
						ibin = ibinradial[iobj2]+bin_radial.n_bin*ibinmain[iobj2];
						legendre(fast_aux[iobj2],leg);
						histo_t weight=get_weight(weight1,&(cat2.weight[dim_weight*iobj2]))/norm2[ibinradial[iobj2]]/norm1;
						for (ill=0;ill<n_ells;ill++) threadcount[ill+ibin*n_ells] += leg[ill]*weight;
					}
				}
				free(fast_aux);
				free(ibinmain);
				free(ibinradial);
			}
			else {
				for (iobj2=0;iobj2<cat2.n_obj;iobj2++) {
					histo_t* pos2=&(cat2.pos[dim_pos*iobj2]);
					histo_t fast_dist_main = get_fast_distance_main(pos1,pos2);
					histo_t fast_dist_aux = get_fast_distance_aux(pos1,pos2);
					//histo_t fast_dist_radial = fast_dist_main*fast_dist_aux;
					histo_t fast_dist_radial = get_fast_distance_radial(pos2);
					if (visit_radial(fast_dist_radial) && visit_main(fast_dist_main) && visit_aux(fast_dist_aux)) {
						ibin = get_bin_index(get_distance_radial(fast_dist_radial),bin_radial)+bin_radial.n_bin*get_bin_index(get_distance_main(fast_dist_main),bin_main);
						legendre(fast_dist_aux,leg);
						histo_t weight=get_weight(weight1,&(cat2.weight[dim_weight*iobj2]));
						for (ill=0;ill<n_ells;ill++) threadcount[ill+ibin*n_ells] += leg[ill]*weight;
					}
				}
			}
		}
#pragma omp critical
		{
			for(ibin=0;ibin<n_bin_tot;ibin++) {
				count[ibin]+=threadcount[ibin];
			}
			free(threadcount);
		}
	} //end omp parallel
}
*/
