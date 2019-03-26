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
//                             3D boxes                              //
/*********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "define.h"
#include "common.h"

#define NUMBER_BOXES 8
#define FRACTION_EXTEND 0.01

static histo_t *bound_min;
static histo_t *bound_max;
static histo_t *l_box;
static size_t n_boxes;
static size_t *n_side;
static size_t n_obj_tot;

static histo_t my_sqrt(histo_t x)
{
#ifdef _FLOAT32
	return sqrtf(x);
#else
	return sqrt(x);
#endif	
}

void cat2box(histo_t cat[],histo_t box[])
{
	//////
	// From cat to box

	if ((corr_type==CORR_SMU)||(corr_type==CORR_SCOS)) {
		size_t idim;
		for(idim=0;idim<dim_box;idim++) box[idim]=cat[idim];
	}
	else if (corr_type==CORR_ANGULAR) {
		box[0]=cos(cat[0]);
		box[1]=cat[1];
	}
}

void box2main(histo_t box[],histo_t main[])
{
	//////
	// From cat to box

	if ((corr_type==CORR_SMU)||(corr_type==CORR_SCOS)) {
		size_t idim;
		for(idim=0;idim<dim_box;idim++) main[idim]=box[idim];
	}
	else if (corr_type==CORR_ANGULAR) {
		main[0]=acos(box[0]);
		main[1]=box[1];
	}
}


void box2pos(histo_t box[],histo_t pos[])
{
	//////
	// From box to pos

	if ((corr_type==CORR_SMU)||(corr_type==CORR_SCOS)) {
		size_t idim;
		for(idim=0;idim<dim_pos;idim++) pos[idim]=box[idim];
	}
	else if (corr_type==CORR_ANGULAR) {
		histo_t sin_theta=my_sqrt(1.-box[0]*box[0]);
		pos[2]=sin_theta*sin(box[1]);
		pos[1]=sin_theta*cos(box[1]);
		pos[0]=box[0];
	}
}


void optimal_nside()
{
	//////
	// Estimates a good candidate for the size
	// of a set of neighbor boxes

	size_t idim=0;
	n_boxes=1;
	size_t max_boxes=(size_t)pow((double)n_obj_tot,1./dim_box);
	
	n_side=(size_t *) malloc(dim_box*sizeof(size_t));
	histo_t* bound_min_main=(histo_t *) malloc(dim_box*sizeof(histo_t));
	histo_t* bound_max_main=(histo_t *) malloc(dim_box*sizeof(histo_t));
	box2main(bound_min,bound_min_main);
	box2main(bound_max,bound_max_main);

	for (idim=0;idim<dim_box;idim++) {
	if (verbose == DEBUG) {
		printf(" - box size %.3f %.3f -> %.3f %.3f; max main %.3f\n",bound_min[idim],bound_max[idim],bound_min_main[idim],bound_max_main[idim],bin_main.max);
	}
		n_side[idim]=(size_t)MIN((NUMBER_BOXES*fabs((bound_max_main[idim]-bound_min_main[idim])/bin_main.max))+1,max_boxes);
		n_boxes*=n_side[idim];
	}

	free(bound_min_main);
	free(bound_max_main);

}

void ind2box(size_t index[],histo_t box[])
{
	//////
	// Returns box position

	size_t idim=0;
	for(idim=0;idim<dim_box;idim++) box[idim]=(histo_t)(index[idim]*l_box[idim]/n_side[idim]+bound_min[idim]);

}

void box2ind(histo_t box[],size_t index[])
{
	//////
	// Returns box index

	size_t idim;
	for(idim=0;idim<dim_box;idim++) index[idim]=(size_t)(floor((box[idim]-bound_min[idim])/l_box[idim]*n_side[idim]));

}

void ind2pos(size_t index[],histo_t pos[])
{
	ind2box(index,pos);
	box2pos(pos,pos);
}


void free_mesh(Mesh mesh)
{
	size_t ibox;
	for(ibox=0;ibox<mesh.n_boxes;ibox++) {
		if(mesh.boxes[ibox].n_obj>0) {
			free(mesh.boxes[ibox].pos);
			free(mesh.boxes[ibox].weight);
			free(mesh.boxes[ibox].bin);
			free(mesh.boxes[ibox].index);
			free(mesh.boxes[ibox].visit_min);
			free(mesh.boxes[ibox].visit_max);
		}
	}

	free(mesh.boxes);
}

void free_params()
{
	free(bound_min);
	free(bound_max);
	free(n_side);
	free(l_box);
}


Box* init_boxes(size_t n_boxes)
{
	size_t ibox;
	Box* boxes = (Box *) malloc(n_boxes*sizeof(Box));
	if (boxes==NULL) error_mem_out();

	for(ibox=0;ibox<n_boxes;ibox++) {
		boxes[ibox].n_obj=0;
		boxes[ibox].pos=NULL;
		boxes[ibox].weight=NULL;
		boxes[ibox].bin=NULL;
		boxes[ibox].index=NULL;
		boxes[ibox].visit_min=NULL;
		boxes[ibox].visit_max=NULL;
	}

	return boxes;
}

void set_catalog_box(Catalog *cat)
{
	size_t iobj;
	size_t n_tot=(*cat).n_obj*dim_box;	
	(*cat).box=(histo_t *) malloc(n_tot*sizeof(histo_t));

	for(iobj=0;iobj<n_tot;iobj+=dim_box) {
		cat2box(&((*cat).pos[iobj]),&((*cat).box[iobj]));
	}
}

void print_boxes()
{
	size_t idim;
	printf(" - there will be (");
	for (idim=0;idim<dim_box;idim++) printf("%zu,",n_side[idim]);
	printf(")=%zu boxes in total\n",n_boxes);
	printf(" - full dimensions are (");
	for (idim=0;idim<dim_box;idim++) printf("%.3f,",l_box[idim]);
	printf(")\n");
	printf(" - box dimensions are (");
	for (idim=0;idim<dim_box;idim++) printf("%.3f,",l_box[idim]/n_side[idim]);
	printf(")\n");
}

void init_params(Catalog *cats,size_t n_cats)
{
	size_t icat,idim,iobj;
	n_obj_tot = 0;
	if ((corr_type==CORR_SMU)||(corr_type==CORR_SCOS)) dim_pos=dim_box;
	if (corr_type==CORR_ANGULAR) dim_pos=dim_box+1;

	bound_min = (histo_t *) malloc(dim_box*sizeof(histo_t));
	bound_max = (histo_t *) malloc(dim_box*sizeof(histo_t));
	l_box = (histo_t *) malloc(dim_box*sizeof(histo_t));

	for (idim=0;idim<dim_box;idim++) {
		bound_min[idim]=cats[0].box[idim];
		bound_max[idim]=cats[0].box[idim];
	}
	for (icat=0;icat<n_cats;icat++) {
		Catalog cat = cats[icat];
		n_obj_tot += cat.n_obj;
		for(iobj=0;iobj<cat.n_obj;iobj++) {
			for (idim=0;idim<dim_box;idim++) {
				bound_min[idim] = MIN(bound_min[idim],cat.box[iobj*dim_box+idim]);
				bound_max[idim] = MAX(bound_max[idim],cat.box[iobj*dim_box+idim]);
			}
		}
	}
	
	for (idim=0;idim<dim_box;idim++) {
		histo_t margin=FRACTION_EXTEND*(bound_max[idim]-bound_min[idim]);
		bound_max[idim]+=margin;
		bound_min[idim]-=margin;
	}
	if (corr_type==CORR_ANGULAR) {
		bound_min[0] = CLAMP(bound_min[0],-1.,1.);
		bound_max[0] = CLAMP(bound_max[0],-1.,1.);
		//bound_max[1] = CLAMP(bound_max[1],0,2*M_PI);
	}
	
	l_box=(histo_t *) malloc(dim_box*sizeof(histo_t));

	for (idim=0;idim<dim_box;idim++) {
		l_box[idim]=bound_max[idim]-bound_min[idim];
	}
	optimal_nside();
	if (verbose == INFO) print_boxes();
}

static void set_visit(Box box,size_t ibox)
{
	size_t idim;
	unravel_index(ibox,n_side,dim_box,box.index);

	if ((corr_type==CORR_SMU)||(corr_type==CORR_SCOS)) {
		for(idim=0;idim<dim_box;idim++) {
			long nbox=(long)(floor(bin_main.max/l_box[idim]*n_side[idim])+1);
			box.visit_min[idim]=(size_t)MAX((long)box.index[idim]-nbox,0); //long for diff
			box.visit_max[idim]=(size_t)MIN((long)box.index[idim]+nbox,n_side[idim]-1);
		}
	}

	else if (corr_type==CORR_ANGULAR) {
		histo_t costheta_min,costheta_max;
		histo_t* thetaphi_high=(histo_t *) malloc(dim_box*sizeof(histo_t));
		histo_t* thetaphi_low=(histo_t *) malloc(dim_box*sizeof(histo_t));
		ind2box(box.index,thetaphi_low);
		for(idim=0;idim<dim_box;idim++) thetaphi_high[idim]=thetaphi_low[idim]+l_box[idim]/n_side[idim];
		box2main(thetaphi_low,thetaphi_low);
		box2main(thetaphi_high,thetaphi_high); //thetaphi_low[0] is theta min (up), thetaphi_high[0] is theta max (down)
	
		if(thetaphi_low[0]>M_PI-bin_main.max) { //Near the South Pole
			costheta_min=-1;
			costheta_max=cos(thetaphi_high[0]-bin_main.max);
			box.visit_min[1]=0;
			box.visit_max[1]=n_side[1]-1;
			//printf("South ");
		}
		else if(thetaphi_high[0]<bin_main.max) { //Near the North Pole
			costheta_min=cos(thetaphi_low[0]+bin_main.max);
			costheta_max=1;
			box.visit_min[1]=0;
			box.visit_max[1]=n_side[1]-1;
			//printf("North ");
		}
		else {
			costheta_min=cos(thetaphi_low[0]+bin_main.max);
    		costheta_max=cos(thetaphi_high[0]-bin_main.max);
			histo_t costheta_bound,dphi;
			histo_t cosalpha=cos(bin_main.max);
			if ((thetaphi_low[0]+thetaphi_high[0])/2.<M_PI/2.) costheta_bound=cos(thetaphi_high[0]); //North
			else costheta_bound=cos(thetaphi_low[0]); //South
			if(costheta_bound>=cosalpha) dphi=M_PI;
      		else dphi=acos(sqrt((cosalpha*cosalpha-costheta_bound*costheta_bound)/(1.-costheta_bound*costheta_bound)));
			//printf("Middle ");
    		if(dphi<M_PI) {
				histo_t phi_min=wrap_phi(thetaphi_low[1]-dphi);
				histo_t phi_max=wrap_phi(thetaphi_high[1]+dphi);
      			box.visit_min[1]=(size_t)(floor((phi_min-bound_min[1])/l_box[1]*n_side[1]));
      			box.visit_max[1]=(size_t)(floor((phi_max-bound_min[1])/l_box[1]*n_side[1]));
				if ((phi_min>bound_max[1]) || (phi_min<bound_min[1])) box.visit_min[1]=0;
				if ((phi_max>bound_max[1]) || (phi_max<bound_min[1])) box.visit_max[1]=n_side[1]-1;
				//printf("%zu %zu %zu  ",box.index[1],box.visit_min[1],box.visit_max[1]);
				//printf("Middle < pi ");
				//box.visit_min[1]=0;
				//box.visit_max[1]=n_side[1]-1;
    		}
    		else {
				box.visit_min[1]=0;
				box.visit_max[1]=n_side[1]-1;
				//printf("Middle > pi ");
			}
		}
		box.visit_min[0]=(size_t)MAX((int)(floor((costheta_min-bound_min[0])/l_box[0]*n_side[0])),0);
		box.visit_max[0]=(size_t)MIN((int)(floor((costheta_max-bound_min[0])/l_box[0]*n_side[0])),n_side[0]-1);

		free(thetaphi_high);
		free(thetaphi_low);
	}

}

Mesh catalog_to_mesh(Catalog cat)
{
	size_t iobj,ibox,idim;
	size_t *index=(size_t *) malloc(dim_box*sizeof(size_t));
	histo_t *pos=(histo_t *) malloc(dim_pos*sizeof(histo_t));
	
	size_t *n_obj=(size_t *)malloc(n_boxes*sizeof(size_t));
	for(ibox=0;ibox<n_boxes;ibox++) n_obj[ibox]=0;
	size_t n_full=0;
	for(iobj=0;iobj<cat.n_obj;iobj++) {
		box2ind(&(cat.box[iobj*dim_box]),index);
		size_t ibox=ravel_index(index,n_side,dim_box);
		if(n_obj[ibox]==0) n_full++;
		n_obj[ibox]++;
	}
	size_t *indices = (size_t *)malloc(n_boxes*sizeof(size_t));
	size_t *full = (size_t *)malloc(n_full*sizeof(size_t));
	n_full=0;
	for(ibox=0;ibox<n_boxes;ibox++) {
		if(n_obj[ibox]>0) {
			indices[ibox] = n_full;
			full[n_full] = ibox;
			n_full++;
		}
	}
	if (verbose == INFO) printf(" - there are %zu objects in %zu out of %zu boxes\n",cat.n_obj,n_full,n_boxes);
	Box *boxes = init_boxes(n_full);
	for(ibox=0;ibox<n_full;ibox++) {
		size_t ifull = full[ibox];
		boxes[ibox].pos=(histo_t *)malloc(dim_pos*n_obj[ifull]*sizeof(histo_t));
		if(boxes[ibox].pos==NULL) error_mem_out();
		boxes[ibox].weight=(histo_t *)malloc(dim_weight*n_obj[ifull]*sizeof(histo_t));
		if(boxes[ibox].weight==NULL) error_mem_out();
		boxes[ibox].bin=(size_t *)malloc(n_obj[ifull]*sizeof(size_t));
		if(boxes[ibox].bin==NULL) error_mem_out();
		boxes[ibox].index=(size_t *)malloc(dim_box*sizeof(size_t));
		if(boxes[ibox].index==NULL) error_mem_out();
		boxes[ibox].visit_min=(size_t *)malloc(dim_box*sizeof(size_t));
		if(boxes[ibox].visit_min==NULL) error_mem_out();
		boxes[ibox].visit_max=(size_t *)malloc(dim_box*sizeof(size_t));
		if(boxes[ibox].visit_max==NULL) error_mem_out();
		boxes[ibox].n_obj=0;
		set_visit(boxes[ibox],ifull);
		//printf("visit %zu %zu %zu   %zu %zu %zu\n",boxes[ibox].index[0],boxes[ibox].visit_min[0],boxes[ibox].visit_max[0],boxes[ibox].index[1],boxes[ibox].visit_min[1],boxes[ibox].visit_max[1]);
		//Get box index
	}
	free(n_obj);
	
	for(iobj=0;iobj<cat.n_obj;iobj++) {
		box2ind(&(cat.box[dim_box*iobj]),index);
		size_t ibox=indices[ravel_index(index,n_side,dim_box)];
		size_t n_obj=boxes[ibox].n_obj;
		box2pos(&(cat.box[dim_box*iobj]),pos);
		for (idim=0;idim<dim_pos;idim++) boxes[ibox].pos[dim_pos*n_obj+idim]=pos[idim];
		histo_t *weight=&(cat.weight[dim_weight*iobj]);
		for (idim=0;idim<dim_weight;idim++) boxes[ibox].weight[dim_weight*n_obj+idim]=weight[idim];
		boxes[ibox].bin[n_obj]=cat.bin[iobj];
		boxes[ibox].n_obj++;
	}
	Mesh mesh;
	mesh.n_boxes = n_full;
	mesh.boxes = boxes;

	free(cat.box);
	free(index);
	free(pos);
	free(indices);
	free(full);
	
	return mesh;
}

void set_meshs(Catalog *cats,Mesh *meshs,size_t n_cats)
{
	size_t icat;
	if (verbose == DEBUG) {
		for (icat=0;icat<n_cats;icat++) write_catalog(cats[icat],"debug_cat1.dat");
	}
	for (icat=0;icat<n_cats;icat++) set_catalog_box(&cats[icat]);
	init_params(cats,n_cats);
	for (icat=0;icat<n_cats;icat++) meshs[icat] = catalog_to_mesh(cats[icat]);
	if (verbose == DEBUG) {
		for (icat=0;icat<n_cats;icat++) write_mesh(meshs[icat],"debug_mesh1.dat");
	}
}

void free_meshs(Mesh *meshs,size_t n_meshs)
{
	size_t imesh;
	for (imesh=0;imesh<n_meshs;imesh++) free_mesh(meshs[imesh]);
	free_params(); 
}
