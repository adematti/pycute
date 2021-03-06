#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "define.h"
#include "common.h"

static histo_t fast_dist_main_min,fast_dist_main_max;
static histo_t fast_dist_aux_min,fast_dist_aux_max;

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

static histo_t power(histo_t base,size_t exp)
{
    histo_t result = 1.;
    while (exp)
    {
        if (exp & 1) result *= base;
        exp /= 2;
        base *= base;
    }
    return result;
}

/*
static float fast_inv_sqrt(float x) {
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
  histo_t weight = 0;
  if (weight_type==WEIGHT_PRODSUM) {
    weight = weight1[0]*weight2[0]*(weight1[1]+weight2[1]);
  }
  else {
    for (idim=0;idim<dim_weight;idim++) {
      weight += weight1[idim]*weight2[idim];
    }
  }
  return weight;
}

static histo_t get_fast_distance_main(histo_t *pos1,histo_t *pos2)
{
  if ((corr_type==CORR_SMU)||(corr_type==CORR_SCOS)) {
    size_t idim;
    histo_t dist=0.;
    histo_t s;
    for (idim=0;idim<dim_pos;idim++){
      s = pos1[idim]-pos2[idim];
      dist += s*s;
    }
    return dist;
  }
  else if (corr_type==CORR_ANGULAR) {
    size_t idim;
    histo_t dist = 0;
    for (idim=0;idim<dim_pos;idim++){
      dist += pos1[idim]*pos2[idim];
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

static histo_t get_fast_distance_losn(histo_t *pos)
{
  size_t idim;
  histo_t dist=0.;
  for (idim=0;idim<dim_pos;idim++) dist+=pos[idim]*pos[idim];
  return dist;
}

static histo_t get_distance_losn(histo_t fast_dist,size_t losn)
{
  if (losn==0) return 1.;
  else if (losn%2==0) return power(fast_dist,losn/2);
  return power(my_sqrt(fast_dist),losn);
}

static histo_t get_fast_distance_aux(histo_t *pos1,histo_t *pos2,LOS los,histo_t *dist_los)
{
  size_t idim;
  histo_t dist=0.,norms=0.,normd=0.;
  if (corr_type==CORR_SMU) {
    histo_t s,d;
    for (idim=0;idim<dim_pos;idim++) {
      if (los.type==LOS_MIDPOINT) d=(pos2[idim]+pos1[idim])/2.;
      else if (los.type==LOS_ENDPOINT) d=pos1[idim];
      else if (los.type==LOS_FIRSTPOINT) d=pos2[idim];
      else d=los.los[idim];
      s=pos1[idim]-pos2[idim];
      norms+=s*s;
      normd+=d*d;
      dist+=s*d;
    }
    if ((los.type==LOS_MIDPOINT)||(los.type==LOS_FIRSTPOINT)) *dist_los=get_distance_losn(normd,los.n);
  }
  else if (corr_type==CORR_SCOS) {
    for (idim=0;idim<dim_pos;idim++){
      norms+=pos1[idim]*pos1[idim];
      normd+=pos2[idim]*pos2[idim];
      dist+=pos1[idim]*pos2[idim];
    }
    *dist_los=1.;
  }
  return dist*my_abs(dist)/(norms*normd); //to get the sign
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

static _Bool visit_bin(size_t bin)
{
  //printf("%.3f  ",fast_dist);
  return ((0<=bin)&&(bin<bin_bin.n_bin));
}

static _Bool visit_box(Box box1,Box box2)
{
  if ((corr_type==CORR_SMU)||(corr_type==CORR_SCOS)) {
    size_t idim;
    for (idim=0;idim<dim_box;idim++) {
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


void cross_2pcf_main(Mesh mesh1,Mesh mesh2,histo_t *meanmain,histo_t *count)
{
  histo_t *htreadmain,*threadcount;
  set_fast_distance_main_limit();
  size_t n_bin_tot=bin_main.n_bin;
  size_t ibin;

  for (ibin=0;ibin<n_bin_tot;ibin++) {
    meanmain[ibin]=0.;
    count[ibin]=0.;
  }

#pragma omp parallel default(none)        \
  shared(mesh1,mesh2,meanmain,count,n_bin_tot,bin_main,dim_pos,dim_weight) private(htreadmain,threadcount)
  {
    size_t n_boxes1 = mesh1.n_boxes;
    Box *boxes1 = mesh1.boxes;
    size_t n_boxes2 = mesh2.n_boxes;
    Box *boxes2 = mesh2.boxes;
    size_t ibin,ibox1;
    htreadmain = (histo_t *) calloc(n_bin_tot,sizeof(histo_t));
    threadcount = (histo_t *) calloc(n_bin_tot,sizeof(histo_t));

#pragma omp for schedule(static)
    for (ibox1=0;ibox1<n_boxes1;ibox1++) {
      Box box1=boxes1[ibox1];
      size_t ibox2;
      for (ibox2=0;ibox2<n_boxes2;ibox2++) {
        Box box2=boxes2[ibox2];
        if (visit_box(box1,box2)) {
          size_t n_obj1=box1.n_obj;
          size_t n_obj2=box2.n_obj;
          size_t iobj1;
          for (iobj1=0;iobj1<n_obj1;iobj1++) {
            histo_t *pos1=&(box1.pos[dim_pos*iobj1]);
            histo_t *weight1=&(box1.weight[dim_weight*iobj1]);
            size_t iobj2;
            for (iobj2=0;iobj2<n_obj2;iobj2++) {
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
      for (ibin=0;ibin<n_bin_tot;ibin++) {
        meanmain[ibin]+=htreadmain[ibin];
        count[ibin]+=threadcount[ibin];
      }
      free(htreadmain);
      free(threadcount);
    }
  } //end omp parallel
  for (ibin=0;ibin<n_bin_tot;ibin++) {
    meanmain[ibin]/=count[ibin];
  }
}


void auto_2pcf_main(Mesh mesh1,histo_t *meanmain,histo_t *count)
{
  histo_t *htreadmain,*threadcount;
  set_fast_distance_main_limit();
  size_t n_bin_tot=bin_main.n_bin;
  size_t ibin;

  for (ibin=0;ibin<n_bin_tot;ibin++) {
    meanmain[ibin]=0.;
    count[ibin]=0.;
  }

#pragma omp parallel default(none)        \
  shared(mesh1,meanmain,count,n_bin_tot,bin_main,dim_pos,dim_weight) private(htreadmain,threadcount)
  {
    size_t n_boxes1 = mesh1.n_boxes;
    Box *boxes1 = mesh1.boxes;
    size_t ibin,ibox1;
    htreadmain = (histo_t *) calloc(n_bin_tot,sizeof(histo_t));
    threadcount = (histo_t *) calloc(n_bin_tot,sizeof(histo_t));

#pragma omp for schedule(static)
    for (ibox1=0;ibox1<n_boxes1;ibox1++) {
      Box box1=boxes1[ibox1];
      size_t n_obj1=box1.n_obj;
      size_t iobj1;
      for (iobj1=0;iobj1<n_obj1;iobj1++) {
        histo_t *pos1=&(box1.pos[dim_pos*iobj1]);
        histo_t *weight1=&(box1.weight[dim_weight*iobj1]);
        size_t iobj2;
        for (iobj2=iobj1+1;iobj2<n_obj1;iobj2++) {
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
      for (ibox2=ibox1+1;ibox2<n_boxes1;ibox2++) {
        Box box2=boxes1[ibox2];
        if (visit_box(box1,box2)) {
          size_t n_obj2=box2.n_obj;
          for (iobj1=0;iobj1<n_obj1;iobj1++) {
            histo_t *pos1=&(box1.pos[dim_pos*iobj1]);
            histo_t *weight1=&(box1.weight[dim_weight*iobj1]);
            size_t iobj2;
            for (iobj2=0;iobj2<n_obj2;iobj2++) {
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
      for (ibin=0;ibin<n_bin_tot;ibin++) {
        meanmain[ibin]+=htreadmain[ibin];
        count[ibin]+=threadcount[ibin];
      }
      free(htreadmain);
      free(threadcount);
    }
  } //end omp parallel
  for (ibin=0;ibin<n_bin_tot;ibin++) {
    meanmain[ibin]/=count[ibin];
  }
}



void cross_2pcf_main_aux(Mesh mesh1,Mesh mesh2,histo_t *meanmain,histo_t *meanaux,histo_t *count,LOS los)
{
  histo_t *htreadmain,*htreadaux,*threadcount;
  set_fast_distance_main_limit();
  set_fast_distance_aux_limit();
  size_t n_bin_tot=bin_main.n_bin*bin_aux.n_bin;
  size_t ibin;

  for (ibin=0;ibin<n_bin_tot;ibin++) {
    meanmain[ibin]=0.;
    meanaux[ibin]=0.;
    count[ibin]=0.;
  }

#pragma omp parallel default(none)        \
  shared(mesh1,mesh2,meanmain,meanaux,count,n_bin_tot,bin_main,bin_aux,dim_pos,dim_weight,los) private(htreadmain,htreadaux,threadcount)
  {
    size_t n_boxes1 = mesh1.n_boxes;
    Box *boxes1 = mesh1.boxes;
    size_t n_boxes2 = mesh2.n_boxes;
    Box *boxes2 = mesh2.boxes;
    size_t ibin,ibox1;
    htreadmain = (histo_t *) calloc(n_bin_tot,sizeof(histo_t));
    htreadaux = (histo_t *) calloc(n_bin_tot,sizeof(histo_t));
    threadcount = (histo_t *) calloc(n_bin_tot,sizeof(histo_t));

    histo_t dist_los=0.;
    if (los.type==LOS_CUSTOM) dist_los=get_distance_losn(get_fast_distance_losn(los.los),los.n);

#pragma omp for schedule(static)
    for (ibox1=0;ibox1<n_boxes1;ibox1++) {
      Box box1=boxes1[ibox1];
      size_t ibox2;
      for (ibox2=0;ibox2<n_boxes2;ibox2++) {
        Box box2=boxes2[ibox2];
        if (visit_box(box1,box2)) {
          size_t n_obj1=box1.n_obj;
          size_t n_obj2=box2.n_obj;
          size_t iobj1;
          for (iobj1=0;iobj1<n_obj1;iobj1++) {
            histo_t *pos1=&(box1.pos[dim_pos*iobj1]);
            histo_t *weight1=&(box1.weight[dim_weight*iobj1]);
            if (los.type==LOS_ENDPOINT) dist_los=get_distance_losn(get_fast_distance_losn(pos1),los.n);
            size_t iobj2;
            for (iobj2=0;iobj2<n_obj2;iobj2++) {
              histo_t *pos2=&(box2.pos[dim_pos*iobj2]);
              histo_t *weight2=&(box2.weight[dim_weight*iobj2]);
              histo_t fast_dist_main=get_fast_distance_main(pos1,pos2);
              if (visit_main(fast_dist_main)) {
                histo_t fast_dist_aux=get_fast_distance_aux(pos1,pos2,los,&dist_los);
                if (visit_aux(fast_dist_aux)) {
                  histo_t dist_main=get_distance_main(fast_dist_main);
                  histo_t dist_aux=get_distance_aux(fast_dist_aux);
                  size_t ibin=get_bin_index(dist_aux,bin_aux)+bin_aux.n_bin*get_bin_index(dist_main,bin_main);
                  histo_t weight=get_weight(weight1,weight2)/dist_los;
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
      for (ibin=0;ibin<n_bin_tot;ibin++) {
        meanmain[ibin]+=htreadmain[ibin];
        meanaux[ibin]+=htreadaux[ibin];
        count[ibin]+=threadcount[ibin];
      }
      free(htreadmain);
      free(htreadaux);
      free(threadcount);
    }
  } //end omp parallel
  for (ibin=0;ibin<n_bin_tot;ibin++) {
    meanmain[ibin]/=count[ibin];
    meanaux[ibin]/=count[ibin];
  }
}

void auto_2pcf_main_aux(Mesh mesh1,histo_t *meanmain,histo_t *meanaux,histo_t *count,LOS los)
{
  histo_t *htreadmain,*htreadaux,*threadcount;
  set_fast_distance_main_limit();
  set_fast_distance_aux_limit();
  size_t n_bin_tot=bin_main.n_bin*bin_aux.n_bin;
  size_t ibin;

  for (ibin=0;ibin<n_bin_tot;ibin++) {
    meanmain[ibin]=0.;
    meanaux[ibin]=0.;
    count[ibin]=0.;
  }

#pragma omp parallel default(none)        \
  shared(mesh1,meanmain,meanaux,count,n_bin_tot,bin_main,bin_aux,dim_pos,dim_weight,los) private(htreadmain,htreadaux,threadcount)
  {
    size_t n_boxes1 = mesh1.n_boxes;
    Box *boxes1 = mesh1.boxes;
    size_t ibin,ibox1;
    htreadmain = (histo_t *) calloc(n_bin_tot,sizeof(histo_t));
    htreadaux = (histo_t *) calloc(n_bin_tot,sizeof(histo_t));
    threadcount = (histo_t *) calloc(n_bin_tot,sizeof(histo_t));

    histo_t dist_los=0.;
    if (los.type==LOS_CUSTOM) dist_los=get_distance_losn(get_fast_distance_losn(los.los),los.n);

#pragma omp for schedule(static)
    for (ibox1=0;ibox1<n_boxes1;ibox1++) {
      Box box1=boxes1[ibox1];
      size_t n_obj1=box1.n_obj;
      size_t iobj1;
      for (iobj1=0;iobj1<n_obj1;iobj1++) {
        histo_t *pos1=&(box1.pos[dim_pos*iobj1]);
        histo_t *weight1=&(box1.weight[dim_weight*iobj1]);
        if (los.type==LOS_ENDPOINT) dist_los=get_distance_losn(get_fast_distance_losn(pos1),los.n);
        size_t iobj2;
        for (iobj2=iobj1+1;iobj2<n_obj1;iobj2++) {
          histo_t *pos2=&(box1.pos[dim_pos*iobj2]);
          histo_t *weight2=&(box1.weight[dim_weight*iobj2]);
          histo_t fast_dist_main=get_fast_distance_main(pos1,pos2);
          if (visit_main(fast_dist_main)) {
            histo_t fast_dist_aux=get_fast_distance_aux(pos1,pos2,los,&dist_los);
            if (visit_aux(fast_dist_aux)) {
              histo_t dist_main=get_distance_main(fast_dist_main);
              histo_t dist_aux=get_distance_aux(fast_dist_aux);
              size_t ibin=get_bin_index(dist_aux,bin_aux)+bin_aux.n_bin*get_bin_index(dist_main,bin_main);
              histo_t weight=get_weight(weight1,weight2)/dist_los;
              htreadmain[ibin]+=weight*dist_main;
              htreadaux[ibin]+=weight*dist_aux;
              threadcount[ibin]+=weight;
            }
          }
        }
      }
      size_t ibox2;
      for (ibox2=ibox1+1;ibox2<n_boxes1;ibox2++) {
        Box box2=boxes1[ibox2];
        if (visit_box(box1,box2)) {
          size_t n_obj2=box2.n_obj;
          for (iobj1=0;iobj1<n_obj1;iobj1++) {
            histo_t *pos1=&(box1.pos[dim_pos*iobj1]);
            histo_t *weight1=&(box1.weight[dim_weight*iobj1]);
            size_t iobj2;
            for (iobj2=0;iobj2<n_obj2;iobj2++) {
              histo_t *pos2=&(box2.pos[dim_pos*iobj2]);
              histo_t *weight2=&(box2.weight[dim_weight*iobj2]);
              histo_t fast_dist_main=get_fast_distance_main(pos1,pos2);
              if (visit_main(fast_dist_main)) {
                histo_t fast_dist_aux=get_fast_distance_aux(pos1,pos2,los,&dist_los);
                if (visit_aux(fast_dist_aux)) {
                  histo_t dist_main=get_distance_main(fast_dist_main);
                  histo_t dist_aux=get_distance_aux(fast_dist_aux);
                  size_t ibin=get_bin_index(dist_aux,bin_aux)+bin_aux.n_bin*get_bin_index(dist_main,bin_main);
                  histo_t weight=get_weight(weight1,weight2)/dist_los;
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
      for (ibin=0;ibin<n_bin_tot;ibin++) {
        meanmain[ibin]+=htreadmain[ibin];
        meanaux[ibin]+=htreadaux[ibin];
        count[ibin]+=threadcount[ibin];
      }
      free(htreadmain);
      free(htreadaux);
      free(threadcount);
    }
  } //end omp parallel
  for (ibin=0;ibin<n_bin_tot;ibin++) {
    meanmain[ibin]/=count[ibin];
    meanaux[ibin]/=count[ibin];
  }
}

void legendre_fast(histo_t fast_dist,histo_t *leg,MULTI_TYPE type) {

  if (type==MULTI_ALL) {
    histo_t dist_aux=get_distance_aux(fast_dist);
    //printf("%.3f %.3f  ",dist_aux,fast_dist);
    //if (dist_aux>0) printf("%.3f  ",dist_aux);
    legendre_all(dist_aux,my_abs(fast_dist),leg);
  }
  else if (type==MULTI_EVEN) {
    //printf("%.3f ",fast_dist);
    legendre_even(my_abs(fast_dist),leg);
  }
  else if (type==MULTI_ODD) {
    histo_t dist_aux=get_distance_aux(fast_dist);
    legendre_odd(dist_aux,my_abs(fast_dist),leg);
  }
}

void cross_2pcf_multi(Mesh mesh1,Mesh mesh2,histo_t *meanmain,histo_t *count,Pole pole,LOS los)
{
  histo_t *threadmeanmain,*threadcount;
  set_fast_distance_main_limit();
  set_fast_distance_aux_limit();
  size_t n_bin_tot=bin_main.n_bin*pole.n_ells;
  size_t ibin;
  histo_t leg[MAX_ELLS];

  for (ibin=0;ibin<n_bin_tot;ibin++) {
    meanmain[ibin]=0;
    count[ibin]=0;
  }

#pragma omp parallel default(none)        \
  shared(mesh1,mesh2,meanmain,count,n_bin_tot,bin_main,bin_aux,dim_pos,dim_weight,pole,los) private(threadmeanmain,threadcount,leg)
  {
    size_t n_boxes1 = mesh1.n_boxes;
    Box *boxes1 = mesh1.boxes;
    size_t n_boxes2 = mesh2.n_boxes;
    Box *boxes2 = mesh2.boxes;
    size_t ibin,ibox1;
    MULTI_TYPE multi_type = pole.type;
    size_t n_ells = pole.n_ells;
    threadmeanmain = (histo_t *) calloc(n_bin_tot,sizeof(histo_t));
    threadcount = (histo_t *) calloc(n_bin_tot,sizeof(histo_t));

    histo_t dist_los=0.;
    if (los.type==LOS_CUSTOM) dist_los=get_distance_losn(get_fast_distance_losn(los.los),los.n);

#pragma omp for schedule(static)
    for (ibox1=0;ibox1<n_boxes1;ibox1++) {
      Box box1=boxes1[ibox1];
      size_t ibox2;
      for (ibox2=0;ibox2<n_boxes2;ibox2++) {
        Box box2=boxes2[ibox2];
        if (visit_box(box1,box2)) {
          size_t n_obj1=box1.n_obj;
          size_t n_obj2=box2.n_obj;
          size_t iobj1;
          for (iobj1=0;iobj1<n_obj1;iobj1++) {
            histo_t *pos1=&(box1.pos[dim_pos*iobj1]);
            histo_t *weight1=&(box1.weight[dim_weight*iobj1]);
            if (los.type==LOS_ENDPOINT) dist_los=get_distance_losn(get_fast_distance_losn(pos1),los.n);
            size_t iobj2;
            for (iobj2=0;iobj2<n_obj2;iobj2++) {
              histo_t *pos2=&(box2.pos[dim_pos*iobj2]);
              histo_t fast_dist_main=get_fast_distance_main(pos1,pos2);
              if (visit_main(fast_dist_main)) {
                histo_t fast_dist_aux=get_fast_distance_aux(pos1,pos2,los,&dist_los);
                if (visit_aux(fast_dist_aux)) {
                  histo_t dist_main=get_distance_main(fast_dist_main);
                  histo_t weight=get_weight(weight1,&(box2.weight[dim_weight*iobj2]))/dist_los;
                  legendre_fast(fast_dist_aux,leg,multi_type);
                  size_t ibin=get_bin_index(dist_main,bin_main)*n_ells;
                  size_t ill;
                  for (ill=0;ill<n_ells;ill++) {
                    histo_t w=weight*leg[ill];
                    threadmeanmain[ill+ibin]+=w*dist_main;
                    threadcount[ill+ibin]+=w;
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


void auto_2pcf_multi(Mesh mesh1,histo_t *meanmain,histo_t *count,Pole pole,LOS los)
{

  histo_t *threadmeanmain,*threadcount;
  set_fast_distance_main_limit();
  set_fast_distance_aux_limit();
  size_t n_bin_tot=bin_main.n_bin*pole.n_ells;
  size_t ibin;
  histo_t leg[MAX_ELLS];

  for (ibin=0;ibin<n_bin_tot;ibin++) {
    meanmain[ibin]=0;
    count[ibin]=0;
  }

#pragma omp parallel default(none)        \
  shared(mesh1,meanmain,count,n_bin_tot,bin_main,bin_aux,dim_pos,dim_weight,pole,los) private(threadmeanmain,threadcount,leg)
  {
    size_t n_boxes1 = mesh1.n_boxes;
    Box *boxes1 = mesh1.boxes;
    size_t ibin,ibox1;
    MULTI_TYPE multi_type = pole.type;
    size_t n_ells = pole.n_ells;
    threadmeanmain = (histo_t *) calloc(n_bin_tot,sizeof(histo_t));
    threadcount = (histo_t *) calloc(n_bin_tot,sizeof(histo_t));

    histo_t dist_los=0.;
    if (los.type==LOS_CUSTOM) dist_los=get_distance_losn(get_fast_distance_losn(los.los),los.n);

#pragma omp for schedule(static)
    for (ibox1=0;ibox1<n_boxes1;ibox1++) {
      Box box1=boxes1[ibox1];
      size_t n_obj1=box1.n_obj;
      size_t iobj1;
      for (iobj1=0;iobj1<n_obj1;iobj1++) {
        histo_t *pos1=&(box1.pos[dim_pos*iobj1]);
        histo_t *weight1=&(box1.weight[dim_weight*iobj1]);
        if (los.type==LOS_ENDPOINT) dist_los=get_distance_losn(get_fast_distance_losn(pos1),los.n);
        size_t iobj2;
        for (iobj2=iobj1+1;iobj2<n_obj1;iobj2++) {
          histo_t *pos2=&(box1.pos[dim_pos*iobj2]);
          histo_t fast_dist_main=get_fast_distance_main(pos1,pos2);
          if (visit_main(fast_dist_main)) {
            histo_t fast_dist_aux=get_fast_distance_aux(pos1,pos2,los,&dist_los);
            if (visit_aux(fast_dist_aux)) {
              histo_t dist_main=get_distance_main(fast_dist_main);
              histo_t weight=get_weight(weight1,&(box1.weight[dim_weight*iobj2]))/dist_los;;
              legendre_fast(fast_dist_aux,leg,multi_type);
              size_t ibin=get_bin_index(dist_main,bin_main)*n_ells;
              size_t ill;
              for (ill=0;ill<n_ells;ill++) {
                histo_t w=weight*leg[ill];
                threadmeanmain[ill+ibin]+=w*dist_main;
                threadcount[ill+ibin]+=w;
              }
            }
          }
        }
      }
      size_t ibox2;
      for (ibox2=ibox1+1;ibox2<n_boxes1;ibox2++) {
        Box box2=boxes1[ibox2];
        if (visit_box(box1,box2)) {
          size_t n_obj2=box2.n_obj;
          for (iobj1=0;iobj1<n_obj1;iobj1++) {
            histo_t *pos1=&(box1.pos[dim_pos*iobj1]);
            histo_t *weight1=&(box1.weight[dim_weight*iobj1]);
            size_t iobj2;
            for (iobj2=0;iobj2<n_obj2;iobj2++) {
              histo_t *pos2=&(box2.pos[dim_pos*iobj2]);
              histo_t fast_dist_main=get_fast_distance_main(pos1,pos2);
              if (visit_main(fast_dist_main)) {
                histo_t fast_dist_aux=get_fast_distance_aux(pos1,pos2,los,&dist_los);
                if (visit_aux(fast_dist_aux)) {
                  histo_t dist_main=get_distance_main(fast_dist_main);
                  histo_t weight=get_weight(weight1,&(box2.weight[dim_weight*iobj2]))/dist_los;;
                  legendre_fast(fast_dist_aux,leg,multi_type);
                  size_t ibin=get_bin_index(dist_main,bin_main)*n_ells;
                  size_t ill;
                  for (ill=0;ill<n_ells;ill++) {
                    histo_t w=weight*leg[ill];
                    threadmeanmain[ill+ibin]+=w*dist_main;
                    threadcount[ill+ibin]+=w;
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

void cross_2pcf_multi_radial_legendre(Mesh mesh1,Mesh mesh2,histo_t *count,Pole *poles,LOS los)
{
  histo_t *threadcount;
  set_fast_distance_main_limit();
  set_fast_distance_aux_limit();
  size_t n_bin_tot=bin_main.n_bin*bin_main.n_bin*poles[0].n_ells*poles[1].n_ells;
  size_t ibin;

  for (ibin=0;ibin<n_bin_tot;ibin++) count[ibin]=0;

#pragma omp parallel default(none)        \
  shared(mesh1,mesh2,count,n_bin_tot,bin_main,bin_aux,dim_pos,dim_weight,poles,los) private(threadcount)
  {
    size_t n_boxes1 = mesh1.n_boxes;
    Box *boxes1 = mesh1.boxes;
    size_t n_boxes2 = mesh2.n_boxes;
    Box *boxes2 = mesh2.boxes;
    size_t ibin,ibox1;
    MULTI_TYPE multi_type1 = poles[0].type;
    MULTI_TYPE multi_type2 = poles[1].type;
    size_t n_ells1 = poles[0].n_ells;
    size_t n_ells2 = poles[1].n_ells;
    threadcount = (histo_t *) calloc(n_bin_tot,sizeof(histo_t));
    histo_t leg1[MAX_ELLS],leg2[MAX_ELLS];
    histo_t dist_los=0.;
    if (los.type==LOS_CUSTOM) dist_los=get_distance_losn(get_fast_distance_losn(los.los),los.n);

#pragma omp for schedule(static)
    for (ibox1=0;ibox1<n_boxes1;ibox1++) {
      Box box1=boxes1[ibox1];
      size_t ibox2;
      for (ibox2=0;ibox2<n_boxes2;ibox2++) {
        Box box2=boxes2[ibox2];
        if (visit_box(box1,box2)) {
          size_t n_obj1=box1.n_obj;
          size_t n_obj2=box2.n_obj;
          size_t iobj1;
          for (iobj1=0;iobj1<n_obj1;iobj1++) {
            histo_t *pos1=&(box1.pos[dim_pos*iobj1]);
            histo_t *weight1=&(box1.weight[dim_weight*iobj1]);
            if (los.type==LOS_ENDPOINT) dist_los=get_distance_losn(get_fast_distance_losn(pos1),los.n);
            histo_t x2=get_fast_distance_losn(pos1);
            histo_t x=get_distance_main(x2);
            size_t iobj2;
            for (iobj2=0;iobj2<n_obj2;iobj2++) {
              histo_t *pos2=&(box2.pos[dim_pos*iobj2]);
              histo_t fast_dist_main=get_fast_distance_main(pos1,pos2);
              if (visit_main(fast_dist_main)) {
                histo_t fast_dist_aux=get_fast_distance_aux(pos1,pos2,los,&dist_los);
                if (visit_aux(fast_dist_aux)) {
                  histo_t dist_main=get_distance_main(fast_dist_main);
                  histo_t weight=get_weight(weight1,&(box2.weight[dim_weight*iobj2]))/dist_los;
                  legendre_fast(fast_dist_aux,leg1,multi_type1);
                  size_t ibin1=get_bin_index(dist_main,bin_main);
                  histo_t xd2=get_fast_distance_losn(pos2);
                  histo_t xd=get_distance_main(xd2);
                  size_t ibin2;
                  for (ibin2=0;ibin2<bin_main.n_bin;ibin2++) {
                    histo_t s=get_bin_mid(ibin2,bin_main);
                    histo_t J=x*s/xd;
                    histo_t mus=(xd2-x2-s*s)/(2.*x*s);
                    //mus=-(x*mus+s)/xd;
                    if (visit_aux(mus*my_abs(mus))) {
                    //if (my_abs(mus)<=1.) {
                      legendre(mus,leg2,multi_type2);
                      size_t ill1,ill2;
                      for (ill1=0;ill1<n_ells1;ill1++) {
                        histo_t tmp = weight*leg1[ill1]/J;
                        for (ill2=0;ill2<n_ells2;ill2++) threadcount[ill2+n_ells2*(ill1+n_ells1*(ibin2+bin_main.n_bin*ibin1))]+=tmp*leg2[ill2];
                      }
                    }
                    //else if (mus<-1.) break;
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
      for (ibin=0;ibin<n_bin_tot;ibin++) count[ibin]+=threadcount[ibin];
      free(threadcount);
    }
  } //end omp parallel

}

void cross_2pcf_multi_angular_legendre(Mesh mesh1,Mesh mesh2,histo_t *count,Pole *poles,LOS los)
{
  histo_t *threadcount;
  set_fast_distance_main_limit();
  set_fast_distance_aux_limit();
  size_t n_bin_tot=bin_main.n_bin*bin_main.n_bin*poles[0].n_ells*poles[1].n_ells;
  size_t ibin;

  for (ibin=0;ibin<n_bin_tot;ibin++) count[ibin]=0;

#pragma omp parallel default(none)        \
  shared(mesh1,mesh2,count,n_bin_tot,bin_main,bin_aux,dim_pos,dim_weight,poles,los) private(threadcount)
  {
    size_t n_boxes1 = mesh1.n_boxes;
    Box *boxes1 = mesh1.boxes;
    size_t n_boxes2 = mesh2.n_boxes;
    Box *boxes2 = mesh2.boxes;
    size_t ibin,ibox1;
    MULTI_TYPE multi_type1 = poles[0].type;
    MULTI_TYPE multi_type2 = poles[1].type;
    size_t n_ells1 = poles[0].n_ells;
    size_t n_ells2 = poles[1].n_ells;
    threadcount = (histo_t *) calloc(n_bin_tot,sizeof(histo_t));
    histo_t leg1[MAX_ELLS],leg2[MAX_ELLS];
    histo_t dist_los=0.;
    if (los.type==LOS_CUSTOM) dist_los=get_distance_losn(get_fast_distance_losn(los.los),los.n);

#pragma omp for schedule(static)
    for (ibox1=0;ibox1<n_boxes1;ibox1++) {
      Box box1=boxes1[ibox1];
      size_t ibox2;
      for (ibox2=0;ibox2<n_boxes2;ibox2++) {
        Box box2=boxes2[ibox2];
        if (visit_box(box1,box2)) {
          size_t n_obj1=box1.n_obj;
          size_t n_obj2=box2.n_obj;
          size_t iobj1;
          for (iobj1=0;iobj1<n_obj1;iobj1++) {
            histo_t *pos1=&(box1.pos[dim_pos*iobj1]);
            histo_t *weight1=&(box1.weight[dim_weight*iobj1]);
            if (los.type==LOS_ENDPOINT) dist_los=get_distance_losn(get_fast_distance_losn(pos1),los.n);
            histo_t x2=get_fast_distance_losn(pos1);
            histo_t x=get_distance_main(x2);
            size_t iobj2;
            for (iobj2=0;iobj2<n_obj2;iobj2++) {
              histo_t *pos2=&(box2.pos[dim_pos*iobj2]);
              histo_t fast_dist_main=get_fast_distance_main(pos1,pos2);
              if (visit_main(fast_dist_main)) {
                histo_t fast_dist_aux=get_fast_distance_aux(pos1,pos2,los,&dist_los);
                if (visit_aux(fast_dist_aux)) {
                  histo_t dist_main=get_distance_main(fast_dist_main);
                  histo_t weight=get_weight(weight1,&(box2.weight[dim_weight*iobj2]))/dist_los;
                  legendre_fast(fast_dist_aux,leg1,multi_type1);
                  size_t ibin1=get_bin_index(dist_main,bin_main);
                  histo_t xd2=get_fast_distance_losn(pos2);
                  histo_t d2=fast_dist_main;
                  histo_t mud2=my_abs(fast_dist_aux);
                  size_t ibin2;
                  for (ibin2=0;ibin2<bin_main.n_bin;ibin2++) {
                    histo_t s=get_bin_mid(ibin2,bin_main);
                    histo_t s2=s*s;
                    histo_t b=2.*(1.-mud2)*d2*x/(xd2*s);
                    histo_t c=-1.+(x2+s2)*(1.-mud2)*d2/(xd2*s2);
                    histo_t det=b*b-4.*c;
                    if (det<0.) continue;
                    det=my_sqrt(det);
                    histo_t solmus[2]={(-b-det)/2.,(-b+det)/2.};
                    size_t imu;
                    for (imu=0;imu<2;imu++) {
                      histo_t mus=solmus[imu];
                      if (visit_aux(mus*my_abs(mus))) {
                        //printf("%zu %zu %.3lf  ",ibin2,imu,mus);
                        histo_t xs2=x2+2.*mus*x*s+s2;
                        histo_t xs=my_sqrt(xs2);
                        histo_t J=my_abs(s/xs-(x+s*mus)*x*s/(xs2*xs));
                        //mus=-(x*mus+s)/xs;
                        //if (visit_aux(mus*my_abs(mus))) {
                        legendre(mus,leg2,multi_type2);
                        size_t ill1,ill2;
                        for (ill1=0;ill1<n_ells1;ill1++) {
                          histo_t tmp = weight*leg1[ill1]/J;
                          for (ill2=0;ill2<n_ells2;ill2++) threadcount[ill2+n_ells2*(ill1+n_ells1*(ibin2+bin_main.n_bin*ibin1))]+=tmp*leg2[ill2];
                        }
                        //}
                      }
                    }
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
      for (ibin=0;ibin<n_bin_tot;ibin++) count[ibin]+=threadcount[ibin];
      free(threadcount);
    }
  } //end omp parallel

}

/*
void cross_2pcf_multi_radial_legendre(Mesh mesh1,Mesh mesh2,histo_t *count,Pole *poles)
{
  histo_t *threadcount;
  set_fast_distance_main_limit();
  set_fast_distance_aux_limit();
  size_t n_bin_tot=bin_main.n_bin*bin_main.n_bin*poles[0].n_ells*poles[1].n_ells;
  size_t ibin;

  for (ibin=0;ibin<n_bin_tot;ibin++) count[ibin]=0;

#pragma omp parallel default(none)        \
  shared(mesh1,mesh2,count,n_bin_tot,bin_main,bin_aux,dim_pos,dim_weight,poles) private(threadcount)
  {
    size_t n_boxes1 = mesh1.n_boxes;
    Box *boxes1 = mesh1.boxes;
    size_t n_boxes2 = mesh2.n_boxes;
    Box *boxes2 = mesh2.boxes;
    size_t ibin,ibox1;
    MULTI_TYPE multi_type1 = poles[0].type;
    MULTI_TYPE multi_type2 = poles[1].type;
    size_t n_ells1 = poles[0].n_ells;
    size_t n_ells2 = poles[1].n_ells;
    threadcount = (histo_t *) calloc(n_bin_tot,sizeof(histo_t));
    for (ibin=0;ibin<n_bin_tot;ibin++) threadcount[ibin]=0;
    histo_t leg1[MAX_ELLS],leg2[MAX_ELLS];

#pragma omp for schedule(static)
    for (ibox1=0;ibox1<n_boxes1;ibox1++) {
      Box box1=boxes1[ibox1];
      size_t ibox2;
      for (ibox2=0;ibox2<n_boxes2;ibox2++) {
        Box box2=boxes2[ibox2];
        if (visit_box(box1,box2)) {
          size_t n_obj1=box1.n_obj;
          size_t n_obj2=box2.n_obj;
          size_t iobj1;
          for (iobj1=0;iobj1<n_obj1;iobj1++) {
            histo_t *pos1=&(box1.pos[dim_pos*iobj1]);
            histo_t *weight1=&(box1.weight[dim_weight*iobj1]);
            size_t iobj2;
            for (iobj2=0;iobj2<n_obj2;iobj2++) {
              histo_t *pos2=&(box2.pos[dim_pos*iobj2]);
              histo_t fast_dist_main=get_fast_distance_main(pos1,pos2);
              if (visit_main(fast_dist_main)) {
                histo_t fast_dist_aux=get_fast_distance_aux(pos1,pos2);
                if (visit_aux(fast_dist_aux)) {
                  histo_t dist_main=get_distance_main(fast_dist_main);
                  histo_t weight=get_weight(weight1,&(box2.weight[dim_weight*iobj2]));
                  legendre_fast(fast_dist_aux,leg1,multi_type1);
                  size_t ibin1=get_bin_index(dist_main,bin_main);
                  long ibin2; //can be <0
                  for (ibin2=bin_main.n_bin-1;ibin2>=0;ibin2--) {
                    histo_t dist_main2=get_bin_mid(ibin2,bin_main);
                    histo_t fast_dist_aux2=fast_dist_main*fast_dist_aux/(dist_main2*dist_main2);
                    if (visit_aux(fast_dist_aux2)) {
                      legendre_fast(fast_dist_aux2,leg2,multi_type2);
                      size_t ill1,ill2;
                      for (ill1=0;ill1<n_ells1;ill1++) {
                        histo_t tmp = weight*leg1[ill1]/dist_main2;
                        for (ill2=0;ill2<n_ells2;ill2++) threadcount[ill2+n_ells2*(ill1+n_ells1*(ibin2+bin_main.n_bin*ibin1))]+=tmp*leg2[ill2];
                      }
                    }
                    else if (my_abs(fast_dist_aux2)>1.) break;
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
      for (ibin=0;ibin<n_bin_tot;ibin++) count[ibin]+=threadcount[ibin];
      free(threadcount);
    }
  } //end omp parallel

}
*/


void cross_3pcf_multi(Mesh* meshs,size_t n_meshs,histo_t *count,Pole *poles,LOS *los)
{
  set_fast_distance_main_limit();
  set_fast_distance_aux_limit();
  size_t ibin,ibin2,ibin3,ill2,ill3;
  size_t n_sec = n_meshs - 1;
  size_t n_bin_sec[2] = {bin_main.n_bin*poles[0].n_ells,bin_main.n_bin*poles[1].n_ells};
  size_t n_bin_tot = (n_sec == 2) ? n_bin_sec[0]*n_bin_sec[1] : n_bin_sec[0]*n_bin_sec[0];
  histo_t leg[MAX_ELLS];
  histo_t* countsec[2];
  //printf("nells %zu\n",n_ells);
  for (ibin=0;ibin<n_bin_tot;ibin++) count[ibin]=0.;
  histo_t* threadcount;

#pragma omp parallel default(none)        \
  shared(count,meshs,n_sec,n_bin_tot,n_bin_sec,bin_main,dim_pos,dim_weight,poles,los) private(countsec,threadcount,leg)
  {
    size_t isec,ibin,ibox1;
    size_t n_boxes1 = meshs[0].n_boxes;
    Box *boxes1 = meshs[0].boxes;
    for (isec=0;isec<n_sec;isec++) countsec[isec] = (histo_t*) malloc(n_bin_sec[isec]*sizeof(histo_t));
    threadcount = (histo_t*) calloc(n_bin_tot,sizeof(histo_t));
    histo_t dist_los[2];
    for (isec=0;isec<n_sec;isec++) {
      if (los[isec].type==LOS_CUSTOM) dist_los[isec]=get_distance_losn(get_fast_distance_losn(los[isec].los),los[isec].n);
    }
#pragma omp for schedule(static)
    for (ibox1=0;ibox1<n_boxes1;ibox1++) {
      Box box1=boxes1[ibox1];
      size_t n_obj1=box1.n_obj;
      size_t iobj1;
      for (iobj1=0;iobj1<n_obj1;iobj1++) {
        histo_t* pos1=&(box1.pos[dim_pos*iobj1]);
        histo_t weight1=box1.weight[dim_weight*iobj1];
        if (weight1==0.) continue;
        size_t isec;
        for (isec=0;isec<n_sec;isec++) {
          histo_t *count2 = countsec[isec];
          size_t ibin2;
          for (ibin2=0;ibin2<n_bin_sec[isec];ibin2++) count2[ibin2]=0.;
        }
        for (isec=0;isec<n_sec;isec++) {
          size_t n_boxes2 = meshs[isec+1].n_boxes;
          Box *boxes2 = meshs[isec+1].boxes;
          histo_t *count2 = countsec[isec];
          MULTI_TYPE multi_type2 = poles[isec].type;
          size_t n_ells2 = poles[isec].n_ells;
          LOS los2 = los[isec];
          histo_t dist_los2=dist_los[isec];
          if (los2.type==LOS_ENDPOINT) dist_los2=get_distance_losn(get_fast_distance_losn(pos1),los2.n);
          size_t ibox2;
          for (ibox2=0;ibox2<n_boxes2;ibox2++) {
            Box box2=boxes2[ibox2];
            if (visit_box(box1,box2)) {
              size_t n_obj2=box2.n_obj;
              size_t iobj2;
              for (iobj2=0;iobj2<n_obj2;iobj2++) {
                histo_t* pos2=&(box2.pos[dim_pos*iobj2]);
                histo_t fast_dist_main = get_fast_distance_main(pos1,pos2);
                if (visit_main(fast_dist_main)) {
                  histo_t fast_dist_aux = get_fast_distance_aux(pos1,pos2,los2,&dist_los2);
                  if (visit_aux(fast_dist_aux)) {
                    size_t ibin = get_bin_index(get_distance_main(fast_dist_main),bin_main);
                    histo_t weight = box2.weight[dim_weight*iobj2]/dist_los2;
                    legendre_fast(fast_dist_aux,leg,multi_type2);
                    size_t ill2;
                    for (ill2=0;ill2<n_ells2;ill2++) count2[ill2+ibin*n_ells2] += leg[ill2]*weight;
                  }
                }
              }
            }
          }
        }
        if (n_sec==1) {
          histo_t *count2 = countsec[0];
          histo_t *count3 = countsec[0];
          size_t n_ells2 = poles[0].n_ells;
          size_t n_ells3 = n_ells2;
          size_t ibin2;
          for (ibin2=0;ibin2<bin_main.n_bin;ibin2++) {
            if (count2[ibin2*n_ells2]!=0.) {
              size_t ibin3;
              for (ibin3=0;ibin3<bin_main.n_bin;ibin3++) {
                if (count3[ibin3*n_ells3]!=0.) {
                  size_t ill2,ill3;
                  for (ill2=0;ill2<n_ells2;ill2++) {
                    for (ill3=ill2;ill3<n_ells3;ill3++) {
                      threadcount[ill3+n_ells3*(ill2+n_ells2*(ibin3+bin_main.n_bin*ibin2))] += weight1*count2[ill2+ibin2*n_ells2]*count3[ill3+ibin3*n_ells3];
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
          size_t n_ells2 = poles[0].n_ells;
          size_t n_ells3 = poles[1].n_ells;
          size_t ibin2;
          for (ibin2=0;ibin2<bin_main.n_bin;ibin2++) {
            if (count2[ibin2*n_ells2]!=0.) {
              size_t ibin3;
              for (ibin3=0;ibin3<bin_main.n_bin;ibin3++) {
                if (count3[ibin3*n_ells3]!=0.) {
                  size_t ill2,ill3;
                  for (ill2=0;ill2<n_ells2;ill2++) {
                    for (ill3=0;ill3<n_ells3;ill3++) {
                      threadcount[ill3+n_ells3*(ill2+n_ells2*(ibin3+bin_main.n_bin*ibin2))] += weight1*count2[ill2+ibin2*n_ells2]*count3[ill3+ibin3*n_ells3];
                    }
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
      for (ibin=0;ibin<n_bin_tot;ibin++) count[ibin]+=threadcount[ibin];
      free(threadcount);
      for (isec=0;isec<n_sec;isec++) free(countsec[isec]);
    }
  } //end omp parallel
  if (n_sec==1) {
    size_t n_ells2 = poles[0].n_ells;
    size_t n_ells3 = n_ells2;
    for (ibin2=0;ibin2<bin_main.n_bin;ibin2++) {
      for (ibin3=0;ibin3<bin_main.n_bin;ibin3++) {
        for (ill2=0;ill2<n_ells2;ill2++) {
          for (ill3=0;ill3<ill2;ill3++) {
            count[ill3+n_ells3*(ill2+n_ells2*(ibin3+bin_main.n_bin*ibin2))] = count[ill2+n_ells2*(ill3+n_ells3*(ibin2+bin_main.n_bin*ibin3))];
          }
        }
      }
    }
  }
}

void cross_3pcf_multi_double_los(Mesh* meshs,size_t n_meshs,histo_t *count,Pole *poles,LOS *los)
{
  set_fast_distance_main_limit();
  set_fast_distance_aux_limit();
  size_t ibin;
  size_t n_sec = n_meshs - 1;
  size_t n_bin_sec[2] = {bin_main.n_bin*poles[0].n_ells,bin_main.n_bin*poles[1].n_ells};
  size_t n_bin_tot = (n_sec == 2) ? n_bin_sec[0]*n_bin_sec[1] : n_bin_sec[0]*n_bin_sec[0];
  histo_t leg[MAX_ELLS];
  histo_t* countsec[2];
  //printf("nells %zu\n",n_ells);
  for (ibin=0;ibin<n_bin_tot;ibin++) count[ibin]=0.;
  histo_t* threadcount;

#pragma omp parallel default(none)        \
  shared(count,meshs,n_sec,n_bin_tot,n_bin_sec,bin_main,dim_pos,dim_weight,poles,los) private(countsec,threadcount,leg)
  {
    size_t isec,ibin,ibox1;
    size_t n_boxes1 = meshs[0].n_boxes;
    Box *boxes1 = meshs[0].boxes;
    for (isec=0;isec<2;isec++) countsec[isec] = (histo_t*) malloc(n_bin_sec[isec]*sizeof(histo_t));
    threadcount = (histo_t*) calloc(n_bin_tot,sizeof(histo_t));
    histo_t dist_los[2];
    for (isec=0;isec<n_sec;isec++) {
      if (los[isec].type==LOS_CUSTOM) dist_los[isec]=get_distance_losn(get_fast_distance_losn(los[isec].los),los[isec].n);
    }
#pragma omp for schedule(static)
    for (ibox1=0;ibox1<n_boxes1;ibox1++) {
      Box box1=boxes1[ibox1];
      size_t n_obj1=box1.n_obj;
      size_t iobj1;
      for (iobj1=0;iobj1<n_obj1;iobj1++) {
        histo_t* pos1=&(box1.pos[dim_pos*iobj1]);
        histo_t weight1=box1.weight[dim_weight*iobj1];
        if (weight1==0.) continue;
        size_t isec;
        for (isec=0;isec<2;isec++) {
          histo_t *count2 = countsec[isec];
          size_t ibin2;
          for (ibin2=0;ibin2<n_bin_sec[isec];ibin2++) count2[ibin2]=0.;
        }
        for (isec=0;isec<n_sec;isec++) {
          size_t n_boxes2 = meshs[isec+1].n_boxes;
          Box *boxes2 = meshs[isec+1].boxes;
          histo_t *count2 = countsec[isec];
          histo_t *count3 = countsec[1];
          MULTI_TYPE multi_type2 = poles[isec].type;
          size_t n_ells2 = poles[isec].n_ells;
          LOS los2 = los[isec];
          histo_t dist_los2=dist_los[isec];
          if (los2.type==LOS_ENDPOINT) dist_los2 = get_distance_losn(get_fast_distance_losn(pos1),los2.n);
          size_t ibox2;
          for (ibox2=0;ibox2<n_boxes2;ibox2++) {
            Box box2=boxes2[ibox2];
            if (visit_box(box1,box2)) {
              size_t n_obj2=box2.n_obj;
              size_t iobj2;
              for (iobj2=0;iobj2<n_obj2;iobj2++) {
                histo_t* pos2=&(box2.pos[dim_pos*iobj2]);
                histo_t fast_dist_main = get_fast_distance_main(pos1,pos2);
                if (visit_main(fast_dist_main)) {
                  histo_t fast_dist_aux1 = get_fast_distance_aux(pos1,pos2,los2,&dist_los2);
                  histo_t fast_dist_aux2 = (isec==0) ? get_fast_distance_aux(pos2,pos1,los2,&dist_los2) : fast_dist_aux1;
                  if (visit_aux(fast_dist_aux1) && visit_aux(fast_dist_aux2)) {
                    size_t ibin = get_bin_index(get_distance_main(fast_dist_main),bin_main);
                    histo_t weight = box2.weight[dim_weight*iobj2]/dist_los2;
                    legendre_fast(fast_dist_aux1,leg,multi_type2);
                    size_t ill2;
                    for (ill2=0;ill2<n_ells2;ill2++) count2[ill2+ibin*n_ells2] += leg[ill2]*weight;
                    if (n_sec==1) {
                      for (ill2=0;ill2<n_ells2;ill2++) count3[ill2+ibin*n_ells2] += leg[ill2]*weight;
                    }
                    if (isec==0) {
                      legendre_fast(fast_dist_aux2,leg,multi_type2);
                      for (ill2=0;ill2<n_ells2;ill2++) count2[ill2+ibin*n_ells2] += leg[ill2]*weight;
                    }
                  }
                }
              }
            }
          }
        }
        histo_t *count2 = countsec[0];
        histo_t *count3 = countsec[1];
        size_t n_ells2 = poles[0].n_ells;
        size_t n_ells3 = poles[1].n_ells;
        size_t ibin2;
        for (ibin2=0;ibin2<bin_main.n_bin;ibin2++) {
          if (count2[ibin2*n_ells2]!=0.) {
            size_t ibin3;
            for (ibin3=0;ibin3<bin_main.n_bin;ibin3++) {
              if (count3[ibin3*n_ells3]!=0.) {
                size_t ill2,ill3;
                for (ill2=0;ill2<n_ells2;ill2++) {
                  for (ill3=0;ill3<n_ells3;ill3++) {
                    threadcount[ill3+n_ells3*(ill2+n_ells2*(ibin3+bin_main.n_bin*ibin2))] += weight1*count2[ill2+ibin2*n_ells2]*count3[ill3+ibin3*n_ells3];
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
      for (ibin=0;ibin<n_bin_tot;ibin++) count[ibin]+=threadcount[ibin];
      free(threadcount);
      for (isec=0;isec<2;isec++) free(countsec[isec]);
    }
  } //end omp parallel
}
/*
void cross_2pcf_multi_binned(Mesh mesh1, Mesh mesh2, histo_t *count,Pole pole,LOS los,_Bool normalize)
{
  set_fast_distance_main_limit();
  set_fast_distance_aux_limit();
  size_t n_bin_tot = bin_main.n_bin*bin_bin.n_bin*pole.n_ells;
  size_t ibin;
  histo_t leg[MAX_ELLS];
  histo_t norm1=1.;
  histo_t *norm2 = (histo_t*) calloc(bin_bin.n_bin,sizeof(histo_t));
  //printf("%zu \n",cat1.n_obj);

  size_t n_boxes1 = mesh1.n_boxes;
  Box* boxes1 = mesh1.boxes;
  size_t n_boxes2 = mesh2.n_boxes;
  Box* boxes2 = mesh2.boxes;

  if (normalize) {
    size_t ibox;
    norm1=0;
    for (ibox=0;ibox<n_boxes1;ibox++) {
      Box box1 = boxes1[ibox];
      size_t n_obj = box1.n_obj;
      size_t iobj;
      for (iobj=0;iobj<n_obj;iobj++) norm1 += box1.weight[dim_weight*iobj];
    }
    for (ibox=0;ibox<n_boxes2;ibox++) {
      Box box2 = boxes2[ibox];
      size_t n_obj = box2.n_obj;
      size_t iobj;
      for (iobj=0;iobj<n_obj;iobj++) {
        size_t bin = box2.bin[iobj];
        if (visit_bin(bin)) norm2[bin] += box2.weight[dim_weight*iobj];
      }
    }
  }

  for (ibin=0;ibin<n_bin_tot;ibin++) count[ibin]=0.;
  histo_t* threadcount;

#pragma omp parallel default(none)        \
  shared(count,boxes1,boxes2,n_boxes1,n_boxes2,n_bin_tot,bin_main,bin_bin,dim_pos,dim_weight,pole,los,normalize,norm1,norm2) private(threadcount,leg)
  {
    size_t ibin,ibox2;
    MULTI_TYPE multi_type = pole.type;
    size_t n_ells = pole.n_ells;
    threadcount = (histo_t*) calloc(n_bin_tot,sizeof(histo_t))
    histo_t dist_los=0.;
    if (los.type==LOS_CUSTOM) dist_los=get_distance_losn(get_fast_distance_losn(los.los),los.n);

#pragma omp for schedule(static)
    for (ibox2=0;ibox2<n_boxes2;ibox2++) {
      Box box2=boxes2[ibox2];
      size_t n_obj2=box2.n_obj;
      size_t iobj2;
      for (iobj2=0;iobj2<n_obj2;iobj2++) {
        histo_t* pos2=&(box2.pos[dim_pos*iobj2]);
        histo_t* weight2=&(box2.weight[dim_weight*iobj2]);
        size_t bin = box2.bin[iobj2];
        if (visit_bin(bin)) {
          size_t ibox1;
          for (ibox1=0;ibox1<n_boxes1;ibox1++) {
            Box box1=boxes1[ibox1];
            if (visit_box(box1,box2)) {
              size_t n_obj1=box1.n_obj;
              size_t iobj1;
              for (iobj1=0;iobj1<n_obj1;iobj1++) {
                histo_t* pos1=&(box1.pos[dim_pos*iobj1]);
                if (los.type==LOS_ENDPOINT) dist_los=get_distance_losn(get_fast_distance_losn(pos1),los.n);
                histo_t fast_dist_main = get_fast_distance_main(pos1,pos2);
                histo_t fast_dist_aux = get_fast_distance_aux(pos1,pos2,los,&dist_los);
                if (visit_main(fast_dist_main) && visit_aux(fast_dist_aux)) {
                  ibin = bin+bin_bin.n_bin*get_bin_index(get_distance_main(fast_dist_main),bin_main);
                  legendre_fast(fast_dist_aux,leg,multi_type);
                  histo_t weight=get_weight(&(box1.weight[dim_weight*iobj1]),weight2)/dist_los;
                  if (normalize) weight /= norm1*norm2[bin];
                  size_t ill;
                  for (ill=0;ill<n_ells;ill++) threadcount[ill+ibin*n_ells] += leg[ill]*weight;
                }
              }
            }
          }
        }
      }
    }
#pragma omp critical
    {
      for (ibin=0;ibin<n_bin_tot;ibin++) count[ibin]+=threadcount[ibin];
      free(threadcount);
    }
  } //end omp parallel
  free(norm2);
}
*/
void cross_2pcf_multi_binned(Mesh mesh1, Mesh mesh2, histo_t *count,Pole pole,LOS los,size_t tobin)
{
  set_fast_distance_main_limit();
  set_fast_distance_aux_limit();
  size_t n_bin_tot = bin_main.n_bin*bin_bin.n_bin*pole.n_ells;
  size_t ibin;
  histo_t leg[MAX_ELLS];

  size_t n_boxes1 = mesh1.n_boxes;
  Box* boxes1 = mesh1.boxes;
  size_t n_boxes2 = mesh2.n_boxes;
  Box* boxes2 = mesh2.boxes;

  for (ibin=0;ibin<n_bin_tot;ibin++) count[ibin]=0.;
  histo_t* threadcount;

#pragma omp parallel default(none)        \
  shared(count,boxes1,boxes2,n_boxes1,n_boxes2,n_bin_tot,bin_main,bin_bin,dim_pos,dim_weight,pole,los,tobin) private(threadcount,leg)
  {
    size_t ibin,ibox1;
    MULTI_TYPE multi_type = pole.type;
    size_t n_ells = pole.n_ells;
    threadcount = (histo_t*) calloc(n_bin_tot,sizeof(histo_t));
    histo_t dist_los=0.;
    if (los.type==LOS_CUSTOM) dist_los=get_distance_losn(get_fast_distance_losn(los.los),los.n);

#pragma omp for schedule(static)
    for (ibox1=0;ibox1<n_boxes1;ibox1++) {
      Box box1=boxes1[ibox1];
      size_t ibox2;
      for (ibox2=0;ibox2<n_boxes2;ibox2++) {
        Box box2=boxes2[ibox2];
        if (visit_box(box1,box2)) {
          size_t n_obj1=box1.n_obj;
          size_t n_obj2=box2.n_obj;
          size_t iobj1;
          size_t bin=0;
          for (iobj1=0;iobj1<n_obj1;iobj1++) {
            histo_t* pos1=&(box1.pos[dim_pos*iobj1]);
            histo_t* weight1=&(box1.weight[dim_weight*iobj1]);
            if (los.type==LOS_ENDPOINT) dist_los=get_distance_losn(get_fast_distance_losn(pos1),los.n);
            if (tobin==1) {
              bin = box1.bin[iobj1];
              if (!visit_bin(bin)) continue;
            }
            size_t iobj2;
            for (iobj2=0;iobj2<n_obj2;iobj2++) {
              if (tobin==2) {
                bin = box2.bin[iobj2];
                if (!visit_bin(bin)) continue;
              }
              histo_t* pos2=&(box2.pos[dim_pos*iobj2]);
              histo_t fast_dist_main = get_fast_distance_main(pos1,pos2);
              histo_t fast_dist_aux = get_fast_distance_aux(pos1,pos2,los,&dist_los);
              if (visit_main(fast_dist_main) && visit_aux(fast_dist_aux)) {
                ibin = bin+bin_bin.n_bin*get_bin_index(get_distance_main(fast_dist_main),bin_main);
                legendre_fast(fast_dist_aux,leg,multi_type);
                histo_t weight=get_weight(weight1,&(box2.weight[dim_weight*iobj2]))/dist_los;
                size_t ill;
                for (ill=0;ill<n_ells;ill++) threadcount[ill+ibin*n_ells] += leg[ill]*weight;
              }
            }
          }
        }
      }
    }
#pragma omp critical
    {
      for (ibin=0;ibin<n_bin_tot;ibin++) count[ibin]+=threadcount[ibin];
      free(threadcount);
    }
  } //end omp parallel
}

void cross_4pcf_multi_binned(Mesh *meshs,histo_t *count,Pole *poles,LOS *los,size_t *tobin)
{
  size_t n_bin_main = bin_main.n_bin;
  size_t n_bin_bin = bin_bin.n_bin;
  size_t n_ells1 = poles[0].n_ells;
  size_t n_ells2 = poles[1].n_ells;
  size_t n_bin_2pcf1 = n_bin_main*n_bin_bin*n_ells1;
  size_t n_bin_2pcf2 = n_bin_main*n_bin_bin*n_ells2;
  size_t n_bin_tot = n_bin_main*n_bin_main*n_ells1*n_ells2;
  histo_t *count1 = (histo_t*) malloc(n_bin_2pcf1*sizeof(histo_t));
  cross_2pcf_multi_binned(meshs[0],meshs[1],count1,poles[0],los[0],tobin[0]);
  histo_t *count2 = (histo_t*) malloc(n_bin_2pcf2*sizeof(histo_t));
  cross_2pcf_multi_binned(meshs[2],meshs[3],count2,poles[1],los[1],tobin[1]-2);

  size_t ibin;
  for (ibin=0;ibin<n_bin_tot;ibin++) count[ibin]=0.;
  histo_t* threadcount;

#pragma omp parallel default(none)        \
  shared(count,count1,count2,n_bin_tot,n_bin_main,n_bin_bin,n_ells1,n_ells2) private(threadcount)
  {
    size_t ibin;
    threadcount = (histo_t*) calloc(n_bin_tot,sizeof(histo_t));

#pragma omp for schedule(static)
    for (ibin=0;ibin<n_bin_bin;ibin++) {
      size_t ibin1;
      for (ibin1=0;ibin1<n_bin_main;ibin1++) {
        size_t ill1;
        for (ill1=0;ill1<n_ells1;ill1++) {
          histo_t tmp1 = count1[ill1+n_ells1*(ibin+ibin1*n_bin_bin)];
          if (tmp1!=0.) {
            size_t ibin2;
            for (ibin2=0;ibin2<n_bin_main;ibin2++) {
              size_t ill2;
              for (ill2=0;ill2<n_ells2;ill2++) {
                threadcount[ill2+n_ells2*(ill1+n_ells1*(ibin2+n_bin_main*ibin1))] += tmp1*count2[ill2+n_ells2*(ibin+ibin2*n_bin_bin)];
              }
            }
          }
        }
      }
    }
#pragma omp critical
    {
      for (ibin=0;ibin<n_bin_tot;ibin++) count[ibin] += threadcount[ibin];
      free(threadcount);
    }
  } //end omp parallel
  free(count1);
  free(count2);
}
