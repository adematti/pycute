#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "define.h"
#include "common.h"

//  Timing variables
#ifdef _HAVE_OMP
#include <omp.h>
static double relbeg,relend,absbeg,absend;
#else //_HAVE_OMP
#include <time.h>
static time_t relbeg,relend,absbeg,absend;
#endif //_HAVE_OMP

///////////////////////////
//General purpose functions

histo_t wrap_phi(histo_t phi)
{
  if(phi<0) return wrap_phi(phi+2*M_PI);
  else if (phi>=2*M_PI) return wrap_phi(phi-2*M_PI);
  return phi;
}

size_t get_bin_index(histo_t x,Bin bin)
{
  if (bin.type==BIN_LIN) return (size_t)((x-bin.min)/bin.step);
  else if (bin.type==BIN_LOG) return (size_t)((log10(x)-bin.log10min)/bin.step);
  else if (bin.type==BIN_CUSTOM) return get_dichotomy_index(x,bin.edges,0,bin.n_bin);
  fprintf(stderr,"Invalid bin type\n");
  exit(1);
}

histo_t get_bin_mid(size_t ibin,Bin bin)
{
  if (bin.type==BIN_LIN) return (ibin+0.5)*bin.step+bin.min;
  else if ((bin.type==BIN_LOG)||(bin.type==BIN_CUSTOM)) {
    if (ibin>=bin.n_bin) return bin.edges[ibin];
    return (bin.edges[ibin]+bin.edges[ibin+1])/2.;
  }
  fprintf(stderr,"Invalid bin type\n");
  exit(1);
}

histo_t get_bin_edge(size_t ibin,Bin bin)
{
  if (bin.type==BIN_LIN) return ibin*bin.step+bin.min;
  else if (bin.type==BIN_LOG) return pow(10.,ibin*bin.step+bin.min);
  else if (bin.type==BIN_CUSTOM) return bin.edges[ibin];
  fprintf(stderr,"Invalid bin type\n");
  exit(1);
}

size_t ravel_index(size_t* ind,size_t* shape,size_t n_dim)
{
  size_t index=ind[0];
  size_t prod=1;
  size_t ii;
  for (ii=1;ii<n_dim;ii++) {
    prod *= shape[ii-1];
    index += prod*ind[ii];
  }
  return index;
}

void unravel_index(size_t ind,size_t* shape,size_t n_dim,size_t* index)
{
  index[0] = ind%shape[0];
  size_t ii;
  for (ii=1;ii<n_dim;ii++) {
    ind = (ind-index[ii-1])/shape[ii-1];
    index[ii] = ind%shape[ii];
  }
}


size_t get_dichotomy_index(histo_t x,histo_t *tab,size_t min,size_t max)
{
  if (min==max-1) return min;
  else {
    size_t ind=(min+max)/2;
    if (x<tab[ind]) return get_dichotomy_index(x,tab,min,ind);
    else if (x>tab[ind]) return get_dichotomy_index(x,tab,ind,max);
    return ind;
  }
}

void legendre_all(histo_t mu,histo_t mu2,histo_t* leg)
{
  histo_t mu3 = mu2*mu;
  histo_t mu4 = mu2*mu2;
  histo_t mu5 = mu4*mu;
  histo_t mu6 = mu4*mu2;
  histo_t mu7 = mu6*mu;
  histo_t mu8 = mu6*mu2;
  //histo_t mu10 = mu8*mu2;
  //histo_t mu12 = mu10*mu2;
  leg[0] = 1.;
  leg[1] = mu;
  leg[2] = 0.5*(3.*mu2-1.);
  leg[3] = 0.5*(5.*mu3-3.*mu);
  leg[4] = 1./8.*(35.*mu4-30.*mu2+3.);
  leg[5] = 1./8.*(63.*mu5-70.*mu3+15.*mu);
  leg[6] = 1./16.*(231.*mu6-315.*mu4+105*mu2-5.);
  leg[7] = 1./16.*(429.*mu7-693.*mu5+315*mu3-35.*mu);
  leg[8] = 1./128.*(6435.*mu8-12012.*mu6+6930.*mu4-1260.*mu2+35.);
  //leg[5] = 1./256.*(46189.*mu10-109395.*mu8+90090.*mu6-30030.*mu4+3465.*mu2-63.);
  //leg[6] = 1./1024.*(676039.*mu12-1939938.*mu10+2078505.*mu8-1021020.*mu6+225225.*mu4-18018.*mu2+231.);
}

void legendre_even(histo_t mu2,histo_t* leg)
{
  histo_t mu4 = mu2*mu2;
  histo_t mu6 = mu4*mu2;
  histo_t mu8 = mu6*mu2;
  histo_t mu10 = mu8*mu2;
  histo_t mu12 = mu10*mu2;
  leg[0] = 1.;
  leg[1] = 0.5*(3.*mu2-1.);
  leg[2] = 1./8.*(35.*mu4-30.*mu2+3.);
  leg[3] = 1./16.*(231.*mu6-315.*mu4+105*mu2-5.);
  leg[4] = 1./128.*(6435.*mu8-12012.*mu6+6930.*mu4-1260.*mu2+35.);
  leg[5] = 1./256.*(46189.*mu10-109395.*mu8+90090.*mu6-30030.*mu4+3465.*mu2-63.);
  leg[6] = 1./1024.*(676039.*mu12-1939938.*mu10+2078505.*mu8-1021020.*mu6+225225.*mu4-18018.*mu2+231.);
}

void legendre_odd(histo_t mu,histo_t mu2,histo_t* leg)
{
  histo_t mu3 = mu2*mu;
  histo_t mu5 = mu3*mu2;
  histo_t mu7 = mu5*mu2;
  histo_t mu9 = mu7*mu2;
  histo_t mu11 = mu9*mu2;
  leg[0] = mu;
  leg[1] = 0.5*(5.*mu3-3.*mu);
  leg[2] = 1./8.*(63.*mu5-70.*mu3+15.*mu);
  leg[3] = 1./16.*(429.*mu7-693.*mu5+315*mu3-35.*mu);
  leg[4] = 1./128.*(12155.*mu9-25740.*mu7+18018*mu5-4620.*mu3+315.*mu);
  leg[5] = 1./256.*(88179.*mu11-230945.*mu9+218790*mu7-90090.*mu5+15015.*mu3-693.*mu);
}

void legendre(histo_t dist,histo_t* leg,MULTI_TYPE type) {

  if (type==MULTI_ALL) {
    legendre_all(dist,dist*dist,leg);
  }
  else if (type==MULTI_EVEN) {
    legendre_even(dist*dist,leg);
  }
  else if (type==MULTI_ODD) {
    legendre_odd(dist,dist*dist,leg);
  }
}

void timer(size_t i)
{
  /////
  // Timing routine
  // timer(0) -> initialize relative clock
  // timer(1) -> read relative clock
  // timer(2) -> read relative clock and initialize it afterwards
  // timer(4) -> initialize absolute clock
  // timer(5) -> read absolute clock
#ifdef _HAVE_OMP
  if(i==0)
    relbeg=omp_get_wtime();
  else if(i==1) {
    relend=omp_get_wtime();
    printf(" - relative time ellapsed %.1f ms\n",1000*(relend-relbeg));
  }
  else if(i==2) {
    relend=omp_get_wtime();
    printf(" - relative time ellapsed %.1f ms\n",1000*(relend-relbeg));
    relbeg=omp_get_wtime();
  }
  else if(i==4)
    absbeg=omp_get_wtime();
  else if(i==5) {
    absend=omp_get_wtime();
    printf(" - total time ellapsed %.1f ms \n",1000*(absend-absbeg));
  }
#else //_HAVE_OMP
  int diff;

  if(i==0)
    relbeg=time(NULL);
  else if(i==1) {
    relend=time(NULL);
    diff=(int)(difftime(relend,relbeg));
    printf(" - relative time ellapsed %02d:%02d:%02d \n",
     diff/3600,(diff/60)%60,diff%60);
  }
  else if(i==2) {
    relend=time(NULL);
    diff=(size_t)(difftime(relend,relbeg));
    printf(" - relative time ellapsed %02d:%02d:%02d \n",
     diff/3600,(diff/60)%60,diff%60);
    relbeg=time(NULL);
  }
  else if(i==4)
    absbeg=time(NULL);
  else if(i==5) {
    absend=time(NULL);
    diff=(size_t)(difftime(absend,absbeg));
    printf(" - total time ellapsed %02d:%02d:%02d \n",
     diff/3600,(diff/60)%60,diff%60);
  }
#endif //_HAVE_OMP
}


void error_mem_out(void)
{
  //////
  // Memory shortage handler
  fprintf(stderr,"CUTE: Out of memory!!\n");
  exit(1);
}

void error_open_file(char *fname)
{
  //////
  // Open error handler
  fprintf(stderr,"CUTE: Could not open file %s \n",fname);
  exit(1);
}

///////////////////////////

void write_mesh(Mesh mesh,char *fn)
{
  //////
  // Writes pixel map into file fn, only used for debugging
  FILE *fr;
  fr=fopen(fn,"w");
  if(fr==NULL) error_open_file(fn);
  size_t ibox;
  for (ibox=0;ibox<mesh.n_boxes;ibox++) {
    if(mesh.boxes[ibox].n_obj>0) {
      size_t iobj;
      for (iobj=0;iobj<mesh.boxes[ibox].n_obj;iobj++) {
        size_t idim;
        for (idim=0;idim<dim_pos;idim++) fprintf(fr,"%f ",mesh.boxes[idim].pos[dim_pos*iobj+idim]);
        for (idim=0;idim<dim_weight;idim++) fprintf(fr,"%f ",mesh.boxes[idim].weight[dim_weight*iobj+idim]);
        fprintf(fr,"\n");
      }
    }
  }
  fclose(fr);
}

void write_catalog(Catalog cat,char *fn)
{
  //////
  // Writes catalog into file fn, only used for debugging
  FILE *fr;
  size_t jj;
  fr=fopen(fn,"w");
  if(fr==NULL) error_open_file(fn);
  for (jj=0;jj<cat.n_obj;jj++) {
    size_t kk;
    fprintf(fr,"%zu ",jj);
    for (kk=0;kk<dim_box;kk++) fprintf(fr,"%f ",cat.pos[dim_box*jj+kk]);
    for (kk=0;kk<dim_weight;kk++) fprintf(fr,"%f ",cat.weight[dim_weight*jj+kk]);
    fprintf(fr,"\n");
  }
  fclose(fr);
}

void write_histo(size_t num_histo,histo_t* hist,char *fn)
{
  //////
  // Writes catalog into file fn, only used for debugging
  FILE *fr;
  size_t ii;
  fr=fopen(fn,"w");
  if(fr==NULL) error_open_file(fn);
  for (ii=0;ii<num_histo;ii++) {
    fprintf(fr,"%f\n",hist[ii]);
  }
  fclose(fr);
}
