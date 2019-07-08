#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "global_vars.h"

#define alpha 1
#define Gamma 0
#define steps 25000
#define additive_steps 1000
#define internal_steps 10
#define printing_step 10
#define h 0.01
#define hh h*0.5
#define h6 h/6
#define max_error 10e-10

const double P[nodes] = {-1, 1, 1, 1, -1, -1, 1, -1};
const int AI[nodes+1] = {0, 3, 5, 8, 11, 13, 16, 18, 20};
const int AV[connections] = {5, 1, 6, 0, 2, 3, 1, 6, 2, 4, 7, 3, 5, 4, 7, 0, 2, 0, 3, 5};
double weights[connections] = {1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03};
//doubling nodes 3-4 capacity -> double weights in positions 5 and 8


void derivs(double *y, double *dydt){

  double sum = 0;

  for (int i=0; i<nodes; i++){
    
    sum = 0;
    for (int j = AI[i]; j < AI[i+1]; j++){
	sum += weights[j] * sin( y[i] - y[(AV[j])] );
      }
    dydt[i] = y[i+nodes];
    dydt[i+nodes] = -alpha*y[i+nodes] - Gamma*y[i] + P[i] - sum;
    
  }
}

void runge_kutta(double* y, double *yt, double *dym, double *dyt, double *dydx){

  for (int t=1; t<=internal_steps; t++){
    
    derivs(yt, dydx);

    for (int i=0; i<2*nodes; i++){
      yt[i] = y[i] + hh*dydx[i];
    }
    derivs(yt, dyt);
  
    for (int i=0; i<2*nodes; i++){
      yt[i] = y[i] + hh*dyt[i];
    }
    derivs(yt, dym);

    for (int i=0; i<2*nodes; i++) {
      yt[i] = y[i] + h*dym[i];
      dym[i] += dyt[i];
    }
    derivs(yt, dyt);

    for (int i=0; i<2*nodes; i++) {
      y[i] = y[i] + h6*(dydx[i]+dyt[i]+2.0*dym[i]);
      
    }
  }
  
}


void stability_check(double* y, double *yt, double *dym, double *dyt, double *dydx){
  
  double theta_save[nodes];
  double error[nodes];
  double sum = 0.;

  FILE* stability;
  stability = fopen("stability", "w");
  
  for (int i=0; i<nodes; i++){
    theta_save[i] = y[i];
    error[i] = 0.;
  }

  runge_kutta(y, yt, dym, dyt, dydx);
  
  for (int i=0; i<nodes; i++){
    error[i] = fabs(theta_save[i] - y[i]);
    sum += error[i];
  }
  if (sum >= max_error) {
    //fprintf(stability, "delta =%16.8e\n", delta);
    fprintf(stability, "error =%16.8e\n", sum);
    fprintf(stability, "Stability not reached\n\n");
    //*unstable = 1;
  }
  else {
    //fprintf(stability, "delta =%16.8e\n", delta);
    fprintf(stability, "error =%16.8e\n", sum);
    fprintf(stability, "Stability reached\n\n");
  }
}


int main(){

  double tstart, tstop, ctime=0;
  struct timespec ts;
  
  double *y = (double*) malloc(2 * nodes * sizeof(double));
  double *yt = (double*) malloc(2 * nodes * sizeof(double));
  double *dym = (double*) malloc(2 * nodes * sizeof(double));
  double *dyt = (double*) malloc(2 * nodes * sizeof(double));
  double *dydx = (double*) malloc(2 * nodes * sizeof(double));
  for (int i=0; i<2*nodes; i++){
    y[i] = 0;
    yt[i] = 0.;
    dym[i] = 0.;
    dyt[i] = 0.;
    dydx[i] = 0.;
  }

  FILE* theta_doc;
  theta_doc = fopen("theta", "w");

  tstart = TCPU_TIME;

  for (int t=1; t<=steps; t+=internal_steps){
    
    runge_kutta(y, yt, dym, dyt, dydx);

    //printing on file
    fprintf(theta_doc, "%16.8f", t*h);
    for (int i=0; i<nodes; i++){
      fprintf(theta_doc, "%16.8e", y[i]/M_PI);
    }
    fprintf(theta_doc, "\n");
    
   }

  stability_check(y, yt, dym, dyt, dydx);

  ctime += TCPU_TIME - tstart;
  printf("%g sec \n", ctime);

  fclose(theta_doc);
  memset(theta_doc, 0, sizeof(*theta_doc));
  free(theta_doc);
  memset(y, 0, sizeof(*y));
  free(y);
  memset(yt, 0, sizeof(*yt));
  free(yt);
  memset(dym, 0, sizeof(*dym));
  free(dym);
  memset(dyt, 0, sizeof(*dyt));
  free(dyt);
  memset(dydx, 0, sizeof(*dydx));
  free(dydx);
  
  return 0;
}



