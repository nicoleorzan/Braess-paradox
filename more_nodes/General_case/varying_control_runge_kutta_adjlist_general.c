#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include "global_vars.h"

#define nodes 8
#define connections 20 //20 or 22
#define alpha 1.0
#define Gamma 0.1
#define Pmax 0.1
#define steps 50000
#define additive_steps 1000
#define internal_steps 10
#define h 0.01
#define hh h*0.5
#define h6 h/6
#define max_error 10e-10

#define delta_step 0.01
#define delta_min 0.0

const double P[nodes] = {-1, 1, 1, 1, -1, -1, 1, -1};
const int AI[nodes+1] = {0, 3, 5, 8, 11, 13, 16, 18, 20};
//const int AI[nodes+1] = {0, 3, 6, 9, 12, 14, 17, 19, 22}; //22 connections
const int AV[connections] = {5, 1, 6, 0, 2, 3, 1, 6, 2, 4, 7, 3, 5, 4, 7, 0, 2, 0, 3, 5};
//const int AV[connections] = {5, 1, 6, 0, 2, 3, 3, 1, 6, 1, 2, 4, 7, 3, 5, 4, 7, 0, 2, 0, 3, 5}; //22 connections
double weights[connections] = {1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03};
//double weights[connections] = {1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03}; //22 connections
//doubling nodes 3-4 capacity -> double weights in positions 5 and 8
double delta = 1.0;


void derivs(double *y, double *dydt){

  double sum = 0;

  for (int i=0; i<nodes; i++){
    
    sum = 0;
    for (int j = AI[i]; j < AI[i+1]; j++){
	sum += weights[j] * sin( y[i] - y[(AV[j])] );
      }
    dydt[i] = y[i+nodes];
    dydt[i+nodes] = -alpha*y[i+nodes] - Gamma*y[i] + P[i] - sum; //Pmax*tanh(delta*y[i]) + P[i] - sum;
    
  }
}

void runge_kutta(double* y, double *yt, double *dym, double *dyt, double *dydx, int rk_steps){

  for (int t=1; t<=rk_steps; t++){
    
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


void stability_check(double* y, double *yt, double *dym, double *dyt, double *dydx, bool *unstable, FILE* control_variation){
  
  double theta_save[nodes];
  double error[nodes];
  double sum = 0.;

  for (int i=0; i<nodes; i++){
    theta_save[i] = y[i];
    error[i] = 0.;
  }
  //fprintf(stdout, "%16.8f\n", y[0]);
  runge_kutta(y, yt, dym, dyt, dydx, additive_steps);
  //fprintf(stdout, "%16.8f\n", y[0]);
  
  for (int i=0; i<nodes; i++){
    error[i] = fabs(theta_save[i] - y[i]);
    sum += error[i];
  }
  if (sum >= max_error) {
    fprintf(control_variation, "delta =%16.8e\n", delta);
    fprintf(control_variation, "error =%16.8e\n", sum);
    fprintf(control_variation, "Stability not reached\n\n");
    *unstable = 1;
  }
  else {
    fprintf(control_variation, "delta =%16.8e\n", delta);
    fprintf(control_variation, "error =%16.8e\n", sum);
    fprintf(control_variation, "Stability reached\n\n");
  }
}


int main(){

  double tstart, tstop, ctime=0;
  struct timespec ts;
  bool unstable = 0;
  
  double *y = (double*) malloc(2 * nodes * sizeof(double));
  double *yt = (double*) malloc(2 * nodes * sizeof(double));
  double *dym = (double*) malloc(2 * nodes * sizeof(double));
  double *dyt = (double*) malloc(2 * nodes * sizeof(double));
  double *dydx = (double*) malloc(2 * nodes * sizeof(double));

  FILE* control_variation;
  control_variation = fopen("control_variation", "w");

  fprintf(control_variation, "Usual network with 8 nodes\n");
  fprintf(control_variation,"control form: Pmax*tanh(delta*theta)\n");
  fprintf(control_variation,"Pmax = %f, variation of delta from %f to %f with step %f\n", Pmax, delta, delta_min, delta_step);
  fprintf(control_variation, "Computation of stability reached/unreached in every variation of delta with comparison after %i and %i steps\n", steps, steps+additive_steps);
  fprintf(control_variation, "Computation is stopped when stability is not reached anymore\n");
  fprintf(control_variation, "Error is computed as sum_i ( |(theta_%i[i] - theta_%i[i]| ) \n", steps, steps+additive_steps);
  fprintf(control_variation, "Stability is not reached if error >= %16.8e\n\n", max_error);

  tstart = TCPU_TIME;

  while (delta > delta_min){

    for (int i=0; i<2*nodes; i++){
      y[i] = 0;
      yt[i] = 0.;
      dym[i] = 0.;
      dyt[i] = 0.;
      dydx[i] = 0.;
    }
  
    for (int t=1; t<=steps; t+=internal_steps){
      runge_kutta(y, yt, dym, dyt, dydx, internal_steps);
    }
    
    stability_check(y, yt, dym, dyt, dydx, &unstable, control_variation);
    
    if (unstable == 1) break;
    delta = delta - delta_step;
  }

  ctime += TCPU_TIME - tstart;
  printf("%g sec \n", ctime);

  fclose(control_variation);
  memset(control_variation, 0, sizeof(*control_variation));
  free(control_variation);
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



