#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "../time_computing.h"

#define nodes 8
#define connections 20 //20 or 22
#define alpha 1.0
#define Gamma 0.1
#define Pmax 0
#define delta 1.0
#define steps 25000
#define additive_steps 1000
#define internal_steps 10
#define printing_step 10
#define h 0.001
#define hh h*0.5
#define h6 h/6
#define max_error 10e-10
#define max_capacity 1.73
#define deltaK 0.1

const double P[nodes] = {-1, 1, 1, 1, -1, -1, 1, -1};
const int AI[nodes+1] = {0, 3, 5, 8, 11, 13, 16, 18, 20};
//const int AI[nodes+1] = {0, 3, 6, 9, 12, 14, 17, 19, 22}; //22 connections
const int AV[connections] = {5, 1, 6, 0, 2, 3, 1, 6, 2, 4, 7, 3, 5, 4, 7, 0, 2, 0, 3, 5};
//const int AV[connections] = {5, 1, 6, 0, 2, 3, 3, 1, 6, 1, 2, 4, 7, 3, 5, 4, 7, 0, 2, 0, 3, 5}; //22 connections
double weights[connections] = {1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03};
//double weights[connections] = {1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03}; //22 connections
//doubling nodes 3-4 capacity -> double weights in positions 5 and 8

void derivs(double *y, double *dydt){

  double sum = 0;

  for (int i=0; i<nodes; i++){
    
    sum = 0;
    for (int j = AI[i]; j < AI[i+1]; j++){
	sum += weights[j] * sin( y[i] - y[(AV[j])] );
      }
    dydt[i] = y[i+nodes];
    dydt[i+nodes] = -alpha*y[i+nodes] - Pmax*tanh(delta*y[i]) + P[i] - sum; //Gamma*y[i] + P[i] - sum;
    
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
    fprintf(stdout, "error =%16.8e\n", sum);
    fprintf(stdout, "Stability not reached\n\n");
  }
  else {
    fprintf(stdout, "error =%16.8e\n", sum);
    fprintf(stdout, "Stability reached\n\n");
  }
}

void printer(double * y, FILE * f){
  fprintf(f, "%16.8e", (y[3]-y[4])/M_PI); //diff nodes 4-5
  fprintf(f, "%16.8e", (y[3]-y[7])/M_PI); //diff nodes 4-8
  fprintf(f, "%16.8e", (y[0]-y[5])/M_PI); //diff nodes 1-6
  fprintf(f, "%16.8e", (y[2]-y[3])/M_PI); //diff nodes 3-4
  fprintf(f, "%16.8e", (y[2]-y[1])/M_PI); //diff nodes 3-2
  fprintf(f, "%16.8e", (y[2]-y[6])/M_PI); //diff nodes 3-7
  fprintf(f, "%16.8e", (y[5]-y[4])/M_PI); //diff nodes 6-5
  fprintf(f, "%16.8e", (y[5]-y[7])/M_PI); //diff nodes 6-8
  fprintf(f, "%16.8e", (y[1]-y[0])/M_PI); //diff nodes 2-1
  fprintf(f, "%16.8e", (y[6]-y[0])/M_PI); //diff nodes 7-1
  fprintf(f, "\n");
}


int main(){

  double tstart, tstop, ctime=0;
  struct timespec ts;
  
  double *y = (double*) malloc(2 * nodes * sizeof(double));
  double *yt = (double*) malloc(2 * nodes * sizeof(double));
  double *dym = (double*) malloc(2 * nodes * sizeof(double));
  double *dyt = (double*) malloc(2 * nodes * sizeof(double));
  double *dydx = (double*) malloc(2 * nodes * sizeof(double));
  int iter = 0;
  double cap = weights[5];

  FILE* capacity_doc;
  capacity_doc = fopen("increasing_capacity", "w");

  tstart = TCPU_TIME;
 
  while (cap < max_capacity){

    fprintf(capacity_doc, "%16.8f", deltaK*iter);
    
    for (int i=0; i<2*nodes; i++){
      y[i] = 0;
      yt[i] = 0.;
      dym[i] = 0.;
      dyt[i] = 0.;
      dydx[i] = 0.;
    }

    for (int t=1; t<=steps; t+=internal_steps){  
      runge_kutta(y, yt, dym, dyt, dydx);
    }
    //stability_check(y, yt, dym, dyt, dydx);

    printer(y, capacity_doc);

    weights[5] += deltaK;
    weights[8] += deltaK;
    iter += 1;
    cap += deltaK;
  }

  ctime += TCPU_TIME - tstart;
  printf("%g sec \n", ctime);

  fclose(capacity_doc);
  memset(capacity_doc, 0, sizeof(*capacity_doc));
  free(capacity_doc);
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



