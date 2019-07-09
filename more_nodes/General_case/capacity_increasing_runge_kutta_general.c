#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "include/time_computing.h"
#include "include/network.h"
#include "include/runge_kutta_varying_control.h"

#define steps 25000
#define additive_steps 1000
#define internal_steps 10
#define max_error 10e-10

#define max_capacity 1.73
#define deltaK 0.01


void stability_check(double* y){
    
  double theta_save[nodes];
  double error[nodes];
  double sum = 0.;
  
  for (int i=0; i<nodes; i++){
    theta_save[i] = y[i];
    error[i] = 0.;
  }

  runge_kutta(y, additive_steps);
  
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
  int iter = 0;
  double cap = weights[5];

  FILE* capacity_doc;
  capacity_doc = fopen("increasing_capacity", "w");

  tstart = TCPU_TIME;

  for (int i=0; i<2*nodes; i++){
    y[i] = 0;
  }

  while (cap < max_capacity){
    fprintf(stdout, "%f\n", cap);
    fprintf(capacity_doc, "%16.8f", deltaK*iter);

    for (int t=1; t<=steps; t+=internal_steps){  
      runge_kutta(y, internal_steps);
    }
    stability_check(y);

    printer(y, capacity_doc);

    weights[5] += deltaK;
    weights[8] += deltaK;
    iter += 1;
    cap += deltaK;
    /*fprintf(stdout, "%16.8e \n", weights[5]);
    for (int i=0; i<2*nodes; i++){
      fprintf(stdout, "%16.8e ", y[i]);
    }
    fprintf(stdout, "\n");*/
  }

  ctime += TCPU_TIME - tstart;
  printf("%g sec \n", ctime);

  fclose(capacity_doc);
  memset(capacity_doc, 0, sizeof(*capacity_doc));
  free(capacity_doc);
  memset(y, 0, sizeof(*y));
  free(y);
    
  return 0;
}



