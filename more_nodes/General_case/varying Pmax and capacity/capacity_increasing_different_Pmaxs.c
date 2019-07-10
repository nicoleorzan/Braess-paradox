#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "../include/time_computing.h"
#include "../include/network_Pmaxs.h"
#include "../include/runge_kutta_varying_control_Pmaxs.h"

#define steps 70000
#define additive_steps 1000
#define internal_steps 10
#define max_error 10e-10

#define max_capacity 2.5//1.73
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

  for (int i=0; i<2*nodes; i++){
    fprintf(f, "%16.8e ",Pmax[i]*tanh(delta[i]*y[i]) );
  }
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
  
  for (int i=0; i< nodes; i++){
    if (i==1 || i==2 || i==3 || i==6){
      Pmax[i] = 0.2;
      delta[i] = 0.1/Pmax[i];
    }
    else Pmax[i] = 0.1;
  }
  
  for (int i=0; i<2*nodes; i++){
    y[i] = 0;
  }

  while (cap <= max_capacity){
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



