#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include "include/time_computing.h"
#include "include/network.h"
#include "include/runge_kutta_varying_control.h"

#define steps 50000
#define additive_steps 1000
#define internal_steps 10
#define max_error 10e-10

#define delta_step 0.01
#define delta_min 0.0

void print_info(FILE *file){
  fprintf(file, "Usual network with 8 nodes\n");
  fprintf(file, "control form: Pmax*tanh(delta*theta)\n");
  fprintf(file, "Pmax = %f, variation of delta from %f to %f with step %f\n", Pmax, delta, delta_min, delta_step);
  fprintf(file, "Computation of stability reached/unreached in every variation of delta with comparison after %i and %i steps\n", steps, steps+additive_steps);
  fprintf(file, "Computation is stopped when stability is not reached anymore\n");
  fprintf(file, "Error is computed as sum_i ( |(theta_%i[i] - theta_%i[i]| ) \n", steps, steps+additive_steps);
  fprintf(file, "Stability is not reached if error >= %16.8e\n\n", max_error);
}

void stability_check(double* y, bool *unstable, FILE* file){
  
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
    fprintf(file, "delta =%16.8e\n", delta);
    fprintf(file, "error =%16.8e\n", sum);
    fprintf(file, "Stability not reached\n\n");
    *unstable = 1;
  }
  else {
    fprintf(file, "delta =%16.8e\n", delta);
    fprintf(file, "error =%16.8e\n", sum);
    fprintf(file, "Stability reached\n\n");
  }
}


int main(){

  double tstart, tstop, ctime=0;
  struct timespec ts;
  bool unstable = 0;
  
  double *y = (double*) malloc(2 * nodes * sizeof(double));

  FILE* file;
  file = fopen("control_variation", "w");

  print_info(file);
 
  tstart = TCPU_TIME;

  while (delta > delta_min){

    for (int i=0; i<2*nodes; i++){
      y[i] = 0;
    }
  
    for (int t=1; t<=steps; t+=internal_steps){
      runge_kutta(y, internal_steps);
    }
    
    stability_check(y, &unstable, file);
    
    if (unstable == 1) break;
    delta = delta - delta_step;
  }

  ctime += TCPU_TIME - tstart;
  printf("%g sec \n", ctime);

  fclose(file);
  memset(file, 0, sizeof(*file));
  free(file);
  memset(y, 0, sizeof(*y));
  free(y);
  
  return 0;
}

