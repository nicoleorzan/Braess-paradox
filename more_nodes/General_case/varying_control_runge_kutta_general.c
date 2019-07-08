#include <stdlib.h>
#include <stdio.h>
#include <string.h>
//#include <math.h>
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

void stability_check(double* y, bool *unstable, FILE* control_variation){
  
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
    }
  
    for (int t=1; t<=steps; t+=internal_steps){
      runge_kutta(y, internal_steps);
    }
    
    stability_check(y, &unstable, control_variation);
    
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
  
  return 0;
}



