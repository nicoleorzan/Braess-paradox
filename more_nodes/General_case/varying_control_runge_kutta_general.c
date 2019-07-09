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
#define max_error 10e-7

#define Pmin 0.0
#define Pstep 0.01

void print_info(FILE *file){
  fprintf(file, "Usual network with 8 nodes\n");
  fprintf(file, "control form: Pmax*tanh(delta*theta)\n");
  fprintf(file, "variation of Pmax from %f to %f with step %f\n", Pmax, Pmin, Pstep);
  fprintf(file, "and modifying delta as 0.1/Pmax to keep slope stable to 0.1 (in approx) \n");
  fprintf(file, "Computation of stability reached/unreached in every variation of Pmax with comparison after %i and %i steps\n", steps, steps+additive_steps);
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
    fprintf(file, "Pmax =%16.8e\n", Pmax);
    fprintf(file, "delta =%16.8e\n", delta);
    fprintf(file, "slope (Pmax*delta) =%16.8e\n", Pmax*delta);
    fprintf(file, "error =%16.8e\n", sum);
    fprintf(file, "Stability not reached\n\n");
    *unstable = 1;
  }
  else {
    fprintf(file, "Pmax =%16.8e\n", Pmax);
    fprintf(file, "delta =%16.8e\n", delta);
    fprintf(file, "slope (Pmax*delta) =%16.8e\n", Pmax*delta);
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
  FILE* f2;
  file = fopen("control_variation", "w");
  print_info(file);
  
  f2 = fopen("Pmax_delta", "w");
  fprintf(f2, "Pmax       delta\n");
   
  tstart = TCPU_TIME;

  Pmax = 1.0;
  delta = 0.1/Pmax;

  while (Pmax >= Pmin){
    for (int i=0; i<2*nodes; i++){
      y[i] = 0;
    }
  
    for (int t=1; t<=steps; t+=internal_steps){
      runge_kutta(y, internal_steps);
    }
    fprintf(f2, "%16.8e %16.8e\n", Pmax, delta);
    
    stability_check(y, &unstable, file);
    
    if (unstable == 1) break;

    Pmax = Pmax - Pstep;
    delta = 0.1/Pmax;
  }

  ctime += TCPU_TIME - tstart;
  printf("%g sec \n", ctime);

  fclose(file);
  memset(file, 0, sizeof(*file));
  free(file);
  fclose(f2);
  memset(f2, 0, sizeof(*f2));
  free(f2);
  memset(y, 0, sizeof(*y));
  free(y);
  
  return 0;
}

