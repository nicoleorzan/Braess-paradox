#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include "../include/time_computing.h"
#include "../include/network.h"
#include "../include/stability_check.h"
#include "../include/runge_kutta.h"

#define steps 190000
#define additive_steps 1000
#define internal_steps 10
//#define max_error 10e-10

#define Pmin 0.0
#define Pstep 0.01
#define Pmax_fixed 0.2
#define delta_min 0.0
#define delta_step 0.1
#define delta_max 4.0

void print_info(FILE *file){
  fprintf(file, "control form: Pmax*tanh(delta*theta)\n");
  fprintf(file, "variation of Pmax from %f to %f with step %f\n", Pmax, Pmin, Pstep);
  fprintf(file, "and modifying Pmax as 0.1/delta to keep slope stable to 0.1 (in approx) \n");
  fprintf(file, "Computation of stability reached/unreached in every variables update with comparison after %i and %i steps\n", steps, steps+additive_steps);
  fprintf(file, "Computation is stopped when stability is not reached anymore\n");
  fprintf(file, "Error is computed as sum_i ( |(theta_%i[i] - theta_%i[i]| ) \n", steps, steps+additive_steps);
  fprintf(file, "Stability is not reached if error >= %16.8e\n\n", max_error);
}

void printer(double * y, FILE * f){
  fprintf(f, "%16.8e ", delta);
  fprintf(f, "%16.8e ", Pmax);
  for (int i=0; i<nodes; i++){
    fprintf(f, "%16.8e ", -Pmax*tanh(delta*y[i]));
  }
  fprintf(f, "\n");
}



int main(){

  double tstart, tstop, ctime=0;
  struct timespec ts;
  bool unstable = 0;
  
  double *y = (double*) malloc(2 * nodes * sizeof(double));
  for (int i=0; i<2*nodes; i++){
    y[i] = 0;
  }

  FILE * f;
  f = fopen("lunedi", "w");
  fprintf(f, "     delta     Pmax    control[0]      control[1]       control[2]      control[3]     control[4]      control[5]       control[6]      control[7]\n");
  
  tstart = TCPU_TIME;

  Pmax = Pmax_fixed;
  delta = delta_max;

  while (delta >= delta_min){
    for (int t=1; t<=steps; t+=internal_steps){
      runge_kutta(y, internal_steps);
    }

    printer(y, f);
    
    stability_check(runge_kutta, y, additive_steps, &unstable);

    if (unstable == 1)  break;

    delta = delta - delta_step;
  }

  ctime += TCPU_TIME - tstart;
  printf("%g sec \n", ctime);

  fclose(f);
  memset(f, 0, sizeof(*f));
  free(f);
  memset(y, 0, sizeof(*y));
  free(y);
  
  return 0;
}

