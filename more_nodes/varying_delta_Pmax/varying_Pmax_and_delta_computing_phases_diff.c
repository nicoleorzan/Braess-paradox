#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include "include/time_computing.h"
#include "include/network.h"
#include "include/stability_check.h"
#include "include/runge_kutta_varying_control.h"

#define steps 110000
#define additive_steps 1000
#define internal_steps 10
//#define max_error 10e-10

#define Pmin 0.0
#define Pstep 0.01

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
  fprintf(f, "%16.8e ", Pmax);
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
  bool unstable = 0;
  
  double *y = (double*) malloc(2 * nodes * sizeof(double));
  for (int i=0; i<2*nodes; i++){
    y[i] = 0;
  }

  FILE * f;
  //file = fopen("control_variation", "w");
  //print_info(file);

  f = fopen("theta", "w");
 
  tstart = TCPU_TIME;


  Pmax = 1.0;
  delta = 0.1/Pmax;

  while (Pmax >= Pmin){
    for (int t=1; t<=steps; t+=internal_steps){
      runge_kutta(y, internal_steps);
    }

    printer(y, f);
    
    stability_check(runge_kutta, y, additive_steps, &unstable);

    if (unstable == 1)  break;

    Pmax = Pmax - Pstep;
    delta = 0.1/Pmax;

  }

  ctime += TCPU_TIME - tstart;
  printf("%g sec \n", ctime);

  //fclose(file);
  //memset(file, 0, sizeof(*file));
  //free(file);
  fclose(f);
  memset(f, 0, sizeof(*f));
  free(f);
  memset(y, 0, sizeof(*y));
  free(y);
  
  return 0;
}

