#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include "include/time_computing.h"
#include "include/network.h"
#include "include/runge_kutta.h"
#include "include/stability_check.h"

#define steps 500000
#define additive_steps 1000
#define internal_steps 10

double delta_step = 0.1;
double delta_max = 4.1;
double Pmax_step = 0.1;
double Pmax_max = 1.1;

int main(){

  double tstart, tstop, ctime=0;
  struct timespec ts;
  
  double *y = (double*) malloc(2 * nodes * sizeof(double));
  for (int i=0; i<2*nodes; i++){
    y[i] = 0;
  }
  bool unstable = 0;
  
  FILE* doc;
  doc = fopen("control_delta", "w");

  tstart = TCPU_TIME;
  
  Pmax = 0;
  fprintf(doc, "delta      Pmax      control[0]      control[1]       control[2]      control[3]       control[4]      control[5]       control[6]      control[7]\n");
  
  while (Pmax <= Pmax_max){

    delta = 0;
    fprintf(stdout, "Pmax = %16.8e \n", Pmax);
    
    while (delta <= delta_max){
      
      fprintf(stdout, "delta = %16.8e \n", delta);
    
      for (int t=1; t<=steps; t+=internal_steps){
	runge_kutta(y, internal_steps);
      }

      fprintf(doc, "%16.8f %16.8f ", delta, Pmax);
      for (int i=0; i<nodes; i++){
	fprintf(doc, "%16.8e", Pmax*tanh(delta*y[i]));
      }
      fprintf(doc, "\n");

      stability_check(runge_kutta, y, additive_steps, &unstable);

      delta = delta + delta_step;
    }

    Pmax = Pmax + Pmax_step;
  }

  ctime += TCPU_TIME - tstart;
  printf("%g sec \n", ctime);

  fclose(doc);
  memset(doc, 0, sizeof(*doc));
  free(doc);
  memset(y, 0, sizeof(*y));
  free(y);
  
  return 0;
}



