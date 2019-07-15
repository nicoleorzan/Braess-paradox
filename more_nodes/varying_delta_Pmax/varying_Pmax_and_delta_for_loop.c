#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include "include/time_computing.h"
#include "include/network.h"
#include "include/runge_kutta.h"
#include "include/stability_check.h"

#define steps 1000000
#define additive_steps 1000
#define internal_steps 10
#define Pmax_num 11
#define delta_num 41

double delta_step = 0.1;
double delta_max = 4.1;
double Pmax_step = 0.1;
double Pmax_max = 1.1;
double Pmax_vals[Pmax_num] = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
double delta_vals[delta_num];


int main(){

  for (int i=0; i<delta_num; i++){
    delta_vals[i] = i*0.1;
    fprintf(stdout, "%16.8f ", delta_vals[i]);
  }
  
  double tstart, tstop, ctime=0;
  struct timespec ts;
  
  double *y = (double*) malloc(2 * nodes * sizeof(double));
  for (int i=0; i<2*nodes; i++){
    y[i] = 0;
  }
  bool unstable = 0;
  
  FILE* doc;
  doc = fopen("prova", "w");

  tstart = TCPU_TIME;
  
  Pmax = 0;
  fprintf(doc, "     delta      Pmax      control[0]      control[1]       control[2]      control[3]       control[4]      control[5]       control[6]      control[7]\n");
  
  for (int i=Pmax_num-1; i>=0; i--){

    Pmax = Pmax_vals[i];
    delta = 0;
    //fprintf(stdout, "Pmax = %16.8e \n", Pmax);

    #pragma omp parallel for
    for (int j=delta_num-1; j>=0; j--){

      delta = delta_vals[j];
      //fprintf(stdout, "delta = %16.8e \n", delta);
    
      for (int t=1; t<=steps; t+=internal_steps){
	runge_kutta(y, internal_steps);
      }

      fprintf(doc, "%16.8f %16.8f ", delta, Pmax);
      for (int i=0; i<nodes; i++){
	fprintf(doc, "%16.8e", Pmax*tanh(delta*y[i]));
      }
      fprintf(doc, "\n");

      stability_check(runge_kutta, y, additive_steps, &unstable);

    }

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



