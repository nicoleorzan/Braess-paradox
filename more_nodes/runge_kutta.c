#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include "../include/time_computing.h"
#include "../include/network.h"
#include "../include/runge_kutta.h"
#include "../include/stability_check.h"

#define steps 2500
#define additive_steps 1000
#define internal_steps 10
#define printing_step 10


int main(){
  delta = 0.1/Pmax;
  //fprintf(stdout, "delta = %16.8e\n", delta);
  //fprintf(stdout, "delta*Pmax = %16.8e\n", delta*Pmax);
  double tstart, tstop, ctime=0;
  struct timespec ts;
  
  double *y = (double*) malloc(2 * nodes * sizeof(double));
  for (int i=0; i<2*nodes; i++){
    y[i] = 0;
  }
  bool unstable = 0;
  
  FILE* theta_doc;
  theta_doc = fopen("theta", "w");

  tstart = TCPU_TIME;

  for (int t=1; t<=steps; t+=internal_steps){
    
    runge_kutta(y, internal_steps);

    //printing on file
    fprintf(theta_doc, "%16.8f", t*h);
    for (int i=0; i<2*nodes; i++){
      fprintf(theta_doc, "%16.8e", y[i]);
    }
    fprintf(theta_doc, "\n");
    
   }

  stability_check(runge_kutta, y, additive_steps, &unstable);

  ctime += TCPU_TIME - tstart;
  printf("%g sec \n", ctime);

  fclose(theta_doc);
  memset(theta_doc, 0, sizeof(*theta_doc));
  free(theta_doc);
  memset(y, 0, sizeof(*y));
  free(y);
  
  return 0;
}



