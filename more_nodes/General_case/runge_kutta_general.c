#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "include/time_computing.h"
#include "include/network.h"
#include "include/runge_kutta.h"

#define steps 25000
#define additive_steps 1000
#define internal_steps 10
#define printing_step 10
#define max_error 10e-7

void stability_check(double* y){
  
  double theta_save[nodes];
  double error[nodes];
  double sum = 0.;
  
  for (int i=0; i<nodes; i++){
    theta_save[i] = y[i];
    error[i] = 0.;
  }

  runge_kutta(y, internal_steps);
  
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


int main(){

  double tstart, tstop, ctime=0;
  struct timespec ts;
  
  double *y = (double*) malloc(2 * nodes * sizeof(double));
  for (int i=0; i<2*nodes; i++){
    y[i] = 0;
  }
  /*for (int i=0; i<2*nodes; i++){
    weights[i] = 1.68;
    }*/
  //weights[5] = 1.68*2;
  //weights[8] = 1.68*2;

  FILE* theta_doc;
  theta_doc = fopen("theta", "w");

  tstart = TCPU_TIME;

  for (int t=1; t<=steps; t+=internal_steps){
    
    runge_kutta(y, internal_steps);

    //printing on file
    fprintf(theta_doc, "%16.8f", t*h);
    /*for (int i=0; i<nodes; i++){
      fprintf(theta_doc, "%16.8e", y[i]/M_PI);
      }*/
    for (int i=0; i<2*nodes; i++){
      fprintf(theta_doc, "%16.8e", y[i]);
    }
    fprintf(theta_doc, "\n");
    
   }

  stability_check(y);

  ctime += TCPU_TIME - tstart;
  printf("%g sec \n", ctime);

  fclose(theta_doc);
  memset(theta_doc, 0, sizeof(*theta_doc));
  free(theta_doc);
  memset(y, 0, sizeof(*y));
  free(y);
  
  return 0;
}



