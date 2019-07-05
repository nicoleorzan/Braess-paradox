#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global_vars.h"
#include <sys/resource.h>
#include <sys/times.h>
#include <time.h>


const double P[nodes] = {-1, 1, 1, 1, -1, -1, 1, -1};

double omega_funct(double omega_old, double theta, double somma, int i){
  return -alpha*omega_old - Gamma*theta + P[i] - somma;
}

void runge_kutta(double* omega, double* theta, int *AI, int *AV, double *weights){
  
  double k1[2], k2[2], k3[2], k4[2];
  double sum = 0;
  
  for (int i=0; i<nodes; i++){

    sum = 0;
    for (int j = AI[i]; j < AI[i+1]; j++){
      sum += weights[j] * sin( theta[i] - theta[(AV[j])] ); 
    }
	
    k1[0] = omega_funct(omega[i], theta[i], sum, i);
    k1[1] = omega[i];
    k2[0] = omega_funct(omega[i]+k1[0]*hh, theta[i]+k1[1]*hh, sum, i);
    k2[1] = omega[i]+k1[1]*hh;
    k3[0] = omega_funct(omega[i]+k2[0]*hh, theta[i]+k2[1]*hh, sum, i);
    k3[1] = omega[i]+k2[1]*hh;
    k4[0] = omega_funct(omega[i]+k3[0]*h, theta[i]+k3[1]*hh, sum, i);
    k4[1] = omega[i]+k3[1]*h;

    omega[i] = omega[i] + h6*(k1[0]+2*k2[0]+2*k3[0]+k4[0]);
    theta[i] = theta[i] + h6*(k1[1]+2*k2[1]+2*k3[1]+k4[1]);
  }

}


int main(){

  double tstart, tstop, ctime = 0;
  struct timespec ts;
  
  int steps = 2500, printing_step = 10;
  double *theta = (double*) malloc(nodes * sizeof(double));
  double *omega = (double*) malloc(nodes * sizeof(double));
  for (int i=0; i<nodes; i++){
    theta[i] = 0;
    omega[i] = 0;
  }
  
  // reading from file

  
  int AI[nodes+1],  AV[connections];
  double weights[connections];
 
  char const * ai = "files_to_read/AI.txt";
  char const * av = "files_to_read/AV.txt";
  char const * ww = "files_to_read/weights.txt";
  
  file_reader(AI, AV, weights, ai, av, ww);

  //doubling nodes 3-4 capacity
  //weights[5] = weights[0] * 2;
  //weights[8] = weights[0] * 2;
    
  FILE *theta_doc;
  theta_doc = fopen("theta", "w");


  tstart = TCPU_TIME;

  for (int t=1; t<=steps; t++){
    if (t % printing_step == 0) fprintf(theta_doc, "%16.8f", t*h);

    runge_kutta(omega, theta, AI, AV, weights);

    if (t % printing_step == 0){
      for (int i=0; i<nodes; i++){
	fprintf(theta_doc, "%16.8e", theta[i]/M_PI);
      }
      fprintf(theta_doc, "\n");
    }
  }

  ctime += TCPU_TIME - tstart;
  printf("%g sec \n", ctime);

  fclose(theta_doc);
  memset(theta_doc, 0, sizeof(*theta_doc));
  free(theta_doc);
  memset(theta, 0, sizeof(*theta));
  memset(omega, 0, sizeof(*omega));
  free(theta);  free(omega);
  
  return 0;
}



 
