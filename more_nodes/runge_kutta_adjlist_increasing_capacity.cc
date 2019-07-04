#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global_vars.h"

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

  int steps = 2500;

  double *theta = (double*) malloc(nodes * sizeof(double));
  double *omega = (double*) malloc(nodes * sizeof(double));

  // reading from file

  int AI[nodes+1],  AV[connections];
  double weights[connections];
  
  char const * ai = "files_to_read/AI.txt";
  char const * av = "files_to_read/AV.txt";
  char const * ww = "files_to_read/weights.txt";
  
  file_reader(AI, AV, weights, ai, av, ww);
  
  // integration using runge-kutta method of 4th order

  FILE *theta_doc;
    
  theta_doc = fopen("theta", "w");
  
  double sum = 0, max_capacity = 1.73, deltaK = 0.1;
  double cap = weights[5];
  int iter = 0, printing_step = 2400;
  
  while (cap < max_capacity){

    fprintf(theta_doc, "%16.8f", deltaK*iter);
    
    for (int i=0; i<nodes; i++){
      theta[i] = 0;
      omega[i] = 0;
    }
   
    for (int t=1; t<=steps; t++){
      runge_kutta(omega, theta, AI, AV, weights);
    }
      
    fprintf(theta_doc, "%16.8e", (theta[3]-theta[4])/M_PI); //diff nodes 4-5
    fprintf(theta_doc, "%16.8e", (theta[3]-theta[7])/M_PI); //diff nodes 4-8
    fprintf(theta_doc, "%16.8e", (theta[0]-theta[5])/M_PI); //diff nodes 1-6
    fprintf(theta_doc, "%16.8e", (theta[2]-theta[3])/M_PI); //diff nodes 3-4
    fprintf(theta_doc, "%16.8e", (theta[2]-theta[1])/M_PI); //diff nodes 3-2
    fprintf(theta_doc, "%16.8e", (theta[2]-theta[6])/M_PI); //diff nodes 3-7
    fprintf(theta_doc, "%16.8e", (theta[5]-theta[4])/M_PI); //diff nodes 6-5
    fprintf(theta_doc, "%16.8e", (theta[5]-theta[7])/M_PI); //diff nodes 6-8
    fprintf(theta_doc, "%16.8e", (theta[1]-theta[0])/M_PI); //diff nodes 2-1
    fprintf(theta_doc, "%16.8e\n", (theta[6]-theta[0])/M_PI); //diff nodes 7-1
       
    weights[5] += deltaK;
    weights[8] += deltaK;
    iter += 1;
    cap += deltaK;
   
  }
   

  fclose(theta_doc);
  memset(theta_doc, 0, sizeof(*theta_doc));
  free(theta_doc);
  memset(theta, 0, sizeof(*theta));
  memset(omega, 0, sizeof(*omega));
  free(theta);  free(omega);
  
  return 0;
}

