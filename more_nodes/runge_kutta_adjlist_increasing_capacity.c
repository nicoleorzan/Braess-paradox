#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global_vars.h"

const double P[nodes] = {-1, 1, 1, 1, -1, -1, 1, -1};
const int AI[nodes+1] = {0, 3, 5, 8, 11, 13, 16, 18, 20};
const int AV[connections] = {5, 1, 6, 0, 2, 3, 1, 6, 2, 4, 7, 3, 5, 4, 7, 0, 2, 0, 3, 5};
double weights[connections] = {1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03};


double omega_funct(double omega_old, double theta, double somma, int i){
  return -alpha*omega_old - Gamma*theta + P[i] - somma;
}

void printer(double * y, FILE * f){
  fprintf(f, "%16.8e", (y[3]-y[4])/M_PI); //diff nodes 4-5
  fprintf(f, "%16.8e", (y[3]-y[7])/M_PI); //diff nodes 4-8
  fprintf(f, "%16.8e", (y[0]-y[5])/M_PI); //diff nodes 1-6
  fprintf(f, "%16.8e", (y[2]-y[3])/M_PI); //diff nodes 3-4
  fprintf(f, "%16.8e", (y[2]-y[1])/M_PI); //diff nodes 3-2
  fprintf(f, "%16.8e", (y[2]-y[6])/M_PI); //diff nodes 3-7
  fprintf(f, "%16.8e", (y[5]-y[4])/M_PI); //diff nodes 6-5
  fprintf(f, "%16.8e", (y[5]-y[7])/M_PI); //diff nodes 6-8
  fprintf(f, "%16.8e", (y[1]-y[0])/M_PI); //diff nodes 2-1
  fprintf(f, "%16.8e\n", (y[6]-y[0])/M_PI); //diff nodes 7-1
}

void runge_kutta(double* omega, double* theta){
  
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

  int steps = 2500, printing_step = 2400;

  double *theta = (double*) malloc(nodes * sizeof(double));
  double *omega = (double*) malloc(nodes * sizeof(double));

  FILE *theta_doc;
  theta_doc = fopen("theta", "w");
  
  double sum = 0, max_capacity = 1.73, deltaK = 0.1;
  double cap = weights[5];
  int iter = 0;

  // iterate increasing capacity
  while (cap < max_capacity){

    fprintf(theta_doc, "%16.8f", deltaK*iter);
    
    for (int i=0; i<nodes; i++){
      theta[i] = 0;
      omega[i] = 0;
    }
   
    for (int t=1; t<=steps; t++){
      runge_kutta(omega, theta);
    }

    printer(theta, theta_doc);
           
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

