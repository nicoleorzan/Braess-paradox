#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include "global_vars.h"

const double P[nodes] = {-1, 1, 1, 1, -1, -1, 1, -1};
const int AI[nodes+1] = {0, 3, 5, 8, 11, 13, 16, 18, 20};
const int AV[connections] = {5, 1, 6, 0, 2, 3, 1, 6, 2, 4, 7, 3, 5, 4, 7, 0, 2, 0, 3, 5};
double weights[connections] = {1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03};
double delta = 1;


double omega_funct(double omega_old, double theta, double somma, int i){
  //return -alpha*omega_old - Pmax*delta*theta + P[i] - somma;
  return -alpha*omega_old - Pmax*tanh(delta*theta) + P[i] - somma;
}


void runge_kutta(double* omega, double* theta, int internal_step){

  double k1[2], k2[2], k3[2], k4[2];
  double sum = 0;
    
  for (int t=0; t<internal_step; t++){
    
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

}


void stability_check(double* omega, double* theta, bool* unstable, FILE* control_variation, double max_error){
  
  double theta_save[nodes];
  double error[nodes];
  double sum = 0.;
  
  for (int i=0; i<nodes; i++){
    theta_save[i] = theta[i];
    error[i] = 0.;
  }

  int additive_steps = 100;
  runge_kutta(omega, theta, additive_steps);
  for (int i=0; i<nodes; i++){
    error[i] = fabs(theta_save[i] - theta[i]);
    sum += error[i];
  }
  if (sum >= max_error) {
    fprintf(control_variation, "delta =%16.8e\n", delta);
    fprintf(control_variation, "error =%16.8e\n", sum);
    fprintf(control_variation, "Stability not reached\n\n");
    *unstable = 1;
  }
  else {
    fprintf(control_variation, "delta =%16.8e\n", delta);
    fprintf(control_variation, "error =%16.8e\n", sum);
    fprintf(control_variation, "Stability reached\n\n");
  }
  
}




int main(){
 
  double tstart, tstop, ctime = 0;
  struct timespec ts;
  
  int steps = 25000, internal_step = 10, printing_step = 10;
  double *theta = (double*) malloc(nodes * sizeof(double));
  double *omega = (double*) malloc(nodes * sizeof(double));
  for (int i=0; i<nodes; i++){
    theta[i] = 0;
    omega[i] = 0;
  }
      
  FILE *control_variation;
  control_variation = fopen("control_variation", "w");
  fprintf(control_variation,"control form: Pmax*tanh(delta*theta)\n");
  
  double delta_step = 0.01, delta_min = 0.5;
  bool unstable = 0;
  double max_error = 10e-10;

  fprintf(control_variation,"Pmax = %f, variation of delta from %f to %f with step %f\n", Pmax, delta, delta_min, delta_step);
  fprintf(control_variation, "Computation of stability reached/unreached in every variation of delta with comparison after 25000 and 25100 steps\n");
  fprintf(control_variation, "Computation is stopped when stability is not reached anymore\n");
  fprintf(control_variation, "Error is computed as sum_i ( |(theta_old[i] - theta_new[i]| ) \n");
  fprintf(control_variation, "Stability is not reached if error >= %16.8e\n\n", max_error);

  tstart = TCPU_TIME;

  while (delta > delta_min){
    
    for (int i=0; i<nodes; i++){
      theta[i] = 0;
      omega[i] = 0;
    }
    
    for (int t=1; t<=steps; t+=internal_step){
      runge_kutta(omega, theta, internal_step);
    }

    stability_check(omega, theta, &unstable, control_variation, max_error);
    if (unstable == 1) break;
    delta = delta - delta_step;
  }

  ctime += TCPU_TIME - tstart;
  printf("%g sec \n", ctime);

  fclose(control_variation);
  memset(control_variation, 0, sizeof(*control_variation));
  free(control_variation);
  memset(theta, 0, sizeof(*theta));
  memset(omega, 0, sizeof(*omega));
  free(theta);  free(omega);
  
  return 0;
}



 
