#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "global_vars.h"

#define alpha 1
#define Gamma 0
#define steps 25000
#define additive_steps 1000
#define internal_steps 10
#define printing_step 10
#define h 0.01
#define hh h*0.5
#define h6 h/6

#define delta_step 0.01
#define delta_min 0.5
#define max_error 10e-10

const double P[nodes] = {-1, 1, 1, 1, -1, -1, 1, -1};
const int AI[nodes+1] = {0, 3, 5, 8, 11, 13, 16, 18, 20};
const int AV[connections] = {5, 1, 6, 0, 2, 3, 1, 6, 2, 4, 7, 3, 5, 4, 7, 0, 2, 0, 3, 5};
double weights[connections] = {1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03};
//doubling nodes 3-4 capacity -> double weights in positions 5 and 8


void derivs(double *y, double *dydt){

  double sum = 0;

  for (int i=0; i<nodes; i++){
    
    sum = 0;
    for (int j = AI[i]; j < AI[i+1]; j++){
	sum += weights[j] * sin( y[i] - y[(AV[j])] );
      }
    dydt[i] = y[i+nodes];
    dydt[i+nodes] = -alpha*y[i+nodes] - Gamma*y[i] + P[i] - sum;
    
  }
}

void runge_kutta(double* y, FILE* theta_doc){

  double yt[nodes*2] = {0}, dym[nodes*2] = {0}, dyt[nodes*2] = {0}, dydx[nodes*2] = {0};
  for (int i=0; i<nodes*1; i++){
    yt[i] = 0.;
    dym[i] = 0.;
    dyt[i] = 0.;
    dydx[i] = 0.;
  }

  for (int t=1; t<=internal_steps; t++){
    //if (t % printing_step == 0) fprintf(theta_doc, "%16.8f", t*h);
    
    derivs(yt, dydx);

    for (int i=0; i<2*nodes; i++){
      yt[i] = y[i] + hh*dydx[i];
    }
    derivs(yt, dyt);
  
    for (int i=0; i<2*nodes; i++){
      yt[i] = y[i] + hh*dyt[i];
    }
    derivs(yt, dym);

    for (int i=0; i<2*nodes; i++) {
      yt[i] = y[i] + h*dym[i];
      dym[i] += dyt[i];
    }
    derivs(yt, dyt);

    for (int i=0; i<2*nodes; i++) {
      y[i] = y[i] + h6*(dydx[i]+dyt[i]+2.0*dym[i]);
      
      /*if (t % printing_step == 0 && (i>=0 && i<nodes)){
	fprintf(theta_doc, "%16.8e", y[i]/M_PI);
	}*/
    }
          
    //if (t % printing_step == 0) fprintf(theta_doc, "\n");
  }
  
}

void stability_check(double* y, bool* unstable, FILE* stability_check){
  
  double theta_save[nodes];
  double error[nodes];
  double sum = 0.;
  
  for (int i=0; i<nodes; i++){
    theta_save[i] = y[i];
    error[i] = 0.;
  }

  runge_kutta(y, additive_steps);
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

  double tstart, tstop, ctime=0;
  struct timespec ts;
  
  double *y = (double*) malloc(2 * nodes * sizeof(double));
  for (int i=0; i<2*nodes; i++){
    y[i] = 0;
  }

  FILE *theta_doc;
  theta_doc = fopen("theta", "w");

  tstart = TCPU_TIME;

  for (int t=1; t<=steps; t+=internal_steps){
    
    fprintf(theta_doc, "%16.8f", t*h);

    runge_kutta(y, theta_doc);
    
    for (int i=0; i<nodes; i++){
      fprintf(theta_doc, "%16.8e", y[i]/M_PI);
    }
    fprintf(theta_doc, "\n");
    
   }

  ctime += TCPU_TIME - tstart;
  printf("%g sec \n", ctime);

  fclose(theta_doc);
  memset(theta_doc, 0, sizeof(*theta_doc));
  free(theta_doc);
  memset(y, 0, sizeof(*y));
  free(y);
  
  return 0;
}



