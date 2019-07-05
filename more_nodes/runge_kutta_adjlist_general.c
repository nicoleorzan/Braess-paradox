#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "global_vars.h"

const double P[nodes] = {-1, 1, 1, 1, -1, -1, 1, -1};
const int AI[nodes+1] = {0, 3, 5, 8, 11, 13, 16, 18, 20};
const int AV[connections] = {5, 1, 6, 0, 2, 3, 1, 6, 2, 4, 7, 3, 5, 4, 7, 0, 2, 0, 3, 5};
double weights[connections] = {1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03};


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



int main(){

  double tstart, tstop, ctime=0;
  struct timespec ts;
  
  int steps = 2500, printing_step = 10;
  double *y = (double*) malloc(2 * nodes * sizeof(double));
  
  for (int i=0; i<2*nodes; i++){
    y[i] = 0;
  }

  //doubling nodes 3-4 capacity
  //weights[5] = weights[0] * 2;
  //weights[8] = weights[0] * 2;


  FILE *theta_doc;
  theta_doc = fopen("theta", "w");

  double yt[nodes*2] = {0}, dym[nodes*2] = {0}, dyt[nodes*2] = {0}, dydx[nodes*2] = {0};
  for (int i=0; i<nodes*1; i++){
    yt[i] = 0.;
    dym[i] = 0.;
    dyt[i] = 0.;
    dydx[i] = 0.;
  }

  tstart = TCPU_TIME;

  for (int t=1; t<=steps; t++){
    if (t % printing_step == 0) fprintf(theta_doc, "%16.8f", t*h);

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
      
      if (t % printing_step == 0 && (i>=0 && i<nodes)){
	fprintf(theta_doc, "%16.8e", y[i]/M_PI);
      }
    }
          
   if (t % printing_step == 0) fprintf(theta_doc, "\n");
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



