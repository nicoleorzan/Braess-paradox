#ifndef RUNGE_KUTTA
#define RUNGE_KUTTA

#include <math.h>

#define h 0.01
#define hh h*0.5
#define h6 h/6

double yt[2*nodes], dydx[2*nodes], dyt[2*nodes], dym[2*nodes];

void derivs(double *y, double *dydt){

  double sum = 0;

  for (int i=0; i<nodes; i++){
    
    sum = 0;
    for (int j = AI[i]; j < AI[i+1]; j++){
	sum += weights[j] * sin( y[i] - y[(AV[j])] );
      }
    dydt[i] = y[i+nodes];
    dydt[i+nodes] = -alpha*y[i+nodes] - Pmax*tanh(delta*y[i]) + P[i] - sum;
    //dydt[i+nodes] = -alpha*y[i+nodes] - Pmax[i]*tanh(delta[i]*y[i]) + P[i] - sum;
    //dydt[i+nodes] = -alpha*y[i+nodes] - Gamma*y[i] + P[i] - sum;  old code
    
  }
}

void runge_kutta(double* y, int internal_steps){

  for (int t=1; t<=internal_steps; t++){
    
    derivs(y, dydx);

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
      
    }
  }
  
}

#endif
