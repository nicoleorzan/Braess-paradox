#include <stdlib.h>
#include <stdio.h>
#include <math.h>

# define nodes 8

void print_mat(int *matrix, int N, int M){
  for (int i=0; i<N; i++){
      for (int j=0; j<M; j++){
	fprintf(stdout, "%d ", matrix[i*M+j]);
      }
      fprintf(stdout, "\n");
  }  
}

double omega_funct(double omega_old, double alpha, double gamma, double theta, double P, double somma){
  return -alpha*omega_old - gamma*theta + P - somma;
}

int main(){

  int steps = 2500;
  
  double fR = 50;
  double wR = 2*M_PI*fR;
  double h = 0.01;
  double *theta = (double*) malloc(nodes * (steps+1) * sizeof(double));
  double *omega = (double*) malloc(nodes * (steps+1) * sizeof(double)); 
  double *gamma = (double*) malloc(nodes * sizeof(double));
  double *alpha = (double*) malloc(nodes * sizeof(double));
  double *somma = (double*) malloc(nodes * sizeof(double));
  int *ADJ = (int*) malloc(nodes * nodes * sizeof(int));
  double *K = (double*) malloc(nodes * nodes * sizeof(double));
  double *P = (double*) malloc(nodes * sizeof(double));

  
  // variables initialization
  
  for (int i=0; i<nodes; i++){
    gamma[i] = 0.1;
    alpha[i] = 1;
    theta[i] = 0;
    omega[i] = 0;
    somma[i] = 0;
    for (int j=0; j<nodes; j++){
      K[i*nodes+j] = 1.03;
    }
  }
  P[0] = -1; //consumers
  P[4] = -1;
  P[7] = -1;
  P[5] = -1;
  
  P[3] = 1; //generators
  P[2] = 1;
  P[6] = 1;
  P[1] = 1;

  // adjacency matrix
  
  ADJ[0+1*nodes] = 1;
  ADJ[0+6*nodes] = 1;
  ADJ[0+5*nodes] = 1;
  ADJ[1+0*nodes] = 1;
  ADJ[1+2*nodes] = 1;
  ADJ[2+1*nodes] = 1;
  ADJ[2+6*nodes] = 1;
  ADJ[2+3*nodes] = 1;
  ADJ[3+2*nodes] = 1;
  ADJ[3+4*nodes] = 1;
  ADJ[3+7*nodes] = 1;
  ADJ[4+3*nodes] = 1;
  ADJ[4+5*nodes] = 1;
  ADJ[5+4*nodes] = 1;
  ADJ[5+7*nodes] = 1;
  ADJ[5+0*nodes] = 1;
  ADJ[6+0*nodes] = 1;
  ADJ[6+2*nodes] = 1;
  ADJ[7+3*nodes] = 1;
  ADJ[7+5*nodes] = 1;

  //adding line 2 - 4 to unstabilize network
  
  //ADJ[1+3*nodes] = 1;
  //ADJ[3+1*nodes] = 1;

  
  // integration using euler method
  
  FILE  *theta_doc;
  theta_doc = fopen("theta", "w");
  
  double k1, k2, k3, k4;

  for (int t=1; t<=steps; t++){
    fprintf(theta_doc, "%16.8f", t*h);
    for (int i=0; i<nodes; i++){  

      somma[i] = 0;
      for (int j=0; j<nodes; j++){
	if (ADJ[i+nodes*j]!=0){
	  somma[i] += 1.03 * sin(theta[i+(t-1)*nodes] - theta[j + (t-1)*nodes]);
	}
      }

      k1 = omega_funct(omega[i+(t-1)*nodes], alpha[i], gamma[i], theta[i+(t-1)*nodes], P[i], somma[i]);
      k2 = omega_funct(omega[i+(t-1)*nodes]+0.5*k1*h, alpha[i], gamma[i], theta[i+(t-1)*nodes], P[i], somma[i]);
      k3 = omega_funct(omega[i+(t-1)*nodes]+0.5*k2*h, alpha[i], gamma[i], theta[i+(t-1)*nodes], P[i], somma[i]);
      k4 = omega_funct(omega[i+(t-1)*nodes]+k3*h, alpha[i], gamma[i], theta[i+(t-1)*nodes], P[i], somma[i]);
      
      omega[i+t*nodes] = omega[i+(t-1)*nodes] + h/6*(k1+2*k2+2*k3+k4);
      theta[i+t*nodes] = theta[i+(t-1)*nodes] + h*(omega[i+(t-1)*nodes]);

      fprintf(theta_doc, "%16.8e", theta[i+t*nodes]/M_PI); 
    }
    fprintf(theta_doc, "\n");

  }
  
  fclose(theta_doc);
  free(theta);  free(gamma);  free(P);  free(alpha);  free(omega); free(somma);  free(K); free(ADJ);
  
  return 0;
}
