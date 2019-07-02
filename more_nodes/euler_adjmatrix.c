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

int main(){

  //int nodes = 8;
  int steps = 2000;
  
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
    gamma[i] = 0;
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
  
  ADJ[1+3*nodes] = 1;
  ADJ[3+1*nodes] = 1;

  print_mat(ADJ, nodes, nodes);

  
  // integration using euler method
  
  FILE  *theta_doc;

  theta_doc = fopen("theta", "w");


  for (int t=1; t<=steps; t++){
    
    fprintf(theta_doc, "%16.8f", t*h);
    
    for (int i=0; i<nodes; i++){  

      somma[i] = 0;
      for (int j=0; j<nodes; j++){
	if (ADJ[i+nodes*j]!=0){
	somma[i] += 1.03 * sin(theta[i+(t-1)*nodes] - theta[j + (t-1)*nodes]); //weights[j]
	}
      }
      
      omega[i+t*nodes] = omega[i+(t-1)*nodes] + h*(-alpha[i]*omega[i+(t-1)*nodes]- gamma[i]*theta[i+(t-1)*nodes] + P[i] - somma[i]);
      theta[i+t*nodes] = theta[i+(t-1)*nodes] + h*(omega[i+(t-1)*nodes]);

      fprintf(theta_doc, "%16.8e", theta[i+t*nodes]/M_PI); 
    }
    fprintf(theta_doc, "\n");

  }
  
  fclose(theta_doc);
  free(theta);  free(gamma);  free(P);  free(alpha);  free(omega); free(somma);  free(K); free(ADJ);
  
  return 0;
}
