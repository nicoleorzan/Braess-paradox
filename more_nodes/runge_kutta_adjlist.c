#include <stdlib.h>
#include <stdio.h>
#include <math.h>

# define nodes 8

void print_mat(double *matrix, int N, int M){
  for (int i=0; i<N; i++){
      for (int j=0; j<M; j++){
	fprintf(stdout, "%f ", matrix[i*M+j]);
      }
      fprintf(stdout, "\n");
  }  
}

double omega_funct(double omega_old, double alpha, double gamma, double theta, double P, double somma){
  return -alpha*omega_old - gamma*theta + P - somma;
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
  double *K = (double*) malloc(nodes * nodes * sizeof(double));
  double *P = (double*) malloc(nodes * sizeof(double));
  double *num_interactions = (double*) malloc(nodes * sizeof(double));

  
  // reading from file
 
  FILE *AII, *AVV, *interact, *W;
  int interactions[nodes], AI[nodes+1],  AV[20], weights[20];
  interact = fopen("files_to_read/interactions.txt", "r");
  AVV = fopen("files_to_read/AV.txt", "r");
  AII = fopen("files_to_read/AI.txt", "r");
  W = fopen("files_to_read/weights.txt", "r");

  if (AII == NULL || AVV == NULL || interact == NULL){
    printf("Error Reading File\n");
    exit (0);
  }
  for (int i = 0; i < nodes; i++){
    fscanf(interact, "%d,", &interactions[i] );
  }
  for (int i = 0; i < nodes+1; i++){
    fscanf(AII, "%d,", &AI[i] );
  }
  for (int i = 0; i < 20; i++){
    fscanf(AVV, "%d,", &AV[i] );
    fscanf(W, "%d,", &weights[i] );
  }
  fclose(interact);
  fclose(AVV);
  fclose(AII);
  fclose(W);

  
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

  
  // integration using runge-kutta method of 4th order
  
  FILE  *theta_doc;
  
  theta_doc = fopen("theta", "w");
  
  double k1, k2, k3, k4;
  
  for (int t=1; t<=steps; t++){
    fprintf(theta_doc, "%16.8f", t*h);
    
    for (int i=0; i<nodes; i++){  

      somma[i] = 0;
      for (int j = AI[i]+1; j<=AI[i+1]; j++){
	 somma[i] += 1.03 * sin(theta[i+(t-1)*nodes] - theta[(AV[j]-1) + (t-1)*nodes]);
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
  free(theta);  free(gamma);  free(P);  free(alpha);  free(omega); free(somma);  free(K);
  
  return 0;
}
