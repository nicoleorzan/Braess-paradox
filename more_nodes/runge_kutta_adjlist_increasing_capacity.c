#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global_vars.h"

# define nodes 8

int connections = 20;
double alpha = 1;
double Gamma = 0.1;
double P[8] = {-1, 1, 1, 1, -1, -1, 1, -1};


void print_mat(double *matrix, int N, int M){
  for (int i=0; i<N; i++){
      for (int j=0; j<M; j++){
	fprintf(stdout, "%f ", matrix[i*M+j]);
      }
      fprintf(stdout, "\n");
  }  
}

double omega_funct(double omega_old, double theta, double somma, int i){
  return -alpha*omega_old - Gamma*theta + P[i] - somma;
}

double theta_funct(double omega){
  return omega;
}


void derivs(double *omega, double *theta, double somma, double *dydt){
  for (int i=0; i<nodes; i++){
    dydt[i] = -alpha*omega[i] - Gamma*theta[i] + P[i] - somma;
  }
}


int main(){

  int steps = 2500;

  double fR = 50;
  double wR = 2*M_PI*fR;
  double h = 0.01;
  double *theta = (double*) malloc(nodes * sizeof(double));
  double *omega = (double*) malloc(nodes * sizeof(double));

  for (int i=0; i<nodes; i++){
    theta[i] = 0;
    omega[i] = 0;
  }
  
  // reading from file
 
  FILE *AII, *AVV, *W;
  int AI[nodes+1],  AV[connections];
  double weights[connections];
  AVV = fopen("files_to_read/AV.txt", "r");
  //AVV = fopen("files_to_read/AV_adding_line.txt", "r");
  AII = fopen("files_to_read/AI.txt", "r");
  //AII = fopen("files_to_read/AI_adding_line.txt", "r");
  W = fopen("files_to_read/weights.txt", "r");
  //W = fopen("files_to_read/weights_adding_line.txt", "r");

  if (AII == NULL || AVV == NULL || W == NULL){
    printf("Error Reading File\n");
    exit(0);
  }
  for (int i = 0; i < connections; i++){
    fscanf(W, "%lf,", &weights[i] );
    fscanf(AVV, "%d,", &AV[i] );
  }
  for (int i = 0; i < nodes+1; i++){
    fscanf(AII, "%d,", &AI[i] );
  }
 
  fclose(AVV);
  fclose(AII);
  fclose(W);
  
  
  // integration using runge-kutta method of 4th order

  FILE *theta_doc;
  double hh = h/0.5, h6 = h/6;
  double sum = 0;
    
  theta_doc = fopen("theta", "w");
  
  double k1[2], k2[2], k3[2], k4[2];

  //double yt[8], dym[8], dyt[8], dydx[8];

  double cap = weights[0];
  int iter = 0;
  double max_capacity = 1.73;
  double deltaK = 0.1;
  int printing_step = 2400;
  
  while (cap < max_capacity){

    fprintf(theta_doc, "%16.8f", deltaK*iter);
    
    for (int i=0; i<nodes; i++){
      theta[i] = 0;
      omega[i] = 0;
    }
   
    for (int t=1; t<=steps; t++){
      //fprintf(stdout, "%d\n", t);
      //fprintf(theta_doc, "%16.8f", t*h);

      for (int i=0; i<nodes; i++){

	sum = 0;
	for (int j = AI[i]; j < AI[i+1]; j++){
	  sum += weights[j] * sin( theta[i] - theta[(AV[j])] ); 
	}
      
	k1[0] = omega_funct(omega[i], theta[i], sum, i);
	k1[1] = theta_funct(omega[i]);
	k2[0] = omega_funct(omega[i]+k1[0]*hh, theta[i], sum, i);
	k2[1] = theta_funct(omega[i]+k1[1]*hh);
	k3[0] = omega_funct(omega[i]+k2[0]*hh,  theta[i], sum, i);
	k3[1] = theta_funct(omega[i]+k2[1]*hh);
	k4[0] = omega_funct(omega[i]+k3[0]*h,  theta[i], sum, i);
	k4[1] = theta_funct(omega[i]+k3[1]*h);

	omega[i] = omega[i] + h6*(k1[0]+2*k2[0]+2*k3[0]+k4[0]);
	theta[i] = theta[i] + h6*(k1[1]+2*k2[1]+2*k3[1]+k4[1]);
      
      
	
	/*if (t % printing_step == 0){
	  fprintf(theta_doc, "%16.8e", theta[i]/M_PI);
	  }
	*/
      }
      if (t == printing_step){
	fprintf(theta_doc, "%16.8e", (theta[3]-theta[4])/M_PI); //diff nodes 4-5
	fprintf(theta_doc, "%16.8e", (theta[3]-theta[7])/M_PI); //diff nodes 4-8
	fprintf(theta_doc, "%16.8e", (theta[0]-theta[5])/M_PI); //diff nodes 1-6
	fprintf(theta_doc, "%16.8e", (theta[2]-theta[3])/M_PI); //diff nodes 3-4
	fprintf(theta_doc, "%16.8e", (theta[2]-theta[1])/M_PI); //diff nodes 3-2
	fprintf(theta_doc, "%16.8e", (theta[2]-theta[6])/M_PI); //diff nodes 3-7
	fprintf(theta_doc, "%16.8e", (theta[5]-theta[4])/M_PI); //diff nodes 6-5
	fprintf(theta_doc, "%16.8e", (theta[5]-theta[7])/M_PI); //diff nodes 6-8
	fprintf(theta_doc, "%16.8e", (theta[1]-theta[0])/M_PI); //diff nodes 2-1
	fprintf(theta_doc, "%16.8e", (theta[6]-theta[0])/M_PI); //diff nodes 7-1
      }
      
      //fprintf(theta_doc, "\n");
    }

    
    
    weights[5] += deltaK;
    weights[8] += deltaK;
    iter+=1;
    cap += deltaK;

    for (int k=0; k<connections; k++){
      fprintf(stdout, "%f ", weights[k]);
    }

    fprintf(theta_doc, "\n");
    fprintf(stdout, "\n");
    
  }
   

  fclose(theta_doc);
  free(theta);  free(omega);
  
  return 0;
}



 /*

derivs(omega, theta, sum, dydx);

for (int i=0; i<nodes; i++){
      yt[i] = omega[i] + hh*dydx[i];
      }
    derivs(yt, theta, somma, dyt);
  
    for (int i=0; i<nodes; i++){
      yt[i] = omega[i] + hh*dyt[i];
    }
    derivs(yt, theta, somma, dym);

    for (int i=0; i<nodes; i++) {
      yt[i] = omega[i] + h*dym[i];
      dym[i] += dyt[i];
    }
    derivs(yt, theta, somma, dyt);

    for (int i=0; i<nodes; i++) {
      omega[i] = omega[i] + h6*(dydx[i]+dyt[i]+2.0*dym[i]);
      theta[i] = theta[i] + h*(omega[i]);
      if (t % printing_step == 0){
	fprintf(theta_doc, "%16.8e", theta[i]/M_PI);
      }
    }
    */
