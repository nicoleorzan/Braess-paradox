#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global_vars.h"
#include <sys/resource.h>
#include <sys/times.h>
#include <time.h>

# define nodes 8

#define TCPU_TIME (clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &ts ), (double)ts.tv_sec +	\
		   (double)ts.tv_nsec * 1e-9)

int connections = 20;
double alpha = 1;
double Gamma = 0.1;
double P[8] = {-1, 1, 1, 1, -1, -1, 1, -1};


double omega_funct(double omega_old, double theta, double somma, int i){
  return -alpha*omega_old - Gamma*theta + P[i] - somma;
}

double theta_funct(double omega){
  return omega;
}


int main(){

  int steps = 2500;

  double tstart, tstop, ctime=0;
  struct timespec ts;
  
  int printing_step = 10;
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

  //doubling nodes 3-4 capacity
  //weights[5] = weights[0] * 2;
  //weights[8] = weights[0] * 2;
  

  
  // integration using runge-kutta method of 4th order

  FILE *theta_doc;
  double hh = h/0.5, h6 = h/6;
  double sum = 0;

  theta_doc = fopen("theta", "w");
  
  double k1[2], k2[2], k3[2], k4[2];


  tstart = TCPU_TIME;

  for (int t=1; t<=steps; t++){
    fprintf(theta_doc, "%16.8f", t*h);

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
      
      
          
      if (t % printing_step == 0){
	fprintf(theta_doc, "%16.8e", theta[i]/M_PI);
      }
    }
    fprintf(theta_doc, "\n");
  }

  ctime += TCPU_TIME - tstart;
  printf("%g sec \n", ctime);

  fclose(theta_doc);
  free(theta);  free(omega);
  
  return 0;
}



 