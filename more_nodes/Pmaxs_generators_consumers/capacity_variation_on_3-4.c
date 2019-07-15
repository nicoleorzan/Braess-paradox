#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include "../include/time_computing.h"
#include "../include/network_generators_consumers.h"
#include "../include/stability_check.h"
#include "../include/runge_kutta_generators_consumers.h"

#define steps 1000000
#define additive_steps 1000
#define internal_steps 10

#define max_capacity 1.73
#define deltaK 0.01


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
  fprintf(f, "%16.8e", (y[6]-y[0])/M_PI); //diff nodes 7-1
  fprintf(f, "\n");
}

void printer_control(double *y, FILE * ff){
  for (int i=0; i<nodes; i++){
    fprintf(ff, "%18.8e ", -Pmax[i]*tanh(delta[i]*y[i]) );
  }
  fprintf(ff, "\n");
}


int main(){

  double tstart, tstop, ctime=0;
  struct timespec ts;
  
  double *y = (double*) malloc(2 * nodes * sizeof(double));
  for (int i=0; i<2*nodes; i++){
    y[i] = 0;
  }
  int iter = 0;
  bool unstable = 0;
  double cap = weights[5];

  FILE* capacity_doc, *control;
  capacity_doc = fopen("tmp_bis", "w");
  control = fopen("control_bis", "w");
  fprintf(capacity_doc, "     deltacapacity     delta45/pi       delta48/pi       delta16/pi       delta34/pi      delta32/pi      delta37/pi     delta65/pi\n     delta68/pi     delta21/pi      delta71/pi\n");
  fprintf(control, "      deltaCapacity      control[0]     control[1]     control[2]     control[3]     control[4]     control[5]     control[6]     control[7]\n");
  
  tstart = TCPU_TIME;

  for (int i=0; i<2*nodes; i++){
    y[i] = 0;
  }

  while (cap < max_capacity){
    fprintf(stdout, "%f\n", cap);
    fprintf(capacity_doc, "%16.8f", deltaK*iter);
    fprintf(control, "%16.8f", deltaK*iter);

    for (int t=1; t<=steps; t+=internal_steps){  
      runge_kutta(y, internal_steps);
    }
    stability_check(runge_kutta, y, additive_steps, &unstable);

    if (unstable == 1) break;

    printer(y, capacity_doc);
    printer_control(y, control);

    weights[5] += deltaK;
    weights[8] += deltaK;
    iter += 1;
    cap += deltaK;

  }

  ctime += TCPU_TIME - tstart;
  printf("%g sec \n", ctime);

  fclose(capacity_doc);
  memset(capacity_doc, 0, sizeof(*capacity_doc));
  free(capacity_doc);
  fclose(control);
  memset(control, 0, sizeof(*control));
  free(control);
  memset(y, 0, sizeof(*y));
  free(y);
    
  return 0;
}



