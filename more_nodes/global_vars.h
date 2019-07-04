#ifndef GLOBAL_VARIABLES
#define GLOBAL_VARIABLES

#include <string.h>


#define TCPU_TIME (clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &ts ), (double)ts.tv_sec +	\
		   (double)ts.tv_nsec * 1e-9)

#define nodes 8
#define connections 20
#define alpha 1
#define Gamma 0.1
#define h 0.01
#define hh h/0.5
#define h6 h/6

extern const double P[8];

const void file_reader(int* AI, int* AV, double* weights, char const * ai_name, char const * av_name, char const * weights_name){

  FILE *AII, *AVV, *W;
  
  AVV = fopen(av_name, "r");
  AII = fopen(ai_name, "r");
  W = fopen(weights_name, "r");

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
  memset(AVV, 0, sizeof(*AVV));
  free(AVV);
  fclose(AII);
  memset(AII, 0, sizeof(*AII));
  free(AII);
  fclose(W);
  memset(W, 0, sizeof(*W));
  free(W);
}


#endif
