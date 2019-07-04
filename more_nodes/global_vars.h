#ifndef GLOBAL_VARIABLES
#define GLOBAL_VARIABLES

#define TCPU_TIME (clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &ts ), (double)ts.tv_sec +	\
		   (double)ts.tv_nsec * 1e-9)
#define nodes 8
#define connections 20
#define alpha 1
#define Gamma 0.1

extern const double P[8];

const void file_reader(int* AI, int* AV, double* weights, char* ai_name, char* av_name, char* weights_name){

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
  fclose(AII);
  fclose(W);
}


#endif
