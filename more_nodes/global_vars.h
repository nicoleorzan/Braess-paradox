#ifndef GLOBAL_VARIABLES
#define GLOBAL_VARIABLES

#define TCPU_TIME (clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &ts ), (double)ts.tv_sec +	\
		   (double)ts.tv_nsec * 1e-9)
#define nodes 8
#define connections 20
#define alpha 1
#define Gamma 0.1

extern const double P[8];

const void file_reader(int* AI, int* AV, double* weights){

  FILE *AII, *AVV, *W;
  
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
}


#endif
