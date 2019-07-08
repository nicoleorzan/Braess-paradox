#ifndef GLOBAL_VARIABLES
#define GLOBAL_VARIABLES

#include <string.h>
#include <sys/resource.h>
#include <sys/times.h>
#include <time.h>


#define TCPU_TIME (clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &ts ), (double)ts.tv_sec +	\
		   (double)ts.tv_nsec * 1e-9)

#define nodes 8
#define connections 20

extern const double P[nodes];
extern const int AI[nodes+1];
extern const int AV[connections];
extern double weights[connections];
extern double delta;




#endif
