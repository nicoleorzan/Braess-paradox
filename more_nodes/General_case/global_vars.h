#ifndef GLOBAL_VARIABLES
#define GLOBAL_VARIABLES

#include <string.h>
#include <sys/resource.h>
#include <sys/times.h>
#include <time.h>


#define TCPU_TIME (clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &ts ), (double)ts.tv_sec +	\
		   (double)ts.tv_nsec * 1e-9)

#endif
