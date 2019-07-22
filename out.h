#ifndef NETWORK 
#define NETWORK
 
#define nodes  7
#define connections  14
#define alpha  1.0
const double P[nodes] =  {1, -1, -1, -1, 1, 1, -1} ;
const int AI[nodes+1] = {0, 3, 5, 6, 8, 10, 12, 14} ;
const int AV[connections] = {4, 6, 7, 4, 7, 5, 0, 1, 3, 6, 0, 5, 0, 1} ;
double weights[connections] = {1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03} ;
double Pmax =  0.1 ;
double delta =  1.0 ;

#endif
