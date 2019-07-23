#ifndef NETWORK 
#define NETWORK
 
#define nodes  8
#define connections  24
#define alpha  1.0
#define Gamma  0.0 

const double P[nodes] =  {1, 1, -1, -1, 1, -1, 1, -1} ;
const int AI[nodes+1] = {0, 4, 7, 11, 12, 15, 17, 21, 24} ;
const int AV[connections] = {1, 5, 6, 7, 0, 2, 6, 1, 4, 6, 7, 4, 2, 3, 6, 0, 7, 0, 1, 2, 4, 0, 2, 5} ;
double weights[connections] = {10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10} ;
double Pmax =  0.1 ;
double delta =  1.0 ;

#endif
