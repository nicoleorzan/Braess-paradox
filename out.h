#ifndef NETWORK 
#define NETWORK
 
#define nodes  8
#define connections  26
#define alpha  1.0
#define Gamma  0

const double P[nodes] =  {1, -1, -1, 1, -1, 1, 1, -1} ;
const int AI[nodes+1] = {0, 2, 6, 8, 13, 15, 18, 21, 26} ;
const int AV[connections] = {5, 6, 2, 3, 5, 7, 1, 3, 1, 2, 4, 6, 7, 3, 7, 0, 1, 7, 0, 3, 7, 1, 3, 4, 5, 6} ;
double weights[connections] = {10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10} ;
//double Pmax =  0.01 ;
//double delta =  0.1/0.01 ;

double Pmax[nodes] = {0.01, 0.1, 0.1, 0.01, 0.1, 0.01, 0.01, 0.1}; // al contrario
double delta[nodes] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

#endif
