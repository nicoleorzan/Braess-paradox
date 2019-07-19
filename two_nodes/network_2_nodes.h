#ifndef NETWORK
#define NETWORK

#define nodes 2
#define connections 2
#define alpha 0.1
#define Gamma 0.1

const double P[nodes] = {1, -1.2};
const int AI[nodes+1] = {0, 1, 2};
const int AV[connections] = {1, 0};
double weights[connections] = {1.5, 1.5};

double delta = 1;
double Pmax = 1.0;

#endif
