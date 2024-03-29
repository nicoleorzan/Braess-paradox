#ifndef NETWORK
#define NETWORK

#define nodes 8
#define connections 20 // 22
#define alpha 1.0
#define Gamma 0


//const double Gamma[nodes] = {0, 0, 0, 0, 0, 0, 0, 0};
const double P[nodes] = {-1, 1, 1, 1, -1, -1, 1, -1};
const int AI[nodes+1] = {0, 3, 5, 8, 11, 13, 16, 18, 20};
//const int AI[nodes+1] = {0, 3, 6, 9, 13, 15, 18, 20, 22}; //22 connections
const int AV[connections] = {5, 1, 6, 0, 2, 3, 1, 6, 2, 4, 7, 3, 5, 4, 7, 0, 2, 0, 3, 5};
//const int AV[connections] = {5, 1, 6, 0, 2, 3, 3, 1, 6, 1, 2, 4, 7, 3, 5, 4, 7, 0, 2, 0, 3, 5}; //22 connections
double weights[connections] = {1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03};
//double weights[connections] = {1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03}; //22 connections
double Pmax = 0.01;
double delta = 1;

#endif
