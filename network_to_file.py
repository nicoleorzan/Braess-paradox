#import matplotlib.pyplot as plt
from networkx import nx
import numpy as np
#import sys
#import random
from random import choice

class Network_to_file:
    
    def __init__(self, Graph, alpha = 1.0, delta = 1.0, Pmax = 0.1, Pvals = None):
        self.G = Graph
        self.alpha = alpha
        self.delta = delta
        self.Pmax = Pmax
        self.Pvals = Pvals
        self.createArrays()
        self.Pvals()
        
    def createArrays(self):
        self.nodes = self.G.number_of_nodes()
        self.AV = [n for i in sorted(self.G.nodes) for n in self.G.neighbors(i)]
        self.connections = len(self.AV)
        
        self.AI = [0]
        tmp = 0
        for i in sorted(G.nodes):
            tmp += len([n for n in self.G.neighbors(i)])
            self.AI.append(tmp)
        
        #self.weights = [ G[i][n][0]['weight'] for i in sorted(G.nodes) for n in G.neighbors(i) ] if multiedge
        self.weights = [ self.G[n][i]['weight'] for i in sorted(self.G.nodes) for n in self.G.neighbors(i) ]
        assert(len(self.weights) == self.connections)
        
    def Pvals(self):
        if self.Pvals == None:
            a = random.sample(G.nodes, int(len(self.G)/2))
            b = self.G.nodes - a
            self.Pvals = []
            for i in range(len(self.G)):
                if i in a:
                    self.Pvals.append(1)
                else:
                    self.Pvals.append(-1)
        
    def PrintOnFile(self):

        with open('out.h', mode='wt', encoding='utf-8') as f:
            print("#ifndef NETWORK \n#define NETWORK\n ", file = f)
            print("#define nodes ", self.nodes, file = f)
            print("#define connections ", self.connections, file = f)
            print("#define alpha ", self.alpha, file = f)
            print("const double P[nodes] = ", "{" + ", ".join([str(x) for x in self.Pvals]) + "}", ";", file = f)
            print("const int AI[nodes+1] =", "{" + ", ".join([str(x) for x in self.AI]) + "}",";", file = f)
            print("const int AV[connections] =", "{" + ", ".join([str(x) for x in self.AV]) + "}", ";", file = f)
            print("double weights[connections] =", "{" + ", ".join([str(x) for x in self.weights]) + "}", ";", file = f)
            print("double Pmax = ", self.Pmax, ";", file = f)
            print("double delta = ", self.delta, ";", file = f)
            print("\n#endif", file = f)
        f.close
        
    def get_generators(self):
        return self.Pvals
        
