from networkx import nx
import numpy as np
#import sys
import random

class Network_to_file:
    
    def __init__(self, G, alpha = 1.0, Gamma = 0., delta = 1.0, Pmax = 0.1, Pvals = None):
        self.G = G
        self.alpha = alpha
        self.Gamma = Gamma
        self.delta = delta
        self.Pmax = Pmax
        self.Pvals = Pvals
        self.createArrays()
        self.create_Pvals()
        
    def createArrays(self):
        self.nodes = self.G.number_of_nodes()
        self.AV = [n for i in sorted(self.G.nodes) for n in self.G.neighbors(i)]
        self.connections = len(self.AV)
        
        self.AI = [0]
        tmp = 0
        for i in sorted(self.G.nodes):
            tmp += len([n for n in self.G.neighbors(i)])
            self.AI.append(tmp)
        
        #self.weights = [ G[i][n][0]['weight'] for i in sorted(G.nodes) for n in G.neighbors(i) ] if multiedge
        self.weights = [ self.G[n][i]['weight'] for i in sorted(self.G.nodes) for n in self.G.neighbors(i) ]
        assert(len(self.weights) == self.connections)
        
    def create_Pvals(self):
        if self.Pvals == None:
            a = random.sample(self.G.nodes, int(len(self.G)/2))
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
            print("#define Gamma ", self.Gamma, file = f)
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
        
