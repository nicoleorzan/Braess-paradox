from networkx import nx
import matplotlib.pyplot as plt
import random 

class Define_Net():
    
    def __init__(self, n, p, w):
        self.n = n
        self.p = p
        self.w = w
        self.create_network()
        self.create_Pvals()
        
    def create_network(self):

        self.G = nx.erdos_renyi_graph(n = self.n, p = self.p)
        self.G.remove_nodes_from(list(nx.isolates(self.G)))
        while nx.is_connected(self.G) == False:
            self.G = nx.erdos_renyi_graph(n = self.n, p = self.p)
            self.G.remove_nodes_from(list(nx.isolates(self.G)))
            
        self.G = nx.relabel_nodes(self.G, dict(zip(self.G.nodes(), range(len(self.G)))))
        for source, target in self.G.edges():
            self.G[source][target]['weight'] = self.w
            
    def create_Pvals(self):
        a = random.sample(self.G.nodes, int(len(self.G)/2))
        b = self.G.nodes - a
        self.Pvals = []
        for i in range(len(self.G)):
            if i in a:
                self.Pvals.append(1)
            else:
                self.Pvals.append(-1)

    def plot_net(self):
        nx.draw(self.G, node_color = self.Pvals, vmin=0, vmax=4, cmap = plt.cm.get_cmap('rainbow'), with_labels = True)
        
    def get_network(self):
        return self.G

