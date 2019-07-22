from networkx import nx

class Define_Net():
    
    def __init__(self, n, p):
        self.n = n
        self.p = p
        self.create_network()
        
    def create_network(self):

        self.G = nx.erdos_renyi_graph(n = self.n, p = self.p)
        self.G.remove_nodes_from(list(nx.isolates(self.G)))
        while nx.is_connected(self.G) == False:
            self.G = nx.erdos_renyi_graph(n = self.n, p = self.p)
            self.G.remove_nodes_from(list(nx.isolates(self.G)))
            
        self.G = nx.relabel_nodes(self.G, dict(zip(self.G.nodes(), range(len(self.G)))))
        for source, target in self.G.edges():
            self.G[source][target]['weight'] = 1.03

    def plot_net(self):
        nx.draw(self.G,with_labels = True)
        
    def get_network(self):
        return self.G

