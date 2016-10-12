# input for cross-testing of test LAR files to a python environment
from larlib import *
import networkx as nx

lines = tuple(open("test/csv/test1.V", 'r'))
V = [list(eval(line)) for line in lines]
lines = tuple(open("test/csv/test1.EV", 'r'))
EV = [list(eval(line)) for line in lines]

VIEW(STRUCT(MKPOLS((V,[[u-1,v-1] for u,v in EV]))))

G = nx.Graph()
G.add_nodes_from(range(1,len(V)))
G.add_edges_from(EV)

bcc = [sub.edges() for sub in nx.biconnected_component_subgraphs(G, copy=True) 
		if len(sub.edges())>1]
VIEW(MKPOL([V,CAT(bcc),1]))
VIEW(EXPLODE(1.2,1.2,1)(MKPOLS((V,[[u-1,v-1] for u,v in CAT(bcc) ]))))
