# coding=utf-8


import networkx as nx
from editdistance import eval as lv
import sys

fileName = sys.argv[1]
kSize = int(sys.argv[2])
kmerCount = {}

print("OPENING FILE")
with open(fileName, 'r') as f:
    loop = True
    ident = ""
    seq = ""
    while(loop):
        try:
            line = next(f)
        except StopIteration:
            loop = False
        else:
            if(line[0] != ">"):
                seq += line[:-1]
            elif(seq != ""):
                for i in range(len(seq) - kSize + 1):
                    km = seq[i:i + kSize]
                    if(km not in kmerCount.keys()):
                        kmerCount[km] = 1
                    else:
                        kmerCount[km] += 1

print("GENERATED KMER COUNT")

# print(list(kmerCount.items()))

G = nx.Graph()
kmers = list([key for key, val in kmerCount.items() if val > 40])
for kmer1 in kmers:
    for kmer2 in kmers:
        if(lv(kmer1, kmer2) < 2):
            G.add_node(kmer1, count=kmerCount[kmer1])
            G.add_node(kmer2, count=kmerCount[kmer2])
            G.add_edge(kmer1, kmer2)

nx.write_graphml(G, "neighbourgGraph.graphml")
print("GENERATED NEIGHBORHOOD GRAPH")
