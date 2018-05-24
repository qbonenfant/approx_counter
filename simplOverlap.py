# coding=utf-8


import networkx as nx
import sys

fileName = sys.argv[1] if(len(sys.argv) >= 2) else "nope"
minExact = int(sys.argv[2]) if(len(sys.argv) >= 3) else 10
kmerCount = {}


def haveOverlap(seq1, seq2):
    minOverlap = min(len(seq1), len(seq2)) - 1
    return(seq1[-minOverlap:] == seq2[:minOverlap])


print("OPENING FILE")
with open(fileName, 'r') as f:
    loop = True
    ident = ""
    seq = ""
    next(f)  # skipping first line
    while(loop):
        try:
            line = next(f)
        except StopIteration:
            loop = False
        else:
            
            noErr, oneErr, twoErr, km = line.rstrip("\n").split("\t")
            kmerCount[km] = (noErr, oneErr, twoErr)

print("GENERATED KMER COUNT")

# print(list(kmerCount.items()))

G = nx.DiGraph()
kmers = list([key for key, val in kmerCount.items()])
for kmer1 in kmers:
    for kmer2 in kmers:
        count1 = kmerCount[kmer1]
        count2 = kmerCount[kmer2]
        G.add_node(kmer1, exact=count1[0],
                   oneErr=count1[1], twoErr=count1[2])
        G.add_node(kmer2, exact=count2[0],
                   oneErr=count2[1], twoErr=count2[2])
        if(haveOverlap(kmer1, kmer2)):
            G.add_edge(kmer1, kmer2)
        elif(haveOverlap(kmer2, kmer1)):
            G.add_edge(kmer2, kmer1)

nx.write_graphml(G, "neighbourgGraph.graphml")

print("GENERATED NEIGHBORHOOD GRAPH")
