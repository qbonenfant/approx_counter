# coding=utf-8


import sys
import networkx as nx

fileName = sys.argv[1] if(len(sys.argv) >= 2) else "nope"
verbose = True
kmerCount = {}
kmerList = []
adapters = {}


def haveOverlap(seq1, seq2):
    minOverlap = min(len(seq1), len(seq2)) - 1
    if(seq1[-minOverlap:] == seq2[:minOverlap]):
        return(seq1 + seq2[minOverlap:])
    else:
        return("")


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
            kmerList.append(km)
            kmerCount[km] = (noErr, oneErr, twoErr)

K = len(kmerList[0])

g = nx.DiGraph()

for km in kmerList[0:1]:
    ov = km
    found = True
    used = []
    while(found):
        found = False
        for km2 in kmerList:
            if(km2 not in used):
                direct = haveOverlap(ov, km2)
                reverse = haveOverlap(km2, ov)

                if(direct != "" and reverse == ""):
                    ov = direct
                    found = True
                    used.append(km2)
                    break
                elif(reverse != "" and direct == ""):
                    ov = reverse
                    found = True
                    used.append(km2)
                    break

print(ov)
for km in [ov[i:i+K] for i in range(len(ov)-K+1)]:
    print(km)
    for km2 in kmerList:
        
        direct = haveOverlap(km, km2)
        reverse = haveOverlap(km2, km)

        if(direct != ""):
            w = kmerCount[km2]
            g.add_node(km2, weight=sum(
                [int(i) for i in w]), zero=w[0], one=w[1], two=w[2])
            g.add_edge(km, km2, direction="RIGHT")
        if(reverse != "" ):
            w = kmerCount[km2]
            g.add_node(km2, weight=sum(
                [int(i) for i in w]), zero=w[0], one=w[1], two=w[2])
            g.add_edge(km2, km, direction="LEFT")

nx.write_graphml(g, 'outputGraph.graphml')
