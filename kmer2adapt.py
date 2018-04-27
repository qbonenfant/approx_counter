# coding=UTF-8
import networkx as nx

# naive version: overlap of most common kmer

G = nx.DiGraph()


def haveOverlap(main, potential, limit):
    overlapped = ""
    lm = len(main)
    lp = len(potential)
    maxOV = min(lm, lp)
    if(maxOV >= limit):
        for i in range(maxOV - 1, limit - 1, -1):
            # checkin if pot is prefix
            if main[lm - i:] == potential[:i]:
                overlapped = main + potential[i:]

            elif potential[lp - i:] == main[:i]:
                overlapped = potential + main[i:]
    return(overlapped)

def graphOverlaps(main, potential, limit, nxGraph):

    lm = len(main)
    lp = len(potential)
    maxOV = min(lm, lp)
    if(maxOV >= limit):
        for i in range(maxOV - 1, limit - 1, -1):
            # checkin if pot is prefix
            if main[lm - i:] == potential[:i]:
                nxGraph.add_node(main)
                nxGraph.add_node(potential)
                nxGraph.add_edge(main, potential, weigth=i)

            elif potential[lp - i:] == main[:i]:
                nxGraph.add_node(main)
                nxGraph.add_node(potential)
                nxGraph.add_edge(potential, main, weigth=i)


kmers = []
with open("comptageFullKmer.txt", "r") as f:
    f.readline()
    for line in f:
        nb, kmer = line.replace("\n", "").split(" ; ")
        kmers.append(kmer)

cutoff = 5

# for kmer1 in kmers:
#     for kmer2 in kmers:
#         graphOverlaps(kmer1, kmer2, cutoff, G)
# nx.write_graphml(G, "exportedGraph.graphml")

longestPattern = []
loop = 0
while(len(kmers)>0 or loop >= 1000):
    newSeq = kmers[0]
    kmers.pop(0)
    found = True
    while found:
        found = False
        currentElem = ""
        for elem in kmers:
            ovl = haveOverlap(newSeq, elem, cutoff)
            if(ovl != ''):
                currentElem = elem
                newSeq = ovl
                found = True
                break
        if(found):
            kmers.pop(kmers.index(currentElem))
    longestPattern.append(newSeq)
    loop += 1
results = list(reversed(sorted(longestPattern, key = lambda x: len(x))))
for elem in results[:10]:
    print(elem)
