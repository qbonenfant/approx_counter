# coding=utf-8

import sys
import networkx as nx

count_file = sys.argv[1]

kmer_count = {}
kmer_list = []
adapters = {}


def get_weight(g,path):
    total = 0
    for node in path:
        total += g.nodes[node]["weight"]
    return(total)

def haveOverlap(seq1, seq2):
    minOverlap = min(len(seq1), len(seq2)) - 1
    if(seq1[-minOverlap:] == seq2[:minOverlap]):
        return(seq1 + seq2[minOverlap:])
    else:
        return("")

def greedy_assembl(g,kmer_list):

    km = kmer_list[0]    
    g.node[km]["path"] = "greedy" if g.node[km]["path"] == "" else "both"
    g.node[km]["position"] = "FIRST"

    ov = km
    found = True
    used = [km]
    while(found):
        found = False
        for km2 in kmer_list:
            if(km2 not in used):
                direct = haveOverlap(ov, km2)
                reverse = haveOverlap(km2, ov)
                if(direct != "" and reverse == ""):
                    ov = direct
                    found = True
                    used.append(km2)
                    g.node[km2]["path"] = "greedy" if g.node[km2]["path"] == "" else "both"
                    g.node[km2]["position"] = "RIGHT"
                    break

                elif(reverse != "" and direct == ""):   
                    ov = reverse
                    found = True
                    used.append(km2)
                    g.node[km2]["path"] = "greedy" if g.node[km2]["path"] == "" else "both"
                    g.node[km2]["position"] = "LEFT"
                    break
    return(ov)


g = nx.DiGraph()

# building graph
with open(count_file, 'r') as f:
    for line in f:
        km, nb = line.rstrip("\n").split("\t")
        kmer_count[km] = int(nb)
        kmer_list.append(km)
        g.add_node(km, weight = int(nb), path = "" )
    

K = len(kmer_list[0])

# searching overlaps
for km in kmer_list:
    for km2 in kmer_list:
        # since i do all the possible combinations
        # There is no need to test both orientation
        direct = haveOverlap(km, km2)
        if(direct != ""):            
            g.add_edge(km, km2)
        
# removing singletons
g.remove_nodes_from( [node for node in g.nodes if  len(list(nx.all_neighbors(g,node))) == 0] )
    
lg_path = nx.dag_longest_path(g)

for node in lg_path:
    nd = g.nodes[node]
    nd['path'] = "main"

print( lg_path[0][:-1] + "".join( el[-1] for el in lg_path ) )
print(greedy_assembl(g,kmer_list))


#heavy path
sources = [n for n in g.nodes if not list(g.predecessors(n))  ]
targets = [n for n in g.nodes if not list(g.successors(n))    ]
# print(sources)
# print(targets)

hv_path = []
w = 0
for source in sources:
    for target in targets:
        try:
            heaviest_path = max((path for path in nx.all_simple_paths(g, source, target)),
                        key=lambda path: get_weight(g,path))
            current_w = get_weight(g,heaviest_path)
            if(current_w > w):
                hv_path = heaviest_path
                w = current_w
        except ValueError:
            pass

for node in hv_path:
    g.nodes[node]["path"] = "heavy" if g.nodes[node]["path"] == "" else "common"

print( hv_path[0][:-1] + "".join( el[-1] for el in hv_path ) )
nx.write_graphml(g, count_file + '.graphml')