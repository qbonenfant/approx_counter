# coding=utf-8

import sys
import networkx as nx


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

def build_graph(count_file):
    """Build the adapter from the kmer count file.
       The way it is done is by building a directed weighted graph
       and searching for the heaviest path.
       I also added the greedy adapter output
    """

    kmer_count = {}
    kmer_list = []

    g = nx.DiGraph()

    # building graph
    with open(count_file, 'r') as f:
        for line in f:
            km, nb = line.rstrip("\n").split("\t")
            kmer_count[km] = int(nb)
            kmer_list.append(km)
            g.add_node(km, weight = int(nb))#, path = "" )
        
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

    # Returning only the biggest connected component
    return( g )



def greedy_assembl(g):
    """Greedy assembly method to compute the adapter
    TODO: modify it so it use directly the graph...
    @param the overlap graph of kmers
    @return the longest debruijn sequence starting by the first kmer
    """
    kmer_dict = g.nodes(data=True)
    kmer_list = list( dict(kmer_dict).keys() )
    kmer_list.sort( key= lambda x: kmer_dict[x]['weight'])
    kmer_list.reverse()

    km = kmer_list[0]
    ov = km
    found = True
    used = [km]
    #annotating graph
    g.node[km]["path"] = (g.node[km]["path"] + ",greedy").lstrip(",")

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
                    g.node[km2]["path"] = (g.node[km2]["path"] + ",greedy").lstrip(",")
                    break

                elif(reverse != "" and direct == ""):   
                    ov = reverse
                    found = True
                    used.append(km2)
                    g.node[km2]["path"] = (g.node[km2]["path"] + ",greedy").lstrip(",")

                    break
    return(ov)

def dag_heaviest_path(G):
    """Returns the heaviest path in a DAG

    Parameters
    ----------
    G : NetworkX DiGraph
        Graph

    Returns
    -------
    path : list
        Heaviest path

    Comment
    -------
    This is a modified version of the dag_longest_path
    using node weight as distance.
    """
    dist = {}  # stores [node, distance] pair
    for node in nx.topological_sort(G):
        # pairs of dist,node for all incoming edges
        pairs = [(dist[v][0] + G.node[v]["weight"], v) for v in G.pred[node]]
        if pairs:
            dist[node] = max(pairs)
        else:
            dist[node] = (0, node)
    node, (length, _) = max(dist.items(), key=lambda x: x[1])
    path = []
    while length > 0:
        path.append(node)
        length, node = dist[node]
    return list(reversed(path))


def heavy_path(g):
    """ Searching the truly heaviest path between all source and target nodes.
        Even if longer path tend to be heavier, very heavy short path can 
        also be selected.
    """
    
    hv_path = dag_heaviest_path(g.subgraph(max(nx.weakly_connected_components(g), key= lambda x: get_weight(g,x))))

    for n in hv_path:
        g.node[n]["path"] = (g.node[n]["path"] + ",heavy").lstrip(",")

    return( hv_path[0][:-1] + "".join( el[-1] for el in hv_path ) )


def longest_path(g):
    lg_path = nx.dag_longest_path(g.subgraph(max(nx.weakly_connected_components(g), key= lambda x: get_weight(g,x))))

    for n in lg_path:
        g.node[n]["path"] = (g.node[n]["path"] + ",long").lstrip(",")

    return(lg_path[0][:-1] + "".join( el[-1] for el in lg_path ))


count_file = sys.argv[1]
  
# Building graph
g = build_graph( count_file )

# preping for anotation
nx.set_node_attributes(g, "", "path" )

# greedy adapters
greedy_adapter = greedy_assembl(g)

# Longest path, catching exception to avoid loops
try:
    long_adapter = longest_path(g)

except:
    print("Could not compute adaper using longest path  method", file = sys.stderr)
    print("The resulting graph probably contains a loop.", file = sys.stderr)
    long_adapter  = ""

# heaviest path adapter  
try:
    heavy_adapter = heavy_path(g)
except:
    print("Could not compute adaper using heaviest path  method", file = sys.stderr)
    print("The resulting graph probably contains a loop.", file = sys.stderr)
    heavy_adapter = ""


#Exporting graph
path = "./adapter_graph.graphml"
nx.write_graphml(g,path)

# printing adapters
print("Greedy")
print(greedy_adapter)
print("Long")
print(long_adapter)
print("Heavy")
print(heavy_adapter)