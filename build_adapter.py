# coding=utf-8

import sys
import networkx as nx
import argparse
import shutil
import os
#This class come from the Porechop tool
# https://github.com/rrwick/Porechop
class MyHelpFormatter(argparse.HelpFormatter):
    """
    This is a custom formatter class for argparse. It allows for some custom formatting,
    in particular for the help texts with multiple options (like bridging mode and verbosity level).
    http://stackoverflow.com/questions/3853722
    """
    def __init__(self, prog):
        terminal_width = shutil.get_terminal_size().columns
        os.environ['COLUMNS'] = str(terminal_width)
        max_help_position = min(max(24, terminal_width // 3), 40)
        super().__init__(prog, max_help_position=max_help_position)

    def _get_help_string(self, action):
        help_text = action.help
        if action.default != argparse.SUPPRESS and 'default' not in help_text.lower() and \
                action.default is not None:
            help_text += ' (default: ' + str(action.default) + ')'
        return help_text


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
    g.nodes[km]["path"] = (g.nodes[km]["path"] + ",greedy").lstrip(",")

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
                    g.nodes[km2]["path"] = (g.nodes[km2]["path"] + ",greedy").lstrip(",")
                    break

                elif(reverse != "" and direct == ""):   
                    ov = reverse
                    found = True
                    used.append(km2)
                    g.nodes[km2]["path"] = (g.nodes[km2]["path"] + ",greedy").lstrip(",")

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
        g.nodes[n]["path"] = (g.nodes[n]["path"] + ",heavy").lstrip(",")

    return( hv_path[0][:-1] + "".join( el[-1] for el in hv_path ) )


def longest_path(g):
    lg_path = nx.dag_longest_path(g.subgraph(max(nx.weakly_connected_components(g), key= lambda x: get_weight(g,x))))

    for n in lg_path:
        g.nodes[n]["path"] = (g.nodes[n]["path"] + ",long").lstrip(",")

    return(lg_path[0][:-1] + "".join( el[-1] for el in lg_path ))



def get_arguments():
    """
    Parse the command line arguments.
    """
    
    parser = argparse.ArgumentParser(description='Adapter builder'
                                                 'Rebuild adapter from kmer counting file',
                                     formatter_class=MyHelpFormatter, add_help=False)

    main_group = parser.add_argument_group('Main options')
    main_group.add_argument('-s', '--start_counter', required=False,
                            help='Counter file for the starting adapter')
    main_group.add_argument('-e', '--end_counter', required=False,
                            help='Counter file for the end adapter')

    main_group.add_argument('-v', '--verbosity', type=int, default=1,
                            help='Level of progress information: 0 = none, 1 = all')

    main_group.add_argument('-g','--export_graph', action='store_true',
                                    help='Export the overlap graph')


    help_args = parser.add_argument_group('Help')
    help_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                           help='Show this help message and exit')
    
    args = parser.parse_args()


    return args


args = get_arguments()

start_file = args.start_counter
end_file = args.end_counter
verbosity = args.verbosity
  

greedy_adapter = {}
long_adapter   = {}
heavy_adapter  = {}



for count_file in [start_file, end_file]:
    if(count_file):
        
        # Building graph        
        g = build_graph( count_file )
        
        # preping for anotation
        nx.set_node_attributes(g, "", "path" )
        
        #################
        # greedy adapters
        greedy_adapter[count_file] = greedy_assembl(g)

        ##############
        # Longest path
        try:
            long_adapter[count_file] = longest_path(g)

        except:
            if(verbosity>0):
                print("Could not compute adapter using longest path  method", file = sys.stderr)
                print("The resulting graph probably contains a loop.", file = sys.stderr)
            long_adapter[count_file]  = ""

        ###############
        # heaviest path
        try:
            heavy_adapter[count_file] = heavy_path(g)
        except:
            if(verbosity>0):
                print("Could not compute adapter using heaviest path  method", file = sys.stderr)
                print("The resulting graph probably contains a loop.", file = sys.stderr)
            heavy_adapter[count_file] = ""

        #Exporting graph
        if(args.export_graph):
            path =  count_file + ".graphml"
            nx.write_graphml(g,path)


# printing adapters

for adapt_name, adapt_dict  in [("greedy",greedy_adapter),("longest_path",long_adapter), ("heaviest_path", heavy_adapter)]:

    print(adapt_name)    
    if(start_file):
        print(adapt_dict[start_file])
    if(end_file):
        print(adapt_dict[end_file])