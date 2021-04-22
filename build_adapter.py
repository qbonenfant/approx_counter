"""
Standalone script for adapter reconstruction from k-mer count file.
"""

import os
import shutil
import argparse
import subprocess
from multiprocessing import cpu_count
from collections import defaultdict as dd
import networkx as nx
import sys


##############################################################################
#                                 CONSTANTS                                  #
##############################################################################

CUT_RATIO = 1.05
METHODS = ["greedy", "heavy"]


##############################################################################
#                              UTILITY FUNCTIONS                             #
##############################################################################


def haveOverlap(seq1, seq2):
    """Check if the sequence 1 is a prefix of sequence 2
    @param first sequence
    @param second sequence
    @return is there an overlap ?
    """
    minOverlap = min(len(seq1), len(seq2)) - 1
    if(seq1[-minOverlap:] == seq2[:minOverlap]):
        return(seq1 + seq2[minOverlap:])
    else:
        return("")


def get_weight(g, path):
    """Compute the weight of a path in the graph
    @param the graph
    @param the path (list of nodes)
    @return the total weight of this path

    """
    total = 0
    for node in path:
        total += g.nodes[node]["weight"]
    return(total)


def concat_path(path):
    """Concat the kmers of a path into a single sequence
    @param The path as a k-mer list
    @return the full sequence.
    """
    return(path[0][:-1] + "".join(el[-1] for el in path))


def check_drop(path, g):
    """Check if there is a frequency drop in the path
    and cut the adapter if necessary.
    This function cut the end of forward adapters.
    @param The path, as a kmer list
    @param The De Bruijn graph.
    @return The adjusted adapter path
    """
    last = 0
    cut = 0
    for i in range(len(path) - 1, 1, -1):
        last = path[i]
        prev = path[i - 1]
        if(g.nodes[prev]["weight"] / float(g.nodes[last]["weight"]) > CUT_RATIO):
            cut -= 1
        else:
            break
    if(cut == 0):
        return(path)
    return(path[:cut])


def check_drop_back(path, g):
    """Check if there is a frequency drop in the path
    and cut the adapter if necessary.
    This function cut the start of reverse adapters.
    @param The path, as a kmer list
    @param The De Bruijn graph.
    @return The adjusted adapter path
    """
    first = 0
    cut = 0
    for i in range(len(path) - 1):
        first = path[i]
        nxt = path[i + 1]
        if(g.nodes[nxt]["weight"] / float(g.nodes[first]["weight"]) > CUT_RATIO):
            cut += 1
        else:
            break
    return(path[cut:])


def print_result(adapters, v=1, print_dest=sys.stdout):
    """Display the result, adjusting format depending on selected verbosity
    @param an adapter dictionnary, using appropriate methods as keys
    @param v, the verbosity level, default to 1
    @param print destination, can be custom, but stdout is default
    """

    global METHODS

    out = print_dest

    if(v > 0):
        print("\n\nINFERRED ADAPTERS:\n",
              file=out)
    else:
        out = sys.stdout

    for meth in METHODS:
        msg = meth
        srt = adapters[meth]["start"]
        end = adapters[meth]["end"]

        if(v >= 1):
            meth += " assembly method"
            meth = meth.capitalize()
            srt = "Start:\t" + srt
            end = "End:\t" + end

        print(msg, file=out)
        print(srt, file=out)
        print(end, file=out)


##############################################################################
#                                 ASSEMBLY                                   #
##############################################################################

def greedy_assembl(g):
    """Greedy assembly method to compute the adapter using the graph as input.
    @param the De Bruijn graph of kmers
    @return the longest debruijn sequence starting by the first kmer
    """
    start = max(g.nodes, key=lambda x: g.nodes[x]["weight"])
    path = [start]

    right_node = start
    left_node = start
    while(left_node or right_node):
        # # DEBUG
        # print("PATH")
        # print(path)

        # forward extension
        if(right_node):
            # # DEBUG
            # print("R node",right_node)
            r_list = [el for el in g.successors(right_node) if el not in path]
            if(not r_list):
                right_node = None
            else:
                right_node = max(r_list, key=lambda x: g.nodes[x]["weight"])
                path.append(right_node)

        # reverse extension
        if(left_node):
            # # DEBUG
            # print("L node", left_node)
            l_list = [el for el in g.predecessors(left_node) if el not in path]

            if(not l_list):
                left_node = None
            else:
                left_node = max(l_list, key=lambda x: g.nodes[x]["weight"])
                path = [left_node] + path
    return(path)

    kmer_dict = g.nodes(data=True)
    kmer_list = list(dict(kmer_dict).keys())
    kmer_list.sort(key=lambda x: kmer_dict[x]['weight'])
    kmer_list.reverse()

    km = kmer_list[0]
    ov = km
    found = True
    used = [km]
    # annotating graph
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
        pairs = [(dist[v][0] + G.nodes[v]["weight"], v) for v in G.pred[node]]
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

    hv_path = dag_heaviest_path(g)

    # annotating graph
    for n in hv_path:
        g.nodes[n]["path"] = (g.nodes[n]["path"] + ",heavy").lstrip(",")

    return(hv_path)


def longest_path(g):
    """ DEPRECATED - This function is not used anymore
        Searching the longest path between all source and target nodes.
    """
    lg_path = nx.dag_longest_path(g)

    # annotating graph
    for n in lg_path:
        g.nodes[n]["path"] = (g.nodes[n]["path"] + ",long").lstrip(",")

    return(lg_path)


def greedy_path(g):
    """Greedy assembly of the adapter from the graph,
    starting from the most frequent
    """

    gd_path = greedy_assembl(g)
    # annotating graph
    for n in gd_path:
        g.nodes[n]["path"] = (g.nodes[n]["path"] + ",greedy").lstrip(",")
    return(gd_path)

##############################################################################
#                              GRAPH BUILDING                                #
##############################################################################


def build_graph(count_file):
    """Build the adapter from the kmer count file.
       The way it is done is by building a directed weighted graph
       and searching for the heaviest path.
       I also added the greedy adapter output
    """

    kmer_count = {}
    kmer_list = []

    g = nx.DiGraph()

    # Avoiding IO error
    try:
        # building graph
        with open(count_file, 'r') as f:
            for line in f:
                km, nb = line.rstrip("\n").split("\t")
                kmer_count[km] = int(nb)
                kmer_list.append(km)
                g.add_node(km, weight=int(nb))  # , path = "" )

    except FileNotFoundError:
        print("\n/!\\ Unable to open k-mer count file:", file=sys.stderr)
        print(count_file, file=sys.stderr)
        print("Either the file was moved, deleted, or filename is invalid.",
              file=sys.stderr)
        print("It is also possible that end adapter ressearch was skipped.",
              file=sys.stderr)
        print("Be sure skip_end / se option is deactivated in adaptFinder.\n",
              file=sys.stderr)

    else:
        # searching overlaps
        for km in kmer_list:
            for km2 in kmer_list:
                # since i do all the possible combinations
                # There is no need to test both orientation
                direct = haveOverlap(km, km2)
                if(direct != ""):
                    g.add_edge(km, km2)

        # removing singletons
        g.remove_nodes_from([node for node in g.nodes if len(
            list(nx.all_neighbors(g, node))) == 0])

        # Returning only the biggest connected component
        g = g.subgraph(max(nx.weakly_connected_components(g),
                           key=lambda x: get_weight(g, x)))
    finally:
        return(g)

##############################################################################
#                              ADPATER BUILDING                              #
##############################################################################


def build_adapter(args):
    """Building adapters from counts using different method:
        - Building a Debruijn graph and searching heaviest path
        - Greedy assembly based on kmer rank (most frequent first)
    @param An args dict containting a kmer count file prefix
    """

    v = args.verbosity
    unable_to_build = False

    adapters = dd(dict)
    adp = []

    # Removing the "start/end" from the filename if specified.
    if("start" in args.input or "end" in args.input):
        out_file_name = ".".join(args.input.split(".")[:-1])
    else:
        out_file_name = args.input

    # building adapters for each ends
    for which_end in ["start", "end"]:
        if(v >= 1):
            print("Assembling " + which_end + " adapters")

        # Building graph
        g = build_graph(out_file_name + "." + which_end)

        # Checking graph is properly build
        if(not nx.is_empty(g)):

            # preping for anotation
            nx.set_node_attributes(g, "", "path")

            # Greedy Adapter
            if(v >= 1):
                print("\tBuilding greedy adapter")

            greedy_p = greedy_path(g)
            cut_greedy_p = []
            if(which_end == "start"):
                cut_greedy_p = check_drop(greedy_p, g)
            elif(which_end == "end"):
                cut_greedy_p = check_drop_back(greedy_p, g)

            adapters["greedy"][which_end] = concat_path(cut_greedy_p)

            # Heavy adapter
            if(v >= 1):
                print("\tBuilding heavy path adapter", file=print_dest)
            try:
                heavy_p = heavy_path(g)
                cut_heavy_p = []
                if(which_end == "start"):
                    cut_heavy_p = check_drop(heavy_p, g)
                elif(which_end == "end"):
                    cut_heavy_p = check_drop_back(heavy_p, g)
                adapters["heavy"][which_end] = concat_path(cut_heavy_p)

            except nx.exception.NetworkXUnfeasible:
                print("\t/!\\ Could not compute " + which_end + " adaper"+
                      "using heaviest path  method",
                      file=sys.stderr)
                print("\t/!\\ The resulting graph probably contains a loop.",
                      file=sys.stderr)
                adapters["heavy"][which_end] = [""]

            # Exporting, if required
            if(args.export_graph is not None):
                if(v >= 1):
                    print("\tExporting assembly graph", file=print_dest)

                base = "_adapter_graph"
                if(".graphml" in args.export_graph):
                    base = "_" + os.path.basename(args.export_graph)
                else:
                    path = os.path.join(os.path.dirname(
                        args.export_graph), which_end + base + ".graphml")

                nx.write_graphml(g, path)
        else:
            unable_to_build = True

    if(not unable_to_build):
        # We just need to print the adapter
        return(adapters)
    else:
        sys.stderr.write("/!\\ERROR: Unable to build adapter. QUIT\n")
        exit(1)


##############################################################################
#                                  MAIN                                      #
##############################################################################


class MyHelpFormatter(argparse.HelpFormatter):

    """
    This is a custom formatter class for argparse. It allows for some custom
    formatting, in particular for the help texts with multiple options
    (like bridging mode and verbosity level).
    http://stackoverflow.com/questions/3853722
    """

    def __init__(self, prog):
        terminal_width = shutil.get_terminal_size().columns
        os.environ['COLUMNS'] = str(terminal_width)
        max_help_position = min(max(24, terminal_width // 3), 40)
        super().__init__(prog, max_help_position=max_help_position)

    def _get_help_string(self, action):
        help_text = action.help
        if action.default != argparse.SUPPRESS and \
           'default' not in help_text.lower() and \
           action.default is not None:

            help_text += ' (default: ' + str(action.default) + ')'
        return help_text


def get_arguments():
    """
    Parse the command line arguments.
    """
    parser = argparse.ArgumentParser(description='Build adapter'
                                     'A tool for building adapters in Oxford '
                                     'Nanopore reads from k-mer count.',
                                     formatter_class=MyHelpFormatter,
                                     add_help=False)

    main_group = parser.add_argument_group('Main options')
    main_group.add_argument('-i', '--input', required=True,
                            help='Path to k-mer count file or file name prefix\
                            (without the ".start" or ".end)".')
    main_group.add_argument('-v', '--verbosity', type=int, default=0,
                            help='Level of info display 0 = minimal, 1 = std,'
                                 '2 = a lot.')

    graph_group = parser.add_argument_group('Graph mangement options')
    graph_group.add_argument('--export_graph', type=str,
                             help='Path to export the graph used for assembly\
                             (.graphml format), if you want to keep it')

    help_args = parser.add_argument_group('Help')
    help_args.add_argument('-h', '--help', action='help',
                           default=argparse.SUPPRESS,
                           help='Show this help message and exit')
    args = parser.parse_args()

    return args


if __name__ == '__main__':

    # Parse args
    args = get_arguments()

    # Build adapt
    adapt = build_adapter(args)

    # Print adapters
    print_result(adapt, v = args.verbosity)
