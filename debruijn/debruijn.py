#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import statistics
import random
import networkx as nx
import matplotlib

random.seed(9001)

__author__ = "Lapeyre Adrien"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Lapeyre Adrien"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Lapeyre Adrien"
__email__ = "lapeyreadr@eisti.eu"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()


def read_fastq(fastq_file):
    """Returns a generator of sequences
      :Parameters:
          fastq_file: Path to the fastq_file
      Returns: generator of sequences
    """
    # Open the file in read mode
    with open(fastq_file, 'r') as my_file:
        lines = my_file.readlines()
        # Starts on the second line and focus on 1 line over 4
        for line in lines[1::4]:
            # Return the line with the sequence and without the \n character
            yield line[:-1]


def cut_kmer(read, kmer_size):
    """Returns a generator of kmers
      :Parameters:
          read: Sequence (string)
          kmer_size: Size of kmers (int)
      Returns: generator of kmers
    """
    # Return each kmer of a sequence
    for i in range(len(read)-kmer_size+1):
        yield read[i:i+kmer_size]


def build_kmer_dict(fastq_file, kmer_size):
    """Builds a dictionnary of kmers
      :Parameters:
          fastq_file: Path to the fastq_file
          kmer_size: Size of kmers (int)
      Returns: dictionnary of kmers and their number of occurences (dict)
    """
    dic = {}
    # Sequence generator
    read1 = read_fastq(fastq_file)
    for seq in read1:
        # Kmer generator
        read2 = cut_kmer(seq, kmer_size)
        # Add each kmer in the dictionnary if it doesn't exist or +1 to his value
        for kmer in read2:
            if kmer not in dic:
                dic[kmer] = 1
            else:
                dic[kmer] += 1
    return dic


def build_graph(kmer_dict):
    """Builds a networkx digraph
      :Parameters:
          kmer_dict: Dictionnary of kmers (dict)
      Returns: digraph (nx.Digraph)
    """
    graph = nx.DiGraph()
    # Build the digraph adding each suffix and prefix kmers as nodes with edges and their weights
    for key, value in kmer_dict.items():
        graph.add_edge(key[:-1], key[1:], weight=value)
    return graph


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    """Removes paths of a list in the graph
      :Parameters:
          graph: (nx.Digraph)
          path_list: list of paths (list)
          delete_entry_node: (bool)
          delete_sink_node: (bool)
      Returns: graph (nx.Digraph)
    """
    for path in path_list:
        # Delete the entry node if true
        if delete_entry_node:
            graph.remove_node(path[0])
        # Delete the sink node if true
        if delete_sink_node:
            graph.remove_node(path[-1])
        # Delete other nodes
        graph.remove_nodes_from(path[1:-1])
    return graph


def std(data):
    """Computes standard deviation of a list of numbers
      :Parameters:
          data: (list)
      Returns: standard deviation (float)
    """
    return statistics.stdev(data)


def select_best_path(graph, path_list, path_length, weight_avg_list,
                     delete_entry_node=False, delete_sink_node=False):
    """Selects the best path of a list and removes others in the graph
      :Parameters:
          graph: (nx.Digraph)
          path_list: list of paths (list)
          path_length: list of paths lengths (list)
          weight_avg_list: list of paths weights (list)
          delete_entry_node: (bool)
          delete_sink_node: (bool)
      Returns: graph (nx.Digraph)
    """
    # Initialize best path to the first one
    best_path_index = 0
    # Search the best path
    for i in range(1, len(path_list)):
        weight1 = weight_avg_list[best_path_index]
        weight2 = weight_avg_list[i]
        length1 = path_length[best_path_index]
        length2 = path_length[i]
        # Select the weightest one
        if weight1 < weight2:
            best_path_index = i
        elif weight1 == weight2:
            # Select the lengthest one
            if length1 < length2:
                best_path_index = i
            elif length1 == length2:
                # Select randomly
                best_path_index = random.choice([best_path_index, i])
    # Remove the best path from paths to remove
    path_list.remove(path_list[best_path_index])
    # Remove other paths from the graph
    graph = remove_paths(graph, path_list, delete_entry_node, delete_sink_node)
    return graph


def path_average_weight(graph, path):
    """Computes the average weight of a path
      :Parameters:
          graph: (nx.Digraph)
          path: list of nodes (list)
      Returns: average weight (int)
    """
    total_weight = 0
    for i in range(len(path)-1):
        total_weight += graph[path[i]][path[i+1]]['weight']
    mean_weight = total_weight/(len(path)-1)
    return mean_weight


def solve_bubble(graph, ancestor_node, descendant_node):
    """Remove a selected bubble in the graph
      :Parameters:
          graph: (nx.Digraph)
          ancestor_node: node
          descendant_node: node
      Returns: graph (nx.Digraph)
    """
    # Compute all possible paths between ancestor and descendant node
    path_list = list(nx.all_simple_paths(graph, source=ancestor_node, target=descendant_node))
    path_length = []
    weight_avg_list = []
    # Compute their length and average weight
    for path in path_list:
        path_length.append(len(path))
        weight_avg_list.append(path_average_weight(graph, path))
    # Remove paths that aren't the best
    graph = select_best_path(graph, path_list, path_length, weight_avg_list)
    return graph


def simplify_bubbles(graph):
    """Remove all bubbles in the graph
      :Parameters:
          graph: (nx.Digraph)
      Returns: graph (nx.Digraph)
    """
    ancestor_nodes = []
    descendant_nodes = []
    for node in graph.nodes:
        # Find predecessors of the node
        predecessors = list(graph.predecessors(node))
        # If the node as more than 1 predecessor (might be a bubble)
        if len(predecessors) > 1:
            for i in range(len(predecessors)-1):
                # Find the lowest common ancestor of two of the predecessors
                ancestor = nx.lowest_common_ancestor(graph, predecessors[i], predecessors[i+1])
                # If an ancestor exists
                if ancestor != None:
                    # Save the two node at the extremity of the bubble
                    ancestor_nodes.append(ancestor)
                    descendant_nodes.append(node)
                    break
    # Delete all the bubles finded in the graph
    for ancestor, descendant in zip(ancestor_nodes, descendant_nodes):
        if ancestor in graph.nodes and descendant in graph.nodes:
            graph = solve_bubble(graph, ancestor, descendant)
    return graph


def solve_entry_tips(graph, starting_nodes):
    pass


def solve_out_tips(graph, ending_nodes):
    pass


def get_starting_nodes(graph):
    """Finds all starting nodes
      :Parameters:
          graph: (nx.Digraph)
      Returns: list of starting nodes (list)
    """
    return [n for n, d in graph.in_degree() if d == 0]


def get_sink_nodes(graph):
    """Finds all sink nodes
      :Parameters:
          graph: (nx.Digraph)
      Returns: list of sink nodes (list)
    """
    return [n for n, d in graph.out_degree() if d == 0]


def get_contigs(graph, starting_nodes, ending_nodes):
    """Finds all contigs
      :Parameters:
          graph: (nx.Digraph)
          starting_nodes: list of starting nodes (list)
          ending_nodes: list of sink nodes (list)
      Returns: tuples of all contigs and their length (list)
    """
    contigs = []
    for start in starting_nodes:
        for end in ending_nodes:
            # Compute all possible paths between start and end node
            paths = nx.all_simple_paths(graph, source=start, target=end)
            for path in paths:
                # Create the contig
                contig = path[0]
                for seq in path[1:]:
                    contig += seq[-1]
                contigs.append((contig, len(contig)))
    return contigs


def save_contigs(contigs_list, output_file):
    """Creates the contigs.fasta file
      :Parameters:
          contigs_list: tuples of all contigs and their length (list)
          output_file: path of the contigs.fasta file
    """
    text = ''
    cpt = 0
    # Prepare the text to fasta format
    for element in contigs_list:
        text += ('>contig_' + str(cpt) + ' len=' + str(element[1]) + '\n'
                 + fill(element[0]) + '\n')
        cpt += 1
    # Write it in the output_file
    with open(output_file, 'w') as my_file:
        my_file.write(text)


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    # Create the digraph
    dic = build_kmer_dict(args.fastq_file, args.kmer_size)
    graph = build_graph(dic)

    # Simplify the graph
    graph = simplify_bubbles(graph)

    # Graph search and creation of the contigs.fasta file
    starting_nodes = get_starting_nodes(graph)
    ending_nodes = get_sink_nodes(graph)
    contigs_list = get_contigs(graph, starting_nodes, ending_nodes)
    save_contigs(contigs_list, args.output_file)


if __name__ == '__main__':
    main()
