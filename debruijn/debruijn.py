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
import networkx as nx
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics

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
    # Open the file in read mode
    with open(fastq_file, 'r') as my_file:
        for line in my_file:
            # Return the line with the sequence and without the \n character
            yield next(my_file)[:-1]
            next(my_file)
            next(my_file)


def cut_kmer(read, kmer_size):
    # Return each kmer of a sequence
    for i in range(len(read)-kmer_size+1):
        yield read[i:i+kmer_size]


def build_kmer_dict(fastq_file, kmer_size):
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
    graph = nx.DiGraph()
    # Build the digraph adding each suffix and prefix kmers as nodes with edges and their weights
    for key, value in kmer_dict.items():
        graph.add_edge(key[:-1], key[1:], weight=value)
    return graph


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    pass


def std(data):
    pass


def select_best_path(graph, path_list, path_length, weight_avg_list,
                     delete_entry_node=False, delete_sink_node=False):
    pass


def path_average_weight(graph, path):
    pass


def solve_bubble(graph, ancestor_node, descendant_node):
    pass


def simplify_bubbles(graph):
    pass


def solve_entry_tips(graph, starting_nodes):
    pass


def solve_out_tips(graph, ending_nodes):
    pass


def get_starting_nodes(graph):
    return [n for n, d in graph.in_degree() if d == 0]


def get_sink_nodes(graph):
    return [n for n, d in graph.out_degree() if d == 0]


def get_contigs(graph, starting_nodes, ending_nodes):
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

    # Graph search and creation of the contigs.fasta file
    starting_nodes = get_starting_nodes(graph)
    ending_nodes = get_sink_nodes(graph)
    contigs_list = get_contigs(graph, starting_nodes, ending_nodes)
    save_contigs(contigs_list, args.output_file)


if __name__ == '__main__':
    main()
