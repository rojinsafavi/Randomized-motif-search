
#Q5
import numpy as np
import sys
from collections import defaultdict

class BruijnGraph(object):
    '''Construct the de Bruijn graph from a list of kmers.
    Attributes:
        edgeSize - the edge length
        edgesList - the sequence
    Methodes:
        nodesDic - dict representation of the graph
        output - Rosalind output format'''

    def __init__(self, edgesList, edgeSize):
        self.edgeSize = edgeSize
        self.edgesList = edgesList

    def nodesDic(self):
        '''Returns a dict of nodes.
        Returns:
            a graph represetation of the nodes in a dict format'''
        nodes = defaultdict(list)
        for kmer in self.edgesList:
            start = kmer[:-1]
            end = kmer[1:]
            if start not in nodes:
                nodes[start] = [end]
            else:
                nodes[start].append(end)
        return nodes

    def output(self):
        '''Rosalind output format.
        Returns:
            STDOUT Rosalind output format'''
        dic = self.nodesDic()
        for key, value in dic.items():
            sys.stdout.write(key + " -> " + ','.join(value) + "\n")


def main():
    '''creates an instance of the BruijnGraph class, sets the arguments, and calls the method'''
    edgesList = []
    kmerLen = ''
    for line in sys.stdin:
        l = line.split("\n")[0]
        kmerLen = len(l)
        edgesList.append(l)


    m = BruijnGraph(edgesList = edgesList, edgeSize = kmerLen)

    m.output()


if __name__ == "__main__":
    main()
