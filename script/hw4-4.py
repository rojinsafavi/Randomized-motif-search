#Q4
import numpy as np
import sys
from collections import defaultdict


class BruijnConstruct(object):
    '''Created DeBruij graph object in the form of an adjacency list.
    Attributes:
        inputSequence - sequence
        edgeSize - edge size
    Methodes:
        stringSplit - split an edge into two nodes
        nodesDic - creates a dictionary of nodes
        output - stdout the Rosalind format'''

    def __init__(self, inputSequence, edgeSize):
        self.edgeSize = edgeSize
        self.inputSequence = inputSequence

    def stringSplit(self):
        '''Split the sequence into edges.
        Returns:
            a list of edges of length edgeSize'''
        edgesList= []
        for i in np.arange(0, len(self.inputSequence) - self.edgeSize + 1):
            edgesList.append(self.inputSequence[i: i + self.edgeSize]) # creating a list of edges
        return edgesList

    def nodesDic(self):
        '''Splitting all esges into two nodes.
        Returns:
            a dict of all nodes and their connection'''
        edges = self.stringSplit() # list of edges
        nodes = defaultdict(list) # dict of nodes
        for kmer in edges:
            start = kmer[:-1] # first part of the kmer ex: ACG from ACGA
            end = kmer[1:] # last part of the kmer ex CGA from ACGA
            if start not in nodes:
                nodes[start] = [end]
            else:
                nodes[start].append(end)
        return nodes

    def output(self):
        '''Rosalind output in an adjacency list format.
        Returns:
            adjacency list of nodes'''
        dic = self.nodesDic()
        for key, value in dic.items():
            sys.stdout.write(key + " -> " + ','.join(value) + "\n")


def main():
    '''Pre-processing the stdin, creating necessary objects, and feed them into BruijnConstruct class.'''
    inputSequence = '' # the sequence
    kmerLen = '' #the edge length
    for index, line in enumerate(sys.stdin):
        l = line.split("\n")[0]
        if index == 0: # first line is the edge length
            kmerLen = int(l)
        elif index == 1: # second line is the seq
            inputSequence += l

    m = BruijnConstruct(inputSequence = inputSequence , edgeSize = kmerLen)

    m.output()

if __name__ == "__main__":
    main()
