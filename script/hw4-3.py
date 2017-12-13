import numpy as np
import itertools
import random
import sys

#Q3
class OverlapGraph(object):

    def __init__(self, kmers, kmerLen):
        '''Construct the overlap graph of a collection of k-mers.
        Attributes:
            kmers - list of kmer compositions.
            kmerLen - the length of the kmer.
        Methodes:
            overlapTest - tests if two kmers have overlap.
            edges - the overlap graph in the form of an adjacency list.'''
        self.kmers = kmers
        self.kmerLen = kmerLen

    #idea was taken from : https://stackoverflow.com/questions/27878067/rosalind-overlap-graphs
    def overlapTest(self,firstKmer, secondKmer): #test if there is an overlap between two motifs
        '''test if there is an overlap between two motifs.
        Returns - True if there is an overlap'''
        return firstKmer[-self.kmerLen + 1:] == secondKmer[:self.kmerLen -1]

    def edges(self):
        '''Connects edges.
        Return - Returns adjacency list of overlaped kmers.
        '''
        for firstKmer, secondKmer in itertools.combinations(self.kmers, 2): #all possible combinations between keys and values
            if self.overlapTest(firstKmer, secondKmer) == True:
                sys.stdout.write(firstKmer + ' -> ' +  secondKmer + '\n')
            elif self.overlapTest(secondKmer, firstKmer) == True:
                sys.stdout.write(secondKmer+ ' -> '  + firstKmer + '\n')
            else:
                pass

def main():
    '''Processing the input stdin, creating required objects for OverlapGraph, and initiation OverlapGraph'''
    kmers = [] #list of all kmers
    for line in sys.stdin:
        s = line.split('\n')[0]
        kmers.append(s)
    kmerLen = len(kmers[0]) # the kmer length

    m = OverlapGraph(kmers= kmers, kmerLen= kmerLen)
    m.edges()

if __name__ == "__main__":
    main()
