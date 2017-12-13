import numpy as np
import sys

#Q1
class KmerComposition(object):
    '''Creates a KmerComposition object with following attributes.
    Attributes:
        sequence - a list containing one sequence.
        kmerLen - the length of the kmer composition.
    Methods:
        randomNums - generats random number based on sequence length WO replacement, stores itin a list.
        kmers - returns a list of kmer compositions
    '''
    def __init__(self, sequence, kmerLen):
        self.sequence = sequence # a seqeunce in a list
        self.kmerLen = kmerLen # the length of the kmer

    def randomNums(self):
        '''Generates random numbers from a list WO replacement.
        Returns: list of randomNums'''
        length = len(self.sequence) #finding the length of the input sequence
        randomNums = np.random.choice(length - self.kmerLen + 1, length - self.kmerLen + 1, replace=False)
        return randomNums # return a list of random numbers

    def kmers(self):
        '''Generates kmer composition of the sequence.
        Returns: kmer composition in stdout format'''
        randomNum = self.randomNums()
        for i in randomNum: #take each number in the randomNum list as a kmer starting point
            sys.stdout.write(self.sequence[i: i + self.kmerLen] + "\n")

def main():
    '''Processes the stdin file, and stores it in seq and kmerLen to be used by Kmercomposition class.'''
    seq = '' # the sequence string
    kmerLen = ''# the length of the kmerComposition given by Rosalind in its first line
    for index, line in enumerate(sys.stdin):
        l = line.split("\n")[0] # removing the "\n"
        if index == 0: # Rosalind has the first line in the text as the kmer size
            kmerLen = int(l)
        elif index == 1: # Rosalind has the second line as the sequence
            seq += l

    m = KmerComposition(sequence= seq, kmerLen = kmerLen) #initiation kmerComposition class
    m.kmers()

if __name__ == "__main__":
    main()
