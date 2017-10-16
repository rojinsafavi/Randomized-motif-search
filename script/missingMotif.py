#Author: Rojin Safavi#
#####################

import sys
import numpy as np
from fastaReader import FastAreader
import itertools
import argparse


class MotifSearch(object):
    '''Returns a MotifSearch object with the following attributes and methods.

    Attributes:
        minMotif -- minimum motif size to evaluate (default 3)
        maxMotif -- maximum motif size to evaluate (default 5)
        cutoff -- Z-score cutoff (default -5)
        inputFile -- STDIN object
    Methods:
        kmerDic -- creates a kmer dic
        inputSeq --
        kmerCount -- counts the number of each motif
        expectedCount -- calculates the expected number for each motif
        kmerProb -- calculates the probability under null
        cutOff -- removes kmers above cutoff threshold
        dataFrame -- sorting the kmer dict
    '''
    def __init__(self, inputFile, minMotif=3, maxMotif=8, cutoff=-5):
        '''MotifSearch constructor, initialize a new instance of MotifSearch.
        '''
        self.minMotif = minMotif  # min len of the mottif, it has to be -3 at least
        self.maxMotif = maxMotif  # max can be atmost 8
        self.cutoff = cutoff  # cutoff default should be -5
        self.inputFile = inputFile  # path of the input file

    def kmerDic(self, alphabet="ATCG"):
        ''' Returns a dictionary of kmers.

        Arguments:
            alphabet -- alphabets to use to create the dict (default = "ATCG")
        Returns:
            a dictionary of kmers
        '''
        dic = {}
        keysList = []
        alphabet = ''.join(sorted(alphabet))
        for motifLength in np.arange(1, self.maxMotif + 1):  # range = 1,2,3,4,5,6,7,8
            kmers = [''.join(p) for p in itertools.product(alphabet, repeat=motifLength)]  # https://github.com/UCSC-nanopore-cgl/nanopore-RNN/blob/master/nanotensor/create_training_data.py
            keysList = keysList + kmers
        for kmer in keysList:
            dic[kmer] = [kmer, len(kmer), 0]  # kmer, kmer_len, actual count (starting from 0)
        return dic  # dic with all possible kmers from 1 through 8 inclusive

    def inputSeq(self):
        '''Concatenate the input fasta with \n between each sequence.

        Returns:
                a list of all sequences in a fasta file
        '''
        fastaInitiate = FastAreader(self.inputFile)  # initiating FastAreader
        all_fasta = list(fastaInitiate.readFasta())
        seq = ''
        for fastaSeq in np.arange(len(all_fasta)):
            seq += all_fasta[fastaSeq][1] + "\n"
        return seq  # concate all fasta files and add \n in between them

    def kmerCount(self):
        '''Counting kmers.

        Returns:
            A dict with the number of each motif
            '''
        dic = self.kmerDic()  # loading the dic
        sequence = self.inputSeq()  # loading the sequence
        sequenceLen = len(sequence)  # findnig the seq len (includes 'n')
        for index in np.arange(0, sequenceLen + 1):
            for kmerLen in np.arange(1, self.maxMotif + 1):
                kmer = sequence[index:index + kmerLen]
                if kmer in dic:  # only counting canonical nucleotides
                    dic[kmer][2] = dic[kmer][2] + 1
        self.sequenceLen = sequenceLen
        return dic

    def expectedCount(self):
        '''Calculates the expected motif count.

        Returns:
            dict with an expected count column
        '''
        dic = self.kmerCount()  # get the dic with counted values
        for kmer, value in dic.items():
            if len(kmer) >= 3:  # ec ATC
                firstPortion = kmer[0:-1]  # ex AT
                secondPortion = kmer[1:]  # ex  TC
                middlePortion = kmer[1:-1]  # ex T
                firstValue = dic[firstPortion][2]  # findin the actual count of AT
                secondValue = dic[secondPortion][2]  # findin the actual count of TC
                middleValue = dic[middlePortion][2]  # findin the actual count of T
                # now, what if the actual count of T is zero? maybe I should set the middle_portion to 1 if there is nothing?
                if middleValue == 0:
                    pass
                else:
                    dic[kmer].append(firstValue*secondValue/middleValue) # return e(kmer)
        return dic

    def kmerProb(self):
        ''' Calculates markov(0) approximation for the expected count (E).

        Returns:
            A dict with kmer probability
        '''
        dic = self.expectedCount()  # loading the dic: kmer, kmer_len, actual count, expected count
        sequenceLen = self.sequenceLen
        for kmer, value in dic.items():
            if len(value) >= 4:
                probability = float(value[3]/sequenceLen)
                if probability == 0:
                    pass
                else:
                    mean = value[3]
                    sd = (sequenceLen*probability*(1-probability))**(1/2)
                    zScore = (value[2]-mean)/sd
                    dic[kmer].append(zScore)
        return dic

    def cutOff(self):
        '''removes rows with zScores above the cutoff.

        Returns:
            a dict with zScores bellow the cutoff'''
        dic = self.kmerProb()
        for key in list(dic):
            if len(dic[key]) <= 4:
                dic.pop(key, None)

        for key in list(dic):
            if len(dic[key]) >= 5:
                if dic[key][4] > self.cutoff:
                    dic.pop(key, None)

        return dic

    def dataFrame(self):
        '''Adds labels and sorts the dictionary and returns STDOUT.

        Returns:
            STDOUT a sorted and labled table
        '''
        finalDic = self.cutOff()
        valueList = []
        for i, k in finalDic.items():
            valueList.append(k)
        valueList.sort(key=lambda x: (x[1],-x[4]), reverse= True)
        for i in valueList:
            del i[1]
        valueList.insert(0, ['Motif', 'actual count', 'expected count', 'Zscore'])
        for i in valueList:
            sys.stdout.write("\t".join(map(str, i)) + '\n')

def main():
    '''creates an instance of the MotifSearch class, sets the arguments, and calls the methods'''
    parser = argparse.ArgumentParser(description="finding underrepesented motifs")
    parser.add_argument('-min', '--minMotif', type=int, default=3, help='the minimum motif size to evaluate')
    parser.add_argument('-max', '--maxMotif', type=int, default=8, help='the maximum motif size to evaluate')
    parser.add_argument('-cut', '--cutoff', type=int, default=-5, help='the Z-score cutoff')
    args = parser.parse_args()
    f = MotifSearch("", args.minMotif, args.maxMotif, args.cutoff)
    f.dataFrame()

if __name__ == "__main__":
    main()
