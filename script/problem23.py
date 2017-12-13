import numpy as np
import itertools
import sys

class RosalindParse(object):
    '''Parsing the Rosalind input.
    Attributes:
        inputFile -  the input rosalind file.
    Methodes:
        rosalindParse - parsing the input'''
    def __init__(self, inputFile):
        self.inputFile = inputFile

    def rosalidParse(self):
        '''parses the input file froom Rosalind.
        Returns:
            observedSeq - the observed sequence
            observedChar - the observed sequence alphabet (list of lists)
            stateSeq - the hidden state path
            states - the states
            '''
        infile = self.inputFile
        readInput = infile.readlines()
        observedSeq = readInput[0].strip()
        observedChar = readInput[2].strip().split('\t')
        stateSeq = readInput[4].strip()
        states = readInput[6].strip().split('\t')
        return observedSeq, observedChar, stateSeq, states

class ParaEstimation(object):
    '''HMM Parameter Estimation.
    Attributes:
        observedSeq - the observed sequence
        observedAlphabet - the observed alphabet
        stateSeq - the hiddenstate sequence
        states - the hidden states
    Methodes:
        transitionBackbone - creates the transition backbone
        emissionBackbone -creates the emission backbone
        possibleCombinations - finds all possible combinations between states
        count - counts the number of substrings in a string
        emission - creates the transition
        transition - creates the emission
    '''
    def __init__(self, observedSeq, observedAlphabet, stateSeq, states):
        self.observedSeq = observedSeq
        self.observedAlphabet = observedAlphabet
        self.stateSeq = stateSeq
        self.states = states

    def transitionBackbone(self): #this is for transition
        '''Constructing the transition backbone, dict of dic.
        Returns:
            the transition backbone, with dict of dict structure'''
        # here I'm creating the transition backbone
        transition = {}
        for state in self.states:
            if state not in transition:
                transition[state] = {}
        for k, v in transition.items():
            for i in self.states:
                transition[k][i] = 0
        return transition

    def emissionBackbone(self):
        '''Construction the emission backbone.
        Returns:
            an empty emission dict of dict'''
        # I'm creating the emission backbone
        emission = {}
        for state in self.states:
            if state not in emission:
                emission[state] = {}
        for k, v in emission.items():
            for i in self.observedAlphabet:
                emission[k][i] = 0
        return emission

    def possibleCombinations(self):
        '''List pf possible combinations.
        Returns:
            A list of possible combination'''
        # finding all possible combinations between states
        comb = [''.join(i) for i in itertools.product(self.states, repeat = 2)]
        return comb

    def emission(self):
        '''Create the emission probability, in dict of dict structure.
        Returns:
            the emission probability'''
        emission = self.emissionBackbone() # getting the backbone
        for i in range(len(self.stateSeq.strip())):
            ob = self.observedSeq[i]
            st = self.stateSeq[i]
            emission[st][ob] += 1

        for k, v in emission.items():
            su = sum(list(v.values()))
            if su != 0: # if the sum is not zero
                for j,p in v.items():
                    emission[k][j] = float(p/su) # if the sum is not zero divide each element by the sum
            else:# if the sum is zero
                for j,p in v.items():
                    emission[k][j] = float(1/len(self.observedAlphabet)) # if the sum is zero just devide the probability by the number of observed alphabets
        return emission

    def count(self,string, substring):
        '''Counting the number of substrings in a string.
        Returns the number of substrings occurences in a string'''
        stringSize = len(string)
        substringSize = len(substring)
        count = 0
        for i in range(0,stringSize-substringSize+1):
            if string[i:i+substringSize] == substring:
                count+=1
        return count

    def transition(self):
        '''Calculating the emission probability, dict of dict.
        Returns:
            the emission probability'''
        comb = self.possibleCombinations() # all possible combinations between states
        transition = self.transitionBackbone() # getting the transition backbone

        stateCount = {state : 0  for state in self.states} # creating a dictionary od states

        for state in self.states:
            stateCount[state] = self.count(self.stateSeq[:-1], state) # adding the count of number of times that each state occured in the observed stateSeq

        for combTrans in comb:
            transition[combTrans[0]][combTrans[1]] += float(self.count(self.stateSeq, combTrans)) # counting the number of times that each transition happend in the stateSeq

        for k,v in transition.items(): # for each key and value in the transition
            if sum(list(v.values())) != 0: # if the sum of values in the v at not zero
                for kk,vv in v.items():
                    transition[k][kk] = vv/stateCount[k] # do the divission to find the probability
            elif sum(list(v.values())) == 0:# if the sum of values in the v at not zero
                for kk,vv in v.items():
                    transition[k][kk] = 1/len(self.states) # just devide by the number of states 
        return transition




class OutputFormat(object):
    '''Parsing rosaling output format.
    Attributes:
        transition - tranistion
        emission - emission
        observedChar - observed characters
    Methods:
        outputParsing - STDOUT the output
     '''
    def __init__(self, transition, emission, observedChar, states):
        self.transition = transition
        self.emission = emission
        self.observedChar = observedChar
        self.states = states

    def outputParsing(self):
        '''STDOUT the output'''
        sys.stdout.write("\t" + '\t'.join(map(str, self.states)) + "\n")
        for key, values in self.transition.items():
            values = list(values.values())
            sys.stdout.write(key)
            for element in values:
                if element == 0:
                    sys.stdout.write('\t' + str(int(element)))
                else:
                    sys.stdout.write('\t' + str(float('%.3g'% element)))
            sys.stdout.write('\n')

        sys.stdout.write('--------' + '\n')

        sys.stdout.write('\t' + '\t'.join(map(str, self.observedChar)) + "\n")
        for key, values in self.emission.items():
            values = list(values.values())
            sys.stdout.write(key)
            for element in values:
                if element == 0:
                    sys.stdout.write('\t' + str(int(element)))
                else:
                    sys.stdout.write('\t' + str(float('%.3g'% element)))
            sys.stdout.write('\n')


def main():
    '''Process the input file, create an object of the class, and create the output'''
    fr = RosalindParse(sys.stdin)
    observedSeq, observedChar, stateSeq, states = fr.rosalidParse()
    x = ParaEstimation(observedSeq = observedSeq, observedAlphabet= observedChar, stateSeq= stateSeq, states= states)
    outputInit = OutputFormat(transition = x.transition(), emission =x.emission(), observedChar = observedChar, states= states)
    outputInit.outputParsing()
if __name__ == "__main__":
    main()
