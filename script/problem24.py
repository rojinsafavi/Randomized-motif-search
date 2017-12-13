from problem23 import ParaEstimation
from problem23 import OutputFormat
from problem18 import ViterbiPath
import sys
import numpy as np


class RosalindParse(object):
    '''Parsing Roslaind input.
    Attributes:
        inputFile
    Methodes:
        rosalindParse - parsing the input'''

    def __init__(self, inputFile):
        self.inputFile = inputFile

    def rosalidParse(self):
        '''parsing the input
        Returns:
            numIteration
            observedSeq
            observedAlphabet
            states
            initTransition - initial transition
            initEmission - initial emission '''
        infile = self.inputFile
        readInput = infile.readlines()
        numIteration = int(readInput[0].strip())
        observedSeq = readInput[2].strip()
        observedAlphabet = readInput[4].strip().split('\t')
        states = readInput[6].strip().split('\t')

        transition = []
        for i in readInput[9:9 + len(states)]:
            transition.append(i.strip().split('\t')[1:])
        transition = np.asanyarray(transition)
        transition = transition.astype(np.float)

        emission = []
        for i in readInput[9 + len(states) + 2:]:
            emission.append(i.strip().split('\t')[1:])
        emission = np.asanyarray(emission)

        emission = emission.astype(np.float)

        initTransition = {} # here I'm just converting the transition matrix to a dict of dict
        for i,j in enumerate(states):
            for k, n in  enumerate(states):
                if j not in initTransition:
                    initTransition[j] = {}
                    initTransition[j][n] = transition[i][k]
                else:
                    initTransition[j][n] = transition[i][k]

        initEmission= {}# here I'm just converting the emission matrix to a dict of dict
        for i,j in enumerate(states):
            for k, n in enumerate(observedAlphabet):
                if j not in initEmission:
                    initEmission[j] = {}
                    initEmission[j][n] = emission[i][k]
                else:
                    initEmission[j][n] = emission[i][k]

        return numIteration, observedSeq, observedAlphabet, states, initTransition, initEmission



class ViterbiLearning(object):
    '''Viterbi Learning
    Attributes:
        numIteration
        observedSeq
        observedAlphabet
        states
        initTransition - initial transition
        initEmission - initial emission
    Methods:
        viterbilearning - outputs the final tranisiton and emission matricis'''
    def __init__(self, numIteration,observedSeq,observedAlphabet,states, initTransition, initEmission ):
        self.numIteration = numIteration
        self.observedSeq = observedSeq
        self.observedAlphabet = observedAlphabet
        self.states = states
        self.initTransition = initTransition
        self.initEmission = initEmission

    def viterbilearning(self):
        ''' A matrix of transition probabilities and a matrix of emission probabilities.
        Returns:
            transition - transition probabilities
            emission - emission probabilities'''
        for i in range(self.numIteration + 1):
            if i == 0:
                transition = self.initTransition # getting the initial transition
                emission = self.initEmission # getting the initial emission
                #here I'm finding the viterbi path based on the initial transition and emission dicts
                vp = ViterbiPath(observation = self.observedSeq , transition = transition, emission= emission, states = self.states, observedChar = self.observedAlphabet)
                transition = vp.transition # the new tranisiton
                emission = vp.emission # the new emission
            else:
                # here I'm using the previouse transition and emission to find the viterbi path
                vp = ViterbiPath(observation = self.observedSeq , transition = transition, emission= emission, states = self.states, observedChar = self.observedAlphabet)
                stateSeq = vp.viterbiPath() # finding the new stateSeq
                # getting the parameter Estimation
                x = ParaEstimation(observedSeq = self.observedSeq , observedAlphabet= self.observedAlphabet, stateSeq= stateSeq, states= self.states)
                emission = x.emission() # new emission
                transition = x.transition() # new transition
        return transition, emission



def main():
    '''Process the input file, create an object of the class, and create the output'''
    fr = RosalindParse(sys.stdin)
    numIteration,observedSeq,observedAlphabet,states, initTransition, initEmission = fr.rosalidParse()

    TE= ViterbiLearning(numIteration,observedSeq,observedAlphabet,states, initTransition, initEmission )
    transition,emission = TE.viterbilearning()
    g = OutputFormat(transition = transition, emission = emission, observedChar = observedAlphabet, states = states)
    g.outputParsing()

if __name__ == "__main__":
    main()
