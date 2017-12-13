import numpy as np
import sys

class PathProb(object):
    '''Finds the probability of a hidden path.
    Attributes:
        sequence - a sequence of states
        states - hidden states
        transitionMatrix - transition matrix'''
    def __init__(self, sequence, states, transitionMatrix):
        self.sequence = sequence
        self.states = states
        self.transitionMatrix = transitionMatrix

    def pathProb(self):
        '''Returns the probability of the path.
        Returns:
            the probability of the path'''
        m = 1
        for i in range(len(self.sequence) - 1):
            stateOne = self.sequence[i]
            stateTwo = self.sequence[i+1]
            indexOne = self.states[stateOne]
            indexTwo = self.states[stateTwo]
            m = m * self.transitionMatrix[indexOne][indexTwo]
        return m/2

def main():
    '''Process the input file, create an object of the class, and create the output'''
    # pre-processing Rosalind input
    reading = sys.stdin.readlines()
    stateSeq = reading[0].strip()
    states = {}
    jojo = reading[2].strip().split('\t')
    for i, j in enumerate(jojo):
        states[j] = i
    transition = []
    for i in reading[5:5 + len(states)]:
        transition.append(i.strip().split('\t')[1:])
    transition = np.asanyarray(transition)
    transition = transition.astype(np.float)
    PP = PathProb(sequence = stateSeq, states = states , transitionMatrix = transition)
    sys.stdout.write(str(PP.pathProb()))
if __name__ == "__main__":
    main()
