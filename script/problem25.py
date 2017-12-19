import sys
import numpy as np

class SoftDecoding(object):
    def __init__(self, observation, transition, emission, states, observedChar ):
        '''Finds the probability Pr(x) that the HMM emits x.
        Attributes:
            observation - a string of observations
            observedChar - alphabets in the observation starting
            states - hidden states
            emission - emission matrix
            transition - transition matrix
        Methods:
            forwardBackward - returns the probability Pr(x) that the HMM emits x.
            softDecoding - The outcome likelihood for pr(x)'''
        self.observation = observation
        self.transition = transition
        self.emission = emission
        self.observedChar = observedChar
        self.states = states

    def forwardBackward(self):
        '''Returns the probability Pr(x) that the HMM emits x.
        Returns:
            forward and backward matrices'''
        forwardProb = [{state : {"currentProb" : (1.0 / len(self.states)) * self.emission[state][self.observation[0]]} for state in self.states}]
        for characterIndex in range(1, len(self.observation)):
            forwardProb.append({}) #for each position append a new dict to the forwardProb
            for state in self.states:
                SUM = 0
                for prev_state in self.states:
                    prob = (forwardProb[characterIndex-1][prev_state]["currentProb"] * self.transition[prev_state][state] * self.emission[state][self.observation[characterIndex]])
                    SUM += prob
                forwardProb[characterIndex][state] = {"currentProb" : SUM}


        backwardProb = [{state : {"currentProb" : 1.0} for state in self.states}]
        for characterIndex in range(len(self.observation)-2, -1,-1):
            backwardProb.append({}) #for each position append a new dict to the backwardProb
            for state in self.states:
                SUM = 0
                for prev_state in self.states:
                    prob = (backwardProb[-2][prev_state]["currentProb"] * self.transition[state][prev_state] * self.emission[prev_state][self.observation[characterIndex + 1]])
                    SUM += prob
                backwardProb[-1][state] = {"currentProb" : SUM}
        return forwardProb, backwardProb[::-1]


    def softDecoding(self):
        forward, backward = self.forwardBackward()
        t = []
        for i, j in enumerate(forward):
            h = backward[i]
            total = 0
            B = []
            for state in self.states:
                total += j[state]['currentProb']*h[state]['currentProb']
            for state in self.states:
                A = j[state]['currentProb']*h[state]['currentProb']/total
                B.append(A)
            t.append(B)
        return t

def main():
    '''Process the input file, create an object of the class, and create the output'''
    # pre-processing Rosalind input
    reading = sys.stdin.readlines()

    observation = reading[0].strip()
    observedChar =reading[2].strip().split()
    states = reading[4].strip().split()

    transition = []
    for i in reading[7:7 + len(states)]:
        transition.append(i.strip().split()[1:])
    transition = np.asanyarray(transition)
    transition = transition.astype(np.float)
    emission = []
    for i in reading[9 + len(states):9 + 2*len(states)]:
        emission.append(i.strip().split()[1:])
    emission = np.asanyarray(emission)
    emission = emission.astype(np.float)


    initTransition = {} # here I'm converting the matrix to dict of dict
    for i,j in enumerate(states):
        for k, n in  enumerate(states):
            if j not in initTransition:
                initTransition[j] = {}
                initTransition[j][n] = transition[i][k]
            else:
                initTransition[j][n] = transition[i][k]

    initEmission= {} # here I'm converting the matrix to dict of dict
    for i,j in enumerate(states):
        for k, n in enumerate(observedChar):
            if j not in initEmission:
                initEmission[j] = {}
                initEmission[j][n] = emission[i][k]
            else:
                initEmission[j][n] = emission[i][k]

    V = SoftDecoding(observation = observation, transition =initTransition, emission = initEmission, states = states  , observedChar =  observedChar)
    sd = V.softDecoding()
    sys.stdout.write('\t'.join(map(str, states)) + "\n")
    for value in sd:
        sys.stdout.write('\t'.join(map(str, [round(i,4) for i in value])) + '\n')

if __name__ == "__main__":
    main()
