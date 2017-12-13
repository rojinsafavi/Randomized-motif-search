import sys
import numpy as np

class LikelihoodProblem(object):
    def __init__(self, observation, transition, emission, states, observedChar ):
        '''Finds the probability Pr(x) that the HMM emits x.
        Attributes:
            observation - a string of observations
            observedChar - alphabets in the observation starting
            states - hidden states
            emission - emission matrix
            transition - transition matrix
        Methods:
            likelihoodProblem - returns the probability Pr(x) that the HMM emits x.'''
        self.observation = observation
        self.transition = transition
        self.emission = emission
        self.observedChar = observedChar
        self.states = states

    def likelihoodProblem(self):
        '''Returns the probability Pr(x) that the HMM emits x.
        Returns:
            the probability'''
        probPrev = [{state : {"currentProb" : (1.0 / len(self.states)) * np.float(self.emission[self.states[state]][self.observedChar[self.observation[0]]])} for state in self.states}]
        for characterIndex in range(1, len(self.observation)): # for each position in the observed, first I calculated the max probability of the path that ends in that position, then I multiply that by the emission
            probPrev.append({}) #for each position append a new dict to the probPrev
            for state in self.states: # A,B,C ...
                SUM = 0
                for prev_state in self.states:
                    prob = (np.float(probPrev[characterIndex-1][prev_state]["currentProb"]) * np.float(self.transition[self.states[prev_state]][self.states[state]]) * np.float(self.emission[self.states[state]][self.observedChar[self.observation[characterIndex]]]))
                    SUM += prob
                probPrev[characterIndex][state] = {"currentProb" : SUM}
        probSum = 0
        for state in self.states:
            probSum += probPrev[-1].get(state, {}).get('currentProb')

        return  probSum




def main():
    '''Process the input file, create an object of the class, and create the output'''
    # pre-processing Rosalind input
    reading = sys.stdin.readlines()
    observation = reading[0].strip()
    observedChar = {}
    for i, j in enumerate(reading[2].strip().split('\t')):
        observedChar[j] = i
    states = {}
    jojo = reading[4].strip().split('\t')
    for i, j in enumerate(jojo):
        states[j] = i
    transition = []
    for i in reading[7:7 + len(states)]:
        transition.append(i.strip().split('\t')[1:])
    transition = np.asanyarray(transition)
    transition = transition.astype(np.float)
    emission = []
    for i in reading[9 + len(states):9 + 2*len(states)]:
        emission.append(i.strip().split('\t')[1:])
    emission = np.asanyarray(emission)
    emission = emission.astype(np.float)

    V = LikelihoodProblem(observation = observation, transition =transition, emission = emission, states = states  , observedChar =  observedChar)
    sys.stdout.write(str(V.likelihoodProblem()) + '\n')
if __name__ == "__main__":
    main()
