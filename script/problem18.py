import sys
import numpy as np

class ViterbiPath(object):
    def __init__(self, observation, transition, emission, states, observedChar ):
        '''Finds the path that maximizes the  probability Pr(x, π) over all possible paths π.
        Attributes:
            observation - a string of observations
            observedChar - alphabets in the observation starting
            states - hidden states
            emission - emission matrix
            transition - transition matrix
        Methods:
            viterbiPath - a path that maximizes the probability Pr(x, π) over all possible paths'''

        self.observation = observation
        self.transition = transition
        self.emission = emission
        self.observedChar = observedChar
        self.states = states

    def viterbiPath(self):

        '''Returns a path that maximizes the probability Pr(x, π) over all possible paths.
        returns:
            A path of states'''
        probPrev = [{state : {"currentProb" : (1.0/len(self.states)) * self.emission[state][self.observation[0]], "prevState" : None} for state in self.states}] # starting a dict with prev state of None and current robability of neing in state A or B as currentProb
        for characterIndex in range(1, len(self.observation)): # for each position in the observed, first I calculated the max probability of the path that ends in that position, then I multiply that by the emission
            probPrev.append({}) #for each position append a new dict to the probPrev
            for state in self.states: # I0, M1,I1, D2 ...
                maxProb = 0.0
                prevState  = ''
                for prev_state in self.states:
                    prob, prev_st = (probPrev[characterIndex-1][prev_state]["currentProb"] * self.transition[prev_state][state] * self.emission[state][self.observation[characterIndex]], prev_state)
                    if prob > maxProb:
                        maxProb = prob
                        prevState = prev_st
                probPrev[characterIndex][state] = {"currentProb" : maxProb, "prevState" : prevState}

        #backtracking
        stateProb = float('-inf')
        L = ''
        statePath = ''
        # first we only look at the end point, and pick the state with the highest probability for the end, and based on that state we record the prev state
        for state in self.states: # from the end, I pick the state with the highest probability
            endProb = probPrev[-1].get(state, {}).get('currentProb') # get the probability
            prevState = probPrev[-1].get(state, {}).get('prevState') # get the prevState
            if endProb >= stateProb: # only if the the probability is higher greated than other states
                stateProb = endProb
                L = prevState
                currentState = state
        statePath =  statePath + currentState + L
        # now we are looking at everything before the end state
        for j in range(-2, -len(self.observation), -1): # we start from -2 since the last node (-1) is resolved
            K =  probPrev[j].get(L, {}).get('prevState') #get the prev
            L = K
            statePath += L
        return statePath[::-1]

def main():
    '''Process the input file, create an object of the class, and create the output'''
    # pre-processing Rosalind input
    reading = sys.stdin.readlines()
    observation = reading[0].strip()
    observedChar =reading[2].strip().split('\t')
    states = reading[4].strip().split('\t')
    transition = []
    for i in reading[7:7 + len(states)]:
        transition.append(i.strip().split('\t')[1:])
    transition = np.asanyarray(transition)
    transition = transition.astype(np.float)
    emission = []
    for i in reading[9 + len(states):]:
        emission.append(i.strip().split('\t')[1:])
    emission = np.asanyarray(emission)
    emission = emission.astype(np.float)


    initTransition = {}
    for i,j in enumerate(states):
        for k, n in  enumerate(states):
            if j not in initTransition:
                initTransition[j] = {}
                initTransition[j][n] = transition[i][k]
            else:
                initTransition[j][n] = transition[i][k]

    initEmission= {}
    for i,j in enumerate(states):
        for k, n in enumerate(observedChar):
            if j not in initEmission:
                initEmission[j] = {}
                initEmission[j][n] = emission[i][k]
            else:
                initEmission[j][n] = emission[i][k]

    V = ViterbiPath(observation = observation, transition =initTransition, emission = initEmission, states = states  , observedChar =  observedChar)
    sys.stdout.write(V.viterbiPath() + '\n')
if __name__ == "__main__":
    main()
