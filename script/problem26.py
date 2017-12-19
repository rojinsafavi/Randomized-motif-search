import sys
import numpy as np
from problem25 import SoftDecoding


class BaumWelch(object):
    ''' Learning the emission and transition using BaumWelch algorithem.
    Attributes:
        states - states
        forwardMatrix - forwardMatrix
        backwardMatrix - backwardMatrix
        observation - observed sequence
        transition - transition
        emission - emission
        nodeRes - node responsiblity metrix
        observedChar - observed characters
    Methodes:
        edgeWeight - finds the edge weights
        forwardSink - find the forward sink values
        edgeResponsibility - finds the edge responsiblity matrix
        redefineParameters - redefins transition and emission parameters
        find - finds the index of a character in a string
    '''

    def __init__(self,forwardMatrix, backwardMatrix, states, observation, transition, emission, nodeRes, observedChar):
        self.states = states
        self.forwardMatrix = forwardMatrix
        self.backwardMatrix = backwardMatrix
        self.observation = observation
        self.transition = transition
        self.emission = emission
        self.nodeRes = nodeRes
        self.observedChar = observedChar

    def edgeWeight(self):
        '''finds the edge weights.
        Returns:
            A list of dictionaries each having the edge probabilities'''
        edgeProb = [{}]
        for characterIndex in range(1, len(self.observation)):
            edgeProb.append({})
            for state in self.states:
                for nextState in self.states:
                    edgeProb[characterIndex][state + nextState] = self.transition[state][nextState] * self.emission[nextState][self.observation[characterIndex]]
        return edgeProb

    def forwardSink(self):
        '''Finding the forwardSink values.
        Returns :
            a list containg the forward sink values'''
        forward, backward = self.forwardMatrix, self.backwardMatrix
        t = []
        for i in range(len(forward)):
            total = 0
            for state in self.states:
                total += forward[i][state]['currentProb']*backward[i][state]['currentProb']
            t.append(total)
        return t

    def edgeResponsibility(self):
        '''Finding the edge responsiblity values.
        Returns:
            edge responsiblity values'''
        fSink = self.forwardSink() # the forward sink value
        edgeRes = [] # edge responsibility list
        edgeProb = self.edgeWeight()[1:] # get the esdge probability values
        for i in range(len(edgeProb)):
            tot= {}
            for edgeName, edgeW in edgeProb[i].items():
                state, nextState = edgeName[0], edgeName[1]
                mult = self.forwardMatrix[i][state]['currentProb'] * self.backwardMatrix[i+1][nextState]['currentProb'] * edgeW/fSink[i]
                tot[state + nextState] = mult
            edgeRes.append(tot)
        return edgeRes

    def find(self, x, ch):
        '''Finds a list of all indexes that a character has occures ina string.
        Returns:
            a list with all the indexes that a character has occured in a string'''
        return [i for i, j in enumerate(x) if j == ch]

    def redefineParameters(self):
        '''Redefining parameters (transition and emission).
        Returns:
            the new emission and transition dicts'''
        edgeRes = self.edgeResponsibility() # get the edge responsibility values
        newTransition = {state : {} for state in self.states} # the new transitin , dict of dict structure
        for edgeValue in edgeRes: # for each edge value in the edgeRes
            for k,v in edgeValue.items(): # for key and value in the edge values
                start = k[0] # the start of the edge
                end = k[1] # the end of the edge
                if end not in newTransition[start]:
                    newTransition[start][end] = v
                else:
                    newTransition[start][end] += v # counting
        # making the probability
        for k,v in newTransition.items():
            s = sum(list(v.values()))
            for i,j in v.items():
                if s != 0:
                    newTransition[k][i] = j/s


        newEmission = {state : {} for state in self.states} # the new emission , dict of dict structure

        for key in newEmission.keys(): # assigning 0 to the new emission states
            newEmission[key] = {state : 0 for state in self.observedChar}

        for i in self.observedChar:
            indices = self.find(self.observation, i) # gives the index

            for n, m in enumerate(self.states):
                newEmission[m][i] += sum(np.asanyarray(self.nodeRes)[indices][:,n])

       # make the probability
        for k,v in newEmission.items():
            s = sum(list(v.values()))
            for i,j in v.items():
                newEmission[k][i] = j/s


        return newTransition, newEmission



def main():
    '''Process the input file, create an object of the class, and create the output'''
    # pre-processing Rosalind input
    reading = sys.stdin.readlines()
    itteration = reading[0].strip()
    observation = reading[2].strip()
    observedChar =reading[4].strip().split()
    states = reading[6].strip().split()

    transition = []
    for i in reading[9:9 + len(states)]:
        transition.append(i.strip().split()[1:])
    if transition[-1] == []: # there is a space in your input file, which I pop it out here
        transition.pop # there is a space in your input file, which I pop it out here
    transition = np.asanyarray(transition)
    transition = transition.astype(np.float)

    emission = []
    for i in reading[11 + len(states):11 + 2*len(states)]:
        emission.append(i.strip().split()[1:])
    if emission[-1] == []: # there is a space in your input file, which I pop it out here
        emission.pop # there is a space in your input file, which I pop it out here
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


    for i in range(int(itteration)):
        V = SoftDecoding(observation = observation, transition =initTransition, emission = initEmission, states = states  , observedChar =  observedChar)
        forward, backward= V.forwardBackward()
        nodeRes = V.softDecoding()
        er = BaumWelch(forwardMatrix = forward, backwardMatrix = backward, states = states, observation= observation, transition = initTransition, emission = initEmission, nodeRes = nodeRes, observedChar = observedChar)
        initTransition, initEmission = er.redefineParameters()

    sys.stdout.write('\t' + '\t'.join(map(str, states)) + "\n")
    for key in states:
        sys.stdout.write(key)
        for insideKey in states:
            value = initTransition[key][insideKey]
            sys.stdout.write( "\t" + '%.3f'%value)
        sys.stdout.write( "\n")
    sys.stdout.write('--------' + '\n')

    sys.stdout.write('\t' + '\t'.join(map(str, observedChar)) + "\n")
    for key in states:
        sys.stdout.write(key)
        for insideKey in observedChar:
            value = initEmission[key][insideKey]
            sys.stdout.write( "\t" + '%.3f'%value)
        sys.stdout.write( "\n")




if __name__ == "__main__":
    main()
