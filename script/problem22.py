import sys
import numpy as np
from problem21 import PseudocountHMMprofile

class SequenceAlignment(object):
    ''''Sequence Alignment with Profile HMM Problem.
    Attributes:
        string - observed Sequence
        theta - theta
        pseu - pseudocount
        alphabet - observed alphabet
        transition - transition
        emission - emission
    Methodes:
        zeroColumns - addes zero states to the transition
        seperateStates - Seperates states
        bounderies - adds the boundries
        viterbiGraph - creates the graph
        backtracking - backtrack to recover the path
    '''
    def __init__(self, string,  theta, pseu, alphabet, transition, emission):
        self.string = string
        self.theta = theta
        self.pseu = pseu
        self.alphabet = alphabet
        self.transition = transition
        self.emission = emission

    def zeroColumns(self):
        '''Adding D0 and M0 to the transition
        Returns:
            transition'''
        self.transition['D0']= {}
        self.transition['M0']= {}
        for i in list(self.transition.keys()):
            self.transition['D0'][i] = 0
            self.transition['M0'][i] = 0
        return self.transition

    def seperateStates(self):
        '''Creating 6 matricis, 3 for backtracking, 3 for graph
        Returns:
            6 matricis, 3 for backtracking, 3 for graph'''
        self.zeroColumns()
        self.rowsNum = len(self.emission)//3
        self.numColumns = len(self.string) + 1

        self.matrixM = np.zeros([self.rowsNum, self.numColumns])
        self.matrixI = np.zeros([self.rowsNum, self.numColumns])
        self.matrixD = np.zeros([self.rowsNum, self.numColumns])

        self.DD = np.zeros([self.rowsNum, self.numColumns], dtype= object) # for backtracking
        self.II = np.zeros([self.rowsNum, self.numColumns], dtype= object) # for backtracking
        self.MM = np.zeros([self.rowsNum, self.numColumns], dtype= object) # for backtracking

        return self.matrixM, self.matrixI, self.matrixD, self.MM, self.II, self.DD

    def bounderies(self):
        '''constructing the graph bounderies
        Returns:
         bounderies of 6 matrices, 3 for backtracking, 3 for graph'''
        self.seperateStates()
        for i in range(1, self.numColumns):
            if i == 1:
                self.matrixI[0][1] = self.transition['S']['I0']*self.emission['I0'][self.string[0]]
                self.II[0][1] = ('S', self.matrixI[0][1])# for backtracking
            else:
                self.matrixI[0][i] = self.matrixI[0][i -1]* self.transition['I0']['I0']*self.emission['I0'][self.string[i-1]]
                self.II[0][1] = ('I0', self.matrixI[0][i])# for backtracking

        for row in range(1, self.rowsNum):
            if row == 1:
                self.matrixD[1][0] = self.transition['S']['D1']
                self.DD[1][0] = ('S', self.matrixD[1][0])# for backtracking
            else:
                self.matrixD[row][0] = self.matrixD[row -1 ][0]* self.transition['D' + str(row-1)]['D' + str(row)]
                self.DD[row][0] = ('D' + str(row-1), self.matrixD[row][0]) #for backtracking
        return self.matrixM, self.matrixI, self.matrixD, self.MM, self.II, self.DD

    def viterbiGraph(self):
        '''constructing the whole graph, and the backtracking matrix
        Returns:
            6 matrices, 3 for backtracking, 3 for graph'''
        self.bounderies()
        for row in range(1, self.rowsNum):
            for column in range(1, self.numColumns):
                if row == 1 and column == 1:
                    self.matrixD[1][1] =self.matrixI[0][1]*self.transition['I0']['D1']
                    self.matrixM[1][1] = self.transition['S']['M1']*self.emission['M1'][self.string[0]]
                    self.matrixI[1][1] = self.matrixD[1][0]*self.transition['D1']['I1']*self.emission['I1'][self.string[0]]

                    self.DD[1][1] =('I0',self.matrixD[1][1])#for backtracking
                    self.MM[1][1] =('S',self.MM[1][1])#for backtracking
                    self.II[1][1] = ('D1',self.II[1][1])#for backtracking

                else:
                    l = (self.matrixD[row-1][column]*self.transition['D' +  str(row -1)]['D' + str(row)], # same column row before
                         self.matrixM[row-1][column]*self.transition['M' +  str(row -1)]['D' + str(row)],
                         self.matrixI[row-1][column]*self.transition['I' +  str(row -1)]['D' + str(row)])
                    n = ['D' +  str(row -1),'M' +  str(row -1), 'I' +  str(row -1)]
                    p = np.argmax(l)
                    q = n[p]
                    self.matrixD[row][column] = l[p]
                    self.DD[row][column]= (q, l[p])#for backtracking


                    l = (self.matrixD[row-1][column -1]*self.transition['D' +  str(row -1)]['M' + str(row)]* self.emission['M' + str(row)][self.string[column -1]], # same column row before
                         self.matrixM[row-1][column -1]*self.transition['M' +  str(row -1)]['M' + str(row)]* self.emission['M' + str(row)][self.string[column -1]],
                         self.matrixI[row-1][column -1]*self.transition['I' +  str(row -1)]['M' + str(row)]* self.emission['M' + str(row)][self.string[column -1]])
                    n = ['D' +  str(row -1), 'M' +  str(row -1), 'I' +  str(row -1)]
                    p = np.argmax(l)
                    q = n[p]
                    self.matrixM[row][column] = l[p]
                    self.MM[row][column] = (q, l[p])#for backtracking


                    l = (self.matrixD[row][column -1]*self.transition['D' +  str(row)]['I' + str(row)]* self.emission['I' + str(row)][self.string[column -1]], # same column row before
                         self.matrixM[row][column -1]*self.transition['M' +  str(row)]['I' + str(row)]* self.emission['I' + str(row)][self.string[column -1]],
                         self.matrixI[row][column -1]*self.transition['I' +  str(row)]['I' + str(row)]* self.emission['I' + str(row)][self.string[column -1]])
                    n = ['D' +  str(row), 'M' +  str(row), 'I' +  str(row)]
                    p = np.argmax(l)
                    q = n[p]
                    self.matrixI[row][column] = l[p]
                    self.II[row][column] = (q, l[p])#for backtracking

        return self.matrixD,self.matrixM, self.matrixI, self.DD, self.MM,self.II

    def backtracking(self):
        '''An optimal hidden path emitting Text in HMM(Alignment,θ,σ).
        Returns:
            hiddem path'''
        self.viterbiGraph()
        orderedList = []
        # first I find the state before the end state
        l = (self.matrixD[self.rowsNum -1][self.numColumns -1]*self.transition['D' +  str(self.rowsNum -1)]['E'],
             self.matrixM[self.rowsNum -1][self.numColumns -1]*self.transition['M' +  str(self.rowsNum -1)]['E'],
             self.matrixI[self.rowsNum -1][self.numColumns -1]*self.transition['I' +  str(self.rowsNum -1)]['E'])
        n = ['D' +  str(self.rowsNum -1), 'M' +  str(self.rowsNum -1), 'I' +  str(self.rowsNum -1)]
        p = np.argmax(l)
        q = n[p]
        lastProb = l[p]
        prevState = q

        orderedList.append(prevState)
        #prevProb = lastProb

        column = len(self.string)
        row = self.rowsNum -1

        # backtracking
        while column != 0:
            if 'M' in prevState:
                if type(self.MM[row][column]) is not tuple:
                    break
                prevState = self.MM[row][column][0]
                if prevState == 'S':
                    break
                row =int(prevState[1:])
                column -=1

            elif 'D' in prevState:
                if type(self.DD[row][column]) is not tuple:
                    break
                prevState = self.DD[row][column][0]
                if prevState == 'S':
                    break
                row =int(prevState[1:])

            elif 'I' in prevState:
                if type(self.II[row][column]) is not tuple:
                    break
                prevState = self.II[row][column][0]
                if prevState == 'S':
                    break
                row =int(prevState[1:])
                column -=1
            orderedList.append(prevState)

        return orderedList

class RosalindParse(object):

    def __init__(self, inputFile):
        self.inputFile = inputFile

    def rosalidParse(self):
        infile = self.inputFile
        readInput = infile.readlines()
        string = readInput[0].strip()
        theta = float(readInput[2].split()[0])
        pseu = float(readInput[2].split()[1].strip())
        observations = [x.strip() for x in readInput[6:]]
        observedChar = readInput[4].strip().split()
        return string, theta, observations, observedChar, pseu

def main():
    fr = RosalindParse(sys.stdin)
    string, theta, observations, observedChar, pseu= fr.rosalidParse()
    h = PseudocountHMMprofile(theta = theta, observations = observations, observedChar = observedChar, pseu = pseu)
    transition = h.transitionProb()
    emission = h.emissionProb()
    j = SequenceAlignment(string= string, pseu= pseu, theta= theta, alphabet =  observedChar, transition= transition, emission= emission)
    b = j.backtracking()
    sys.stdout.write(' '.join(map(str, b[::-1])) + '\n')

if __name__ == "__main__":
    main()
