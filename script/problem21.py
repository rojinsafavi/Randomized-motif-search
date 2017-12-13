import numpy as np
import sys
from problem20 import OutputFormat

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
            theta - theta value
            observation - the observed alignment (list of lists)
            observedChar - the observed characters
            pseu - pseudocount'''
        infile = self.inputFile
        readInput = infile.readlines()
        theta = float(readInput[0].split('\t')[0])
        observations = [x.strip() for x in readInput[4:]]
        observedChar = readInput[2].strip().split('\t')
        pseu = float(readInput[0].split('\t')[1].strip())
        return theta, observations, observedChar, pseu

class PseudocountHMMprofile(object):
    def __init__(self, theta, observations, observedChar, pseu):
        '''Construct a Profile HMM.
        Attributes:
        Methodes:
        '''
        self.theta = theta
        self.observations = observations
        self.observedChar = observedChar
        self.pseu = pseu

    def mpMatrix(self):
        '''Create a multiple alignment matrix.
        Input:
            a list of all alignments
        Returns:
            A matrix of multiple aligned sequence
        output example:
            [['B' 'E' 'A' 'A' 'C' '-' 'D' 'D' 'D']
            ['-' 'E' 'C' 'A' 'C' 'E' 'D' 'D' 'D']
            ['B' 'E' 'C' 'B' 'A' '-' 'B' 'D' 'D']]'''
        multipleAlignment= []
        for elemetns in self.observations:
            multipleAlignment.append(list(elemetns))
        multipleAlignment = np.asanyarray(multipleAlignment)
        self.multipleAlignment = multipleAlignment
        return multipleAlignment

    def belowTheta(self):
        '''Finding columns that the ratio of (-)/the number of sequences are above the threshold.
        Returns:
            column indecis that are above the threshold ( in list format).'''
        multipleAlignment = self.mpMatrix()
        insertColumns = []
        for i in range(multipleAlignment.shape[1]):
            x = list(multipleAlignment[:, i]).count('-')/len(multipleAlignment)
            if x >= self.theta:
                insertColumns.append(i)
        self.insertColumns = insertColumns
        return insertColumns


    def emissionBackbone(self):
        '''Constructing the emission backbone, dict of dic.
        Returns:
            the emission backbone, with dict of dict structure'''
        #here I'm creating the emission backbone
        insertColumns = self.belowTheta()
        emission = {}
        emission['S'] = {}
        emission['I0'] = {}
        for i in range(1, len(self.observations[0])-len(self.insertColumns) +1 ):
            emission['M' + str(i)] = {}
            emission['D' + str(i)] = {}
            emission['I' + str(i)] = {}
        emission['E'] = {}
        for i,j in emission.items():
            for state in self.observedChar:
                if 'I' in i or 'M' in i:
                    emission[i][state] = self.pseu
                else:
                    emission[i][state] = 0
        return emission



    def emissionProb(self):
        '''Calculating the emission probability, dict of dict.
        Returns:
            the emission probability'''
        emission = self.emissionBackbone()
        e = {}
        for j in range(self.multipleAlignment.shape[0]):
            index = 1
            for i in range(self.multipleAlignment.shape[1]):
                if i not in self.insertColumns:
                    if self.multipleAlignment[j][i] != '-':
                        if 'M' + str(index) not in e:
                            e['M' + str(index)] = [self.multipleAlignment[j][i]]
                            index += 1
                        else:
                            e['M' + str(index)].append(self.multipleAlignment[j][i])
                            index += 1
                    elif self.multipleAlignment[j][i] == '-':
                        if 'D' + str(index) not in e:
                            e['D' + str(index)] = [self.multipleAlignment[j][i]]
                            index += 1
                        else:
                            e['D' + str(index)].append(self.multipleAlignment[j][i])
                            index += 1

                else:
                    if self.multipleAlignment[j][i] != '-':
                        if 'I' + str(index -1) not in e:
                            e['I' + str(index -1)] = [self.multipleAlignment[j][i]]
                        else:
                            e['I' + str(index -1)].append(self.multipleAlignment[j][i])

        for k,v in e.items():
            for i in v:
                if i != '-':
                    emission[k][i] +=1


        for i,j in emission.items():
            if 'D' not in i and 'S' not in i and 'E' not in i:
                realSum = round(sum(list(j.values())) - len(self.observedChar)*self.pseu, 4)
                if 'I' in i or 'M' in i:
                    for k, n in j.items():
                        realCount = np.round(n - self.pseu, 4)
                        if realSum == 0:
                            emission[i][k] = float(1.0/len(self.observedChar))
                        else:
                            emission[i][k] = float((realCount/realSum)*(1-len(self.observedChar)*(self.pseu))+ self.pseu)

        return emission

    def transitionStates(self):
        '''Finding the transition states for wach sequence.
        Returns:
            for each seqeuence, it returns the a list of all the transitions.
            output example:
                if we have 3 sequence, we'll get something like this:
                [['M1', 'M2', 'M3', 'I3', 'M4', 'M5', 'M6', 'M7'],
                 ['D1', 'M2', 'M3', 'I3', 'M4', 'I4', 'M5', 'M6', 'M7'],
                  ['M1', 'M2', 'M3', 'I3', 'M4', 'M5', 'M6', 'M7']'''

        insertColumns = self.belowTheta()
        observedCharList = []
        for j in range(self.multipleAlignment.shape[0]):
            temp = []
            index = 1
            for i in range(self.multipleAlignment.shape[1]):
                if i in insertColumns:
                    if self.multipleAlignment[j][i] != '-':
                        temp.append('I' + str(index -1))
                else:
                    if self.multipleAlignment[j][i] != '-':
                        temp.append('M' + str(index))
                        index += 1
                    elif self.multipleAlignment[j][i] == '-':
                        temp.append('D' + str(index))
                        index += 1
            observedCharList.append(temp)
        self.observedCharList = observedCharList
        return observedCharList

    def transitionBackbone(self):
        '''Construction the transition backbone.
        Returns:
            an empty transition dict of dict'''
        observedCharList = self.transitionStates()
        transition = {}
        transition['S'] = {}
        transition['I0'] = {}
        for i in range(1, len(self.observations[0])-len(self.insertColumns) +1 ):
            transition['M' + str(i)] = {}
            transition['D' + str(i)] = {}
            transition['I' + str(i)] = {}
        transition['E'] = {}

        for k,v in transition.items():
            for kk in transition.keys():
                transition[k][kk] = 0

        for k in list(transition.keys()):
            if k == 'S':
                transition['S']['I0'] = self.pseu
                transition['S']['M1'] = self.pseu
                transition['S']['D1'] = self.pseu
            elif k == 'E':
                transition[list(transition.keys())[-2]][k] = self.pseu
                transition[list(transition.keys())[-3]][k]= self.pseu
                transition[list(transition.keys())[-4]][k]= self.pseu
            else:
                if "I" + str(int(k[1:])) in transition[k]:
                    transition[k]["I" + str(int(k[1:]))] = self.pseu
                if "M" + str(int(k[1:]) + 1) in transition[k]:
                    transition[k]["M" + str(int(k[1:]) + 1)] = self.pseu
                if "D" + str(int(k[1:]) + 1) in transition[k]:
                    transition[k]["D" + str(int(k[1:]) + 1)] = self.pseu
        self.transition = transition
        return transition

    def transitionProb(self):
        '''Create the transition probability, in dict of dict structure.
        Returns:
            the transition probability'''
        transition = self.transitionBackbone()
        for row in self.observedCharList:
            prev = ''
            for index,item in enumerate(row):
                if index == 0:
                    if item != "-":
                        transition['S'][item] += 1
                        prev = item
                    else:
                        transition['I0'][item] += 1
                        prev = item
                else:
                    self.transition[prev][item] += 1
                    prev = item
                    if index == len(row) - 1:
                        self.transition[prev]['E'] += 1
        for k,v in transition.items():
            numNonZeroz = np.count_nonzero(list(v.values())) #counting the number of nonzero elements
            realSum = sum(transition[k].values()) - numNonZeroz*(self.pseu)
            if realSum > 0:
                for j,ll in v.items():
                    if ll != 0 :
                        realCount = ll -self.pseu
                        if realCount != 0:
                            transition[k][j] = float((realCount/realSum)*(1-numNonZeroz*(self.pseu))+ self.pseu)

            if realSum == 0:
                for j,ll in v.items():
                    if ll != 0:
                        transition[k][j] = float(1.0/numNonZeroz)

        return transition

def main():
    fr = RosalindParse(sys.stdin)
    theta, observations, observedChar, pseu = fr.rosalidParse()
    h = PseudocountHMMprofile(theta, observations, observedChar, pseu)
    outputInit = OutputFormat(transition = h.transitionProb(), emission = h.emissionProb(), observedChar = observedChar)
    outputInit.outputParsing()
if __name__ == "__main__":
    main()
