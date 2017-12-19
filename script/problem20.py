import numpy as np
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
            theta - theta value
            observation - the observed alignment (list of lists)
            observedChar - the observed characters
            '''
        infile = self.inputFile
        readInput = infile.readlines()
        theta = float(readInput[0].strip()) #extracting theta
        observations = [x.strip() for x in readInput[4:]] # the observed sequence
        observedChar = readInput[2].strip().split()# the observed characters
        return theta, observations, observedChar

class HmmProfile(object):
    '''Construct a Profile HMM.
    Attributes:
        theta - theta value
        observations - the observed sequence
        observedChar - the observed characters
    Methodes:
        mpMatrix - Creates a matrix of the multiple alignment
        belowTheta - Finding the columns that are to be treated as inserts
        transitionStates - Finds the transition states for each sequence
        transitionBackbone - Creates the transition backbone
        transitionProb - Creates the transition probability  matrix
        emissionBackbone - creates the emission backbone
        emissionProb - Creates the emission probability matrix
    '''
    def __init__(self, theta, observations, observedChar):
        self.theta = theta
        self.observations = observations
        self.observedChar = observedChar

    def mpMatrix(self):
        '''Creates a multiple alignment matrix.
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
            column indices that are above the threshold ( in list format).'''
        multipleAlignment = self.mpMatrix()# get the alignment matrix
        insertColumns = [] # the  indices of columns that should be considered as insertion
        for i in range(multipleAlignment.shape[1]):
            x = list(multipleAlignment[:, i]).count('-')/len(multipleAlignment)
            if x >= self.theta:
                insertColumns.append(i)
        self.insertColumns = insertColumns
        return insertColumns

    def transitionStates(self):
        '''Finding the transition states for each sequence.
        Returns:
            for each seqeuence, it returns the a list of all the transitions.
            output example:
                if we have 3 sequence, we'll get something like this:
                [['M1', 'M2', 'M3', 'I3', 'M4', 'M5', 'M6', 'M7'],
                 ['D1', 'M2', 'M3', 'I3', 'M4', 'I4', 'M5', 'M6', 'M7'],
                  ['M1', 'M2', 'M3', 'I3', 'M4', 'M5', 'M6', 'M7']'''
        insertColumns = self.belowTheta() # a list of column indices that must be considered as insertion
        stateList = [] # a list containing the state observations for each sequence
        for alignment in range(self.multipleAlignment.shape[0]):
            temp = []
            index = 1
            for columnIndex in range(self.multipleAlignment.shape[1]):
                if columnIndex in insertColumns: # if the column is in the insertColumn list
                    if self.multipleAlignment[alignment][columnIndex] != '-': # if we have '-'
                        temp.append('I' + str(index -1)) # we are going to get insetion
                else: # if it is not in the inserColumn list
                    if self.multipleAlignment[alignment][columnIndex] != '-': # if we see '-'
                        temp.append('M' + str(index)) # we are going to have match
                        index += 1
                    elif self.multipleAlignment[alignment][columnIndex] == '-': # if we don't have '-'
                        temp.append('D' + str(index)) # we will have deletion
                        index += 1
            stateList.append(temp)
        self.stateList = stateList
        return stateList

    def transitionBackbone(self):
        '''Construction the transition backbone.
        Returns:
            an empty transition dict, with structure of dict of dict'''
        stateList = self.transitionStates() # get the list of transition states ( list of lists, each sub list coresponding to one of the sequences in order)
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
                transition[k][kk] = 0 # assigning 0 to dic of dict values
        self.transition = transition
        return transition

    def transitionProb(self):
        '''Create the transition probability, in dict of dict structure.
        Returns:
            the transition probability'''
        transition = self.transitionBackbone() #get the backbone
        for row in self.stateList: # self.stateList is explained in transitionStates docstring
            prev = '' # here I'm just counting
            for index, item in enumerate(row):
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
        #here I'm creating the probability
        for k,v in transition.items():
            s = sum(transition[k].values()) # sum of the values
            if s !=0: # if the sum is not zero ( to avoid devision by zero)
                for j,ll in v.items(): #for each key and value in v
                    if ll != 0: # if the value is not zero , because the rosalind does not want to see 0.0
                            transition[k][j] = float(ll/s) # calculating the probability
        return transition

    def emissionBackbone(self):
        '''Constructing the emission backbone, dict of dic.
        Returns:
            the emission backbone, with dict of dict structure'''
        insertColumns = self.belowTheta() # a list of column indices that must be considered as insertion
        emission = {}
        emission['S'] = {}
        emission['I0'] = {}
        for i in range(1, len(self.observations[0])-len(self.insertColumns) +1 ):
            emission['M' + str(i)] = {}
            emission['D' + str(i)] = {}
            emission['I' + str(i)] = {}
        emission['E'] = {}
        #here I'm just assigning 0
        for i,j in emission.items():
            for state in self.observedChar:
                emission[i][state] = 0
        return emission

    def emissionProb(self):
        '''Calculating the emission probability, with structure of dict of dict.
        Returns:
            the emission probability'''
        emission = self.emissionBackbone() # getting the emission backbone
        #here I'm just counting
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

        #making the probability
        for wq,er in emission.items():
            ss = sum(er.values())
            if ss != 0:
                for j,ll in er.items():
                    if ll != 0:
                        emission[wq][j] = float(ll/ss)
        return emission

    def columnLabel(self):
        '''Constructing the emission backbone, dict of dic.
        Returns:
            the emission backbone, with dict of dict structure'''
        insertColumns = self.belowTheta() # a list of column indices that must be considered as insertion
        x = []
        x.append('S')
        x.append('I0')
        for i in range(1, len(self.observations[0])-len(self.insertColumns) +1 ):
            x.append('M' + str(i))
            x.append('D' + str(i))
            x.append('I' + str(i))
        x.append('E')
        return x

class OutputFormat(object):
    '''Parsing rosaling output format.
    Attributes:
        transition - tranistion
        emission - emission
        observedChar - observed characters
    Methods:
        outputParsing - STDOUT the output
     '''
    def __init__(self, transition, emission,observedChar, columns):
        self.transition = transition
        self.emission = emission
        self.observedChar = observedChar
        self.columns = columns

    def outputParsing(self):
        '''STDOUT the output'''
        sys.stdout.write("\t" + '\t'.join(map(str, self.columns)) + "\n")
        for key in self.columns:
            sys.stdout.write(key)
            for insideKey in self.columns:
                element = self.transition[key][insideKey]
                if element == 0:
                    sys.stdout.write('\t' + str(element))
                else:
                    sys.stdout.write('\t' + str(float('%.3g'% element)))
            sys.stdout.write('\n')

        sys.stdout.write('--------' + '\n')

        sys.stdout.write("\t" + '\t'.join(map(str, self.observedChar)) + "\n")

        for key in self.columns:
            sys.stdout.write(key)
            for insideKey in self.observedChar:
                element = self.emission[key][insideKey]
                if element == 0:
                    sys.stdout.write('\t' + str(element))
                else:
                    sys.stdout.write('\t' + str(float('%.3g'% element)))
            sys.stdout.write('\n')


def main():
    '''Intitiation HmmProfile class'''
    fr = RosalindParse(sys.stdin)
    theta, observations, observedChar = fr.rosalidParse()
    h = HmmProfile(theta, observations, observedChar)
    outputInit = OutputFormat(transition = h.transitionProb(), emission = h.emissionProb(), observedChar = observedChar, columns= h.columnLabel())
    outputInit.outputParsing()
if __name__ == "__main__":
    main()
