#Q5
import numpy as np
import sys
from collections import defaultdict

class BruijnGraph(object):
    '''Construct the de Bruijn graph from a list of kmers.
    Attributes:
        edgeSize - the edge length
        inputSequence - the sequence
    Methodes:
        nodesDic - dict representation of the graph
        output - Rosalind output format'''

    def __init__(self, inputSequence, edgeSize):
        self.edgeSize = edgeSize
        self.inputSequence = inputSequence

    def nodesDic(self):
        '''Returns a dict of nodes.
        Returns:
            a graph represetation of the nodes in a dict format'''
        nodes = defaultdict(list)
        for kmer in self.inputSequence:
            start = kmer[:-1]
            end = kmer[1:]
            if start not in nodes:
                nodes[start] = [end]
            else:
                nodes[start].append(end)
        return nodes


class EulerianPath(object):
    '''Finds an Eulerian path in a graph.
    Attributes:
        nodesDic - a directed graph in dict format
        startNode - the start node if it exists
        lastNode - last node of it exists
    Methodes:
        edgesFreq - the difference between the input and output edges from all nodes.
        endNodes - finds the start and end nodes.
        edgesNum - finds the total number of edges
        count - counts the number of outgoing edges given a node
        eulerianPath - finds an Eulerian path
        outputFormat - Rosalind output format of the eulerianPath
        connectEdges - connecting the end nodes with an edge
    '''
    def __init__(self, nodesDic):
        self.nodesDic = nodesDic
        self.startNode = []
        self.lastNode = []

    def edgesFreq(self):
        '''Generates frequency dicts.
        Returns:
            two dicts'''
        L = list()
        inFreq = {} #a dictionary containg the frequency of input edges for all node
        outFreq = defaultdict(int) #a dictionary containg the frequency of output edges for all node
        for key,value in self.nodesDic.items():
            valueLen = len(value)
            inFreq[key] = valueLen
            L = L + value
        for i in L:
            outFreq[i] += 1
        inOut = {key: inFreq[key] - outFreq.get(key, 0) for key in inFreq.keys()} #difference betweein the number of outgouing edges and in goin edges
        outIn = {key: outFreq[key] - inFreq.get(key, 0) for key in outFreq.keys()} #difference betweein the number of ingoing edges and in outgoing edges
        return inOut, outIn

    def endNodes(self):
        '''Find the end nodes.
        Returns:
            end nodes'''
        ends = self.edgesFreq()
        s = ends[0]
        e = ends[1]
        for key,value in s.items():
            if value == 1:
                self.startNode = key #the start node
        startNode = self.startNode
        for key,value in e.items():
            if value == 1:
                self.lastNode = key #the start node
        lastNode= self.lastNode
        return startNode, lastNode

    def connectEdges(self):
        '''Connecting end nodes if we do not have a cycle.
        Returns:
            modified dict withan extra edge'''
        startNode,  lastNode = self.endNodes()
        if startNode != []:
            self.nodesDic[lastNode].append(startNode)

    def edgesNum(self):
        '''Count the total number of edges.
        Returns:
            the total number of edges'''
        totalEdges = 0
        for k,v in self.nodesDic.items():
            totalEdges += len(v)
        return totalEdges

    def count(self,nodesDic,node): #countinf the number of untouched edges for a given node
        '''Counts the number of edges from a node.
        Argumnets:
            nodesDic - a directed graph of nodes in dict format.
            node - the value of the node
        Returns:
            the number of remaind edges from that node
            '''
        remaindEdges = len(nodesDic[node])
        return remaindEdges

    #idea was taken from http://www.graph-magics.com/articles/euler.php
    #witten in collaboration with Julia Philipp
    def eulerianPath(self):
        '''Finds a eulerian path.
        Returns:
            a eulerian path in a list format.'''
        self.connectEdges()
        totalEdges = self.edgesNum() # finding the total number of edges
        stack = [] # stack of the walking path
        location = '' # current location
        circuit = [] # the path that we are confident about
        nodesDic = self.nodesDic # directed graph of nodes in the dict format
        if self.startNode != []:
            location = self.startNode
        else:
            location = random.choice(list(nodesDic.keys()))
        while all(value == [] for value in self.nodesDic.values()) == False: # As long as any all edges have not been used
            while self.count(nodesDic, location) != 0: # counts the number of edges coming out of this location, as long as it is not 0 ( meaning we have not made a loop)
                stack.append(location) # give the location to the stack
                newLocation = nodesDic[location].pop() # go to a new location
                location = newLocation # set the new location to location
            while self.count(nodesDic,location) == 0: #if we have made a loop
                circuit.append(location) # append the location to circuit
                if stack != []: #if the stack is not empty
                    location = stack.pop() #pop from the stack, and put that as a new location
                else:
                    break
        circuit = circuit[::-1] # reverse the list
        if self.startNode != []:
            del circuit[-1] # deleting the last nuc in the sequence in case that we extract a real start point from out data
        return circuit


class SequenceReconstruct(object):
    '''Returns a SequenceReconstruct object to construct a sequence.
    Attributes:
        inputFile - a list of kmers.
    Methodes:
        sequence - returns the sequence.
        '''

    def __init__(self, inputFile):
        self.inputFile = inputFile

    def sequence(self):
        '''assembels the sequences from the kmer.
        Returns:
            a sequence'''
        s = ''
        for i, j in enumerate(self.inputFile):
            if i == 0:
                s = s + j
            else:
                s = s + j[-1]
        return s

class StringReconstruction(object):
    def __init__(self, inputSequence, kmerLen):
        self.inputSequence = inputSequence
        self.kmerLen = kmerLen

    def f(self):
        graph = BruijnGraph(edgeSize= self.kmerLen, inputSequence= self.inputSequence)
        dic = graph.nodesDic() # making a dictionary of nodes
        path = EulerianPath(nodesDic= dic) # intitation EulerianPath class
        eulerPath = path.eulerianPath() # find euler path
        seq = SequenceReconstruct(inputFile = eulerPath) # initiate SequenceReconstruct class
        s = seq.sequence() # get the sequence
        return sys.stdout.write(s) # stdout


def main():
    inputSequence = []
    kmerLen = ''
    for index, line in enumerate(sys.stdin):
        l = line.split("\n")[0]
        if index == 0:
            kmerLen = int(l)
        else:
            inputSequence.append(l)


    m = StringReconstruction(inputSequence = inputSequence, kmerLen = kmerLen)
    m.f()


if __name__ == "__main__":
    main()
