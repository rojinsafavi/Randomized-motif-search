#Q6
import numpy as np
import random
from collections import defaultdict
import sys


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
        return circuit

    def outputFormat(self):
        '''Rosaling stdout format.
        Returns:
            Rosaling stdout format'''
        c = self.eulerianPath()
        if self.startNode != []:
            del c[-1]
        for i,j in enumerate(c):
            if i == 0:
                sys.stdout.write(str(j))
            else:
                sys.stdout.write('->' + str(j))

def main():
    '''Preprocessing the input stdin file, create necessary attributes, initiate the class'''
    dic = defaultdict(list)
    for line in sys.stdin.readlines():
        tmp = line[:-1].split(" -> ")
        k = np.int(tmp[0])
        v = [np.int(i) for i in tmp[1].split(",")]
        dic[k] += v

    d = EulerianPath(nodesDic= dic)
    d.outputFormat()
if __name__ == "__main__":
    main()
