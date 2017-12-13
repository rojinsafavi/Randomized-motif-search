from collections import defaultdict
import sys

class DagPath(object):
    '''Find a longest path between two nodes in an edge-weighted DAG.
    Attributes:
        source _ the start node
        sink  _ the end node
        setDic _ input dictionary with adjacent nodes as tuple, and the edge as value
    Methodes:
        startNodes _ finds nodes with no input
        topologicalSorting _ sorts a graph topologicaly
        reverse _ switchs the incoming and outcoming nodes
        longestPath _ finds the acore of the longest path
        dagPath _ finds the path in DAG
        '''

    def __init__(self, source, sink, setDic):
        self.source = source
        self.sink = sink
        self.setDic = setDic
        self.topologicalSort = []

    def startNodes(self):# possible start nodes are nodes with no incoming edge
        '''Returns a list of nodes with no incoming edge.
        Returns:
            a list of nodes with no incoming edge'''
        candidates = list({pair[0] for pair in self.setDic} - {pair[1] for pair in self.setDic})
        return candidates

    def topologicalSorting(self):
        '''Returns a list of topological sorted nodes.
        Returns:
            A topological ordered list.'''
        setNodes = set(self.setDic.keys()) # a set with all nodes in the setDic
        self.topologicalSort = [] #ordered set
        candidates =  self.startNodes() # a list candidates for starting the topological sorting( nodes with no incoming edges)

        while candidates != []:
            a = candidates.pop() #take the last element of the candidate list, an removing that element from the candidate list
            self.topologicalSort.append(a) # append the elment to the self.topologicalSort

            k = [] # subset of the setDic containing all nodes connected to the candidate node
            for i in setNodes:
                if i[0] == a: # first node is a
                    k.append(i[1]) #making a list of all

            for b in k:
                setNodes.remove((a,b)) # removing the edge from the setDic
                if b not in {edge[1] for edge in setNodes}:
                    candidates.append(b)
        if setNodes != set():
            print('the input setDic is not a dag')
        else:
            return self.topologicalSort

    def reverse(self, subOrdered):
        '''switchs the incoming and outcoming nodes.
        Input :
            a topological list
        Returns:
            a dict with (outcoming, incoming): edge'''
        incomingNodes = defaultdict(list)
        for pairs in self.setDic.keys():
            if pairs[0] in subOrdered:
                if pairs[1] not in incomingNodes:
                    incomingNodes[pairs[1]] = [pairs[0]]
                else:
                    incomingNodes[pairs[1]].append(pairs[0])
        return incomingNodes

    def longestPath(self):
        '''Returns the score of the longest path.
        Returns:
            Returns the score of the longest path'''
        topOrder = self.topologicalSorting()
        topOrder = topOrder[topOrder.index(self.source):topOrder.index(self.sink) + 1] # only keeping the topological ordered list from the source to sink node
        reversedDic = self.reverse(topOrder) # reversing the input dict from (incoming, outcomin): edge to (outcoming, incoming) = edge

        tempGraph = {}
        for i in topOrder:
            tempGraph[i] = float('-inf')
        tempGraph[self.source] = 0

        scoreDic = defaultdict(list)
        scoreDic[self.source].append(0)
        for node in topOrder:
            if node != self.source:
                for i in reversedDic[node]: # itereating over the incoming nodes to this node
                    temp =( tempGraph[i] + self.setDic[(i, node)])
                    if temp > tempGraph[node]:
                        tempGraph[node] = temp
                        scoreDic[node] = [temp] # dictionary of scores (node:score)
        return tempGraph[self.sink], scoreDic, reversedDic


    def dagPath (self):
        '''Returns the path from source to sink.
        Returns:
            Path from source to sink that coresponds to calculated path score'''
        final, scoreDic, reverseDic = self.longestPath()
        end = self.sink
        p = [end]
        score = float('-inf')
        while score != 0: # as long as the score is not 0
            endScore = scoreDic[end][0] # take the score of the node
            prevNode = reverseDic[end] # find the incoming nodes to this node
            for k in prevNode: # itterate over them, and pick the one that satisfies nodeScore - edgeScore = previous Node Score
                if k in scoreDic.keys():
                    if (endScore - self.setDic[(k, end)]) == scoreDic[k][0] :
                        p.append(k)
                        score = scoreDic[k][0]
            end = p[-1]
        return final, p[::-1]

def main():
    '''Process the input file, create an object of the class, and create the output'''

    setDic = {}# here I'm pre-processing rosalind input format
    for j, i in enumerate(sys.stdin):
        if j ==0:
            source = int(i.strip())
        elif j == 1:
            sink = int(i.strip())
        else:
            total = i.split("\n")[0].split(':')
            nodes = total[0].split('->')
            value = int(total[1])
            k = []
            for node in nodes:
                k.append(int(node))
            setDic[tuple(k)] = value

    if source == sink :
        sys.stdout.write(str(0) + "\n" + str(source))
    else:
        dagInit = DagPath(setDic= setDic, source= source, sink= sink)

        pathScore, finalPath = dagInit.dagPath()

        sys.stdout.write(str(pathScore) + "\n") #here I'm pre-processing rosalind output format
        for index, node in enumerate(finalPath):
            if index ==0:
                sys.stdout.write(str(node))
            else:
                sys.stdout.write('->' + str(node))

if __name__ == "__main__":
    main()
