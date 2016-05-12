
import numpy as np

timeStepNum = 5
nx = 25
ny = 2
nz = 2

dims = (nz, ny, nx)

numOfNodes = nx*ny*nz

someList = np.zeros([25])

print(someList)

someList = np.arange(0, 25)

print(someList)

someList = np.reshape(someList, (5, 5))

print(someList)












class Node(object):
    def __init__(self, name):
        self.name = int(name)
        self.coordIndx = np.unravel_index(self.name, dims)
    def getName(self):
        return self.name
    def getCoordIndx(self):
        return self.coordIndx
    def __eq__(self, other):
        return self.name == other.name
    def __hash__(self):
        return self.name.__hash__()

class Rock(object):
    def __init__(self):
        pass

class Fluid(object):
    def __init__(self):
        pass



nodes = []

for i in range(numOfNodes):
    nodes.append(Node(i))

matrixOfIndices = np.arange(0, numOfNodes).reshape((nz, ny, nx))

for zIndx in range(nz):
    for yIndx in range(ny):
        for xIndx in range(nx):
            index = (zIndx, yIndx, xIndx)
            #print(matrixOfIndices[index])
            nodes[matrixOfIndices[index]].setIndex(index)
            #Node(matrixOfIndices[(zIndx, yIndx, xIndx)]).setIndex((zIndx, yIndx, xIndx))
            #print(Node(matrixOfIndices[(zIndx, yIndx, xIndx)]) in nodes)
            #if(zIndx == 0) and (yIndx == 0) and (xIndx == 2):
                #print(matrixOfIndices[(zIndx, yIndx, xIndx)])
                #print(Node(matrixOfIndices[(zIndx, yIndx, xIndx)]) in nodes)

































