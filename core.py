
import numpy as np
import differentiator

class Node(object):
    def __init__(self, nodeIndx, dims):
        self.nodeIndx = int(nodeIndx)
        self.coordIndx = np.unravel_index(self.nodeIndx, dims)
        self.src = None
        self.initBoundaryNode(dims)
        
    def initBoundaryNode(self, dims):
        self.boundaryWRTx = np.array([False, {'before':None, 'after':None}])
        self.boundaryWRTy = np.array([False, {'before':None, 'after':None}])
        self.boundaryWRTz = np.array([False, {'before':None, 'after':None}])
        
        if(self.coordIndx[2] == 0):
            self.boundaryWRTx[0] = True
            self.boundaryWRTx[1]['before'] = BoundaryCondition()
        elif(self.coordIndx[1] == 0):
            self.boundaryWRTy[0] = True
            self.boundaryWRTy[1]['before'] = BoundaryCondition()
        elif(self.coordIndx[0] == 0):
            self.boundaryWRTz[0] = True
            self.boundaryWRTz[1]['before'] = BoundaryCondition()
        
        if(self.coordIndx[2] == dims[2]-1):
            self.boundaryWRTx[0] = True
            self.boundaryWRTx[1]['after'] = BoundaryCondition()
        elif(self.coordIndx[1] == dims[1]-1):
            self.boundaryWRTy[0] = True
            self.boundaryWRTy[1]['after'] = BoundaryCondition()
        elif(self.coordIndx[0] == dims[0]-1):
            self.boundaryWRTz[0] = True
            self.boundaryWRTz[1]['after'] = BoundaryCondition()
    
    def setSrc(self):
        pass


class Grid(object):
    def __init__(self, dims):
        nz, ny, nx = dims
        self.dims = dims
        self.numOfNodes = nz*ny*nx
        self.nodes = self.initNodes()
    
    def initNodes(self):
        nodes = np.zeros([self.numOfNodes], dtype='O')
        for i in range(self.numOfNodes):
            nodes[i] = Node(i, self.dims)
        return nodes
    
    def initGridblockGeometry(self, resDim):
        self.deltaX = resDim[2]/self.dims[2]
        self.deltaY = resDim[1]/self.dims[1]
        self.deltaZ = resDim[0]/self.dims[0]
        self.Vb = self.deltaX * self.deltaY * self.deltaZ
        

class Rock(object):
    def __init__(self, name):
        self.name = name
        pass
    
    def poro(self, pres):
        return 1
    
    def compress(self, pres):
        return 1

class Fluid(object):
    def __init__(self, name):
        self.name = name
        pass
    
    def rho(self, pres):
        return 1
    
    def mu(self, pres):
        pass
    
    def compress(self, pres):
        return 1

class Reservoir(object):
    def __init__(self, grid, fluid, rock, resDim):
        self.grid = grid
        self.fluid = fluid
        self.rock = rock
        self.resDim = resDim
        self.arrayOfIndices = np.arange(self.grid.numOfNodes).reshape(self.grid.dims)
        self.grid.initGridblockGeometry(self.resDim)
    
    def setInitPressure(self, initPressure):
        self.initPressure = np.full(self.grid.dims, initPressure, dtype='float64')
    
    def addBoundaryCondition(self, bc, **kwargs):
        if not all([key in ["x", "y", "z"] and kwargs[key] in ['after', 'before'] for key in kwargs]):
            raise ValueError("Key argument must be one of 'x', 'y', or 'z' and the argument must be either 'after' or 'before'")
        
        nz, ny, nx = self.grid.dims
        
        for key in kwargs:
            if(key=='x'):
                if(kwargs[key] == 'before'):
                    for nodes in self.grid.nodes[self.arrayOfIndices[:, :, 0]]:
                        for node in nodes:
                            node.boundaryWRTx[1]['before'] = bc
                elif(kwargs[key] == 'after'):
                    for nodes in self.grid.nodes[self.arrayOfIndices[:, :, nx-1]]:
                        for node in nodes:
                            node.boundaryWRTx[1]['after'] = bc
            elif(key == 'y'):
                if(kwargs[key] == 'before'):
                    for nodes in self.grid.nodes[self.arrayOfIndices[:, 0, :]]:
                        for node in nodes:
                            node.boundaryWRTy[1]['before'] = bc
                elif(kwargs[key] == 'after'):
                    for nodes in self.grid.nodes[self.arrayOfIndices[:, ny-1, :]]:
                        for node in nodes:
                            node.boundaryWRTy[1]['after'] = bc
            elif(key == 'z'):
                if(kwargs[key] == 'before'):
                    for nodes in self.grid.nodes[self.arrayOfIndices[0, :, :]]:
                        for node in nodes:
                            node.boundaryWRTz[1]['before'] = bc
                elif(kwargs[key] == 'after'):
                    for nodes in self.grid.nodes[self.arrayOfIndices[nz-1, :, :]]:
                        for node in nodes:
                            node.boundaryWRTz[1]['after'] = bc
        
        


class BoundaryCondition(object):
    def __init__(self, bcType="n", value=0):
        if(bcType not in ['n', 'd']):
            raise ValueError("Wrong boundary condition specification!")
        self.bcType = bcType
        self.value = value
        pass




def runOneStepSimulation():
    # This is a function to run a one time-step simulation
    pass































