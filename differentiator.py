import numpy as np
from scipy import linalg
#import example

def oneStepDifferentiator(res, presBefore, dt):
    # this function should be designed such that it returns
    # a system of linear equations
    # with its known solution for each linear equation
    
    # sle stands for system of linear equations :P
    sle = np.zeros([res.grid.numOfNodes, res.grid.numOfNodes], dtype='float64')
    known = np.zeros(res.grid.numOfNodes, dtype='float64')
    
    for indx in range(res.grid.numOfNodes):
        known[indx] += knownRHS(indx, res, presBefore, dt)
        
        sle[indx] += differentialInTime(indx, res, presBefore, dt) \
        - differentialInX(indx, res, presBefore) \
        - differentialInY(indx, res, presBefore) \
        - differentialInZ(indx, res, presBefore) 
    
    return sle, known

def knownRHS(nodeIndx, res, presBefore, dt):
    coordIndx = res.grid.nodes[nodeIndx].coordIndx
    coordIndxBefore = coordIndx
    coordIndxBefore[0] -= 1
    try:
        nodeIndxBefore = np.ravel_multi_index(coordIndxBefore, res.grid.dims)
    except ValueError:
        nodeIndxBefore = None
    
    coordIndxAfter = coordIndx
    coordIndxAfter[0] += 1
    try:
        nodeIndxAfter = np.ravel_multi_index(coordIndxAfter, res.grid.dims)
    except ValueError:
        nodeIndxBefore = None
    
    knownTerm = res.grid.nodes[nodeIndx].src/res.grid.Vb \
                + res.fluid.rho(presBefore[coordIndx])*res.fluid.poro(presBefore[coordIndx])*totalCompressibility(res, presBefore[coordIndx])/dt * presBefore[coordIndx] \
                + (transmissibility(coordIndx, coordIndxAfter, res, presBefore)-transmissibility(coordIndx, coordIndxBefore, res, presBefore))*res.fluid.rho(presBefore[coordIndx])*32.1740485*res.grid.deltaZ
    
    return knownTerm
    


def differentialInTime(nodeIndx, res, presBefore, dt):
    # In a linear equation with some unknowns variables,
    # the left-hand side is usually arranged to contain the unknown terms
    # while the right-hand side contains the constant terms.
    # Hence, the name linEq and knownRHS
    coordIndx = res.grid.nodes[nodeIndx].coordIndx
    linEq = np.zeros(res.grid.numOfNodes, dtype='float64')
    
    linEq[nodeIndx] += res.fluid.rho(presBefore[coordIndx])*res.fluid.poro(presBefore[coordIndx])*totalCompressibility(res, presBefore[coordIndx])/dt
    
    return linEq

def differentialInX(nodeIndx, res, presBefore, linEq):
    coordIndx = res.grid.nodes[nodeIndx].coordIndx
    
    coordIndxBefore = coordIndx
    coordIndxBefore[2] -= 1
    try:
        nodeIndxBefore = np.ravel_multi_index(coordIndxBefore, res.grid.dims)
    except ValueError:
        nodeIndxBefore = None
        
    coordIndxAfter = coordIndx
    coordIndxAfter[2] += 1
    try:
        nodeIndxAfter = np.ravel_multi_index(coordIndxAfter, res.grid.dims)
    except ValueError:
        nodeIndxBefore = None
    
    
    linEq = np.zeros(res.grid.numOfNodes, dtype='float64')
    linEq[nodeIndxBefore] += 
    
    
    pass

def differentialInY(nodeIndx, res, presBefore, linEq):
    coordIndx = res.grid.nodes[nodeIndx].coordIndx
    coordIndxBefore = coordIndx
    coordIndxBefore[1] -= 1
    coordIndxAfter = coordIndx
    coordIndxAfter[1] += 1
    
    linEq = np.zeros(res.grid.numOfNodes, dtype='float64')
    
    
    pass

def differentialInZ(nodeIndx, res, presBefore, linEq):
    coordIndx = res.grid.nodes[nodeIndx].coordIndx
    coordIndxBefore = coordIndx
    coordIndxBefore[0] -= 1
    coordIndxAfter = coordIndx
    coordIndxAfter[0] += 1
    
    linEq = np.zeros(res.grid.numOfNodes, dtype='float64')
    
    
    pass




def totalCompressibility(res, pres):
    return res.rock.poro(pres)*res.fluid.rho(pres)*(res.fluid.compress(pres)+res.rock.compress(pres))


def transmissibility(coordIndx, wrtCoord, res, presBefore):
    pass




