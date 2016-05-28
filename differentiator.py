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
        
        linEqWRTt = differentialInTime(indx, res, presBefore, dt)
        linEqWRTx, knownWRTx = differentialInX(indx, res, presBefore)
        linEqWRTy, knownWRTy = differentialInY(indx, res, presBefore)
        linEqWRTz, knownWRTz = differentialInZ(indx, res, presBefore)
        
        known[indx] += knownWRTx + knownWRTy + knownWRTz
        sle[indx] += linEqWRTt - linEqWRTx - linEqWRTy - linEqWRTz
        
    
    return sle, known

def knownRHS(nodeIndx, res, presBefore, dt):
    coordIndx = res.grid.nodes[nodeIndx].coordIndx
    coordIndxBefore = coordIndx
    coordIndxAfter = coordIndx
    
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



















def differentialInX(nodeIndx, res, presBefore):
    # determine the coordIndx that corresponds to nodeIndx
    coordIndx = res.grid.nodes[nodeIndx].coordIndx
    
    # coordIndxAfter and coordIndxBefore correspond to the coordinate that interact
    # with coordIndx
    coordIndxAfter, coordIndxBefore = coordIndx, coordIndx
    coordIndxAfter[2] += 1
    coordIndxBefore[2] -= 1
    
    linEq = np.zeros(res.grid.numOfNodes, dtype='float64')
    known = 0.0
    
    deltaLen = res.grid.deltaX
    
    boundaryPresCri = res.grid.nodes[nodeIndx].boundaryWRTx[1]
    
    # check if there's a 'before' boundary condition w.r.t. coordIndx
    bc = boundaryPresCri['before']
    if(bc != None):
        if(bc.bcType == 'd'):
            known -= transmissibility(coordIndx, coordIndxBefore, res, presBefore)*bc.value
        elif(bc.bcType == 'n'):
            linEq[nodeIndx] += transmissibility(coordIndx, coordIndxBefore, res, presBefore)
            known -= -1*bc.value*deltaLen
    else:
        nodeIndxBefore = np.ravel_multi_index(coordIndxBefore, res.grid.dims)
        linEq[nodeIndxBefore] += transmissibility(coordIndx, coordIndxBefore, res, presBefore)
    
    # check if there's an 'after' boundary condition w.r.t. coordIndx
    bc = boundaryPresCri['after']
    if(bc != None):
        if(bc.bcType == 'd'):
            known -= transmissibility(coordIndx, coordIndxAfter, res, presBefore)*bc.value
        elif(bc.bcType == 'n'):
            linEq[nodeIndx] += transmissibility(coordIndx, coordIndxAfter, res, presBefore)
            known -= 1*bc.value*deltaLen
    else:
        nodeIndxAfter = np.ravel_multi_index(coordIndxAfter, res.grid.dims)
        linEq[nodeIndxAfter] += transmissibility(coordIndx, coordIndxAfter, res, presBefore)
    
    
    # finally, with confidence, we perform calculation for nodeIndx (aka coordIndx)
    linEq[nodeIndx] += -1*(transmissibility(coordIndx, coordIndxBefore, res, presBefore) + transmissibility(coordIndx, coordIndxAfter, res, presBefore))
    
    return linEq, known




















def differentialInY(nodeIndx, res, presBefore):
    # determine the coordIndx that corresponds to nodeIndx
    coordIndx = res.grid.nodes[nodeIndx].coordIndx
    
    # coordIndxAfter and coordIndxBefore correspond to the coordinate that interact
    # with coordIndx
    coordIndxAfter, coordIndxBefore = coordIndx, coordIndx
    coordIndxAfter[1] += 1
    coordIndxBefore[1] -= 1
    
    linEq = np.zeros(res.grid.numOfNodes, dtype='float64')
    known = 0.0
    
    deltaLen = res.grid.deltaY
    
    boundaryPresCri = res.grid.nodes[nodeIndx].boundaryWRTy[1]
    
    # check if there's a 'before' boundary condition w.r.t. coordIndx
    bc = boundaryPresCri['before']
    if(bc != None):
        if(bc.bcType == 'd'):
            known -= transmissibility(coordIndx, coordIndxBefore, res, presBefore)*bc.value
        elif(bc.bcType == 'n'):
            linEq[nodeIndx] += transmissibility(coordIndx, coordIndxBefore, res, presBefore)
            known -= -1*bc.value*deltaLen
    else:
        nodeIndxBefore = np.ravel_multi_index(coordIndxBefore, res.grid.dims)
        linEq[nodeIndxBefore] += transmissibility(coordIndx, coordIndxBefore, res, presBefore)
    
    # check if there's an 'after' boundary condition w.r.t. coordIndx
    bc = boundaryPresCri['after']
    if(bc != None):
        if(bc.bcType == 'd'):
            known -= transmissibility(coordIndx, coordIndxAfter, res, presBefore)*bc.value
        elif(bc.bcType == 'n'):
            linEq[nodeIndx] += transmissibility(coordIndx, coordIndxAfter, res, presBefore)
            known -= 1*bc.value*deltaLen
    else:
        nodeIndxAfter = np.ravel_multi_index(coordIndxAfter, res.grid.dims)
        linEq[nodeIndxAfter] += transmissibility(coordIndx, coordIndxAfter, res, presBefore)
    
    
    # finally, with confidence, we perform calculation for nodeIndx (aka coordIndx)
    linEq[nodeIndx] += -1*(transmissibility(coordIndx, coordIndxBefore, res, presBefore) + transmissibility(coordIndx, coordIndxAfter, res, presBefore))
    
    return linEq, known
























def differentialInZ(nodeIndx, res, presBefore):
    # determine the coordIndx that corresponds to nodeIndx
    coordIndx = res.grid.nodes[nodeIndx].coordIndx
    
    # coordIndxAfter and coordIndxBefore correspond to the coordinate that interact
    # with coordIndx
    coordIndxAfter, coordIndxBefore = coordIndx, coordIndx
    coordIndxAfter[0] += 1
    coordIndxBefore[0] -= 1
    
    linEq = np.zeros(res.grid.numOfNodes, dtype='float64')
    known = 0.0
    
    deltaLen = res.grid.deltaZ
    
    boundaryPresCri = res.grid.nodes[nodeIndx].boundaryWRTy[1]
    
    # check if there's a 'before' boundary condition w.r.t. coordIndx
    bc = boundaryPresCri['before']
    if(bc != None):
        if(bc.bcType == 'd'):
            known -= transmissibility(coordIndx, coordIndxBefore, res, presBefore)*bc.value
        elif(bc.bcType == 'n'):
            linEq[nodeIndx] += transmissibility(coordIndx, coordIndxBefore, res, presBefore)
            known -= -1*bc.value*deltaLen
    else:
        nodeIndxBefore = np.ravel_multi_index(coordIndxBefore, res.grid.dims)
        linEq[nodeIndxBefore] += transmissibility(coordIndx, coordIndxBefore, res, presBefore)
    
    # check if there's an 'after' boundary condition w.r.t. coordIndx
    bc = boundaryPresCri['after']
    if(bc != None):
        if(bc.bcType == 'd'):
            known -= transmissibility(coordIndx, coordIndxAfter, res, presBefore)*bc.value
        elif(bc.bcType == 'n'):
            linEq[nodeIndx] += transmissibility(coordIndx, coordIndxAfter, res, presBefore)
            known -= 1*bc.value*deltaLen
    else:
        nodeIndxAfter = np.ravel_multi_index(coordIndxAfter, res.grid.dims)
        linEq[nodeIndxAfter] += transmissibility(coordIndx, coordIndxAfter, res, presBefore)
    
    
    # finally, with confidence, we perform calculation for nodeIndx (aka coordIndx)
    linEq[nodeIndx] += -1*(transmissibility(coordIndx, coordIndxBefore, res, presBefore) + transmissibility(coordIndx, coordIndxAfter, res, presBefore))
    
    return linEq, known





















def totalCompressibility(res, pres):
    return res.rock.poro(pres)*res.fluid.rho(pres)*(res.fluid.compress(pres)+res.rock.compress(pres))




















def transmissibility(coordIndx, wrtCoord, res, presBefore):
    def getBoundaryPres(pres, bc, direction, deltaLen):
        if(bc.bcType == 'n'):
            if(direction=='before'):
                return pres - (bc.value*deltaLen)
            elif(direction=='after'):
                return pres + (bc.value*deltaLen)
        elif(bc.bcType == 'd'):
            return bc.value
    
    
    
    # roughly speaking, wr is used to determine which 3D direction this transmissibility term is calculated
    # forward is used to determine if it's the forward transmissibility (x+) or not (x-)
    wr = None
    forward = True
    for i in range(len(coordIndx)):
        differ = coordIndx[i]-wrtCoord[i]
        if(differ != 0):
            wr = i
            if(differ > 0):
                forward = False
            break
    
    # if wr = 0, we're looking at the transmissibility with respect to z
    # if wr = 1, it is with respect to y
    # if wr = 2, it is with respect to x
    
    # we need to complete variables rho, perm, area, mu, Vb, and deltaLen
    # then return the transmissibility
    
    deltaLen = None
    area = None
    hasBC = None
    boundaryPresCri = None
    
    if(wr == 0):
        deltaLen = res.grid.deltaZ
        area = res.grid.deltaX * res.grid.deltaY
        hasBC = res.grid.nodes[np.ravel_multi_index(coordIndx, res.grid.dims)].boundaryWRTz[0]
        if(hasBC):
            boundaryPresCri = res.grid.nodes[np.ravel_multi_index(coordIndx, res.grid.dims)].boundaryWRTz[1]
    elif(wr == 1):
        deltaLen = res.grid.deltaY
        area = res.grid.deltaX * res.grid.deltaZ
        hasBC = res.grid.nodes[np.ravel_multi_index(coordIndx, res.grid.dims)].boundaryWRTy[0]
        if(hasBC):
            boundaryPresCri = res.grid.nodes[np.ravel_multi_index(coordIndx, res.grid.dims)].boundaryWRTy[1]
    elif(wr == 2):
        deltaLen = res.grid.deltaX
        area = res.grid.deltaY * res.grid.deltaZ
        hasBC = res.grid.nodes[np.ravel_multi_index(coordIndx, res.grid.dims)].boundaryWRTx[0]
        if(hasBC):
            boundaryPresCri = res.grid.nodes[np.ravel_multi_index(coordIndx, res.grid.dims)].boundaryWRTx[1]
    
    Vb = res.grid.Vb
    perm = res.rock.perm
    
    wrtPres = None
    
    if(hasBC):
        if(forward):
            wrtPres = getBoundaryPres(presBefore[coordIndx], boundaryPresCri['after'], 'after', deltaLen)
        else:
            wrtPres = getBoundaryPres(presBefore[coordIndx], boundaryPresCri['before'], 'before', deltaLen)
    else:
        if(forward):
            wrtPres = presBefore[wrtCoord]
        else:
            wrtPres = presBefore[wrtCoord]
    
    presAvg = (presBefore[coordIndx] + wrtPres)/2
    
    rho = res.fluid.rho(presAvg)
    mu = res.fluid.mu(presAvg)
    
    transmiss = rho*perm*area/(mu*Vb*deltaLen)
    return transmiss

#def getBoundaryPres(pres, bc, direction, deltaLen):
#    if(bc.bcType == 'n'):
#        if(direction=='before'):
#            return pres - (bc.value*deltaLen)
#        elif(direction=='after'):
#            return pres + (bc.value*deltaLen)
#    elif(bc.bcType == 'd'):
#        return bc.value


