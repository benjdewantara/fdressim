import numpy as np


def generateMFile(resultFilename, gridDim, resDim, timeLevel):
    # MFile refers to typical MATLAB file extension (.m file)
    
    MFile = open(resultFilename+"-visualize-t%s-usingMRST.m" %(timeLevel), "w")
    resultFile = open(resultFilename+"-results.txt", "r")
    
    specString = """
[nx,ny,nz] = deal( %s, %s, %s);
[Lx,Ly,Lz] = deal( %s, %s, %s);
G = cartGrid([nx, ny, nz], [Lx, Ly, Lz ]);
G = computeGeometry(G);


""" %(gridDim[2], gridDim[1], gridDim[0], resDim[2], resDim[1], resDim[0])
    
    pressureData = np.loadtxt(resultFile)
    
    pressureArr = pressureData[timeLevel].reshape(gridDim)
    newPressureArr = np.zeros(gridDim, dtype='float').reshape(gridDim)
    
    for zIndx in range(gridDim[0]):
        newPressureArr[zIndx] = pressureArr[gridDim[0]-zIndx-1]
    
    newPressureArrStr = generatePresArr(newPressureArr.reshape(gridDim[0]*gridDim[1]*gridDim[2]))
    
    otherSpecString = """

pressure = reshape(pressure, [], 1);


show = true(G.cells.num,1);

plotCellData(G, pressure, show, ...
'EdgeColor','k' );

view(3), camproj perspective, colorbar
"""
    
    print(specString+newPressureArrStr+otherSpecString, file=MFile)
    
    resultFile.close()
    MFile.close()


def generatePresArr(pressureDist):
    pressureArrStr = "pressure = ["
    
    for indx in range(len(pressureDist)):
        pressureArrStr += str(pressureDist[indx])
        if(indx < len(pressureDist)-1):
            pressureArrStr += ", "
    
    pressureArrStr += "];"
    
    return pressureArrStr

if(__name__ == '__main__'):
    generateMFile('example-4', (1, 51, 51), (75, 5000, 5000), 13)


