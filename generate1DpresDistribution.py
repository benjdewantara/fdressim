import numpy as np

def generatePresDist(resultFilename, gridDim, timeLevel, yIndx=None, xIndx=None):
    # This only works for a 2D-flow problem!!
    resultFile = open(resultFilename+"-results.txt", "r")
    oneDResultFile = open(resultFilename+"-results-y%s-x%s.txt" %(yIndx, xIndx), "a")
    
    pressureData = np.loadtxt(resultFile)
    
    pressureDataAtOneTimelevel = pressureData[timeLevel].reshape(gridDim)
    
    oneDimPresDist = None
    
    
    if(yIndx != None):
        oneDimPresDist = pressureDataAtOneTimelevel[0, yIndx, :]
        pass
    elif(xIndx != None):
        oneDimPresDist = pressureDataAtOneTimelevel[0, :, xIndx]
        pass
    
    for elm in oneDimPresDist:
        print(str(elm)+" ", file=oneDResultFile, end="")
    
    print("", file=oneDResultFile)
    
    oneDResultFile.close()
    resultFile.close()




if(__name__ == '__main__'):
    nTime = 14
    for tIndx in range(nTime):
        generatePresDist('example-5', (1, 101, 101), tIndx, yIndx=50)
    
    print("1D pressure distribution has been generated")
