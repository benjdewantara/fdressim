
def printArrayToFile(f, nparray):
    for elm in nparray:
        print("%.4f" %(elm), end=' ', file=f)
    print("\n", end='', file=f)


import numpy as np


readfile = open("example-4-results.txt", 'r')
workfile = open("workfile.txt", "w")

dims = (5, 51, 51)

pressureData = np.loadtxt(readfile)

pressureDist = pressureData[23].reshape(dims)

newPressureDist = np.zeros(dims, dtype=float)

for zIndx in range(dims[0]):
    newPressureDist[zIndx] = pressureDist[dims[0]-zIndx-1]


print(pressureDist)


printThis = "["

for i in range(1000):
    printThis += "3000"
    if(i < 1000-1):
        printThis += ", "

printThis += "];"

#print(newPressureDist.reshape(dims[0]*dims[1]*dims[2]), file=workfile)

#np.savetxt(workfile, newPressureDist.reshape(dims[0]*dims[1]*dims[2]))

printArrayToFile(workfile, newPressureDist.reshape(dims[0]*dims[1]*dims[2]))


workfile.close()
readfile.close()

