import numpy as np


readfile = open("example-2-results.txt", 'r')
workfile = open("workfile.txt", "w")

dims = (1, 10, 20)

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

print(newPressureDist.reshape(dims[0]*dims[1]*dims[2]).__str__(), file=workfile)







workfile.close()
readfile.close()

