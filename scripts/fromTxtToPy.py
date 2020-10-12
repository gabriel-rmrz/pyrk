import os
import sys

path = 'nanoAOD_filesPaths/'



filePath = path + sys.argv[1]


fileRead = open(filePath,'r')
print("Working with file ",filePath)

lines = fileRead.readlines()
fileName = filePath.strip('.txt')
print(fileName)
fileName = fileName.split('/')[-1]
print(fileName)


fileWrite = open(filePath.strip('.txt')+'.py', 'w')

print(fileWrite)
fileWrite.write(fileName.split('_')[0] + '_files = [')
for line in lines:
    linemod = line.split('trivcat')[1]
    linemod2 = linemod.strip('\n')
    fileWrite.write("'root://cms-xrd-global.cern.ch/" + linemod2 + "',\n")
    
fileWrite.write(']')

fileRead.close()
fileWrite.close()
