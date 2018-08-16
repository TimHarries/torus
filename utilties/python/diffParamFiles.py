#!/usr/local/bin/python3
#
# Compare two Torus parameter files and report differences
# Run with: python3 diffParamFiles.py file1 file2
#
# Author: D. Acreman, July 2018

import sys
import torusParam

if len(sys.argv) != 3:
    print("Two arguments are required which are the names of the files to be compared")
    raise RuntimeError 

filename1=str(sys.argv[1])
filename2=str(sys.argv[2])
print("Comparing " + filename1 + " and " + filename2)

# Read parameter files
file1=torusParam.parameterList(file=filename1)
file2=torusParam.parameterList(file=filename2)

def compareParams(firstFile, secondFile):
    firstFileOnly=[]
    
    for p in firstFile.parameters:
        if not p in secondFile.parameters:
            firstFileOnly.append(p)

    if firstFileOnly== []:
        print("All parameters in " + str(firstFile.filename) + " occur in " + str(secondFile.filename))
    else:
        print("Parameters only in " + firstFile.filename + str(firstFileOnly) )
            
compareParams(file1,file2)
compareParams(file2,file1)

# Report differences in parameter values
num_diffs=0
for p in file1.parameters:
    if p in file2.parameters:
        if file1.parameters[p] != file2.parameters[p]:
            num_diffs += 1
            print(p + " differs: " + str(file1.parameters[p]) + str(file2.parameters[p]) )

if num_diffs == 0:
    print("No differences in parameter settings found")
elif num_diffs == 1:
    print(str(num_diffs) + " difference in parameter settings found")
else:
    print(str(num_diffs) + " differences in parameter settings found")
