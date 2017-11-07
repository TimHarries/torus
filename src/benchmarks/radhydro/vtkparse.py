from numpy import *

def getvtk (filename):
    "reads in an unstructured vtk file (raw format not compressed xml format) and returns a dictionary with the cell positions and their associated scalars and vectors"
    values={}
    gotpoints=False
    gotcells=False


    with open(filename) as f:
        line=f.readline()
        while line:
            if not(gotpoints) and line[:6]=='POINTS':
                n=int(line.split()[1])
                pos=empty([3,n], dtype=float32)
                line=f.readline()
                try:
                    i=0
                    while True:
                        pos[:,i]=[float(a) for a in line.split()]
                        i+=1
                        line=f.readline()
                except ValueError:
                    gotpoints=True
            elif not(gotcells) and line[:5]=='CELLS':
                n=int(line.split()[1])
                cellpos=empty([3,n], dtype=float32)
                line=f.readline()
                try:
                    i=0
                    while True:
                        _,a,b,c,d=[int(x) for x in line.split()]
                        cellpos[:,i]=(pos[:,a]+pos[:,b]+pos[:,c]+pos[:,d])/4.0
                        i+=1
                        line=f.readline()
                except ValueError:
                    values['position']=cellpos
                    gotcells=True
            elif gotpoints and gotcells and line[:7]=='SCALARS':
                scalar=line.split()[1]
                values[scalar]=[]
                line=f.readline() #unused line for type of lookup
                try:
                    line=f.readline()
                    s=float(line)
                    values[scalar].append(s)
                    line=f.readline()
                    while True :
                        s=float(line)
                        values[scalar].append(s)
                        line=f.readline()
                except ValueError :
                    values[scalar]=array(values[scalar], dtype=float32)
            elif gotpoints and gotcells and line[:7]=='VECTORS':
                vector=line.split()[1]
                values[vector]=[]
                line=f.readline() #unused line for type of lookup
                try:
                    line=f.readline()
                    v=[float(x) for x in line.split()]
                    values[vector].append(v)
                    line=f.readline()
                    while True :
                        v=[float(x) for x in line.split()]
                        values[vector].append(v)
                        line=f.readline()
                except ValueError :
                    values[vector]=array(values[vector], dtype=float32)
            else :
                line=f.readline()
            
    return values
