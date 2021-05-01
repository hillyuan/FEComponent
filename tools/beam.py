# -*- coding: utf-8 -*-
"""
Mesher specified for 

Mar. 2021
"""

import math
import numpy as np

### 入出力定義 ###

## 出力ファイル名のヘッダー。
FILEOUT = "beam.vtk"
fo=open(FILEOUT, 'w')

## initial ##
n_node = 0
xyz = np.array([])
ndleft = np.array([], dtype=int)
ndright = np.array([], dtype=int)
n_element = 0
elements = np.array([], dtype=int)
surface = np.array([], dtype=int)
elethick = np.array([], dtype=float)

L= 100.0
W = 5.0
d = 0.25
ddL = d
dL = divmod( L, d )
nL = int(dL[0])
if( dL[1]>0.0 ):
    ddL = L/nL
ddW = d
dW = divmod( W, d )
nW = int(dW[0])
if( dW[1]>0.0 ):
    ddW = W/nW

nd = -1;
nc = 0
for i in range(0,nL+1):
    x = ddL*i
    if( i==nL ):
        nc1 = nd+1
    for j in range(0,nW+1):
        y = j*ddW
        xyz = np.append(xyz, [x, y])
        nd +=1
        if( i==0 ):
            ndleft = np.append(ndleft, nd)
        if( i==nL ):
            ndright = np.append(ndright, nd)


nd = 0
ne = 0
for i in range(0,nL):
    for j in range(0,nW):
        elements = np.append(elements, [nd,nd+nW+1,nd+2+nW,nd+1])
        cy = 0.25*( xyz[2*nd+1] + xyz[2*(nd+nW+1)+1] + xyz[2*(nd+2+nW)+1] + xyz[2*(nd+1)+1] )
        r = 2.0*math.sqrt( 0.25*W*W - (cy-0.5*W)*(cy-0.5*W) )
        elethick = np.append(elethick, r)
        nd += 1
        ne += 1
    nd += 1
    surface = np.append(surface, [ne-1,2])

## Output ##
fo.write("# vtk DataFile Version 4.0\n")
fo.write("Cook's membrane problem\n")
fo.write("ASCII\n")
fo.write("DATASET UNSTRUCTURED_GRID\n")
n_node = int(len(xyz)/2)
coord = xyz.reshape(n_node,2)
fo.write("POINTS "+str(n_node)+ " double\n")
for i in range(0,n_node):
    fo.write(str(coord[i,0])+' '+str(coord[i,1]) +' 0.0\n')
    
n_element = int(len(elements)/4)
quad = elements.reshape(n_element,4)
fo.write("CELLS "+str(n_element)+ " " + str(5*n_element) + "\n")
for i in range(0,n_element):
    fo.write('4 '+str(quad[i,0])+' '+str(quad[i,1]) + ' '+str(quad[i,2])+' '+str(quad[i,3]) + '\n')
fo.write("CELL_TYPES "+str(n_element) + "\n")
for i in range(0,n_element):
    fo.write('9\n')
fo.write("FIELD NODESET 4\n")
fo.write("NLEFT 1 "+str(len(ndleft)) +" int\n")
for i in range(0,len(ndleft)):
    fo.write(str(ndleft[i])+'\n')
fo.write("NRIGHT 1 "+str(len(ndright)) +" int\n")
for i in range(0,len(ndright)):
    fo.write(str(ndright[i])+'\n')
fo.write("NC0 1 1 int\n")
fo.write('0\n')
fo.write("NC1 1 1 int\n")
fo.write(str(nc1)+'\n')
fo.write("FIELD EDGESET 1\n")
ns = int(len(surface)/2)
fo.write("S1 2 "+str(ns) +" int\n")
for i in range(0,ns):
    fo.write(str(surface[2*i])+' '+str(surface[2*i+1])+'\n')
fo.write("CELL_DATA "+str(n_element) + "\n")
fo.write("SCALARS cell_thickness float 1\n")
fo.write("LOOKUP_TABLE default\n")
for i in range(0,n_element):
    fo.write(str(elethick[i]) +'\n')

fo.close()