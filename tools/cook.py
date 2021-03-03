# -*- coding: utf-8 -*-
"""
Mesher specified for 

Mar. 2021
"""

import sys
import numpy as np

### 入出力定義 ###

## 出力ファイル名のヘッダー。
FILEOUT = "cook.vtk"
fo=open(FILEOUT, 'w')

## initial ##
n_node = 0
xyz = np.array([])
ndgroup = np.array([], dtype=int)
n_element = 0
elements = np.array([], dtype=int)
surface = np.array([], dtype=int)

LY1 = 44.0
LY2 = 16.0
LX1 = 48.0

nL = 32
if len(sys.argv) == 2:
    nL = int( sys.argv[1] )
    
for i in range(0,nL+1):
    y0 = LY1/nL * i
    y1 = LY1 + LY2/nL *i
    for j in range(0,nL+1):
        x = LX1/nL *j
        y = y0 + (y1-y0)/nL * j
        xyz = np.append(xyz, [x, y])

nd = 0
ne = 0
ndgroup = np.append(ndgroup, nd)
for i in range(0,nL):
    for j in range(0,nL):
        elements = np.append(elements, [nd,nd+1,nd+nL+2,nd+1+nL])
        nd += 1
        ne += 1
    nd += 1
    ndgroup = np.append(ndgroup, nd)
    surface = np.append(surface, [ne-1,1])

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
fo.write("FIELD NODESET 1\n")
fo.write("N1 1 "+str(len(ndgroup)) +" int\n")
for i in range(0,len(ndgroup)):
    fo.write(str(ndgroup[i])+'\n')
fo.write("FIELD EDGESET 1\n")
ns = int(len(surface)/2)
fo.write("S1 2 "+str(ns) +" int\n")
for i in range(0,ns):
    fo.write(str(surface[2*i])+' '+str(surface[2*i+1])+'\n')

fo.close()