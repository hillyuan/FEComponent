# -*- coding: utf-8 -*-
"""
Mesher specified for 

Jan. 2021
"""

import os
import yaml
import numpy as np
import math
import matplotlib.pyplot   as plt


### 入出力定義 ###

## 入力パラメータファイル名。
FILEIN     = ".\data\parameters.yaml"  
## 出力ファイル名のヘッダー。
FILEOUT = "rollsystem.vtk"
fo=open(FILEOUT, 'w')

f=open(FILEIN, 'r')
data = yaml.safe_load(f)
f.close()
print(data)

for key, value in data.items():
    if key == "MeshSize":
        meshsize = value
    elif key == "Backup Roll":
        for k2, v2 in value.items():
            if k2 == "Db":
                backupDb = v2
            elif k2 == "Lb":
                backupLb = v2
            elif k2 == "db":
                backupdb = v2
            elif k2 == "lb":
                backuplb = v2
            elif k2 == "lf":
                backuplf = v2
            elif k2 == "chamfer":
                chamfer = v2
            print (key,k2,v2)
    elif key == "Intermediate Roll":
        for k2, v2 in value.items():
            if k2 == "Di":
                interDi = v2
            elif k2 == "Li":
                interLi = v2
            elif k2 == "di":
                interdi = v2
            elif k2 == "li":
                interli = v2
            elif k2 == "offset":
                offset = v2
            print (key,k2,v2)
    elif key == "Work Roll":
        for k2, v2 in value.items():
            if k2 == "Dw":
                workDw = v2
            elif k2 == "Lw":
                workLw = v2
            elif k2 == "dw":
                workdw = v2
            elif k2 == "lw":
                worklw = v2
            print (key,k2,v2)
            

## initial ##
n_node = 0
xyz = np.array([])
n_element = 0
elements = np.array([], dtype=int)
zw = 0.0
zi = zw + 0.5*workDw + 0.5*interDi
zb = zi + 0.5*interDi + 0.5*backupDb
print(zw,zi,zb)
xw = -0.5*workLw - worklw
xi = -0.5*interLi + offset - interli
xb = -0.5*backupLb - backuplb
print(xw,xi,xb)

## Backup roll ##
msl = meshsize
d0 = divmod( backupdb, msl )
n0 = int(d0[0])
print( d0 )
if( d0[1]>0.0 ):
    msl = backupdb/n0

msd = meshsize
d1 = divmod( backuplf, msd )
n1 = int(d1[0])
print(d1)
if( d1[1]>0.0 ):
    msd = backuplf/n1
    
msf = meshsize
d2 = divmod( backuplb-backuplf, msf )
n2 = int(d2[0])
print(d2)
if( d2[1]>0.0 ):
    msf = (backuplb-backuplf)/n2
    
for i in range(0,n1+1):
    x = xb +i*msd
    for j in range(0,n0+1):
        xyz = np.append(xyz, [x, zb+0.5*backupdb-msl*j])
for i in range(1,n2):
    x = xb+backuplf +i*msf
    for j in range(0,n0+1):
        xyz = np.append(xyz, [x, zb+0.5*backupdb-msl*j])
        
for i in range(0,n1+n2-1):
    nd0 = i*(n0+1)
    nd1 = nd0 + n0+1
    for j in range(0,n0):
        elements = np.append(elements, [nd0+j,nd1+j,nd1+1+j,nd0+1+j])
minnd = nd0+n0+1
maxnd = nd1+n0
print("maxnd", maxnd, minnd)
        
## Intermediate Roll ##
z0 = zb + 0.5*backupDb
x0 = xb + backuplb

msd = meshsize
dd = divmod( chamfer[0], msd )
nd = int(dd[0])
if( dd[1]>0.0 ):
    msd = chamfer[0]/n1
print("Lw_offset", nd, msd, nd*msd)

# chamfer
msc = meshsize
dd = divmod( 0.5*(backupDb-backupdb)-chamfer[1], msl )
nc = int(dd[0])
if( dd[1]>0.0 ):
    msc = ( 0.5*(backupDb-backupdb) -chamfer[1] )/nc
dz = chamfer[1]/chamfer[0]*msd
print("Dc", nc, msc, nc*msc, dz)

# element link Db and dd
nd0 = minnd
nd1 = maxnd+nc+1
for j in range(0,n0):
    elements = np.append(elements, [nd0+j,nd1+j,nd1+1+j,nd0+1+j])

for i in range(0,nd+1):
    x = x0 +i*msd
    xyz = np.append(xyz, [x, z0-(chamfer[1]-dz*i)])
    for j in range(1,nc+1):
        xyz = np.append(xyz, [x, z0-msc*j-chamfer[1]])
    for j in range(1,n0+1):
        xyz = np.append(xyz, [x, zb+0.5*backupdb-msl*j])
    for j in range(1,nc):
        xyz = np.append(xyz, [x, zb-0.5*backupdb-msc*j])
    xyz = np.append(xyz, [x, zb-0.5*backupDb-dz*(i-nd)])
    
for i in range(0,nd):
    nd0 = maxnd + 1 + (2*nc+n0+1)*i
    nd1 = nd0 + 2*nc+n0 + 1
    print( "nd", nd0, nd1)
    for j in range(0,2*nc+n0):
        elements = np.append(elements, [nd0+j,nd1+j,nd1+1+j,nd0+1+j])

## Output ##
fo.write("# vtk DataFile Version 4.0\n")
fo.write("3-Roller System\n")
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

fo.close()