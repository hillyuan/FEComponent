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

class IntermediateRoll :
    D1 = 0.0
    D2 = 0.0
    D3 = 0.0
    L1 = 0.0
    L21 = 0.0
    L22 = 0.0
    L3 = 0.0
    offset = 0.0
    Lf = 0.0
    z0 = 0.0
    gnd0 = 0
    pxb = 0.0
    pxw = 0.0
    
    px0 = 0.0
    px1 = 0.0
    
    def generate(self):
        global meshsize, xyz, elements
        print("Generating mesh begin with:", self.gnd0, self.z0)
        self.px0 = 0.5*self.L1 - self.offset
        self.px1 = 0.5*self.L1 + self.offset
        
        x0 = -0.5*self.L1 + self.offset - self.L21 -self.L3
        
        # division along radius direction
        dd3 = meshsize
        d0 = divmod( self.D3, dd3 )
        nd3 = int(d0[0])
        if( d0[1]>0.0 ):
            dd3 = self.D3/nd3
        print("nd3",nd3, dd3, nd3*dd3)
        
        # division along load edge
        ddf = meshsize
        d0 = divmod( self.Lf, ddf )
        ndf = int(d0[0])
        if( d0[1]>0.0 ):
            ddf = self.Lf/ndf
        print(ndf, ddf, ndf*ddf)
        # division along no-load edge
        dds = meshsize
        d0 = divmod( self.L3-self.Lf, dds )
        nds = int(d0[0])
        if( d0[1]>0.0 ):
            dds = (self.L3-self.Lf)/nds
        print(nds, dds, nds*dds)
        
        cnt_xyz = self.gnd0
        for i in range(0,ndf+1):
            x = x0 +i*ddf
            for j in range(0,nd3+1):
                xyz = np.append(xyz, [x, self.z0+0.5*self.D3-dd3*j])
                cnt_xyz += 1
        for i in range(1,nds):
            x = x0 +self.Lf+i*dds
            for j in range(0,nd3+1):
                xyz = np.append(xyz, [x, self.z0+0.5*self.D3-dd3*j])
                cnt_xyz += 1
        cnt_temp = cnt_xyz
                
        for i in range(0,ndf+nds-1):
            nd0 = self.gnd0 + i*(nd3+1)
            nd1 = nd0 + nd3+1
            for j in range(0,nd3):
                elements = np.append(elements, [nd0+j,nd1+j,nd1+1+j,nd0+1+j])
                
        # division along radius direction of linker
        dd21 = meshsize
        d0 = divmod( 0.5*(self.D2-self.D3), dd21 )
        nd21 = int(d0[0])
        if( d0[1]>0.0 ):
            dd21 = 0.5*(self.D2-self.D3)/nd21
        print("p21",nd21, dd21, nd21*dd21)
        
        # division along horizontal direction of linker
        ddlf = meshsize
        d0 = divmod( self.L21, ddlf )
        ndlf = int(d0[0])
        if( d0[1]>0.0 ):
            ddlf = self.L21/ndlf
        print("ndlf",ndlf, ddlf, ndlf*ddlf)
        
        # L3 to L21
        nd0 = cnt_temp - nd3 -1
        nd1 = cnt_temp + nd21
        print(nd0,nd1)
        for j in range(0,nd3):
            elements = np.append(elements, [nd0+j,nd1+j,nd1+1+j,nd0+1+j])
        
        for i in range(0,ndlf):
            x = x0 +self.L3+i*ddlf
            for j in range(0,nd21+1):
                xyz = np.append(xyz, [x, self.z0+0.5*self.D2-dd21*j])
                cnt_xyz += 1
            for j in range(1,nd3+1):
                xyz = np.append(xyz, [x, self.z0+0.5*self.D3-dd3*j])
                cnt_xyz += 1
            for j in range(1,nd21+1):
                xyz = np.append(xyz, [x, self.z0-0.5*self.D3-dd21*j])
                cnt_xyz += 1
        
        # save above z-volue        
        z_value = np.array([])
        for j in range(0,nd21+1):
            z_value = np.append(z_value, self.z0+0.5*self.D2-dd21*j)
        for j in range(1,nd3+1):
            z_value = np.append(z_value, self.z0+0.5*self.D3-dd3*j)
        for j in range(1,nd21+1):
            z_value = np.append(z_value, self.z0-0.5*self.D3-dd21*j)
        
        for i in range(0,ndlf-1):
            nd0 = cnt_temp + i*(2*nd21+nd3+1) 
            nd1 = nd0 + 2*nd21+nd3 + 1
            for j in range(0,2*nd21+nd3):
                elements = np.append(elements, [nd0+j,nd1+j,nd1+1+j,nd0+1+j])
                
        cnt_temp = cnt_xyz
                
        # division along radius direction of L1
        dd1 = meshsize
        d0 = divmod( 0.5*(self.D1-self.D2), dd1 )
        ndd1 = int(d0[0])
        if( d0[1]>0.0 ):
            dd1 = 0.5*(self.D1-self.D2)/ndd1
        print("p1",ndd1, dd1, ndd1*dd1)
        
        # division along horizontal direction of L1
        dl1 = meshsize
        d0 = divmod( self.L1, dl1 )
        ndl1 = int(d0[0])
        if( d0[1]>0.0 ):
            dl1 = self.L1/ndl1
        print("dl1",ndl1, dl1, ndl1*dl1)
        
        # L21 to L1
        nd0 = cnt_temp - 2*nd21 -nd3 -1
        nd1 = cnt_temp + ndd1
        for j in range(0,nd3+2*nd21):
            elements = np.append(elements, [nd0+j,nd1+j,nd1+1+j,nd0+1+j])

        nz = len(z_value)
        for i in range(0,ndl1+1):
            x = x0 +self.L3 +self.L21 +i*ddlf
            for j in range(0,ndd1+1):
                xyz = np.append(xyz, [x, self.z0+0.5*self.D1-dd1*j])
                cnt_xyz += 1
            for j in range(1,nz):
                xyz = np.append(xyz, [x, z_value[j]])
                cnt_xyz += 1
            for j in range(1,ndd1+1):
                xyz = np.append(xyz, [x, z_value[nz-1]-dd1*j])
                cnt_xyz += 1
            
        for i in range(0,ndl1):
            nd0 =  cnt_temp + i*(2*ndd1+nz)
            nd1 = nd0 + 2*ndd1+nz
            for j in range(0,2*ndd1+nz-1):
                elements = np.append(elements, [nd0+j,nd1+j,nd1+1+j,nd0+1+j])


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

iRoll = IntermediateRoll()

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
            if k2 == "D1":
                iRoll.D1 = v2
            elif k2 == "D2":
                iRoll.D2 = v2
            elif k2 == "D3":
                iRoll.D3 = v2
            elif k2 == "L1":
                iRoll.L1 = v2
            elif k2 == "L21":
                iRoll.L21 = v2
            elif k2 == "L22":
                iRoll.L22 = v2
            elif k2 == "L3":
                iRoll.L3 = v2
            elif k2 == "Lf":
                iRoll.Lf = v2
            elif k2 == "offset":
                iRoll.offset = v2
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
            
iRoll.pxb = 0.5*backupLb
iRoll.pxw = 0.5*workLw

## initial ##
n_node = 0
xyz = np.array([])
n_element = 0
elements = np.array([], dtype=int)
zw = 0.0
iRoll.z0 = zw + 0.5*workDw + 0.5*iRoll.D1
zb = iRoll.z0 + 0.5*iRoll.D1 + 0.5*backupDb
print(zw,iRoll.z0,zb)
xw = -0.5*workLw - worklw
#xi = -0.5*interLi + offset - interli
xb = -0.5*backupLb - backuplb
#print(xw,xi,xb)

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
    
## main part
msl = meshsize
dd = divmod( backupLb-2.0*chamfer[0], msl )
nl = int(dd[0])
if( dd[1]>0.0 ):
    msl = ( backupLb-2.0*chamfer[0] )/nc
print("Dl", nl, msl, nl*msl)

for i in range(0,nd+1):
    x = x0 +backuplb +i*msl
    for j in range(1,nc+1):
        xyz = np.append(xyz, [x, z0-msc*j-chamfer[1]])
    
for i in range(0,nd):
    nd0 = maxnd + 1 + (2*nc+n0+1)*i
    nd1 = nd0 + 2*nc+n0 + 1
    print( "nd", nd0, nd1)
    for j in range(0,2*nc+n0):
        elements = np.append(elements, [nd0+j,nd1+j,nd1+1+j,nd0+1+j])

iRoll.gnd0 = 1227     
iRoll.generate()

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