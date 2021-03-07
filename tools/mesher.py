# -*- coding: utf-8 -*-
"""
Mesher specified for 

Jan. 2021
"""

import os
import yaml
import numpy as np
import math

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
    z0 = 0.0    # z-coordiante of center
    gnd0 = 0    # start node number
    pxb = 0.0   # x position of L1 - chamfer of backup roll 
    pxw = 0.0   # x position of L1 of working roll 
    
    xwr = np.array([])   # output: x coordinates of work roll
    nwr = np.array([],dtype=int)   # output: node number of work roll
    xbr = np.array([])   # output: x coordinates of support roll
    nbr = np.array([],dtype=int)   # output: node number of support roll
    
    edload = np.array([],dtype=int) # output: loading edge group
    
    n_nd = 0    # ouput: number of nodes 
    
    def generate(self):
        global meshsize, xyz, elements
        print("Generating mesh Intermediate Roll begin with:", self.gnd0, self.z0)
        
        x0 = -0.5*self.L1 + self.offset - self.L21 -self.L3
        
        # division along D3 radius direction
        dd3 = meshsize
        d0 = divmod( self.D3, dd3 )
        nd3 = int(d0[0])
        if( d0[1]>0.0 ):
            dd3 = self.D3/nd3
        #print("nd3",nd3, dd3, nd3*dd3)
        
        # division along L31 load edge
        ddf = meshsize
        d0 = divmod( self.Lf, ddf )
        ndf = int(d0[0])
        if( d0[1]>0.0 ):
            ddf = self.Lf/ndf
        #print(ndf, ddf, ndf*ddf)
        # division along no-load edge
        dds = meshsize
        d0 = divmod( self.L3-self.Lf, dds )
        nds = int(d0[0])
        if( d0[1]>0.5*dds ):
            nds = nds+1
        if( d0[1]>0.0 ):
            dds = (self.L3-self.Lf)/nds
        #print(nds, dds, nds*dds)
        
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
                elements = np.append(elements, [nd0+j,nd0+1+j,nd1+1+j,nd1+j])
            if( i<ndf ):
                self.edload = np.append(self.edload, [int(len(elements)/4)-1, 1] )
                
        # division along radius direction of L21
        dd21 = meshsize
        d0 = divmod( 0.5*(self.D2-self.D3), dd21 )
        nd21 = int(d0[0])
        if( d0[1]>0.0 ):
            dd21 = 0.5*(self.D2-self.D3)/nd21
        #print("radius of D2:",nd21, dd21, nd21*dd21)
        
        # division along horizontal direction of L21
        ddlf = meshsize
        d0 = divmod( self.L21, ddlf )
        ndlf = int(d0[0])
        if( d0[1]>0.0 ):
            ddlf = self.L21/ndlf
        #print("L21",ndlf, ddlf, ndlf*ddlf)
        
        # L3 to L21
        nd0 = cnt_temp - nd3 -1
        nd1 = cnt_temp + nd21
        for j in range(0,nd3):
            elements = np.append(elements, [nd0+j,nd0+1+j,nd1+1+j,nd1+j])
        
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
                elements = np.append(elements, [nd0+j,nd0+1+j,nd1+1+j,nd1+j])
                
        cnt_temp = cnt_xyz
                
        # division along radius direction of L1
        dd1 = meshsize
        d0 = divmod( 0.5*(self.D1-self.D2), dd1 )
        ndd1 = int(d0[0])
        if( d0[1]>0.0 ):
            dd1 = 0.5*(self.D1-self.D2)/ndd1
        print("p1",ndd1, dd1, ndd1*dd1)
        
        
        # Elements of L21 to L1
        nd0 = cnt_temp - 2*nd21 -nd3 -1
        nd1 = cnt_temp + ndd1
        for j in range(0,nd3+2*nd21):
            elements = np.append(elements, [nd0+j,nd0+1+j,nd1+1+j,nd1+j])
            
        # x coordinates of L1 division
        xmr = np.array([])
        xmr = np.append(xmr, -0.5*self.L1 + self.offset)  # start point
        xmr = np.append(xmr,-self.pxb)                    # chamfer
        xmr = np.append(xmr,0.0)                          # center
        xmr = np.append(xmr,self.pxb)                     # chamfer
        xmr = np.append(xmr,self.pxw)                     # end point of woring roll
        xmr = np.append(xmr, 0.5*self.L1 + self.offset)   # end point
        print("points along L1 of inter roll:",xmr)
        
        # save x-volue 
        x_value = np.array([])
        xs = -0.5*self.L1 + self.offset
        for i in range(0,len(xmr)-1):
            ms = meshsize
            d0 = divmod( xmr[i+1] - xmr[i], ms )
            nm = int(d0[0])
            if( d0[1]>0.0 ):
                ms = (xmr[i+1] - xmr[i])/nm
            print("division:",i, nm, ms, nm*ms)
            for j in range(0,nm):
                cx = xs+j*ms
                x_value = np.append(x_value, cx)
                if( i>0 and i<3 ):
                    self.xbr = np.append(self.xbr, cx)
                if( i<4 ):
                    self.xwr = np.append(self.xwr, cx)
            xs = cx + ms
        x_value = np.append(x_value, xs)
        self.xbr = np.append(self.xbr, self.pxb)
        self.xwr = np.append(self.xwr, self.pxw)
        print( "x_value", x_value )
        print( "xbr", self.xbr )
        print( "xwr", self.xwr )
        
        nz = len(z_value)
        for i in range(0,len(x_value)):
            x = x_value[i]
            if( x>=xmr[1] and x<=xmr[3] ):
                self.nbr = np.append(self.nbr, cnt_xyz)
            for j in range(0,ndd1+1):
                xyz = np.append(xyz, [x, self.z0+0.5*self.D1-dd1*j])
                cnt_xyz += 1
            for j in range(1,nz):
                xyz = np.append(xyz, [x, z_value[j]])
                cnt_xyz += 1
            for j in range(1,ndd1+1):
                xyz = np.append(xyz, [x, z_value[nz-1]-dd1*j])
                cnt_xyz += 1
            if( x<=xmr[4] ):
                self.nwr = np.append(self.nwr, cnt_xyz-1)
                
        print( "nbr", self.nbr )
        print( "nwr", self.nwr )

        for i in range(0,len(x_value)-1):
            nd0 =  cnt_temp + i*(2*ndd1+nz)
            nd1 = nd0 + 2*ndd1+nz
            for j in range(0,2*ndd1+nz-1):
                elements = np.append(elements, [nd0+j,nd0+1+j,nd1+1+j,nd1+j])
                
        cnt_temp = cnt_xyz
                
        # division along L22
        dl4 = meshsize
        d0 = divmod(self.L22, dl4 )
        ndl4 = int(d0[0])
        if( d0[1]>0.0 ):
            dl4 = self.L22/ndl4
        print("L22",ndl4, dl4, ndl4*dl4)
        
        for i in range(1,ndl4+1):
            x = 0.5*self.L1 + self.offset + i*dl4
            for j in range(0,nz):
                xyz = np.append(xyz, [x, z_value[j]])
                cnt_xyz += 1
        cxpos = x  # current x position
        
        # first row of L22
        nd0 = cnt_temp - nd3 -2*nd21 -ndd1-1
        nd1 = cnt_temp
        for j in range(0,nz-1):
            elements = np.append(elements, [nd0+j,nd0+1+j,nd1+1+j,nd1+j])

        #other L22
        for i in range(1,ndl4):
            nd0 = cnt_temp 
            nd1 = nd0 + 2*nd21+nd3 + 1
            for j in range(0,nz-1):
                elements = np.append(elements, [nd0+j,nd0+1+j,nd1+1+j,nd1+j])
                
        cnt_temp = cnt_xyz
        
        # division along L32 no load edge
        for i in range(1,nds+1):
            x = cxpos + i*dds
            for j in range(0,nd3+1):
                xyz = np.append(xyz, [x, self.z0+0.5*self.D3-dd3*j])
                cnt_xyz += 1
                
        cxpos = x  # current x position
                
        # first row of L32
        nd0 = cnt_temp - nd3 - nd21-1
        nd1 = cnt_temp
        for j in range(0,nd3):
            elements = np.append(elements, [nd0+j,nd0+1+j,nd1+1+j,nd1+j])
            
        # other row of L32
        for i in range(1,nds):
            nd0 = cnt_temp + (i-1)*(nd3+1)
            nd1 = nd0 + nd3+1
            for j in range(0,nd3):
                elements = np.append(elements, [nd0+j,nd0+1+j,nd1+1+j,nd1+j])
                
        cnt_temp = cnt_xyz
                
        # division along L32 load edge
        for i in range(1,ndf+1):
            x = cxpos +i*ddf
            for j in range(0,nd3+1):
                xyz = np.append(xyz, [x, self.z0+0.5*self.D3-dd3*j])
                cnt_xyz += 1
                
        for i in range(0,ndf):
            nd0 = cnt_temp + (i-1)*(nd3+1)
            nd1 = nd0 + nd3+1
            for j in range(0,nd3):
                elements = np.append(elements, [nd0+j,nd0+1+j,nd1+1+j,nd1+j])
            self.edload = np.append(self.edload, [int(len(elements)/4)-1, 1] )
        print("edload", self.edload)
                
        self.n_nd = cnt_xyz - self.gnd0
                
class BackupRoll :
    D1 = 0.0
    D2 = 0.0
    L1 = 0.0
    L2 = 0.0
    Lf = 0.0                       # Length of constrained edge
    chamfer = np.array([])
    z0 = 0.0                       # z-coordiante of center
    gnd0 = 0                       # start node number 
    n_nd = 0                       # ouput: number of nodes 
    
    xbr = np.array([])             # input: x coordinates of support roll
    nbr = np.array([],dtype=int)   # input: node number of support roll
    
    ndfix = np.array([],dtype=int) # output: fixed node group
    
    def generate(self):
        global meshsize, xyz, elements
        print("Generating mesh of backup roll begin with:", self.gnd0, self.z0)
        
        ## L2 ##
        mslf = meshsize
        d0 = divmod( self.Lf, mslf )
        nlf = int(d0[0])
        if( d0[1]>0.0 ):
            mslf = self.Lf/nlf
            
        msl = meshsize
        d0 = divmod( self.L2-self.Lf, msl )
        nl = int(d0[0])
        if( d0[1]>0.0 ):
            msl = (self.L2-self.Lf)/nl
        
        ## D2 ##
        msd = meshsize
        d0 = divmod( self.D2, msd )
        nd = int(d0[0])
        if( d0[1]>0.0 ):
            msd = self.D2/nd
       # print("D2:", nd, msd, nd*msd)

        cnt_xyz = self.gnd0
        
        x0 = -0.5*self.L1 - self.L2
        for i in range(0,nlf+1):
            x = x0 +i*mslf
            self.ndfix = np.append(self.ndfix, cnt_xyz)
            for j in range(0,nd+1):
                xyz = np.append(xyz, [x, self.z0+0.5*self.D2-msd*j])
                cnt_xyz += 1
            self.ndfix = np.append(self.ndfix, cnt_xyz-1)
        for i in range(1,nl):
            x = x0 + self.Lf + i*msl
            for j in range(0,nd+1):
                xyz = np.append(xyz, [x, self.z0+0.5*self.D2-msd*j])
                cnt_xyz += 1
                
        cnt_temp = cnt_xyz
                
        for i in range(0,nlf+nl-1):
            nd0 = self.gnd0 + i*(nd+1)
            nd1 = nd0 + nd+1
            for j in range(0,nd):
                elements = np.append(elements, [nd0+j,nd0+1+j,nd1+1+j,nd1+j])
              
        # D1-D2
        msc = meshsize
        d0 = divmod( 0.5*(self.D1-self.D2), msc )
        nc = int(d0[0])
        if( d0[1]>0.0 ):
            msc = 0.5*(self.D1-self.D2)/nc
       # print("D2 to D1", nc, msc, nc*msc)

        # L chamfer        
        mscf = meshsize
        d0 = divmod( self.chamfer[0], mscf )
        ncf = int(d0[0])
        if( d0[1]>0.0 ):
            mscf = self.chamfer[0]/ncf
        dz = self.chamfer[1]/self.chamfer[0]*mscf
       # print("L chamfer", ncf, mscf, ncf*mscf,dz)
        
        cx = x0 + self.L2
        for i in range(0,ncf+1):
            x = cx +i*mscf
            zt = self.z0 + 0.5*self.D1 - self.chamfer[1] + dz*i
            zt1 = self.z0 - 0.5*self.D1 + self.chamfer[1] - dz*i
            xyz = np.append(xyz, [x, zt])
            cnt_xyz += 1
            for j in range(1,nc+1):
                xyz = np.append(xyz, [x, self.z0 + 0.5*self.D1 - msc*j])
                cnt_xyz += 1
            for j in range(1,nd+1):
                xyz = np.append(xyz, [x, self.z0 + 0.5*self.D2 - msd*j])
                cnt_xyz += 1
            for j in range(1,nc):
                xyz = np.append(xyz, [x, self.z0 - 0.5*self.D2 - msc*j])
                cnt_xyz += 1
            xyz = np.append(xyz, [x, zt1])
            cnt_xyz += 1
            
        # Element D2->D1
        nd0 = cnt_temp -nd -1
        nd1 = cnt_temp + nc
        for j in range(0,nd):
            elements = np.append(elements, [nd0+j,nd0+1+j,nd1+1+j,nd1+j])
            
        # Element D1 to first chamfer
        for i in range(0,ncf):
            nd0 = cnt_temp + i*(2*nc+nd+1)
            nd1 = nd0 + 2*nc+nd+1
            for j in range(0,2*nc+nd):
               elements = np.append(elements, [nd0+j,nd0+1+j,nd1+1+j,nd1+j])
        ce = len(elements)
        elements[ce-2] = self.nbr[0]
        xyz = np.delete(xyz, [2*cnt_xyz - 1, 2*cnt_xyz-2])
        cnt_xyz -= 1
        cnt_temp = cnt_xyz
               
        # Node: left chamfer to right chamfer
        nx = len(self.xbr)
        for i in range(1,nx):
            for j in range(0,nc+1):
                xyz = np.append(xyz, [self.xbr[i], self.z0 + 0.5*self.D1 - msc*j])
                cnt_xyz += 1
                #print(self.xbr[i], self.z0 + 0.5*self.D1 - msc*j)
            for j in range(1,nd+1):
                xyz = np.append(xyz, [self.xbr[i], self.z0 + 0.5*self.D2 - msd*j])
                #print(self.xbr[i], self.z0 + 0.5*self.D2 - msd*j)
                cnt_xyz += 1
            for j in range(1,nc):
                xyz = np.append(xyz, [self.xbr[i], self.z0 - 0.5*self.D2 - msc*j])
                cnt_xyz += 1
                
        # Element: left chamfer to right chamfer
        for i in range(0,nx-1):
            nd0 = cnt_temp + (i-1)*(2*nc+nd)
            nd1 = nd0 + 2*nc+nd
            for j in range(0,2*nc+nd-1):
               elements = np.append(elements, [nd0+j,nd0+1+j,nd1+1+j,nd1+j])
            nd0 = nd0 + 2*nc+nd-2
            nd1 = nd1 + 2*nc+nd-2
            elements = np.append(elements, [nd0,self.nbr[i],self.nbr[i+1],nd1])
            
        cnt_temp = cnt_xyz
            
        # right chamfer to right edge
        cx = 0.5*self.L1- self.chamfer[0]
        for i in range(1,ncf+1):
            x = cx +i*mscf
            zt = self.z0 + 0.5*self.D1 - self.chamfer[1] + dz*(ncf-i)
            zt1 = self.z0 - 0.5*self.D1 + self.chamfer[1] - dz*(ncf-i)
            xyz = np.append(xyz, [x, zt])
            cnt_xyz += 1
            for j in range(1,nc+1):
                xyz = np.append(xyz, [x, self.z0 + 0.5*self.D1 - msc*j])
                cnt_xyz += 1
            for j in range(1,nd+1):
                xyz = np.append(xyz, [x, self.z0 + 0.5*self.D2 - msd*j])
                cnt_xyz += 1
            for j in range(1,nc):
                xyz = np.append(xyz, [x, self.z0 - 0.5*self.D2 - msc*j])
                cnt_xyz += 1
            xyz = np.append(xyz, [x, zt1])
            cnt_xyz += 1
            
        # Element::right chamfer to right edge
        nd0 = cnt_temp - (2*nc+nd)
        nd1 = nd0 + 2*nc+nd
        print(nd0,nd1)
        for j in range(0,2*nc+nd-1):
            elements = np.append(elements, [nd0+j,nd0+1+j,nd1+1+j,nd1+j])
        nd0 = nd0 + 2*nc+nd-1
        nd1 = nd1 + 2*nc+nd-1
        elements = np.append(elements, [nd0,self.nbr[nx-1],nd1+1,nd1])

        for i in range(0,ncf-1):
            nd0 = cnt_temp + i*(2*nc+nd+1)
            nd1 = nd0 + 2*nc+nd+1
            for j in range(0,2*nc+nd):
               elements = np.append(elements, [nd0+j,nd0+1+j,nd1+1+j,nd1+j])
               
        cnt_temp = cnt_xyz
        
        # Nodes along left Lf
        for i in range(1,nl):
            x = 0.5*self.L1 + i*msl
            for j in range(0,nd+1):
                xyz = np.append(xyz, [x, self.z0+0.5*self.D2-msd*j])
                cnt_xyz += 1
                
        # Nodes to left End
        for i in range(0,nlf+1):
            x = 0.5*self.L1 + (self.L2-self.Lf) +i*mslf
            self.ndfix = np.append(self.ndfix, cnt_xyz)
            for j in range(0,nd+1):
                xyz = np.append(xyz, [x, self.z0+0.5*self.D2-msd*j])
                cnt_xyz += 1
            self.ndfix = np.append(self.ndfix, cnt_xyz-1)
        print("ndfix",self.ndfix)
               
        # Element: L1 to left L2
        nd0 = cnt_temp - (nc+nd+1)
        nd1 = cnt_temp
        for j in range(0,nd):
            elements = np.append(elements, [nd0+j,nd0+1+j,nd1+1+j,nd1+j])
            
        # Element to left end
        for i in range(0,nlf+nl-1):
            nd0 = cnt_temp +i*(nd+1)
            nd1 = nd0 + nd+1
            for j in range(0,nd):
               elements = np.append(elements, [nd0+j,nd0+1+j,nd1+1+j,nd1+j])
               
        self.n_nd = cnt_xyz - self.gnd0
               
class WorkingRoll :
    D1 = 0.0
    D2 = 0.0
    D3 = 0.0
    L1 = 0.0
    L2 = 0.0
    L3 = 0.0
    Lf = 0.0                       # Length of constrained edge
    z0 = 0.0                       # z-coordiante of center
    gnd0 = 0                       # start node number 
    n_nd = 0                       # ouput: number of nodes 
    xb = 0.0                       # start position of inter roll
    
    xbr = np.array([])             # input: x coordinates with inter roll
    nbr = np.array([],dtype=int)   # input: node number with inter roll
    
    def generate(self):
        global meshsize, xyz, elements
        print("Generating mesh of working roll begin with:", self.gnd0, self.z0)
        
        x0 = -0.5*self.L1 - self.L2 -self.L3
        
        # division along D3 radius direction
        dd3 = meshsize
        d0 = divmod( self.D3, dd3 )
        nd3 = int(d0[0])
        if( d0[1]>0.0 ):
            dd3 = self.D3/nd3
        print("nd3",nd3, dd3, nd3*dd3)
        
        # division along L31 load edge
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
                
        # division along radius direction of L21
        dd21 = meshsize
        d0 = divmod( 0.5*(self.D2-self.D3), dd21 )
        nd21 = int(d0[0])
        if( d0[1]>0.0 ):
            dd21 = 0.5*(self.D2-self.D3)/nd21
        print("radius of D2:",nd21, dd21, nd21*dd21)
        
        # division along horizontal direction of L2
        ddlf = meshsize
        d0 = divmod( self.L2, ddlf )
        ndlf = int(d0[0])
        if( d0[1]>0.0 ):
            ddlf = self.L2/ndlf
        print("L2",ndlf, ddlf, ndlf*ddlf)
        
        # L3 to L2
        nd0 = cnt_temp - nd3 -1
        nd1 = cnt_temp + nd21
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
        print("D1-D2",ndd1, dd1, ndd1*dd1)
        
        # division to edge of Intermediate roll
        dl11 = meshsize
        d0 = divmod( 0.5*self.L1-self.xb, dl11 )
        ndl11 = int(d0[0])
        if( d0[1]>0.0 ):
            dl11 = (0.5*self.L1-self.xb)/ndl11
        print("L11",ndl11, dl11, ndl11*dl11)
        print("z_value:",z_value)     
        nz = len(z_value)
        for i in range(0,ndl11):
            x = -0.5*self.L1 + i*dl11
            for j in range(0,ndd1+1):
                xyz = np.append(xyz, [x, self.z0+0.5*self.D1-dd1*j])
                cnt_xyz += 1
            for j in range(1,nz):
                xyz = np.append(xyz, [x, z_value[j]])
                cnt_xyz += 1
            for j in range(1,ndd1+1):
                xyz = np.append(xyz, [x, z_value[nz-1]-dd1*j])
                cnt_xyz += 1
        
        # Elements of L21 to L1
        nd0 = cnt_temp - 2*nd21 -nd3 -1
        nd1 = cnt_temp + ndd1
        for j in range(0,nd3+2*nd21):
            elements = np.append(elements, [nd0+j,nd1+j,nd1+1+j,nd0+1+j])
            
        # Element in L11
        for i in range(0,ndl11-1):
            nd0 = cnt_temp + i*(2*ndd1+nz) 
            nd1 = nd0 + 2*ndd1+nz
            for j in range(0,2*ndd1+nz-1):
                elements = np.append(elements, [nd0+j,nd1+j,nd1+1+j,nd0+1+j])
                
        cnt_temp = cnt_xyz
                
        # Nodes to left edge
        print(cnt_xyz)
        print(self.nbr)
        print(self.xbr)
        for i in range(0,len(self.nbr)):
            for j in range(1,ndd1):
                xyz = np.append(xyz, [self.xbr[i], self.z0+0.5*self.D1-dd1*j])
                cnt_xyz += 1
            for j in range(0,nz):
                xyz = np.append(xyz, [self.xbr[i], z_value[j]])
                cnt_xyz += 1
            for j in range(1,ndd1+1):
                xyz = np.append(xyz, [self.xbr[i], z_value[nz-1]-dd1*j])
                cnt_xyz += 1
         
        # first row of L12
        nd0 = cnt_temp - (2*ndd1+nz)
        nd1 = nd0 + 2*ndd1+nz
        elements = np.append(elements, [nd0,self.nbr[0],nd1,nd0+1])
        for i in range(0,len(self.nbr)-1):
            nd0 = cnt_temp + i*(2*ndd1+nz-1) 
            nd1 = nd0 + (2*ndd1+nz-1)
            elements = np.append(elements, [self.nbr[i],self.nbr[i+1],nd1,nd0])
        
        # ramins row of L12
        for i in range(0,len(self.nbr)):
            nd0 = cnt_temp - (2*ndd1+nz) + 1 + i*(2*ndd1+nz-1) 
            nd1 = nd0 + (2*ndd1+nz-1)
            for j in range(0,2*ndd1+nz-2):
                elements = np.append(elements, [nd0+j,nd1+j,nd1+1+j,nd0+1+j])
                
        cnt_temp = cnt_xyz
        
        # L22
        for i in range(1,ndlf+1):
            x = 0.5*self.L1 + i*ddlf
            for j in range(0,nz):
                xyz = np.append(xyz, [x, z_value[j]])
                cnt_xyz += 1
                
        # Element: L1 to L22
        nd0 = cnt_temp - (2*ndd1+nz-1)
        nd1 = cnt_temp
        for j in range(0,nz-1):
            elements = np.append(elements, [nd0+j,nd1+j,nd1+1+j,nd0+1+j])
        
        #Element: L22
        for i in range(1,ndlf):
            nd0 = cnt_temp + (i-1)*nz 
            nd1 = nd0 + nz
            for j in range(0,nz-1):
                elements = np.append(elements, [nd0+j,nd1+j,nd1+1+j,nd0+1+j])
                
        cnt_temp = cnt_xyz
                
        #L32
        for i in range(1,nds):
            x = 0.5*self.L1 + self.L2 + i*dds
            for j in range(0,nd3+1):
                xyz = np.append(xyz, [x, self.z0+0.5*self.D3-dd3*j])
                print(x, self.z0+0.5*self.D3-dd3*j)
                cnt_xyz += 1
        for i in range(0,ndf+1):
            x = 0.5*self.L1 + self.L2 + self.L3- self.Lf +i*ddf
            for j in range(0,nd3+1):
                xyz = np.append(xyz, [x, self.z0+0.5*self.D3-dd3*j])
                print(x, self.z0+0.5*self.D3-dd3*j)
                cnt_xyz += 1
                
        # Element: L32-L32f
        nd0 = cnt_temp - (nd3+nd21+1)
        nd1 = cnt_temp
        print( cnt_temp, nd0,nd1)
        for j in range(0,nd3):
            elements = np.append(elements, [nd0+j,nd1+j,nd1+1+j,nd0+1+j])
            
        # Element: L32f
        for i in range(1,ndf+nds):
            nd0 = cnt_temp + (i-1)*(nd3+1)
            nd1 = nd0 + nd3+1
            print(nd0,nd1)
            for j in range(0,nd3):
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
bRoll = BackupRoll()
wRoll = WorkingRoll()

for key, value in data.items():
    if key == "MeshSize":
        meshsize = value
    elif key == "Backup Roll":
        for k2, v2 in value.items():
            if k2 == "D1":
                bRoll.D1 = v2
            elif k2 == "L1":
                bRoll.L1 = v2
            elif k2 == "D2":
                bRoll.D2 = v2
            elif k2 == "L2":
                bRoll.L2 = v2
            elif k2 == "Lf":
                bRoll.Lf = v2
            elif k2 == "chamfer":
                bRoll.chamfer = v2
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
            if k2 == "D1":
                wRoll.D1 = v2
            elif k2 == "L1":
                wRoll.L1 = v2
            elif k2 == "D2":
                wRoll.D2 = v2
            elif k2 == "L2":
                wRoll.L2 = v2
            elif k2 == "D3":
                wRoll.D3 = v2
            elif k2 == "L3":
                wRoll.L3 = v2
            elif k2 == "Lf":
                wRoll.Lf = v2
            print (key,k2,v2)
            

## initial ##
n_node = 0
xyz = np.array([])
n_element = 0
elements = np.array([], dtype=int)

zw = 0.0   # z-coordinate of working roll

iRoll.gnd0 = 0  
iRoll.pxb = 0.5*bRoll.L1 - bRoll.chamfer[0]
iRoll.pxw = 0.5*wRoll.L1
iRoll.z0 = zw + 0.5*wRoll.D1 + 0.5*iRoll.D1
iRoll.generate()

bRoll.gnd0 = iRoll.n_nd
bRoll.xbr = iRoll.xbr
bRoll.nbr = iRoll.nbr
bRoll.z0 = iRoll.z0 + 0.5*iRoll.D1 + 0.5*bRoll.D1
bRoll.generate()

wRoll.gnd0 = iRoll.n_nd + bRoll.n_nd
wRoll.z0 = zw
wRoll.xb = 0.5*iRoll.L1 - iRoll.offset
wRoll.xbr = iRoll.xwr
wRoll.nbr = iRoll.nwr
wRoll.generate()

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