#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 13:49:58 2021

@author: Sergio Lucarini
"""

import gmsh
import numpy as np






def locate(npart,rad,packing,micro_type):
    # vol_tot=L[0]*L[1]*L[2]
  if micro_type=='l': return None,None
  if micro_type=='r':
    h=np.array([[-L[0],-L[1],-L[2]],[-L[0],-L[1],0],[-L[0],-L[1],L[2]],\
        [-L[0],0,-L[2]],[-L[0],0,0],[-L[0],0,L[2]],\
        [-L[0],L[1],-L[2]],[-L[0],L[1],0],[-L[0],L[1],L[2]],\
        [0,-L[1],-L[2]],[0,-L[1],0],[0,-L[1],L[2]],\
        [0,0,-L[2]],[0,0,0],[0,0,L[2]],\
        [0,L[1],-L[2]],[0,L[1],0],[0,L[1],L[2]],\
        [L[0],-L[1],-L[2]],[L[0],-L[1],0],[L[0],-L[1],L[2]],\
        [L[0],0,-L[2]],[L[0],0,0],[L[0],0,L[2]],\
        [L[0],L[1],-L[2]],[L[0],L[1],0],[L[0],L[1],L[2]]])
    images=int(3**3)
    stop=True;trials=0
    while stop:
        trials+=1
        stop=False
        rr=np.repeat(rad,npart*images)
        collocation=np.zeros([npart*images])
        packings=np.zeros([npart*images])
        xc=np.zeros([npart*images,3])
        import random
        np.random.seed(random.randint(0,4294967296))
        for i in range(npart):
            itercol=0
            ies=np.where(collocation==2)[0]
            rres=rr[ies]
            xces=xc[ies]
            while np.all(collocation[images*i:images*(i+1)]==0):
                x=np.random.rand(3)*L
                # xi=np.random.randint(0,xc0[:,0].size)
                # x=xc0[xi]
                diff=x-xces
                distances=np.sqrt(np.einsum('xi,xi->x',diff,diff))
                if not np.any(distances<rr[i*images+images//2]/packing+rres):
                    collocation[images*i:images*(i+1)]=2
                    packings[images*i:images*(i+1)]=packing
                    xc[images*i:images*(i+1)]=(x+h).reshape(images,3)
                    break
                # packing+=0.001
                itercol+=1
        #        break ####artifi
                if itercol>100000:
                    print('restarting',i)
                    stop=True
                    break
            if stop: break
            if (i+1)%100==0:
                print('grain ',i+1,' / ',npart,'     ||   iters ',itercol,'pack',packing)
    
    #            break
    print('grain ',i+1,' / ',npart,'     ||   iters ',itercol,'pack',packing)
  if micro_type=='s':
    xc=np.array([[0.5,0.5,0.5]])
    rr=np.repeat(rad,len(xc))
  if micro_type=='f':
    xc=np.array([[0.5,0.5,0.],
              [0.,0.5,0.5],
              [0.5,0.,0.5],
              [1.,0.5,0.5],
              [0.5,0.5,1.],
              [0.5,1.,0.5]])
    rr=np.repeat(rad,len(xc))
  if micro_type=='b':
    xc=np.array([[0.5,0.5,0.5],
              [0.,0.,0.],
              [0.,0.,1.],
              [0.,1.,0.],
              [0.,1.,1.],
              [1.,0.,0.],
              [1.,0.,1.],
              [1.,1.,0.],
              [1.,1.,1.]])
    rr=np.repeat(rad,len(xc))
  return xc,rr

def select(xc,rr,micro_type):
    if micro_type=='l': return None
    orig=np.array([0.,0.,0.])
    idex=np.array([],dtype='int')
    for i in range(len(xc)):
        # dist1=np.abs(xc[i]-orig)
        # dist2=np.abs(xc[i]-L)
        
        if np.all(orig<xc[i]) and np.all(xc[i]<L):
            idex=np.append(idex,int(i))
            print(xc[i],'inside',rr[i])
            continue
            
        elif np.all(orig[[0,1]]<xc[i,[0,1]]) and np.all(xc[i,[0,1]]<L[[0,1]]):
            if np.abs(xc[i,2]-orig[2])<rr[i] or np.abs(xc[i,2]-L[2])<rr[i]:
                print(xc[i],'border',rr[i])
                idex=np.append(idex,int(i))
                continue
        elif np.all(orig[[1,2]]<xc[i,[1,2]]) and np.all(xc[i,[1,2]]<L[[1,2]]):
            if np.abs(xc[i,0]-orig[0])<rr[i] or np.abs(xc[i,0]-L[0])<rr[i]:
                print(xc[i],'border',rr[i])
                idex=np.append(idex,int(i))
                continue
        elif np.all(orig[[0,2]]<xc[i,[0,2]]) and np.all(xc[i,[0,2]]<L[[0,2]]):
            if np.abs(xc[i,1]-orig[1])<rr[i] or np.abs(xc[i,1]-L[1])<rr[i]:
                print(xc[i],'border',rr[i])
                idex=np.append(idex,int(i))
                continue
           
        elif orig[0]<xc[i,0] and xc[i,0]<L[0]:
            if (xc[i,1]-orig[1])**2+(xc[i,2]-orig[2])**2<rr[i]**2 or\
                (xc[i,1]-orig[1])**2+(xc[i,2]-L[2])**2<rr[i]**2 or\
                (xc[i,1]-L[1])**2+(xc[i,2]-orig[2])**2<rr[i]**2 or\
                (xc[i,1]-L[1])**2+(xc[i,2]-L[2])**2<rr[i]**2:
                print(xc[i],'border2',rr[i])
                idex=np.append(idex,int(i))
                continue
        elif orig[1]<xc[i,1] and xc[i,1]<L[1]:
            if (xc[i,0]-orig[0])**2+(xc[i,2]-orig[2])**2<rr[i]**2 or\
                (xc[i,0]-orig[0])**2+(xc[i,2]-L[2])**2<rr[i]**2 or\
                (xc[i,0]-L[0])**2+(xc[i,2]-orig[2])**2<rr[i]**2 or\
                (xc[i,0]-L[0])**2+(xc[i,2]-L[2])**2<rr[i]**2:
                print(xc[i],'border2',rr[i])
                idex=np.append(idex,int(i))
                continue
        elif orig[2]<xc[i,2] and xc[i,2]<L[2]:
            if (xc[i,1]-orig[1])**2+(xc[i,0]-orig[0])**2<rr[i]**2 or\
                (xc[i,1]-orig[1])**2+(xc[i,0]-L[0])**2<rr[i]**2 or\
                (xc[i,1]-L[1])**2+(xc[i,0]-orig[0])**2<rr[i]**2 or\
                (xc[i,1]-L[1])**2+(xc[i,0]-L[0])**2<rr[i]**2:
                print(xc[i],'border2',rr[i])
                idex=np.append(idex,int(i))
                continue
    
        elif (xc[i,0]-orig[0])**2+(xc[i,1]-orig[1])**2+(xc[i,2]-orig[2])**2<rr[i]**2 or\
            (xc[i,0]-orig[0])**2+(xc[i,1]-orig[1])**2+(xc[i,2]-L[2])**2<rr[i]**2 or\
            (xc[i,0]-orig[0])**2+(xc[i,1]-L[1])**2+(xc[i,2]-orig[2])**2<rr[i]**2 or\
            (xc[i,0]-orig[0])**2+(xc[i,1]-L[1])**2+(xc[i,2]-L[2])**2<rr[i]**2 or\
            (xc[i,0]-L[0])**2+(xc[i,1]-orig[1])**2+(xc[i,2]-orig[2])**2<rr[i]**2 or\
            (xc[i,0]-L[0])**2+(xc[i,1]-orig[1])**2+(xc[i,2]-L[2])**2<rr[i]**2 or\
            (xc[i,0]-L[0])**2+(xc[i,1]-L[1])**2+(xc[i,2]-orig[2])**2<rr[i]**2 or\
            (xc[i,0]-L[0])**2+(xc[i,1]-L[1])**2+(xc[i,2]-L[2])**2<rr[i]**2:
                print(xc[i],'border3',rr[i])
                idex=np.append(idex,int(i))
                continue

    cc=np.array([xc[idex,0],xc[idex,1],xc[idex,2],rr[idex]]).transpose()
    return cc
    
def get_curves(L,cc, micro_type):
    if micro_type=='s': cc2=cc
    if micro_type=='b': cc2=cc[0:2]
    if micro_type=='f': cc2=cc[0:3]
    pp=np.array([[0,0,1,-0.5],
                 [0,1,0,-0.5],
                 [1,0,0,-0.5]])
                 # [0,1,1,-1],
                 # [1,0,1,-1],
                 # [1,1,0,-1],
                 # [1,1,1,-1.5]])
    
    gmsh.initialize()
    
    gmsh.model.add("trial")
    iv=0
    if micro_type=='r' or micro_type=='l': 
        cc2=cc
        gmsh.model.occ.addBox(0, 0, 0, L[0], L[1], L[2])
        gmsh.model.occ.synchronize()
        for i in range(len(cc)):
            gmsh.model.occ.addSphere(cc[i,0],cc[i,1],cc[i,2],cc[i,3])
            gmsh.model.occ.synchronize()
        inde1=[]
        for i in range(len(cc)):
            gmsh.model.occ.intersect([(3,i+2)],[(3,1)],len(cc)+i+2,removeObject=True,removeTool=False) 
            inde1.append((3,len(cc)+i+2))
            gmsh.model.occ.synchronize()   
        gmsh.model.occ.cut([(3,1)],inde1,2*len(cc)+2, removeObject = True, removeTool = False) 
        gmsh.model.occ.synchronize()  
    else:
     list1=[]
     for ii in range(8):
        gmsh.model.occ.addBox(0, 0, 0, L[0]/2, L[1]/2, L[2]/2);
        iv+=1;ibox=iv;gmsh.model.occ.synchronize()
        inde1=[]
        for i in range(len(cc2)):
            gmsh.model.occ.addSphere(cc2[i,0],cc2[i,1],cc2[i,2],cc2[i,3])
            iv+=1;gmsh.model.occ.synchronize()
            inde1.append((3,iv))
        list1.append([])
        for inde11 in inde1:
            # print(inde11)
            if inde11==inde1[-1]: gmsh.model.occ.intersect([inde11],[(3,ibox)],iv+1,removeObject=True,removeTool=True) 
            else: gmsh.model.occ.intersect([inde11],[(3,ibox)],iv+1,removeObject=True,removeTool=False) 
            iv+=1;gmsh.model.occ.synchronize()
        # gmsh.model.occ.cut([(3,ibox)],inde1,iv+1, removeObject = True, removeTool = False) 
        # iv+=1;gmsh.model.occ.synchronize()  
            list1[-1].append((3,iv))
     # print(list1)
     # gmsh.model.occ.symmetrize([(3,22)],pp[0,0],pp[0,1],pp[0,2],pp[0,3])
     for ii in [1,3,5,7]:
        # gmsh.model.occ.symmetrize([list1[ii][0]],pp[0,0],pp[0,1],pp[0,2],pp[0,3])
        # iv+=1;gmsh.model.occ.synchronize()
        # list1[ii][0]=(3,iv)
        gmsh.model.occ.symmetrize(list1[ii],pp[0,0],pp[0,1],pp[0,2],pp[0,3])
        gmsh.model.occ.synchronize()
        # list1[ii][1]=(3,iv)
     for ii in [2,3,6,7]:
        # gmsh.model.occ.symmetrize([list1[ii][0]],pp[1,0],pp[1,1],pp[1,2],pp[1,3])
        # iv+=1;gmsh.model.occ.synchronize()
        # list1[ii][0]=(3,iv)
        gmsh.model.occ.symmetrize(list1[ii],pp[1,0],pp[1,1],pp[1,2],pp[1,3])
        gmsh.model.occ.synchronize()
        # list1[ii][1]=(3,iv)
     for ii in [4,5,6,7]:
        # gmsh.model.occ.symmetrize([list1[ii][0]],pp[2,0],pp[2,1],pp[2,2],pp[2,3])
        # iv+=1;gmsh.model.occ.synchronize()
        # list1[ii][0]=(3,iv)
        gmsh.model.occ.symmetrize(list1[ii],pp[2,0],pp[2,1],pp[2,2],pp[2,3])
        gmsh.model.occ.synchronize()
        # list1[ii][1]=(3,iv)
     listb=[]
     for iix in range(2):
      for iiy in range(2):
       for iiz in range(2):
        gmsh.model.occ.addBox(0+iix*L[0]/2,0+iiy*L[1]/2,0+iiz*L[2]/2,L[0]/2,L[1]/2,L[2]/2);
        gmsh.model.occ.synchronize()
        iv+=1;listb.append((3,iv))
     for ii in range(8):
        gmsh.model.occ.cut([listb[ii]],list1[ii],iv+1, removeObject = True, removeTool = False) 
        gmsh.model.occ.synchronize()
        iv+=1;
        
        
    gmsh.model.occ.removeAllDuplicates()
    gmsh.model.occ.synchronize()
    # for ii in range(8):
    #     print(ii)
    #     for i in range(len(cc2)):
    #         gmsh.model.occ.addSphere(cc2[i,0],cc2[i,1],cc2[i,2],cc2[i,3])
    #         # list1.append((3,int(i*8+ii+2)))
    #         # gmsh.model.occ.synchronize()

    #     gmsh.model.occ.synchronize()
    #     inde1=[]
    #     print(cc2)
    #     for i in range(len(cc2)):
    #         gmsh.model.occ.intersect([(3,2*len(cc2)*ii+i+2)],[(3,1)],2*len(cc2)*ii+len(cc2)+i+2,removeObject=True,removeTool=False) 
    #         inde1.append((3,2*len(cc2)*ii+len(cc2)+i+2))
            
    #     gmsh.model.occ.synchronize()   
    #     list1.append(inde1)
    # aainde=2*len(cc2)*ii+len(cc2)+i+2
    # gmsh.model.occ.cut([(3,1)],inde1,2*len(cc2)+2, removeObject = True, removeTool = True) 
    # gmsh.model.occ.synchronize()  
    
    
    # for ii in range(3):
    #     gmsh.model.occ.symmetrize(list1[ii],pp[ii,0],pp[ii,1],pp[ii,2],pp[ii,3])
    # # gmsh.model.occ.symmetrize([(3,3)],1,0,0,-0.5)
    # gmsh.model.occ.synchronize()  
 
    # gmsh.model.occ.addBox(0, 0, 0, L[0], L[1], L[2])
    # gmsh.model.occ.synchronize()   
    # gmsh.model.occ.synchronize()  
 
    # nsurf=np.array(gmsh.model.getEntities(2))[:-6,1]
    ncurv=np.array(gmsh.model.getEntities(1))[:,1]
    eps = 1e-6;xlist=np.array([]);
    xlist = np.append(xlist,np.array(gmsh.model.occ.getEntitiesInBoundingBox(-eps,-eps,-eps,L[0]+eps,+eps,+eps,1)).reshape(-1)[1::1])
    xlist = np.append(xlist,np.array(gmsh.model.occ.getEntitiesInBoundingBox(-eps,-eps,L[2]-eps,L[0]+eps,+eps,L[2]+eps,1)).reshape(-1)[1::1])
    xlist = np.append(xlist,np.array(gmsh.model.occ.getEntitiesInBoundingBox(-eps,L[1]-eps,-eps,L[0]+eps,L[1]+eps,+eps,1)).reshape(-1)[1::1])
    xlist = np.append(xlist,np.array(gmsh.model.occ.getEntitiesInBoundingBox(-eps,L[1]-eps,L[2]-eps,L[0]+eps,L[1]+eps,L[2]+eps,1)).reshape(-1)[1::1])
    xlist = np.append(xlist,np.array(gmsh.model.occ.getEntitiesInBoundingBox(-eps,-eps,-eps,+eps,L[1]+eps,+eps,1)).reshape(-1)[1::1])
    xlist = np.append(xlist,np.array(gmsh.model.occ.getEntitiesInBoundingBox(-eps,-eps,L[2]-eps,+eps,L[1]+eps,L[2]+eps,1)).reshape(-1)[1::1])
    xlist = np.append(xlist,np.array(gmsh.model.occ.getEntitiesInBoundingBox(L[0]-eps,-eps,-eps,L[0]+eps,L[1]+eps,+eps,1)).reshape(-1)[1::1])
    xlist = np.append(xlist,np.array(gmsh.model.occ.getEntitiesInBoundingBox(L[0]-eps,-eps,L[2]-eps,L[0]+eps,L[1]+eps,L[2]+eps,1)).reshape(-1)[1::1])
    xlist = np.append(xlist,np.array(gmsh.model.occ.getEntitiesInBoundingBox(-eps,-eps,-eps,+eps,+eps,L[2]+eps,1)).reshape(-1)[1::1])
    xlist = np.append(xlist,np.array(gmsh.model.occ.getEntitiesInBoundingBox(L[0]-eps,-eps,-eps,L[0]+eps,+eps,L[2]+eps,1)).reshape(-1)[1::1])
    xlist = np.append(xlist,np.array(gmsh.model.occ.getEntitiesInBoundingBox(-eps,L[1]-eps,-eps,+eps,L[1]+eps,L[2]+eps,1)).reshape(-1)[1::1])
    xlist = np.append(xlist,np.array(gmsh.model.occ.getEntitiesInBoundingBox(L[0]-eps,L[1]-eps,-eps,L[0]+eps,L[1]+eps,L[2]+eps,1)).reshape(-1)[1::1])
    
    xlist  = np.append(xlist,np.array(gmsh.model.occ.getEntitiesInBoundingBox(-eps,L[1]/2-eps,L[2]/2-eps,L[0]+eps,L[1]/2+eps,L[2]/2+eps,1)).reshape(-1)[1::1])
    xlist  = np.append(xlist,np.array(gmsh.model.occ.getEntitiesInBoundingBox(L[0]/2-eps,-eps,L[2]/2-eps,L[0]/2+eps,L[1]+eps,L[2]/2+eps,1)).reshape(-1)[1::1])
    xlist  = np.append(xlist,np.array(gmsh.model.occ.getEntitiesInBoundingBox(L[0]/2-eps,L[1]/2-eps,-eps,L[0]/2+eps,L[1]/2+eps,L[2]+eps,1)).reshape(-1)[1::1])
    
    xlist  = np.append(xlist,np.array(gmsh.model.occ.getEntitiesInBoundingBox(-eps,L[1]/2-eps,-eps,+eps,L[1]/2+eps,L[2]+eps,1)).reshape(-1)[1::1])
    xlist  = np.append(xlist,np.array(gmsh.model.occ.getEntitiesInBoundingBox(-eps,-eps,L[2]/2-eps,+eps,L[1]+eps,L[2]/2+eps,1)).reshape(-1)[1::1])
    xlist  = np.append(xlist,np.array(gmsh.model.occ.getEntitiesInBoundingBox(L[0]-eps,L[1]/2-eps,-eps,L[0]+eps,L[1]/2+eps,L[2]+eps,1)).reshape(-1)[1::1])
    xlist  = np.append(xlist,np.array(gmsh.model.occ.getEntitiesInBoundingBox(L[0]-eps,-eps,L[2]/2-eps,L[0]+eps,L[1]+eps,L[2]/2+eps,1)).reshape(-1)[1::1])
    
    xlist  = np.append(xlist,np.array(gmsh.model.occ.getEntitiesInBoundingBox(-eps,-eps,L[2]/2-eps,L[0]+eps,+eps,L[2]/2+eps,1)).reshape(-1)[1::1])
    xlist  = np.append(xlist,np.array(gmsh.model.occ.getEntitiesInBoundingBox(L[0]/2-eps,-eps,-eps,L[0]/2+eps,+eps,L[2]+eps,1)).reshape(-1)[1::1])
    xlist  = np.append(xlist,np.array(gmsh.model.occ.getEntitiesInBoundingBox(-eps,L[1]-eps,-eps,L[0]+eps,L[1]+eps,L[2]/2+eps,1)).reshape(-1)[1::1])
    xlist  = np.append(xlist,np.array(gmsh.model.occ.getEntitiesInBoundingBox(L[0]/2-eps,L[1]-eps,-eps,L[0]/2+eps,L[1]+eps,L[2]+eps,1)).reshape(-1)[1::1])
    
    xlist  = np.append(xlist,np.array(gmsh.model.occ.getEntitiesInBoundingBox(-eps,L[1]/2-eps,-eps,L[0]+eps,L[1]/2+eps,+eps,1)).reshape(-1)[1::1])
    xlist  = np.append(xlist,np.array(gmsh.model.occ.getEntitiesInBoundingBox(L[0]/2-eps,-eps,-eps,L[0]/2+eps,L[1]+eps,+eps,1)).reshape(-1)[1::1])
    xlist  = np.append(xlist,np.array(gmsh.model.occ.getEntitiesInBoundingBox(-eps,L[1]/2-eps,L[2]-eps,L[0]+eps,L[1]/2+eps,L[2]+eps,1)).reshape(-1)[1::1])
    xlist  = np.append(xlist,np.array(gmsh.model.occ.getEntitiesInBoundingBox(L[0]/2-eps,-eps,L[2]-eps,L[0]/2+eps,L[1]+eps,L[2]+eps,1)).reshape(-1)[1::1])
    
    
    
    gmsh.model.occ.synchronize()  
    # gmsh.write('aa.geo_unrolled')
    gmsh.finalize()
    # print(ncurv)
    return ncurv,xlist

def write_files(filename,maxsize,micro_type,clen,L,cc,ncurv,xlist):
    order=1
    if micro_type=='s': cc2=cc
    if micro_type=='b': cc2=cc[0:2]
    if micro_type=='f': cc2=cc[0:3]

    Lattice = open (filename+".geo",'w')
    iv=1
    list1=[]
    Lattice.write( "SetFactory('OpenCASCADE');" + "\n" )
    if micro_type=='r' or micro_type=='l': 
        cc2=cc
        Lattice.write("Box(1) = {0,0,0,"+str(L[0])+","+str(L[1])+","+str(L[2])+"};"+ "\n" )
        for i in range(len(cc)):
            Lattice.write("Sphere("+str(i+2)+") = {"+ str(cc[i,0])+","+ str(cc[i,1])+","+ str(cc[i,2])+","+\
                          str(cc[i,3])+",-Pi/2,Pi/2,2*Pi};"+"\n")
        for i in range(len(cc)):
            Lattice.write("BooleanIntersection("+str(i+len(cc)+2)+\
             ") = {Volume{"+str(i+2)+"};Delete;}{Volume{1};};\n")
        Lattice.write("BooleanDifference("+str(2*len(cc)+2)+") = {Volume{1};Delete;}{")
        for i in range(len(cc)):
            Lattice.write("Volume{"+str(i+len(cc)+2)+"};")
        Lattice.write("};\n")
        for ii in range(8): list1.append([])
        for i in range(len(cc)): list1[-1].append(i+len(cc)+2)
        listc=[2*len(cc)+2]

    else:
     for ii in range(8):
        Lattice.write( "Box("+str(iv)+") = {0,0,0,"+str(L[0]/2)+","+str(L[1]/2)+","+str(L[2]/2)+"};"+ "\n" )
        ibox=iv
        iv+=1
        inde1=[]
        for i in range(len(cc2)):
            Lattice.write("Sphere("+str(iv)+") = {"+ str(cc2[i,0])+","+ str(cc2[i,1])+","+ str(cc2[i,2])+","+\
                          str(cc2[i,3])+",-Pi/2,Pi/2,2*Pi};"+"\n")
            inde1.append(iv)
            iv+=1
            #Lattice.write("Characteristic Length{ PointsOf{ Curve{"+str(12+3*i+1)+"}; } } = "+str(clen)+";"+"\n")
            #Lattice.write("Characteristic Length{ PointsOf{ Curve{"+str(12+3*i+2)+"}; } } = "+str(clen)+";"+"\n")
            #Lattice.write("Characteristic Length{ PointsOf{ Curve{"+str(12+3*i+3)+"}; } } = "+str(clen)+";"+"\n")
        list1.append([])
        for inde11 in inde1:
            # print(inde11)
            Lattice.write("BooleanIntersection("+str(iv)+") = {")
            Lattice.write("Volume{"+str(inde11)+"};Delete;")
            if inde11==inde1[-1]: Lattice.write("}{Volume{"+str(ibox)+"};Delete;};\n")
            else: Lattice.write("}{Volume{"+str(ibox)+"};};\n")
            list1[-1].append(iv)
            iv+=1
     pp=np.array([[0,0,1,-0.5],
                 [0,1,0,-0.5],
                 [1,0,0,-0.5]])
     for ii in [1,3,5,7]:
        for list11 in list1[ii]:
            Lattice.write("Symmetry{"+str(pp[0,0])+","+str(pp[0,1])+","+str(pp[0,2])+","+str(pp[0,3])+\
                      "}{ Volume{"+str(list11)+"};}\n")

     for ii in [2,3,6,7]:
        for list11 in list1[ii]:
            Lattice.write("Symmetry{"+str(pp[1,0])+","+str(pp[1,1])+","+str(pp[1,2])+","+str(pp[1,3])+\
                      "}{ Volume{"+str(list11)+"};}\n")
        
     for ii in [4,5,6,7]:
        for list11 in list1[ii]:
            Lattice.write("Symmetry{"+str(pp[2,0])+","+str(pp[2,1])+","+str(pp[2,2])+","+str(pp[2,3])+\
                      "}{ Volume{"+str(list11)+"};}\n")
        
     listb=[]
     for iix in range(2):
      for iiy in range(2):
       for iiz in range(2):
        Lattice.write("Box("+str(iv)+") = {"+str(iix*L[0]/2)+","+str(iiy*L[1]/2)+","+str(iiz*L[2]/2)+","+str(L[0]/2)+","+str(L[1]/2)+","+str(L[2]/2)+"};"+ "\n") 
        listb.append(iv)
        iv+=1;

     listc=[]
     for ii in range(8):
        Lattice.write("BooleanDifference("+str(iv)+") = {")
        Lattice.write("Volume{"+str(listb[ii])+"};Delete;}{")
        for list11 in list1[ii]:
            Lattice.write("Volume{"+str(list11)+"};")
        Lattice.write("};\n")
        listc.append(iv)
        iv+=1;
    
    Lattice.write("Coherence;\n")

    # Lattice.write("BooleanDifference("+str(iv)+") = {")
    # for listb1 in listb:
    #     Lattice.write("Volume{"+str(listb1)+"};Delete;")
    # Lattice.write("}{")
    # for list11 in list1:
    #     Lattice.write("Volume{"+str(list11)+"};;")
    # Lattice.write("};\n")
    
    # for i in nsurf:
    #     Lattice.write("Characteristic Length{ PointsOf{ Surface{"+str(i)+"}; } } = "+str(clen)+";"+"\n")
    
    for i in ncurv:
        # if not i in np.array(x00).reshape(-1)[1::1] and\
        # not i in np.array(x01).reshape(-1)[1::1] and\
        # not i in np.array(x10).reshape(-1)[1::1] and\
        # not i in np.array(x11).reshape(-1)[1::1] and\
        # not i in np.array(y00).reshape(-1)[1::1] and\
        # not i in np.array(y01).reshape(-1)[1::1] and\
        # not i in np.array(y10).reshape(-1)[1::1] and\
        # not i in np.array(y11).reshape(-1)[1::1] and\
        # not i in np.array(z00).reshape(-1)[1::1] and\
        # not i in np.array(z01).reshape(-1)[1::1] and\
        # not i in np.array(z10).reshape(-1)[1::1] and\
        # not i in np.array(z11).reshape(-1)[1::1] and\
        # not i in np.array(xm).reshape(-1)[1::1] and\
        # not i in np.array(ym).reshape(-1)[1::1] and\
        # not i in np.array(zm).reshape(-1)[1::1]:
        if not i in xlist:
              Lattice.write("Characteristic Length{ PointsOf{ Curve{"+str(i)+"}; } } = "+str(clen)+";"+"\n")
    # for i in ncurv:
    #     Lattice.write("Characteristic Length{ PointsOf{ Curve{"+str(i)+"}; } } = "+str(clen)+";"+"\n")
             
    
    
    Lattice.write("Physical Volume(\"0\") = {")
    for listc1 in listc:
        Lattice.write(str(listc1))
        if listc1 == listc[-1]: Lattice.write("")
        else: Lattice.write(",")
    Lattice.write("};\n")
    Lattice.write("Physical Volume(\"1\") = {")
    for ii in range(8):
      for list11 in list1[ii]:
        Lattice.write(str(list11))
        if list11 == list1[ii][-1] and ii==7: Lattice.write("")
        else: Lattice.write(",")
    Lattice.write("};\n")
    Lattice.write ("e=1e-6;" +"\n" +
      "bottomx() = Surface In BoundingBox{-e,-e,-e, e,"+str(L[1])+"+e,"+str(L[2])+"+e};"+"\n" +
      "bottomy() = Surface In BoundingBox{-e,-e,-e, "+str(L[0])+"+e,e,"+str(L[2])+"+e};"+"\n" +
      "bottomz() = Surface In BoundingBox{-e,-e,-e, "+str(L[0])+"+e,"+str(L[1])+"+e,e};"+"\n" +
      "topx() = Surface In BoundingBox{"+str(L[0])+"-e,-e,-e,"+str(L[0])+" +e,"+str(L[1])+"+e,"+str(L[2])+"+e};"+"\n" +
      "topy() = Surface In BoundingBox{-e,"+str(L[1])+"-e,-e, "+str(L[0])+"+e,"+str(L[1])+"+e,"+str(L[2])+"+e};"+"\n" +
      "topz() = Surface In BoundingBox{-e,-e,"+str(L[2])+"-e, "+str(L[0])+"+e,"+str(L[1])+"+e,"+str(L[2])+"+e};"+"\n" +
      "Periodic Surface{topx()} = {bottomx()} Translate{"+str(L[0])+", 0, 0};"+"\n" +
      "Periodic Surface{topy()} = {bottomy()} Translate{0,"+str(L[1])+" ,0};"+"\n" +
      "Periodic Surface{topz()} = {bottomz()} Translate{0,0,"+str(L[2])+"};")
    Lattice.write("Coherence;\n")
    Lattice.close()
    
    
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.open(filename+'.geo')
    gmsh.option.setNumber("Mesh.ElementOrder", order)#### order quadratic
    # gmsh.option.setNumber("Mesh.ElementOrder", 1)
    gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", maxsize)
    #gmsh.model.mesh.setSize(0.1)
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.mesh.generate(3)
    # gmsh.model.mesh.removeDuplicateNodes()
    
    gmsh.plugin.setNumber("MeshVolume", "Dimension", 3)
    gmsh.plugin.setNumber("MeshVolume", "PhysicalGroup", 2)
    gmsh.plugin.run("MeshVolume")
    
    gmsh.write(filename+'.msh')
    gmsh.write(filename+'.geo.opt')
    gmsh.write(filename+'.inp')
    gmsh.finalize()
    
    
    
    
    import meshio
    #msh = meshio.read(filename+'.msh')
    #meshio.write(filename+'.xdmf',msh)
    
    
    ## original mesh file
    msh = meshio.read(filename+'.msh')
    ## physical surface & volume data
    for key in msh.cell_data_dict["gmsh:physical"].keys():
        if order==1:
            if key == "tetra":
                tetra_data = msh.cell_data_dict["gmsh:physical"][key]
        if order==2:
            if key == "tetra10":
                tetra_data = msh.cell_data_dict["gmsh:physical"][key]
    ## cell data
    tetra_cells = np.array([None])
    for cell in msh.cells:
        if order==1:
            if cell.type == "tetra":
                if tetra_cells.all() == None:
                    tetra_cells = cell.data
                else:
                    tetra_cells = np.concatenate((tetra_cells,cell.data))
        if order==2:
            if cell.type == "tetra10":
                if tetra_cells.all() == None:
                    tetra_cells = cell.data
                else:
                    tetra_cells = np.concatenate((tetra_cells,cell.data))
    ## put them together
    if order==1:
        tetra_mesh = meshio.Mesh(points=msh.points, 
                                 cells={"tetra": tetra_cells},
                                 cell_data={"phase":[tetra_data]})
    if order==2:
        tetra_mesh = meshio.Mesh(points=msh.points, 
                                 cells={"tetra10": tetra_cells},
                                 cell_data={"phase":[tetra_data]})
    ## output
    meshio.write(filename+'.xdmf', tetra_mesh)
        # meshio.write("mf.xdmf", triangle_mesh)
    
    # convertMesh(filename+'.msh',filename+'.xdmf')
    return









# filenames=['s','b','f']
filenames=['r']
# filenamesv=['10','20','30']
filenamesv=['10']
# vfs=[[0.1,0.2,0.3],[0.1,0.2,0.3],[0.1,0.2,0.3],[0.1,0.2,0.3]]
vfs=[[0.1]]
# nparts=[1,2,3,20]
nparts=[20]
# packs=[0.5,0.72,0.92] # maximums rand dist
# packs=[0.82,0.82,0.82]
packs=[0.82]





for istru in range(len(filenames)):
    for ivf in range(len(vfs[istru])):
        filename=filenames[istru]+filenamesv[ivf]
        vf=vfs[istru][ivf]
        npart=nparts[istru]
        rad=(vf/npart*3/4/np.pi)**(1./3.)
        packing=packs[ivf]

        L=np.array([1.,1.,1.])
        
        clen=rad/2
        # 0.065 more or less good for r10 / 0.12 for f30
        # 0.08 more or less good for r20 / 0.14 for f30
        # 0.1 more or less good for r30 / 0.15 for f30
        # 0.2 more or less good rest
        maxsize=0.2 # maximum element size
        # 0.2 more or less good all      
        
        xc,rr=locate(npart,rad,packing,filenames[istru])
        cc=select(xc,rr,filenames[istru])
        if filenames[istru]=='l': cc=np.loadtxt(filename+'cc.txt')
        else: np.savetxt(filename+'cc.txt', cc)
        print('here',filename)
        
        
        ncurv,xlist=get_curves(L,cc,filenames[istru])
        write_files(filename, maxsize,filenames[istru], clen, L, cc, ncurv, xlist)




# filename='b30'
# vf=0.3
# npart=2
# rad=(vf/npart*3/4/np.pi)**(1./3.)
# packing=0.82 # 1-tocuh 0-infinite distance
# # 0.45 for r10
# # 0.65 for r20
# # 0.82 for r30
# micro_type='bcc' #sc fcc bcc random load
# L=np.array([1.,1.,1.])
# clen=0.15 # characteristic length for spheres
# # 0.065 more or less good for r10 / 0.12 for f30
# # 0.08 more or less good for r20 / 0.14 for f30
# # 0.1 more or less good for r30 / 0.15 for f30
# # 0.2 more or less good rest
# maxsize=0.2 # maximum element size
# # 0.2 more or less good all














