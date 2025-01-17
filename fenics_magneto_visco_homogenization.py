#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    fenics magnetomechanical homogenization: 3D full-field magnetomechanical homogenization fenics
     code: finite strains, viscosity, incompressibility and mixed stress/flux controlled
    Copyright (C) 2021  Sergio Lucarini

    If using this code for research or industrial purposes, please cite:
    S. Lucarini, M.A. Moreno-Mateos, K. Danas, D. Garcia-Gonzalez
    Insights into the viscohyperelastic response of soft magnetorheological
    elastomers: competition of macrostructural versus microstructural players
    International Journal of Solids and Structures, Under review, 2022
"""

import numpy as np
import fenics as fn
from aux_homog_nl import importMesh,PeriodicBoundary,post,adjust_xdmf


# this is the input of the simulations


# tiemr is a list of steps like abaqus, where each row of the timer is a verctor of:
# [initial time increment, final time of step, minimum time increment, maximum time increment]
timer=[]
t_ramp=0.15;
timer.append([t_ramp/10,t_ramp,t_ramp/100,2*t_ramp/10]) #ramp of h, remove if
timer.append([t_ramp/10,t_ramp+100*t_ramp,t_ramp/200,20*t_ramp]) #hold

# final magnetic ramp
b_end=0.015

# name of the output folder
filename='trial'

# filename of the input xdmf mesh
model='s30.xdmf'
    
# for mechanical stress control, a tensor where comopnents =1 are stress controlled / =0 strain controlled
control=[[1, 0, 0], [1, 1, 0], [0,0,0]]

# fenics expression for the imposition of the macroscopic deformation gradient if control is 0 in that component
fmacro=fn.Expression((("0","0","0"),\
              ("0","0"," 0 "),\
              ("0","0","0")), degree=0, t = 0, t_ramp=t_ramp)
              
# fenics expression for the imposition of the macroscopic Eulerian magnetic density flux              
b_eul=fn.Expression((("0.00","0.00","t<=t_ramp ? "+str(b_end/t_ramp)+\
             "*t : "+str(b_end))), degree=0, t = 0, t_ramp=t_ramp)

            
# materials properties mechanical: [phase1:(shear modulus, bulk modulus = 0, shear modulus of the viscous part, viscosity parameter tau), phase2:(the same)]
me = [(1.034e-3,0,1.515e-3,0.217e-3),(81.78e3,0,0,1)]
# materials properties magnetic: [phase1:(relative magnetic permeability, magnetic saturation), phase2:(the same)]
ma = [(1,0),(30,2)]

# activate and deactivate magnetic field, incompresivility, and viscosity
magnetoactive=True
incompressibility=True
viscosity=True


#run the simulation: (output folder,input xdmf mesh,control tensor,macroscopic deformation gradient,macroscopic Eulerian magnetic density flux,
#  activation of the magnetic field, activation of the incompressibility condition, activation of the viscous part,
#  mechanical properties, magnetic properities)
 
#tolerances and optimization
maxiter_nw=7;tol_nw=5e-3;tol_nw_forces=1e-5
fn.parameters["form_compiler"]["optimize"] = True
fn.parameters["form_compiler"]["cpp_optimize"] = True

# set steps [initial dtime, total time end, dtime min, dtime max]
timer=np.array(timer)

# initialize incerment of time
DT_macro = fn.Expression("dt", dt=0, degree=0)

# boundary conditions # control P or F -> 0-F / 1-P
control_mech_macro=np.array(control)
control_mech_index=np.zeros([3,3]);
control_mech_index[np.where(control_mech_macro==1)[0],\
 np.where(control_mech_macro==1)[1]]=np.arange(np.sum(control_mech_macro))

# boundary conditions # imposing F macro
F_macro=fmacro

# boundary conditions # imposing P macro
P_macro=fn.Expression((("0","0","0"),\
 ("0","0","0"),("0","0","0")), degree=0, t = 0)

# boundary conditions # control B or H -> 0-H / 1-B
control_mag_macro=np.array([1, 1, 1])
if not magnetoactive: control_mag_macro=np.array([0, 0, 0])
control_mag_index=np.zeros([3]);
control_mag_index[np.where(control_mag_macro==1)]=\
 np.arange(np.sum(control_mag_macro))

# boundary conditions # imposing H macro
H_macro=fn.Expression((("0","0","0.0001")), degree=0, t = 0) 

# boundary conditions # imposing B macro
B_macro=b_eul

# mesh reading and preparation
mesh,materials,mvc = importMesh(model)
domains = fn.MeshFunction("size_t", mesh,dim=3)
domains.set_all(0)
domains.array()[:]=materials.array()[:]-1
vertices = np.array([[0, 0.,0],[1, 0.,0],[1, 1.,0],[0, 1.,0],\
                     [0, 0.,1],[1, 0.,1],[1, 1.,1],[0, 1.,1]])
dx = fn.Measure('dx')(domain=mesh,subdomain_data=domains)
dx = dx(metadata={'quadrature_degree': 2})

#material parameters ( shear bulk shearvisco visco )
material_elas = me

#material magnetic param ( mu_r m_s )
material_magnet = ma

# derivatives  of the energy functional Psi, constituvie relations

# dPsi/dF
def P(F,p, i, H,FV):
    shear,bulk,shearv,visco = material_elas[i];mu,ms = material_magnet[i];chi=mu-1
    Finv=fn.inv(F);FinvT=Finv.T;J=fn.det(F)
    FE = fn.dot(F,fn.inv(FV))
    I5=fn.dot(fn.dot(H,Finv),fn.dot(H,Finv));I5r=fn.sqrt(I5)
    dI5dF=-2*fn.outer(fn.dot(FinvT,H),fn.dot(fn.dot(H,Finv),FinvT))
    return shear*F-p*J*FinvT+\
      fn.dot(shearv*FE,fn.inv(FV))-\
      J*dI5dF*4e-1*np.pi/2-\
      J*I5*4e-1*np.pi/2*FinvT-\
      J*4e-1*np.pi*ms/2/(I5r+fn.DOLFIN_EPS)*(-1+2/(fn.exp(-2*chi/(ms+fn.DOLFIN_EPS)*I5r)+1))*dI5dF-\
      J*4e-1*np.pi*FinvT*(ms**2/(chi+fn.DOLFIN_EPS)*fn.ln(0.5*(fn.exp(chi/(ms+fn.DOLFIN_EPS)*I5r)+fn.exp(-chi/(ms+fn.DOLFIN_EPS)*I5r))))

# dPsi/dF_mech
def P0(F,p, i,FV):
    shear,bulk,shearv,visco = material_elas[i]
    Finv=fn.inv(F);FinvT=Finv.T;J=fn.det(F)
    FE = fn.dot(F,fn.inv(FV))
    return shear*F-p*J*FinvT+\
      fn.dot(shearv*FE,fn.inv(FV))

# dPsi/dH  
def B(H, i, F):
    mu,ms = material_magnet[i];chi=mu-1
    Finv=fn.inv(F);FinvT=Finv.T;J=fn.det(F)
    I5r=fn.sqrt(fn.dot(fn.dot(H,Finv),fn.dot(H,Finv)))
    dI5dH=2*fn.dot(fn.dot(FinvT,H),FinvT)
    return J*4e-1*np.pi*(1/2+ms/2/(I5r+fn.DOLFIN_EPS)*(-1+2/(fn.exp(-2*chi/(ms+fn.DOLFIN_EPS)*I5r)+1)))*dI5dH

# viscous flow rule
def compare(F, i, FV,DT):
    shear,bulk,shearv,visco = material_elas[i]
    FE = fn.dot(F,fn.inv(FV));J=fn.det(F)
    # Viscous flow
    # Pnm1 = fn.dot(shearv*FE-shearv*fn.inv(FE).T,fn.inv(FV))
    Pnm1 = fn.dot(shearv*FE,fn.inv(FV))
    sigman = fn.dot(Pnm1,F.T)/J
    sigmaDev = sigman - 1./3.*fn.tr(sigman)*fn.Identity(3)
    DV = 1/2**0.5/visco*sigmaDev
    CompareFV = DT*fn.dot(fn.dot(fn.inv(FE),DV),F)
    return fn.as_vector((CompareFV[0,0], CompareFV[1,1], CompareFV[2,2], CompareFV[0,1]+\
     CompareFV[1,0],CompareFV[1,2]+CompareFV[2,1], CompareFV[0,2]+CompareFV[2,0])) 

#defining functional space
all_elements=[]
# displacement finite element
Vu = fn.VectorElement('CG', mesh.ufl_cell(),2, dim=3);all_elements.append(Vu)
# magnetic potential finite element
if magnetoactive: Vv=fn.FiniteElement('CG', mesh.ufl_cell(), 2);all_elements.append(Vv)
# presure potential finite element
if incompressibility: Vpr=fn.FiniteElement('DG', mesh.ufl_cell(), 0);all_elements.append(Vpr)
# presure potential finite element
if viscosity: VFV=fn.VectorElement('DG', mesh.ufl_cell(), 1,dim=6);all_elements.append(VFV)
# lagrange multiplyer constant element
if np.sum(control_mech_macro)>0: Vlambda = fn.VectorElement('R',\
  mesh.ufl_cell(), 0, int(np.sum(control_mech_macro)));all_elements.append(Vlambda)
# lagrange multiplyer constant element
if magnetoactive and np.sum(control_mag_macro)>0: Vlambdav=fn.VectorElement('R',\
  mesh.ufl_cell(), 0, int(np.sum(control_mag_macro)));all_elements.append(Vlambdav)
 
# mixed formulation with periodic pbcs
V = fn.FunctionSpace(mesh, fn.MixedElement(all_elements),\
 constrained_domain=PeriodicBoundary(vertices))
  
# test, trials and solutions
W_ = fn.TestFunction(V);all_test=fn.split(W_);
dW = fn.TrialFunction(V);
W_sol = fn.Function(V);all_func=fn.split(W_sol);W_sol_t = fn.Function(V);
deltaW_sol = fn.Function(V);diffW = fn.Function(V);diffW2 = fn.Function(V); #solution field in previous iteration
iiii=0;u_=all_test[iiii];u_sol=all_func[iiii];
if magnetoactive: iiii+=1;v_=all_test[iiii];v_sol=all_func[iiii]; 
else: v_sol=fn.Constant(0)
if incompressibility: iiii+=1;pr_=all_test[iiii];pr_sol=all_func[iiii];
else: pr_sol=fn.Constant(0)
if viscosity: iiii+=1;FV_=all_test[iiii];FV_sol=all_func[iiii];ifv=1*iiii
else: FV_sol=fn.Constant((0,0,0,0,0,0))
if np.sum(control_mech_macro)>0:
    iiii+=1;lambda_=all_test[iiii];lambdas=all_func[iiii]
if magnetoactive and np.sum(control_mag_macro)>0: 
    iiii+=1;lambdav_=all_test[iiii];lambdavs=all_func[iiii]

# get the dofs of each field
dofindex=[]
for ai in range(iiii+1):
    dofindex.append(V.sub(ai).dofmap().dofs())
    
# state variable 
all_func_t=fn.split(W_sol_t);
if viscosity: FV_sol_t=all_func_t[ifv] 
else: FV_sol_t=fn.Constant((0,0,0,0,0,0))


# form vectors and tensors for the Lagrange miltipliers in the variational
lambda_tensor = fn.as_tensor([[ lambdas[int(control_mech_index[i,j])]\
 if control_mech_macro[i,j]==1 else 0 for i in range(3)] for j in range(3) ])
lambda__tensor= fn.as_tensor([[ lambda_[int(control_mech_index[i,j])]\
 if control_mech_macro[i,j]==1 else 0 for i in range(3)] for j in range(3) ])
lambdav_vector = fn.as_vector([lambdavs[int(control_mag_index[i])]\
  if control_mag_macro[i]==1 else 0 for i in range(3)])
lambdav__vector = fn.as_vector([lambdav_[int(control_mag_index[i])]\
  if control_mag_macro[i]==1 else 0 for i in range(3)])

# form vectors and tensors for the solution fields
FVV_sol=fn.as_tensor(((FV_sol[0], FV_sol[3]/2, FV_sol[5]/2),(FV_sol[3]/2, FV_sol[1],\
                        FV_sol[4]/2),(FV_sol[5]/2, FV_sol[4]/2, FV_sol[2])))+ fn.Identity(3)
FF_sol=F_macro+lambda_tensor+fn.Identity(3)+fn.grad(u_sol)
H_macro2=H_macro
B_macro2=B_macro
if magnetoactive: 
    FF_macro=F_macro+lambda_tensor+fn.Identity(3)
    HH_sol=H_macro2+lambdav_vector-fn.grad(v_sol)
    B_macro2=fn.det(FF_macro)*fn.dot(fn.inv(FF_macro),B_macro)
else: HH_sol=fn.Constant((0,0,0))

#weak form 
form =  0    
fint=[]

#mechanical energy including F macro and F due to P macro
if magnetoactive: 
    Ptrial_0=P(FF_sol,pr_sol,0,HH_sol,FVV_sol)
    Ptrial_1=P(FF_sol,pr_sol,1,HH_sol,FVV_sol)
else:
    Ptrial_0=P0(FF_sol,pr_sol,0,FVV_sol)
    Ptrial_1=P0(FF_sol,pr_sol,1,FVV_sol)
Ptrial0_0=P0(FF_sol,pr_sol,0,FVV_sol)
Ptrial0_1=P0(FF_sol,pr_sol,1,FVV_sol)
mech_form0 = fn.inner(Ptrial_0,fn.grad(u_))*dx(0)
mech_form1 = fn.inner(Ptrial_1,fn.grad(u_))*dx(1)
form += mech_form0 + mech_form1 
fint.append(mech_form0 + mech_form1)

#magnetic energy including H macro and H dute to B macro
if magnetoactive: 
    Btrial_0=B(HH_sol,0,FF_sol)
    Btrial_1=B(HH_sol,1,FF_sol)
    mag_form0 = -fn.dot(Btrial_0,fn.grad(v_))*dx(0)
    mag_form1 = -fn.dot(Btrial_1,fn.grad(v_))*dx(1)
    form += mag_form0 + mag_form1
    fint.append(mag_form0 + mag_form1)

#incompressibility condition
if incompressibility:
    incompr_form0 = - pr_*(fn.det(FF_sol) - 1)*dx(0)
    incompr_form1 = - pr_*(fn.det(FF_sol) - 1)*dx(1)
    form += incompr_form0 + incompr_form1
    fint.append(incompr_form0 + incompr_form1)

# viscous residual
if viscosity:
    Rvisco0 = - 1e-3*fn.dot(FV_sol-FV_sol_t-compare(FF_sol,0,FVV_sol,DT_macro),FV_)*dx(0)
    Rvisco1 = - 1e-3*fn.dot(FV_sol-FV_sol_t-compare(FF_sol,1,FVV_sol,DT_macro),FV_)*dx(1)
    form += Rvisco0 + Rvisco1
    fint.append(Rvisco0 + Rvisco1)
    
#imposition of macroscopic P for 2 subdomains
if np.sum(control_mech_macro)>0:
    therm_form0 = - fn.inner(lambda__tensor,(Ptrial0_0-P_macro))*dx(0)
    therm_form1 = - fn.inner(lambda__tensor,(Ptrial0_1-P_macro))*dx(1)
    form += therm_form0 + therm_form1
    fint.append(therm_form0 + therm_form1)
    
#imposition of macroscopic B for 2 subdomains
if magnetoactive and np.sum(control_mag_macro)>0:
    therm_formv0 = - fn.inner(lambdav__vector,(Btrial_0-B_macro2))*dx(0)
    therm_formv1 = - fn.inner(lambdav__vector,(Btrial_1-B_macro2))*dx(1)
    form += therm_formv0 + therm_formv1
    fint.append(therm_formv0 + therm_formv1)

# definiton of the global Jacobian
Jac = fn.derivative(form, W_sol, dW)

# fix origin point to 0 displacement and 0 potential to the origin
bcss=[]
def bnd_func(x,on_boundary):
    return fn.near(x[0],0) and fn.near(x[1],0) and fn.near(x[2],0)
bcsu=fn.DirichletBC(V.sub(0),fn.Constant((0.0+fn.DOLFIN_EPS,0.0+fn.DOLFIN_EPS,\
 0.0+fn.DOLFIN_EPS)),bnd_func,method="pointwise"); bcss.append(bcsu)
if magnetoactive: bcsv=fn.DirichletBC(V.sub(1),fn.Constant(0.0+fn.DOLFIN_EPS),\
                 bnd_func,method="pointwise"); bcss.append(bcsv)

# initialize variables and postprocessing
H22=np.array([0]);ttt=np.array([0]);F22=np.array([0]);t=0;dtflag=False;F12=np.array([0])
diffW.vector()[:]=0;P22=np.array([0]);B22=np.array([0]);P12=np.array([0])
file_results=fn.XDMFFile(filename+'/unihshear_out.xdmf')
post(file_results,t,mesh,domains,u_sol,pr_sol,lambda_tensor,\
  F_macro,P0,dx,v_sol,lambdav_vector,H_macro2,B,FVV_sol,FF_sol,HH_sol)
qq_t=np.zeros(iiii+1);qqave=np.zeros(iiii+1);qqave[:]=1e-2;qq_t[:]=1e-2;
qm_t=np.zeros(iiii+1);qmave=np.zeros(iiii+1);qmave[:]=1e-2;qm_t[:]=1e-2;
allinactive=np.zeros(iiii+1,dtype='bool');lowflux=np.zeros(iiii+1,dtype='bool')

# solving sequence
for kstep in range(len(timer)):

  #initialize residuals
  qmcounter=np.zeros(iiii+1);qqcounter=np.zeros(iiii+1)
  qq_t[np.logical_or(qq_t<1e-2,qm_t<1e-2)]=1e-2;qm_t[np.logical_or(qq_t<1e-2,qm_t<1e-2)]=1e-2;
  
  # set timer
  dtime=timer[kstep,0]
  t_end=timer[kstep,1]
  dtold=1*dtime;
  dtmin=timer[kstep,2]
  dtmax=timer[kstep,3]
  diffW2.vector()[:]=0
  
  #time loop
  inc=0;
  while t<t_end-1e-8:
    inc+=1
    
    #update time and bcs
    if t+dtime>t_end: dtime=t_end+fn.DOLFIN_EPS-t
    t+=dtime
    F_macro.t = t;P_macro.t = t;B_macro.t = t;H_macro.t = t;DT_macro.dt=dtime

    #coupled solver
    iter_nw = 0
    deltaW_sol.vector()[:]=1;diffW.vector()[:]=1;error_nwall1=np.ones(iiii+1);error_nwall2=np.ones(iiii+1); #init errors
    while iter_nw<maxiter_nw:
        iter_nw += 1
        
        # assemble linearized system
        A, b = fn.assemble_system(Jac,-form, bcss)
        
        # get residual
        rres=b.get_local()
        
        #check tolerances
        error_nw_old1=1*error_nwall1;error_nw_old2=1*error_nwall2;
        error_nwall1=np.ones(iiii+1);error_nwall2=np.ones(iiii+1);
        allinactive[:]=False;lowflux[:]=False;
        converged1=np.zeros(iiii+1,dtype='bool');converged2=np.zeros(iiii+1,dtype='bool');
        # for each field
        for ai in range(iiii+1):
            # get internal forces and compute relative residuals similar to abaqus documentation for nonlinear problems 
            q0=np.abs(fn.assemble(fint[ai]).get_local()[dofindex[ai]])
            qm = np.max(q0)
            if qm>tol_nw_forces*qm_t[ai]:
                qmave[ai]=qm_t[ai]*(qmcounter[ai])/(qmcounter[ai]+1)+qm*1/(qmcounter[ai]+1)
            else: 
                allinactive[ai]=True
                qmave[ai]=qm_t[ai]

            if qm>0.1*qm_t[ai]:
                    qq=np.mean(q0[q0>=tol_nw_forces*qm_t[ai]])
            else:
                    qq=np.mean(q0)
            Dumax=np.linalg.norm(diffW.vector()[dofindex[ai]][:], ord=np.Inf)
            rmax=np.linalg.norm(rres[dofindex[ai]], ord=np.Inf)
            dumax=np.linalg.norm(deltaW_sol.vector()[dofindex[ai]][:], ord=np.Inf)
            
            if qq>tol_nw_forces*qq_t[ai]:
                qqave[ai]=qq_t[ai]*(qqcounter[ai])/(qqcounter[ai]+1)+qq*1/(qqcounter[ai]+1)
                if rmax<1e-8*qqave[ai]:  converged1[ai]=True; converged2[ai]=True;continue
                error_nwall1[ai]=rmax/qqave[ai]
                error_nwall2[ai]=dumax/Dumax
                if np.abs(error_nwall1[ai])<tol_nw: 
                    converged1[ai]=True
                if np.abs(error_nwall2[ai])<tol_nw or Dumax<1e-8*0.2: 
                    converged2[ai]=True
                continue
            else: 
                lowflux[ai]=True
                qqave[ai]=qq_t[ai]
                if rmax<1e-8*qqave[ai]:  converged1[ai]=True; converged2[ai]=True;continue
                error_nwall1[ai]=tol_nw/tol_nw_forces*rmax/qm_t[ai]
                error_nwall2[ai]=tol_nw/1e-3*dumax/Dumax
                if np.abs(error_nwall1[ai])<tol_nw: 
                    converged1[ai]=True
                    converged2[ai]=True
                    continue
                else:
                    if np.abs(error_nwall2[ai])<tol_nw or Dumax<1e-8*0.2: 
                        converged1[ai]=True
                        converged2[ai]=True
                        continue
                    else:
                        continue



        # check if converged and stop iterating if so
        print(iter_nw-1,'Coupled:',error_nwall1,error_nwall2) 
        print(converged1,converged2,np.logical_or(converged1,converged2))
        if np.any(np.isnan(error_nwall1)) or np.any(np.isnan(error_nwall2)): converged1[:]=False;converged2[:]=False;iter_nw=maxiter_nw;break
        if np.all(np.logical_and(converged1,converged2)) and iter_nw>2: break
        
        # if maximum iterations reached stop
        if iter_nw==maxiter_nw: break

        # solve the linear system and apply corrections for the next iteration
        try: fn.solve(A, deltaW_sol.vector(), b,"mumps")
        except: converged1[:]=False;converged2[:]=False;iter_nw=maxiter_nw;break
        W_sol.vector()[:] = W_sol.vector() + deltaW_sol.vector()
        diffW.vector()[:]=W_sol.vector()-W_sol_t.vector()

        
    # resotre to time t if not converged 
    if iter_nw==maxiter_nw and not np.all(np.logical_and(converged1,converged2)):
        t-=dtime;inc-=1
        dtime=dtime/2.
        qqave[:]=qq_t[:]
        qmave[:]=qm_t[:]
        print('No convergency, decrease t newdt',t,dtime)
        W_sol.vector()[:]=W_sol_t.vector()
        if dtime<dtmin: dtflag=True;break
        continue


 
    # postprocessing if converged
    print('Converged with iters=',iter_nw,t,dtime)
    #get some macroscopic quantities
    F12_t=fn.assemble((FF_sol-fn.Identity(3))[1,2]*dx)
    F22_t=fn.assemble((FF_sol-fn.Identity(3))[2,2]*dx)
    H22_t=fn.assemble(fn.dot(HH_sol,fn.inv(FF_sol))[2]*dx)
    P22_t=fn.assemble(P0(FF_sol,pr_sol,0,FVV_sol)[2,2]*dx(0)+\
                      P0(FF_sol,pr_sol,1,FVV_sol)[2,2]*dx(1))
    P12_t=fn.assemble(P0(FF_sol,pr_sol,0,FVV_sol)[1,2]*dx(0)+\
                      P0(FF_sol,pr_sol,1,FVV_sol)[1,2]*dx(1))
    B22_t=fn.assemble((fn.dot(FF_sol,B(HH_sol,0,FF_sol))/fn.det(FF_sol))[2]*dx(0)+\
                       (fn.dot(FF_sol,B(HH_sol,1,FF_sol))/fn.det(FF_sol))[2]*dx(1))
    H22=np.append(H22,H22_t)
    F12=np.append(F12,F12_t)
    F22=np.append(F22,F22_t)
    ttt=np.append(ttt,t)
    P12=np.append(P12,P12_t)
    P22=np.append(P22,P22_t)
    B22=np.append(B22,B22_t)
    # update output files
    post(file_results,t,mesh,domains,u_sol,pr_sol,lambda_tensor,\
      F_macro,P0,dx,v_sol,lambdav_vector,H_macro2,B,FVV_sol,FF_sol,HH_sol)
    np.savetxt(filename+'/unihshear_data.txt',np.array([ttt,H22,F22,P22,B22,F12,P12]).transpose())
    adjust_xdmf(filename+'/unihshear_out.xdmf')

    #update solution and review residuals
    W_sol_t.vector()[:]=W_sol.vector()
    diffW2.vector()[:]=diffW.vector()
    diffW.vector()[:]=0
    qq_t[:]=qqave[:]
    qm_t[:]=qmave[:]
    qqcounter[lowflux==False]+=1
    qmcounter[allinactive==False]+=1

    # adjust dt if converged easily
    if(iter_nw<=5):
        dtime=dtime*1.5
        if dtime>dtmax: dtime=1*dtmax
        print('Increase dt',dtime)
        
  # stop simulation if minimum time increment requires
  if dtflag: break
