import itertools
import numpy as np
import qutip
import qutip.states
import matplotlib.pyplot as plt
from qutip import *
from polyconfiregion import *
#basis elements
paulibasis=[sigmax().full(), sigmay().full(), sigmaz().full()] #w/o id
pauli=[qeye(2).full(),sigmax().full(), sigmay().full(), sigmaz().full()] # w/ id
GellMann3=[np.matrix('1. 0 0; 0 1. 0;0 0 1.')*np.sqrt(2/3.),np.matrix('0 1. 0; 1. 0 0;0 0 0'),
           np.matrix('0 -1j 0; 1j 0. 0;0 0 0'),np.matrix('1. 0. 0; 0 -1 0;0 0 0'),np.matrix('0 0 1.; 0 0 0;1. 0 0'),
           np.matrix('0 0 -1j; 0. 0 0;1j 0 0'),np.matrix('0 0 0; 0 0 1.;0 1. 0'),np.matrix('0 0. 0; 0 0 -1j;0 1j 0'),
           np.matrix('1 0. 0; 0 1 0;0 0 -2.')/np.sqrt(3)]

#measurements
#single qubit Pauli mmt (MUB-POVM)
mub = [
         np.array([[.5, .5],[.5, .5]])/3,     # X, +1 outcome
         np.array([[.5, -.5],[-.5, .5]])/3,   # X, -1 outcome
         
            np.array([[.5, -.5j],[.5j, .5]])/3,  # Y, +1 outcome
            np.array([[.5, .5j],[-.5j, .5]])/3,  # Y, -1 outcome
    
               np.array([[1.,0],[0,0]])/3,          # Z, +1 outcome
               np.array([[0,0],[0,1.]])/3,         # Z, -1 outcome
       ]



#sic-povm d=2
normal1 = np.array([0, 0, 1.])
normal2 = np.array([2*np.sqrt(2)/3., 0, -1./3])
normal3 = np.array([-np.sqrt(2)/3., np.sqrt(2./3), -1./3])
normal4 = np.array([-np.sqrt(2)/3., -np.sqrt(2./3), -1./3])
sic2=[]
for i in [normal1,normal2,normal3,normal4]:
    sic2.append(r2state(i,paulibasis,2)/2.)
#sic-povm d=3
v1=np.array([1.,0.,0.])
v2=np.array([1./2,np.sqrt(3)*1j/2.,0.])
v3=np.array([1./2,-np.sqrt(3)*1j/2.,0.])
v4=np.array([1./2,1./2,1./np.sqrt(2)])
v5=np.array([1./2,1./2,np.exp(2j*np.pi/3)/np.sqrt(2)])
v6=np.array([1./2,1./2,np.exp(-2j*np.pi/3)/np.sqrt(2)])
v7=np.array([1./2,-1./2,1./np.sqrt(2)])
v8=np.array([1./2,-1./2,np.exp(2j*np.pi/3)/np.sqrt(2)])
v9=np.array([1./2,-1./2,np.exp(-2j*np.pi/3)/np.sqrt(2)])
sic3=[]
for i in [v1,v2,v3,v4,v5,v6,v7,v8,v9]:
    sic3.append(qutip.states.ket2dm(qutip.Qobj(i/np.sqrt(3))).full())









