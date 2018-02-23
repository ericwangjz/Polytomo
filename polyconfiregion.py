
import itertools 
import numpy as np
import qutip.states
from IPython.display import display, Markdown
import polytope.polytope as pc
import matplotlib.pyplot as plt
from qutip import *
from scipy.special import erfinv


def state2r(rho,basis,d):
    r=[]
    for j in basis:
        r.append((d/2.)/np.sqrt(d*(d-1)/2)*np.trace(np.matrix(rho.full())*np.matrix(j)))
    return np.array(r).astype(np.double)
def r2state(r,basis,d):
    state=0.
    for i in range(len(r)):
        state+=r[i]*basis[i]
    state=(np.sqrt(d*(d-1.)/2.)*state+np.eye(d))/d
    return state

'''

def planeintersectpara(pl,d1,d2):
    p1=pc.Polytope(pl,np.array(d1))
    p2=pc.Polytope(-pl,-np.array(d2))
    intersect=pc.intersect(p1,p2)
    interval=intersect.bounding_box
    return intersect.chebXc,intersect.chebR,intersect,interval
 '''

def mmtrep(basis,mmt,d):
    mmtrep=[]
    for i in mmt:
        r=[]
        for j in basis:
            r.append(((d/(2*np.trace(i)))/np.sqrt(d*(d-1)/2.))*np.trace(i*j))
        mmtrep.append(np.array(r))
    return mmtrep
#define general su(d) basis:
def basis(n,pauli):
    kk=list(itertools.product(pauli,repeat=n))
    productbasis=[]
    for i in kk:
        a=np.array([1])
        for j in i:
            a=np.kron(a,j)
        productbasis.append(np.matrix(a)/np.sqrt(2**n/2))
    del productbasis[0]
    return productbasis
    
def allmmts(mmt,n,nested=True):
    mmtcomb=list(itertools.product(mmt,repeat=n))
    mmts=[]
    indices=list(itertools.product(range(2),repeat=n))
    if nested==True:
        for k in mmtcomb: 
            for index in indices:
                indx=list(reversed(index))
                a=np.array([1])
                for i in range(n):
                    a=np.kron(a,k[i][indx[i]])
                mmts.append(a)
    else:
        for i in mmtcomb:
            a=np.array([1])
            for j in i:
               a=np.kron(a,j)
            mmts.append(a)
    return mmts

#errorbar
def probit(eps):
    return -np.sqrt(2)*erfinv(eps-1.)

def errorbar(eps,n,data,m,d):
    upperrorbars=[]
    lowerrorbars=[]
    uper=[]
    lwer=[]
    z=probit(eps)
    for i in data:
            p=(i*n+z**2/2.)/(n+z**2)
            uper.append(p-i+z*np.sqrt(p*(1.-p)/(n+z**2)))
            lwer.append(p-i-z*np.sqrt(p*(1.-p)/(n+z**2)))
    for i in range(len(data)):
            upperrorbars.append((m*d*(data[i]+uper[i])-1.)/(d-1.))
            lowerrorbars.append((m*d*(data[i]+lwer[i])-1.)/(-d+1.))
    return np.append(upperrorbars,lowerrorbars)

def errorbar_cp(eps,n,upper,lower,data,m,d):
    upperrorbars=[]
    lowerrorbars=[]
    for i in range(len(data)):
            upperrorbars.append((m*d*(data[i]+upper[i])-1.)/(d-1.))
            lowerrorbars.append((m*d*(data[i]+lower[i])-1.)/(-d+1.))
    return np.append(upperrorbars,lowerrorbars)


def polytopeCR(errorbars,data,basis,mmt,d):
    pl=np.array(mmtrep(basis,mmt,d)).astype(np.double)
    planenorm=np.append(pl,-pl,axis=0)
    return pc.Polytope(planenorm,errorbars)


def boundingboxsize(polyto):
	return polyto.bounding_box[1]-polyto.bounding_box[0]


