import logging
logging.basicConfig()
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
    #determinant=min(np.linalg.eig(state)[0])
    #if determinant<0.:
    #    print 'unphysical',determinant
    #else:
    #   print 'physical',determinant
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

def errorbar(eps,N,data,mmts,d):
    upperrorbars=[]
    lowerrorbars=[]
    uper=[]
    lwer=[]
    z=probit(eps)
    assert len(data)==len(mmts) #sanity check
    l=len(data)
    for i in data:
            p=(i*N+z**2/2.)/(N+z**2)
            uper.append(p-i+z*np.sqrt(p*(1.-p)/(N+z**2)))
            lwer.append(p-i-z*np.sqrt(p*(1.-p)/(N+z**2)))
    for i in range(l):  
            upperrorbars.append((d/(np.trace(mmts[i]).astype(np.double))*(data[i]+uper[i])-1.)/(d-1.))
            lowerrorbars.append((d/(np.trace(mmts[i]).astype(np.double))*(data[i]+lwer[i])-1.)/(-d+1.))
    return  np.append(upperrorbars,lowerrorbars)#np.array(upperrorbars)



def errorbar_cp(eps,n,upper,lower,data,mmts,d):
    upperrorbars=[]
    lowerrorbars=[]
    for i in range(len(data)):
            upperrorbars.append((d/(np.trace(mmts[i]).astype(np.double))*(data[i]+upper[i])-1.)/(d-1.))
            lowerrorbars.append((d/(np.trace(mmts[i]).astype(np.double))*(data[i]+lower[i])-1.)/(-d+1.))
    return np.append(upperrorbars,lowerrorbars)


def singlepolytopeCR(errorbars,basis,mmt,d,eps):
    pl=np.array(mmtrep(basis,mmt,d)).astype(np.double)
    planenorm=np.append(pl,-pl,axis=0)
    #planenorm=pl
    return pc.Polytope(planenorm,errorbars),eps

def polytomerge(polytopes):
    eps=0.
    mergedA=polytopes[0][0].A
    mergedb=polytopes[0][0].b
    for i in polytopes[1:]:
        eps+=i[1]
        mergedA = np.vstack([mergedA, i[0].A])
        mergedb = np.hstack([mergedb, i[0].b])
    return pc.Polytope(mergedA, mergedb),eps

def polytopeCR(data,basis,mmts,N,d,eps,onePOVM=True):
    if onePOVM == True:
        errorbars=errorbar(eps,N,data,mmts,d)
        polyto,error=singlepolytopeCR(errorbars,basis,mmts,d,eps)
    else:
        polytopes=[]
        assert len(data)==len(mmts) #sanity check
        l=len(data)
        for i in range(l):
            errorbars=errorbar(eps/l,N[i],data[i],mmts[i],d)
            polytopes.append(singlepolytopeCR(errorbars,basis,mmts[i],d,eps/l))
        polyto,error=polytomerge(polytopes)
    return polyto


def boundingboxsize(polyto,n,d,pauli):
    boundingbox={}
    basiselements=list(itertools.product(range(len(pauli)),repeat=n))
    width=polyto.bounding_box[1]-polyto.bounding_box[0]
    del basiselements[0]
    for i in range(d**2-1):
        boundingbox['basis'+str(basiselements[i])]=width[i]
    return boundingbox



