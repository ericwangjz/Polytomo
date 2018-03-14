
import itertools 
import numpy as np
import qutip.states
from IPython.display import display, Markdown
import polytope.polytope as pc
import matplotlib.pyplot as plt
from qutip import *
from scipy.special import erfinv
import math
from scipy.optimize import brentq

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
            traceprod=np.trace(np.matrix(j)*np.matrix(i))
            r.append(((d/(2*np.trace(i)))/np.sqrt(d*(d-1)/2.))*traceprod)
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

def errorbar(eps,N,data,mmts,d,method='CP'):
    assert len(data)==len(mmts) #sanity check
    l=len(data)
    upperrorbars=[]
    lowerrorbars=[]
    uper=[]
    lwer=[]
    if method=='AC':
        z=probit(eps)
        for i in data:
                p=(i*N+z**2/2.)/(N+z**2)
                offset=z*np.sqrt(p*(1.-p)/(N+z**2))
                if offset+p<1.:
                    uper.append(offset+p-i)
                else:
                    uper.append(1.-i)
                if p-offset>0.:
                    lwer.append(p-i-offset)
                else:
                    lwer.append(-i)
    if method=='CP':
        for i in data:
            if i!=0. and i!=1.:
                f=lambda x:i*math.log(i/(i+x))+(1.-i)*math.log((1.-i)/(1.-i-x))-math.log(2.*l/eps)/N
                a=brentq(f,2e-15,1.-i-2e-15)
                uper.append(a)#(m*d*(i+a[0])-1.)/(d-1.))
                b=brentq(f,-i+2e-15,-2e-15)
                lwer.append(b)#(m*d*(i+b[0])-1.)/(1.-d))
            if i==0.:
                g=lambda x:math.log(1/(1.-x))-math.log(2.*l/eps)/N
                a=brentq(g,2e-15,1.-i-2e-15)
                uper.append(a)
                lwer.append(0.)
            elif i==1.:
                h=lambda x:math.log(1/(1.+x))-math.log(2.*l/eps)/N
                a=brentq(h,-i+2e-15,-2e-15)
                uper.append(0.)
                lwer.append(a)
    for i in range(l):
        physicalupperbound=min(data[i]+uper[i],max(np.linalg.eig(mmts[i])[0].real)) #enforce physicality
        upperrorbars.append((d/(np.trace(mmts[i]).astype(np.double))*physicalupperbound-1.)/(d-1.))
        lowerrorbars.append((d/(np.trace(mmts[i]).astype(np.double))*(data[i]+lwer[i])-1.)/(-d+1.))
    return  np.append(upperrorbars,lowerrorbars)#np.array(upperrorbars)



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

def polytopeCR(data,basis,mmts,N,d,eps,onePOVM=True,interval='CP'):
    method=interval
    if onePOVM == True:
        errorbars=errorbar(eps,N,data,mmts,d,method)
        polyto,error=singlepolytopeCR(errorbars,basis,mmts,d,eps)
    else:
        polytopes=[]
        assert len(data)==len(mmts) #sanity check
        l=len(data)
        for i in range(l):
            errorbars=errorbar(eps/l,N[i],data[i],mmts[i],d,method)
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



