import itertools 
import numpy as np
import qutip
import qutip.states
import tomographer
import tomographer.tools.densedm
import tomographer.tools.densedm.mle
import tomographer.querrorbars
import tomographer.jpyutil
from IPython.display import display, Markdown
import polytope.polytope as pc
import matplotlib.pyplot as plt
from qutip import *
from scipy.special import erfinv
from sampling import *
from polyconfiregion import *



#measurements

#single qubit Pauli mmt (MUB-POVM)

paulibasis=[sigmax().full(), sigmay().full(), sigmaz().full()]
GellMann3=[np.matrix('1. 0 0; 0 1. 0;0 0 1.')*np.sqrt(2/3.),np.matrix('0 1. 0; 1. 0 0;0 0 0'),
np.matrix('0 -1j 0; 1j 0. 0;0 0 0'),np.matrix('1. 0. 0; 0 -1 0;0 0 0'),np.matrix('0 0 1.; 0 0 0;1. 0 0'),
np.matrix('0 0 -1j; 0. 0 0;1j 0 0'),np.matrix('0 0 0; 0 0 1.;0 1. 0'),np.matrix('0 0. 0; 0 0 -1j;0 1j 0'),
np.matrix('1 0. 0; 0 1 0;0 0 -2.')/np.sqrt(3)]


mub = [ [
        np.array([[.5, .5],[.5, .5]])/3,     # X, +1 outcome
        np.array([[.5, -.5],[-.5, .5]])/3,   # X, -1 outcome
   ],[
        np.array([[.5, -.5j],[.5j, .5]])/3,  # Y, +1 outcome
        np.array([[.5, .5j],[-.5j, .5]])/3,  # Y, -1 outcome
   ],[
        np.array([[1.,0],[0,0]])/3,          # Z, +1 outcome
        np.array([[0,0],[0,1.]])/3,  ]       # Z, -1 outcome
    ]

pauli=[qeye(2).full(),sigmax().full(), sigmay().full(), sigmaz().full()]


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


#parameters
n=2
eps=0.001
d=2**n
N=10000
basis=basis(n,pauli)
mmts=allmmts(sic2,n,False)
m=1/np.trace(mmts[0]).astype(np.double)

print len(basis),len(mmts),m #sanity check


#data for bell-pair
data=np.array([0.0046, 0.0807, 0.0841, 0.0829, 0.0806, 0.0571, 0.014, 0.0972, \
0.0794, 0.1037, 0.0576, 0.0119, 0.0801, 0.0139, 0.0972, 0.055])
'''
upper=np.array([0.00380018, 0.0129861, 0.0132163, 0.0131358, 0.0129793, 0.0111843, \
0.00604063, 0.0140484, 0.0128964, 0.0144324, 0.0112269, 0.0056331, \
0.0129448, 0.00602198, 0.0140484, 0.0110029])
lower=np.array([-0.00243735, -0.0118272, -0.0120667, -0.0119829, -0.0118201, \
-0.00996026, -0.00469918, -0.012935, -0.0117339, -0.013337, \
-0.0100042, -0.00428624, -0.0117843, -0.00468027, -0.012935, \
-0.00977308])'''
ref=qutip.states.ket2dm(qutip.Qobj(np.array([0,1,1j,0.]/np.sqrt(2))))
ref.dims=[[2, 2], [2, 2]]

print mle(data,N,mmts,d),ref

errorbars=errorbar(eps,N,data,m,d)

cr=polytopeCR(errorbars,data,basis,mmts,d)

print 'polytope done'
print cr.chebXc,cr.chebR,cr.volume
print boundingboxsize(cr)




sample=sampling(d,mmts,data,N)
print'sampling done'
print len(sample)



print figureofinterest_bell(cr,sample,basis,ref,d)


