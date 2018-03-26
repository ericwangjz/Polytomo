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
import matplotlib.pyplot as plt
from qutip import *
from scipy.special import erfinv
from sampling import *
from polyconfiregion import *
from parameters import *

class _Ns: pass # dummy object acts as namespace
def simulate_measurements(rho, Mk, num_samples_per_setting):

    
    num_settings = len(Mk)
    
    num_total_effects = np.sum([ len(Mk[k]) for k in range(len(Mk)) ])
    
    dim = rho.dims[0][0]
    
    # prepare the list of POVM effects, and simulate the measurements
    Emn = []
    Nm = []
    
    idx = 0
    for k in range(num_settings):
        
        # get some random numbers for this measurement setting
        x = np.random.rand(num_samples_per_setting)
        proboffset = 0

        
        for i in range(len(Mk[k])):
            
            Ek = qutip.Qobj(Mk[k][i])
            p = qutip.expect(Ek, rho) # = trace(Ek * rho)
            Emn.append( Ek.data.toarray().astype(dtype=complex) )
            Nm.append( np.count_nonzero( (proboffset <= x) & (x < proboffset+p) ) )
            
            proboffset += p
        
        # sanity check
        assert np.abs(proboffset-1.0) < 1e-6
    
    d = _Ns()
    d.Emn = Emn
    d.Nm = np.array(Nm)
    return d


#parameters
eps=0.01 #confidence level is 1.-eps
n=2       #number of qubits
#d=2**n    #dimension
N=20000  #number of shots
#basis=basis(n,pauli)  #basis of real space R^n where the polytope confidence region is represented.
#mmts=allmmts(sic2,n,False)  #'allmmts' generates all the POVM elements from single qubit measurements,
# nested=False indicates no rearrangement of the POVM elements is needed here.

#target state
#rho_sim=qutip.states.ket2dm(qutip.Qobj(np.array([0,1,1j,0.]/np.sqrt(2))))
#dims indicates the quantum state is bipartite so to calculate negativity.
Nscaling={}
error=[]
fulldimerror=[]

'''
for i in range(1,3):
    d=3**i
    bases=basis(i,GellMann3)
    mmts=allmmts(sic3,i,False)
    values=[]
    for j in range(10):
        dd = simulate_measurements(rand_dm(d), [mmts], N)
        data=dd.Nm/float(N)
        cr=polytopeCR(data,bases,mmts,N,d,eps,True,'CP')
        if pc.is_fulldim(cr)==False:
            fulldimerror.append((N,j))
        sample=sampling(d,mmts,data,N,20000)
        count=0.
        print d,j
        for k in sample:
            state=state2r(Qobj(np.matrix(k)*np.matrix(k).H),bases,d)
            if pc.is_inside(cr,state)==False:
                count+=1.
        print count
        if count !=0.:
            values.append(eps/(count/len(sample)))
        else:
            error.append(('error at',N,j))
        print values
    Nscaling[str(d)]=values
    print Nscaling
print Nscaling



dimes=[2*2*3,2*3*3]
mmtss=[[],[]]
basises=[[],[]]
for i in sic2:
    for j in sic3:
        mmtss[0].append(np.kron(i,j))
for i in pauli:
    for j in GellMann3:
        basises[0].append(np.matrix(np.kron(i,j))/np.sqrt(2))
del basises[0][0]

for i in allmmts(sic2,2,False):
    for j in sic3:
        mmtss[0].append(np.kron(i,j))
for i in pauli:
    for k in pauli:
        kk=np.kron(i,k)/np.sqrt(2)
        for j in GellMann3:
            basises[0].append(np.matrix(np.kron(kk,j))/np.sqrt(2))
del basises[0][0]

for i in sic2:
    for j in allmmts(sic3,2,False):
        mmtss[1].append(np.kron(i,j))
for i in pauli:
    for j in GellMann3:
        kk=np.kron(i,j)/np.sqrt(2)
        for k in GellMann3:
            basises[1].append(np.matrix(np.kron(kk,k))/np.sqrt(2))
del basises[1][0]

for i in range(2):
    d=dimes[i]
    bases=basises[i]
    mmts=mmtss[i]
    values=[]
    for j in range(10):
        dd = simulate_measurements(rand_dm(d), [mmts], N)
        data=dd.Nm/(np.array(N).astype(np.double))
        cr=polytopeCR(data,bases,mmts,N,d,eps,True,'CP')
        if pc.is_fulldim(cr)==False:
            fulldimerror.append((N,j))
        sample=sampling(d,mmts,data,N,20000)
        count=0.
        print d,j
        for k in sample:
            state=state2r(Qobj(np.matrix(k)*np.matrix(k).H),bases,d)
            if pc.is_inside(cr,state)==False:
                count+=1.
        print count
        if count !=0.:
            values.append(eps/(count/len(sample)))
        else:
            error.append(('error at',N,j))
        print values
    Nscaling[str(d)]=values
    print Nscaling
print Nscaling


for i in range(1,5):
    d=2**i
    bases=basis(i,pauli)  #basis of real space R^n where the polytope confidence region is represented.
    mmts=allmmts(sic2,i,False)
    values=[]
    for j in range(5):
        dd = simulate_measurements(rand_dm(d), [mmts], N)
        data=dd.Nm/(np.array(N).astype(np.double))
        cr=polytopeCR(data,bases,mmts,N,d,eps,True,'AC')
        if pc.is_fulldim(cr)==False:
            fulldimerror.append((N,j))
        sample=sampling(d,mmts,data,N,20000)
        count=0.
        print d,j
        for k in sample:
            state=state2r(Qobj(np.matrix(k)*np.matrix(k).H),bases,d)
            if pc.is_inside(cr,state)==False:
                count+=1.
        print count
        if count !=0.:
            values.append(eps/(count/len(sample)))
        else:
            error.append(('error at',N,j))
        print values,Nscaling
    Nscaling[str(d)]=values
    print Nscaling
print Nscaling



print error
print fulldimerror

 

for i in range(2,8):
    N=1000*2**i
    values=[]
    for j in range(10):
        dd = simulate_measurements(rand_dm(d), [mmts], N)
        data=dd.Nm/(np.array(N).astype(np.double))
        cr=polytopeCR(data,basis,mmts,N,d,eps,True,'CP')
        if pc.is_fulldim(cr)==False:
            fulldimerror.append((N,j))
        sample=sampling(d,mmts,data,N,10000)
        count=0.
        print i,j
        for k in sample:
            state=state2r(Qobj(np.matrix(k)*np.matrix(k).H),basis,d)
            if pc.is_inside(cr,state)==False:
                count+=1.
        print count
        if count !=0.:
            values.append(eps/(count/len(sample)))
        else:
            error.append(('error at',N,j))
        print values
    Nscaling[str(N)]=values
    print Nscaling'''

'''
data=[]
for i in mmts:
    dd = simulate_measurements(rho_sim, [[i,np.eye(d)-i]], N)
    data.append(dd.Nm[1]/N)
'''


