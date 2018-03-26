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



#Here we give some examples to demostrate how to use our method to yield confidence region and error bars in QST.




'''
    
    # bell state QST example with data simulatde from single qubit sic2 measurements
    
#parameters
eps=0.001 #confidence level is 1.-eps
n=2       #number of qubits
d=2**n    #dimension
N=10000   #number of shots
basis=basis(n,pauli)  #basis of real space R^n where the polytope confidence region is represented.
mmts=allmmts(sic2,n,False)  #'allmmts' generates all the POVM elements from single qubit measurements,
                            # nested=False indicates no rearrangement of the POVM elements is needed here.

#target state
ref=qutip.states.ket2dm(qutip.Qobj(np.array([0,1,1j,0.]/np.sqrt(2))))
ref.dims=[[2, 2], [2, 2]] #dims indicates the quantum state is bipartite so to calculate negativity.

#data for bell-pair
data=np.array([0.0046, 0.0807, 0.0841, 0.0829, 0.0806, 0.0571, 0.014, 0.0972, \
0.0794, 0.1037, 0.0576, 0.0119, 0.0801, 0.0139, 0.0972, 0.055])
'''


    #GHZ4 QST with data contributed by J.H., raw data from average readout is reduced to single shot readout format for this example.

#parameters
eps=0.01
n=4
d=2**n
N=np.full((1, 256), 1000.)[0] #1000 shots for each of 256 measurements
basis=basis(n,pauli)

#projector gggg with ZZZZ measurement
Mgggg=np.zeros((16,16))
Mgggg[15,15]=1.

#generate all measurement projectors

rot=[ry(0),rx(np.pi/2),ry(np.pi/2),rx(-np.pi)]#4 single qubit measurement directions

def allrotations(rot,n):
    kk=list(itertools.product(rot,repeat=n))
    rotations=[]
    for i in kk:
        a=np.array([1])
        for j in i:
            a=np.kron(a,j.full())
        rotations.append(a)
    return rotations

mmts=[] # on 4 qubits we have 256 combinations of measurment directions, and each is binary-outcome projective measurement.
for i in allrotations(rot,n):
    E=np.array(np.matrix(np.dot(i,Mgggg))*np.matrix(i).H)
    mmts.append([E,np.array(np.eye(d))-E])
mmts=np.array(mmts)

#target state GHZ4
GHZ4=np.zeros((16,16))
GHZ4[0][0],GHZ4[0][15],GHZ4[15][0],GHZ4[15][15]=0.5,0.5,0.5,0.5
ref=Qobj(GHZ4)
#raw counts
rawdata = [0.292532, 0.234211, 0.114336, 0.0438945, 0.121145, 0.138389, \
               0.0425058, 0.0611099, 0.133628, 0.0851014, 0.0716454, 0.0376756, \
               0.0206641, 0.0432682, 0.0568809, 0.0483958, 0.181202, 0.142092, \
               0.075136, 0.0232614, 0.0851456, 0.0984391, 0.038348, 0.0249861, \
               0.0979488, 0.0744694, 0.0327593, 0.0281836, 0, 0.0173787, 0.024647, \
               0.0218541, 0.158443, 0.126502, 0.0506197, 0.0177575, 0.0820284, \
               0.10667, 0.0401649, 0.0514757, 0.0548502, 0.057118, 0.0271862, \
               0.0392239, 0.00804861, 0.0338878, 0.0498848, 0.0337192, 0, 0.0177866, \
               0.00137944, 0.00105243, 0.0110952, 0.0121762, 0.00817515, 0.00877429, \
               0.015751, 0.0266023, 0.0248605, 0.013478, 0.00427652, 0.0112655, \
               0.0311058, 0.0195332, 0.174914, 0.14432, 0.0668748, 0.0242261, \
               0.0736019, 0.0718272, 0.0174085, 0.0363661, 0.0803654, 0.0538751, \
               0.0425813, 0.0169116, 0.0103555, 0.0291196, 0.0311825, 0.0203856, \
               0.0802728, 0.0699345, 0.03476, 0.017505, 0.066943, 0.117461, \
               0.100073, 0.074547, 0.015284, 0.0591221, 0.032619, 0.0926716, \
               0.0381729, 0.0746415, 0.145636, 0.133959, 0.0578487, 0.0364706, \
               0.0204115, 0.0137642, 0.0104791, 0.0451087, 0.0418586, 0.0815989, \
               0.0472619, 0.0183495, 0.0903417, 0.0811906, 0.0449509, 0.100372, \
               0.171499, 0.145606, 0.0103187, 0.0210738, 0., 0.00583616, 0.0493935, \
               0.0698696, 0.136903, 0.0990162, 0.038865, 0.115471, 0.13227, \
               0.117544, 0.0656657, 0.123136, 0.206127, 0.190925, 0.191936, \
               0.153785, 0.0807358, 0.0370085, 0.078236, 0.102192, 0.0266069, \
               0.0378574, 0.0796009, 0.0691021, 0.0334089, 0.0267876, 0.0109755, \
               0.0281196, 0.0419588, 0.0300273, 0.0730362, 0.0596704, 0.0316473, \
               0.0111006, 0.0205293, 0.0505562, 0.0226231, 0.0572812, 0.0555574, \
               0.0267446, 0.0604946, 0.0471528, 0.0345579, 0.0600498, 0.130209, \
               0.0954661, 0.111391, 0.0950624, 0.0417133, 0.0222937, 0.0682222, \
               0.0314955, 0.0359251, 0.0301169, 0.0889993, 0.0628046, 0.0834015, \
               0.0225935, 0.0163621, 0.0482229, 0.0665892, 0.0493484, 0, 0.00973482, \
               0.00153566, 0., 0.0296796, 0.0264542, 0.0825183, 0.0384618, \
               0.0274792, 0.0487861, 0.0751672, 0.0589827, 0.0500777, 0.0516663, \
               0.132913, 0.117011, 0.0125592, 0.00592565, 0.00244603, 0., 0, 0, 0, \
               0, 0.00248995, 0.000690256, 0.00171179, 0., 0., 0.0121102, \
               0.00646147, 0.00681893, 0.0188521, 0.00987664, 0.00931304, 0, \
               0.0479011, 0.0678369, 0.108076, 0.0780356, 0.0367367, 0.0723673, \
               0.106802, 0.0811194, 0.0579062, 0.10295, 0.189735, 0.166155, \
               0.0167237, 0.0188275, 0.0132671, 0.00703162, 0.046196, 0.0568548, \
               0.0943446, 0.0646387, 0.0566357, 0.0954569, 0.127821, 0.0937569, \
               0.0633961, 0.11502, 0.195524, 0.17758, 0.0388139, 0.0295614, \
               0.00749195, 0.00787807, 0.100157, 0.123629, 0.225818, 0.145807, \
               0.0914927, 0.169168, 0.238117, 0.175406, 0.113877, 0.20382, 0.397459, \
           0.30561];
data=[]
for i in rawdata:
    data.append([i,1.-i]) #data format should fit the POVM elements, i.e. binary counts in this case.
data=np.array(data)





#execution

print 'constructing the polytope confidence region......'
#print errorbar(eps,N,data,mmts,d)
cr=polytopeCR(data,basis,mmts,N,d,eps,False)
#onePOVM=False indicates we have a combinations of POVM's in the experiments, which in the GHZ4 example is 256 projective measurments. The default is True as in the bell pair example.

print pc.is_fulldim(cr)#check whether the polytope is full-dimensional.

print 'polytope done'

#return the bounding box indicating the error distribution in different axes (as a dictionary), corresponding to the basis chosen (tensor power of Pauli basis is used in our examples).
print boundingboxsize(cr,n,d,pauli)
#return the state corresponding to the Chebyshev center of the polytop,(not necessarily physical)
print r2state(cr.chebXc,basis,d)
#return the MLE estimate, you can set the max numeber of iterations. Default is 10000
#reshape the array to fit into the MLE function
s=data.size
data=data.reshape(s,)
mmts=mmts.reshape(1,s,d,d)[0]
mlestate = mle(data,N[0],mmts,d,40000)
print mlestate


#draw confidence intervals w.r.t. certain figure of merit from the polytope confidence region.

sample=sampling(d,mmts,data,N[0],4000) # sample from the polytope
print'sampling done'

print confidenceinterval(cr,sample,basis,ref,d,'fidelity','trdist') #output figure of merit


