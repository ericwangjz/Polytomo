
import itertools 
import numpy as np
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
from polyconfiregion import *

def mle(data,N,mmts,d,iterations=10000):
    k=iterations
    llh = tomographer.densedm.IndepMeasLLH(tomographer.densedm.DMTypes(d))
    llh.setMeas(mmts, (data*N).astype(int))
    (rho_MLE, dd) = tomographer.tools.densedm.mle.find_mle(llh, solver_opts={'max_iters':k})
    return rho_MLE

def store_sample_figofmerit(T):
    samples.append(T)
    return 0
samples=[]
def sampling(d,mmts,data,N,samplesize=2000):
    #r = None # global variable

    with tomographer.jpyutil.RandWalkProgressBar() as prg:
        r = tomographer.tomorun.tomorun(
                                    # the dimension of the quantum system
                                    dim=d,
                                    # the tomography data
                                    Emn=mmts,
                                    Nm=data*N,
                                    # Random Walk parameters: step size, sweep size, number of thermalization sweeps, number of live sweeps
                                    mhrw_params=tomographer.MHRWParams(0.01,100,2000,samplesize),
                                     fig_of_merit=store_sample_figofmerit,rng_base_seed=None,
                                    progress_fn=prg.progress_fn
                                    )
        prg.displayFinalInfo(r['final_report_runs'])
    return samples

def figureofinterest_bell(polyto,sample,basis,ref,d):
    fide=[]
    trdist=[]
    negativity=[]
    for i in sample:
        rho=Qobj(np.matrix(i)*np.matrix(i).H,dims=[[2, 2], [2, 2]])
        if pc.is_inside(polyto,state2r(rho,basis,d))==True:
            fide.append(fidelity(ref,rho))
            trdist.append(tracedist(ref,rho))
            negativity.append(partial_transpose(rho,[1,0]).norm()/2.-0.5)
    result={'fidelity':[min(fide),max(fide)],'negativity':[min(negativity),max(negativity)],'tracedistance':[min(trdist),max(trdist)]}
    return result
def figureofinterest2(polyto,sample,basis,ref,d):
    fide=[]
    trdist=[]
    rela_ent=[]
    for i in sample:
        rho=Qobj(np.matrix(i)*np.matrix(i).H)
        if pc.is_inside(polyto,state2r(rho,basis,d))==True:
            fide.append(fidelity(ref,rho))
            trdist.append(tracedist(ref,rho))
            rela_ent.append(entropy._entropy_relative(ref,Qobj(np.matrix(i)*np.matrix(i).H)))
    result={'fidelity':[min(fide),max(fide)],'tracedistance':[min(trdist),max(trdist)],'rela_entropy':[min(rela_ent),max(rela_ent)]}
    return result