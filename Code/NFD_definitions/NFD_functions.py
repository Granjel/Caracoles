"""
author: J. W. Spaak
Functions to compute NFD for the annual plant model
"""

import numpy as np
from NFD_definitions.numerical_NFD import NFD_model, InputError

def apm(N, a, lamb, s, g, pr = False):
    # log transformed annual plant model with exponential growth
    F = lamb*np.exp(-a.dot(np.log(g*N+1)))
    if pr:
        print(np.round(N,2),np.round(np.log((1-g)*s + g*F),2))
    return np.log((1-g)*s + g*F)


def compute_NFD(A,lam,s,g):
    # compute  NFD for the annual plant model
    
    sub_equi = np.zeros(A.shape) # equilibrium densities
    r_i = np.zeros(len(A)) # invasion growth rates
    
    # first compute sub-equilibrium densities and invasion growth rates
    for i in range(len(A)):
        ind = np.arange(len(A))[i != np.arange(len(A))]
        pre_equi = -np.ones(len(A)-1)
        
        # remove species with negative equi densities until all are positive
        while min(pre_equi)<0:
            # compute equilibrium
            pre_equi = np.linalg.solve(A[ind[:,np.newaxis],ind],
                            -np.log((1-(1-g[ind])*s[ind])/g[ind]/lam[ind]))
            pre_equi = (np.exp(pre_equi)-1)/g[ind]
            try:
                if min(pre_equi)<0: #is one species negative?
                    ind = ind[pre_equi>min(pre_equi)]
            except ValueError: # no species can survive in monoculutre
                raise InputError
        
        sub_equi[i,ind] = pre_equi
        if not np.all(np.isfinite(sub_equi)):
            raise InputError
        # check whether growth is close to zero
        growth_res = apm(sub_equi[i], A, lam, s , g)
        if (np.abs(growth_res[ind])>1e-5).any():
            raise # equilibrium not computed correctly
            
        # save invasion growth rate
        r_i[i] = growth_res [i] 
        
    # compute NFD
    pars = {"r_i":r_i, "N_star":sub_equi, "f": lambda N: apm(N,A,lam,s,g),
            "f0": apm(0,A,lam,s,g)}
    # experimental = True, to overwride equilibrium search
    pars = NFD_model(lambda N: apm(N, A, lam, s,g,False), len(A), 
                     pars = pars, experimental=True)
        
    # some species have c = 0, beacuse they're not present
    # replace c with two species c
    for i in range(len(A)):
        for j in range(len(A)):
            if j<i:
                continue
            c = pars["c"][[i,j],[j,i]]
            ind = np.array([i,j])
            if ( (c==0) | (c==np.inf)).all():
                pars_new = {}
                try:
                    pars_2 = NFD_model(apm, args = (A[ind[:,np.newaxis],ind],
                                    lam[ind], s[ind], g[ind]), pars = pars_new)
                    pars["c"][[i,j],[j,i]] = pars_2["c"][[0,1],[1,0]]
                except InputError:
                    continue
    
    # recompute ND and FD with new computed c
    pars["fc"] = np.empty(len(A)) # compute no-niche growth rate
    for i in range(len(A)):
        N = np.zeros(len(A))
        N[i] = np.nansum((pars["N_star"]*pars["c"])[i])
        pars["fc"][i] = apm(N, A, lam, s, g)[i]
    
    pars["ND"] = (pars["r_i"]-pars["fc"])/(pars["f0"]-pars["fc"])
    pars["FD"] = pars["fc"]/pars["f0"]
    
    # compute_equi of species with positive invasion growth rates
    equi_all = np.full(g.shape, np.nan)
    ind = np.arange(len(A))[pars["r_i"]>0]
    equi = np.linalg.solve(A[ind[:,np.newaxis],ind],
                            -np.log((1-(1-g[ind])*s[ind])/g[ind]/lam[ind]))
    equi = (np.exp(equi)-1)/g[ind]
    equi_all[ind] = equi
    
    # check whether invasion growth rate correctly predicts equilibrium
    pars["surv"] = np.full(len(A), False, dtype = bool)
    pars["surv"][ind] = equi>0
    return pars, equi_all