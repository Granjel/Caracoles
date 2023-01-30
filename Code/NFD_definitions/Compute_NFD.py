"""
author: J. W. Spaak
load empirical data and compute NFD
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from NFD_functions import compute_NFD, InputError
from timeit import default_timer as timer

###############################################################################
# load empirical data
optim = "bobyqa" # the optimizer chosen
data_alpha = pd.read_csv("data/alpha_RK_{}.csv".format(optim),
                         sep = ";", decimal = ",")
data_alpha_se = pd.read_csv("data/alpha_standard_error_RK_{}.csv".format(optim)
                , sep = ";", decimal = ",")

data_lamb = pd.read_csv("data/lambda_RK_{}.csv".format(optim),
                        sep = ";", decimal = ",")
data_lamb_se = pd.read_csv("data/lambda_standard_error_RK_{}.csv".format(optim),
                        sep = ";", decimal = ",")

species_rates = pd.read_csv("data/plant_species_traits.csv",
                        sep = ";", decimal = ",")

max_N_mono = pd.read_csv("upper_bound_eq_abundances.csv", sep = ";",
                         index_col = 0)

years = sorted(set(data_alpha.year))
plots = sorted(set(data_alpha["plot"]))
species = data_alpha.keys()[3:] # remove year, plot and focal
species_rates = species_rates[[sp in species for
                               sp in species_rates["species.code"]]]
n_spec = len(species)
cases = ["{}_{}".format(year, plot) for year in years for plot in plots]

# fill in germination and survival rates
g_fix, s_fix, N_max = np.empty((3,len(species)))

# ensure that all files have the same order
for spec in species:
    ind = species_rates["species.code"] == spec
    g_fix[species == spec] = species_rates["germination.rate"][ind]
    s_fix[species == spec] = species_rates["seed.survival"][ind]
    N_max[species == spec] = max_N_mono.upper_bound_N_star[spec]

# dictionaries that contain alpha and lambda values per plot and year
alpha = {case: pd.DataFrame(index = species, columns = species, dtype = float)
            for case in cases}
alpha_se = {case: pd.DataFrame(index = species, columns = species,
                               dtype = float)
            for case in cases}
lamb = {case: pd.DataFrame(index = species, columns = ["lambda"],dtype = float)
            for case in cases}
lamb_se = {case: pd.DataFrame(index = species, columns = ["lambda"],
                              dtype = float)
            for case in cases}

# fill alpha values
for i, r in data_alpha.iterrows():
    alpha["{}_{}".format(r.year, r["plot"])].loc[r.focal] = r[species]
for i, r in data_alpha_se.iterrows():
    alpha_se["{}_{}".format(r.year, r["plot"])].loc[r.focal] = r[species]

# fill lambda values   
for i, row in data_lamb.iterrows():
    lamb["{}_{}".format(row.year, row["plot"])].loc[row.sp] = row["lambda"]
for i, row in data_lamb_se.iterrows():
    lamb_se["{}_{}".format(row.year, row["plot"])].loc[row.sp] = row["lambda"]


ND_cols = ["ND_{}".format(sp) for sp in species]
FD_cols = ["FD_{}".format(sp) for sp in species]
equi_cols = ["equi_{}".format(sp) for sp in species]
state_cols = ["rand_state_{}".format(i) for i in range(4)]
cols = ["case"] + state_cols + ND_cols + FD_cols + equi_cols
ND_cols = np.array(ND_cols)
FD_cols = np.array(FD_cols)
equi_cols = np.array(equi_cols)
       
iters = 100 # number of iterations per plot
NFD_dataframe = pd.DataFrame(columns = cols, index = range(len(cases)*iters))

# compute niche and fitness differences
counter = -1
for case in cases:
    for i in range(iters):
        counter += 1
        start = timer()
        NFD_dataframe.loc[counter, state_cols] = np.random.get_state()[1:]
        NFD_dataframe.loc[counter, "case"] = case
        A = alpha[case].values
        lam = lamb[case]["lambda"].values
        
        # randomize values
        A = A + np.random.uniform(-1,1, A.shape)*alpha_se[case].values
        lam = lam + np.random.uniform(-1,1, lam.shape)*lamb_se[case]["lambda"]
        lam = lam.values
        
        # remove species with non-finite equilibria densities
        N_star_mono = (np.exp(-1/np.diag(A)*
                              np.log((1-(1-g_fix)*s_fix)/g_fix/lam))-1)/g_fix
        
    
        finite = (np.isfinite(N_star_mono) &
                  (N_star_mono>0) &
                  (N_star_mono<N_max))
        N_star_mono = N_star_mono[finite]
        case_species = lamb[case].index[finite]
        
        ind = np.arange(len(A))[finite]
        if len(ind)<2:
            continue
        A = A[ind,ind[:,np.newaxis]]
        lam = lam[ind]
        s = s_fix[ind]
        g = g_fix[ind]
        try:
            pars, equi = compute_NFD(A,lam, s, g)
        except InputError:
            continue
        
        NFD_dataframe.loc[counter, ND_cols[finite]] = pars["ND"]
        NFD_dataframe.loc[counter, FD_cols[finite]] = 1-1/(1-pars["FD"])
        NFD_dataframe.loc[counter, equi_cols[finite]] = equi
        

        print(case, sum(finite), timer()-start)

NFD_dataframe.to_csv("NFD_uncertanty.csv", index = False)