# The same as the LCT_data_analysis notebook, but in a script...
import sys

import os
import estimate, sim
import matplotlib.pyplot as plt
import numpy as np
from betamix import sample_paths, BetaMixture
from sim import sim_and_fit
from scipy.stats import beta
from math import exp, log
import pickle

# Accept seed value from command-line arguments
dir_input = sys.argv[1]
dir_output = sys.argv[2]
seed = int(sys.argv[3])

# Load variables from the file
with open(os.path.join(dir_input, "amy.environment.pkl"), 'rb') as f:
    environment = pickle.load(f)

# Access variables from the loaded environment
s_hat = environment['s_hat']
log10_lam = environment['log10_lam']
n = environment['n']
k = environment['k']
Ne = environment['Ne']
paths = environment['paths']

# Simulate and fit the model with the specified seed
res = sim_and_fit(
    {"s": s_hat}, 
    seed=seed, 
    lam=10 ** log10_lam, 
    n=n, 
    k=k, 
    Ne=Ne, 
    af=paths[seed], 
    em_iterations=10
)

# Save the resulting selection coefficients to a file
output_path = os.path.join(dir_output, f"amy.s_hat.{seed}.txt")
np.savetxt(output_path, res["s_hat"])
