# From https://limix.github.io/CellRegMap/usage.html

from numpy import ones
from numpy.random import RandomState
import sys
import os
import cProfile, pstats, io

# Get the directory of the current script
current_script_path = os.path.dirname(os.path.abspath(__file__))

# Get the path to the 'cellremap' directory
cellremap_speed_path = os.path.join(current_script_path, '..', '..', 'CellRegMap')

# Add the 'cellremap_Speed' directory to the Python path
sys.path.append(cellremap_speed_path)

from cellregmap import run_association, run_interaction, estimate_betas

random = RandomState(10)
n = 30                               # number of samples (cells)
p = 5                                # number of individuals
k = 4                                # number of contexts
y = random.randn(n, 1)               # outcome vector (expression phenotype, one gene only)
C = random.randn(n, k)               # context matrix (cells by contexts/factors)
W = ones((n, 1))                     # intercept (covariate matrix)
hK = random.randn(n, p)              # decomposition of kinship matrix (K = hK @ hK.T)
g = 1.0 * (random.rand(n, 1) < 0.2)  # SNP vector

# ## Association test
# pv0 = run_association(y=y, G=g, W=W, E=C, hK=hK)[0]

# ## Interaction test
# pv1 = run_interaction(y=y, G=g, W=W, E=C, hK=hK)[0]

# turn profiling on
pr = cProfile.Profile()
pr.enable()

betas = estimate_betas(y=y, G=g, W=W, E=C, hK=hK)

pr.disable()
s = io.StringIO()
sortby = 'cumulative'
ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
ps.print_stats()

beta_G = betas[0]                         # persistent effect (scalar)
beta_GxC = betas[1][0]                    # GxC effects (vector)

# save profiling results
with open("estimate_3after.log", "w") as f:
    f.write(s.getvalue())

# save both effects to log file
with open("effects_3after.log", "a") as f:
    f.write(f'persistent genetic effect (betaG): {beta_G}\n')
    f.write(f'cell-level effect sizes due to GxC (betaGxC): {beta_GxC}\n')


# ## Estimate betas
# betas = estimate_betas(y=y, G=g, W=W, E=C, hK=hK)
