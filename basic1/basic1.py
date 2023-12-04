# From https://limix.github.io/CellRegMap/usage.html

from numpy import ones
from numpy.random import RandomState
import sys
import os

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

## Association test
pv0 = run_association(y=y, G=g, W=W, E=C, hK=hK)[0]

## Interaction test
pv1 = run_interaction(y=y, G=g, W=W, E=C, hK=hK)[0]

## Estimate betas
betas = estimate_betas(y=y, G=g, W=W, E=C, hK=hK)
