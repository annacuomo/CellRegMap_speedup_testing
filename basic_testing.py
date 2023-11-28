# From https://limix.github.io/CellRegMap/usage.html

from numpy import ones
from numpy.random import RandomState

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
print(f'Association test p-value: {pv0}')

## Interaction test
pv = run_interaction(y=y, G=g, W=W, E=C, hK=hK)[0]
print(f'Interaction test p-value: {pv}')

# Effect sizes estimation
betas = estimate_betas(y=y, G=g, W=W, E=C, hK=hK)
beta_G = betas[0]                         # persistent effect (scalar)
beta_GxC = betas[1][0]                    # GxC effects (vector)

print(f'persistent genetic effect (betaG): {beta_G}')
print(f'cell-level effect sizes due to GxC (betaGxC): {beta_GxC}')
