from pandas_plink import read_plink1_bin
import pandas as pd
from numpy import array, split, cumsum, zeros, append, ones
import matplotlib.pyplot as plt
import xarray as xr
from limix.qc import quantile_gaussianize
from cellregmap import run_interaction
import itertools
import cProfile, pstats, io
import numpy as np

plink_file = "../test_data/n.indep_100_n.cell_1.bed"
G = read_plink1_bin(plink_file)

txt_file = "../test_data/seed_1_100_nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_Poisson.txt"
txt_df = pd.read_csv(txt_file, sep="\t")

pheno_df = txt_df.iloc[:,5:]

contexts_df = txt_df[["X1","X2","pf1","pf2"]]

txt_df['cell'] = ["cell_" + str(val) for val in txt_df.index.values]
smf_df = txt_df[["IND_ID","cell"]]

def get_groups_from_smf(smf_df):
    n_samples = smf_df.shape[0]
    donors = smf_df['IND_ID'].unique()
    n_donors = len(donors)
    n_cells = array([],dtype=int)
    for donor in donors:
        n_cells = append(n_cells, array(smf_df[smf_df['IND_ID']==donor].shape[0], dtype=int))
    groups = split(range(n_samples), cumsum(n_cells))[:-1]
    return groups

import itertools
def get_block_hK_from_groups(groups):
    n_samples = len(list(itertools.chain.from_iterable(groups)))
    hM = zeros((n_samples, len(groups)))
    for i, idx in enumerate(groups):
        hM[idx, i] = 1.0
    return hM

groups = get_groups_from_smf(smf_df)
hK = get_block_hK_from_groups(groups)

K = hK @ hK.T
plt.matshow(K)

C = xr.DataArray(contexts_df.values, dims=('cell', 'pc'), coords={'cell': contexts_df.index.values, 'pc': contexts_df.columns.values})
C = C.sel(cell=smf_df.index.values)
C = quantile_gaussianize(C)

y = xr.DataArray(pheno_df.values, dims=["cell", "gene"], coords={"cell": pheno_df.index.values, "gene": pheno_df.columns.values})
y = y.sel(gene="gene_1")
y = y.values.reshape(y.shape[0],1)

# LD pruning helper functions
def is_correlated(variant, pruned_variants, threshold=0.2):
    for pruned_variant in pruned_variants:
        # Check for constant variants
        if np.std(variant) == 0 or np.std(pruned_variant) == 0:
            continue
        correlation = np.corrcoef(variant, pruned_variant)[0, 1] ** 2
        if correlation >= threshold:
            return True
    return False

def ld_prune(G, num_variants=100, variance_threshold=1e-5):
    # Start with the first variant
    selected_indices = [0]
    pruned_variants = [G[:, 0].values]

    # Sequentially check each variant
    for i in range(1, G.shape[1]):
        if len(selected_indices) >= num_variants:
            break

        candidate_variant = G[:, i].values
        if not is_correlated(candidate_variant, pruned_variants):
            pruned_variants.append(candidate_variant)
            selected_indices.append(i)

    return selected_indices

selected_variants = ld_prune(G)
G_sel = G[:, selected_variants]
G_sel = xr.DataArray(G_sel.values, dims=('sample', 'variant'), coords={'sample': G_sel['sample'].values, 'variant': G_sel['variant'].values})
G_sel = G_sel.sel(sample=smf_df['IND_ID'].values)
GG = G_sel.values

W = ones((10000,1))

# turn profiling on
pr = cProfile.Profile()
pr.enable()

pv = run_interaction(y=y, G=GG, W=W, E=C.values[:,0:10], hK=hK)[0]

pr.disable()
s = io.StringIO()
sortby = 'cumulative'
ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
ps.print_stats()

# save profiling results to ld_100_inter_before.log
with open("ld_500_inter_after.log", "w") as f:
    f.write(s.getvalue())

# save pv to 100_pv_inter_before.log up to 8 decimal places
np.savetxt("500_pv_inter_after.log", pv, fmt="%.8f")
