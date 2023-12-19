from pandas_plink import read_plink1_bin
import pandas as pd
from numpy import array, split, cumsum, zeros, append, ones
import matplotlib.pyplot as plt
import xarray as xr
from limix.qc import quantile_gaussianize
from cellregmap import run_association
import itertools
import cProfile, pstats, io

plink_file = "../test_data/n.indep_100_n.cell_1.bed"
G = read_plink1_bin(plink_file)

txt_file = "../test_data/seed_1_100_nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_Poisson.txt"
txt_df = pd.read_csv(txt_file, sep="\t")
# get the data for a1 in IND_ID
txt_df_a1 = txt_df[txt_df["IND_ID"] == "a1"]

pheno_df_a1 = txt_df_a1.iloc[:,5:]

contexts_df_a1 = txt_df_a1[["X1","X2","pf1","pf2"]]

txt_df_a1['cell'] = ["cell_" + str(val) for val in txt_df_a1.index.values]
smf_df = txt_df_a1[["IND_ID","cell"]]

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

C = xr.DataArray(contexts_df_a1.values, dims=('cell', 'pc'), coords={'cell': contexts_df_a1.index.values, 'pc': contexts_df_a1.columns.values})
C = C.sel(cell=smf_df.index.values)
C = quantile_gaussianize(C)

y = xr.DataArray(pheno_df_a1.values, dims=["cell", "gene"], coords={"cell": pheno_df_a1.index.values, "gene": pheno_df_a1.columns.values})
y = y.sel(gene="gene_1")
y = y.values.reshape(y.shape[0],1)

G_sel = xr.DataArray(G.values, dims=('sample', 'variant'), coords={'sample': G['sample'].values, 'variant': G['variant'].values})
G_sel = G_sel.sel(sample=smf_df['IND_ID'].values)
GG = G_sel.values

W = ones((100,1))

# turn profiling on
pr = cProfile.Profile()
pr.enable()

pv = run_association(y=y, G=GG, W=W, E=C.values[:,0:10], hK=hK)[0]

pr.disable()
s = io.StringIO()
sortby = 'cumulative'
ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
ps.print_stats()
with open("a1_new.log", 'w') as f:
    print(s.getvalue(), file=f)

pv1 = pd.DataFrame({"sample":G_sel.sample.values,
               "pv":pv,
               "variant":G_sel.variant.values})
pv1.to_csv("a1.csv")
