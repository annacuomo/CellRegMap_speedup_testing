## Test data:

* ```seed_1_100_nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_Poisson.txt``` contains a few different pieces of information:
  * the first column is the individual ID -- based on this you can construct a block diagonal kinship file following [these instructions](https://github.com/annacuomo/CellRegMap_analyses/blob/main/endodiff/preprocessing/block_diagonal_K.ipynb)
  * columns 2-5 include covariates, or contexts -- use these as your context file
  * the remaining columns represent the simulated single-cell expression file for 100 genes
* ```n.indep_100_n.cell_1.bed```, ```n.indep_100_n.cell_1.bim``` and ```n.indep_100_n.cell_1.fam``` are simulated genotype files, in the [plink](https://en.wikipedia.org/wiki/PLINK_(genetic_tool-set)) format

## How to load the data

### Genotypes

Open genotype files using the ```pandas_plink``` module:

```python
from pandas_plink import read_plink1_bin
plink_file = "n.indep_100_n.cell_1.bed"
G = read_plink1_bin(plink_file)
```

### Text file containing all other info

Load using pandas:

```python
import pandas as pd

txt_file = "seed_1_100_nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_Poisson.txt"
txt_df = pd.read_csv(txt_file, sep="\t")
```

### Gene expression

```python
pheno_df = txt_df.iloc[:,5:]
```

### Contexts

```python
contexts_df = txt_df[["X1","X2","pf1","pf2"]]
```

#### Kinship

First, extract cells and donor ids:

```python
# sample mapping file mapping cells to individuals
txt_df['cell'] = ["cell_" + str(val) for val in txt_df.index.values]
smf_df = txt_df[["IND_ID","cell"]]
```

build useful functions:

```python
from numpy import array, split, cumsum, zeros, append
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
```

and apply

```python
# indices for each group of cells (group=individual)
groups = get_groups_from_smf(smf_df)
hK = get_block_hK_from_groups(groups)
```

plot to check that full kinship matrix is positive and symmetric:

```python
import matplotlib.pyplot as plt
K = hK @ hK.T
plt.matshow(K)
```

Then follow e.g., https://github.com/annacuomo/CellRegMap_analyses/blob/main/endodiff/usage/notebooks/run_association_chr19_test.ipynb

