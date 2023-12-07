## Test data:

* ```seed_1_100_nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_Poisson.txt``` contains a few different pieces of information:
  * the first column is the individual ID -- based on this you can construct a block diagonal kinship file following [these instructions](https://github.com/annacuomo/CellRegMap_analyses/blob/main/endodiff/preprocessing/block_diagonal_K.ipynb)
  * columns 2-5 include covariates, or contexts -- use these as your context file
  * the remaining columns represent the simulated single-cell expression file for 100 genes
* ```n.indep_100_n.cell_1.bed```, ```n.indep_100_n.cell_1.bim``` and ```n.indep_100_n.cell_1.fam``` are simulated genotype files, in the [plink](https://en.wikipedia.org/wiki/PLINK_(genetic_tool-set)) format

## How to load the data

### Genotypes

```python
from pandas_plink import read_plink1_bin
plink_file = "n.indep_100_n.cell_1.bed"
G = read_plink1_bin(plink_file)
```

### Kinship


### Contexts


### Gene expression

