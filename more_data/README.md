# Simulated Dataset Testing
This folder contains the scripts and profiling data logs for association test and interaction test

## Association Test
### Random Variants Selection
* [random_selection_association.py](./random_selection_association.py): Randomly select variants from the Genotype (G) file and perform association test. The number of variants to be selected can be adjusted by changing the second parameter in the line below:
      <pre>selected_variants = np.random.choice(G.variant.size, <b>100</b>, replace=False)</pre>
### Linkage Disequilibrium (LD) Pruning Variants Selection
* [LDpruning_association.py](./LDpruning_association.py): Perform LD pruning on the Genotype (G) file and perform association test. Note: the variant size can be adjusted by changing the `num_variants` parameter in the `ld_prune` function:
      <pre>def ld_prune(G, <b>num_variants=100</b>, variance_threshold=1e-5)</pre>
## Interaction Test
* [LDpruning_ineraction.py](./LDpruning_interaction.py): Perform LD pruning on the Genotype (G) file and perform interaction test. Note: the variant size can be adjusted in the same way as above.

## Note
All the scripts above are designed to test one signle gene at a time. The default gene is specified in the line below:
```python
y = y.sel(gene="gene_1")
```
 If you want to test different ones, you can simply change the gene name in the bracket.