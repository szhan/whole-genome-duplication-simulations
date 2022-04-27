# Whole Genome Duplication Simulations

Simulate ploidy shifts along a phylogenetic tree using the BiSSE model.

## Introduction

This script is intended to be run with ``Rscript`` via commandline, so that it can be easily parallelized and run on HPCs.

The following input BiSSE parameters are required:
- [ ] ``spec0`` - Speciation rate when the state is 0.
- [ ] ``spec1`` - Speciation rate when the state is 1.
- [ ] ``ext0`` - Extinction rate when the state is 0.
- [ ] ``ext1`` - Extinction rate when the state is 1.
- [ ] ``q01`` - Rate of transition from 0 to 1.
- [ ] ``q10`` - Rate of transition from 1 to 0.

Additionally, the following input parameters are required:
- [ ] ``ntaxa`` - Number of tip taxa (termination criterion).
- [ ] ``nreps`` - Number of simulation replicates to run in parallel.
- [ ] ``out_file`` - Path to the output file.

BiSSE simulations are run with the root state set to 0.


## Dependencies
- [ ] ``dplyr``
- [ ] ``diversitree``
- [ ] ``optparse``
- [ ] ``parallel``
- [ ] ``stringr``


## Examples

The following command sets off a BiSSE simulation until a phylogeny with 10,000 tip taxa is generated.
```
Rscript --no-save simulate_histories.R --spec0=1 --spec1=0.4 --ext0=0.5 --ext1=0.50 --q01=0.5 --q10=0.65 --ntaxa=10000 --nreps=100 --out_file=s1_e3_q4.csv
```
