# DE analyses of scRNA-seq data with plate effects

To run the simulation and analysis code:

1. Download the count table from http://www.ebi.ac.uk/teichmann-srv/espresso, and store it in a (newly created) `reference/ESpresso` subdirectory.
2. Run `reference/submitter.sh` to construct a simulation function based on real data.
3. Run `simulations/submitter.sh` to perform the simulations for type I error control and power.
This assumes you have an LFS system, otherwise remove the `bsub` preamble in front of each job to run it in a standard fashion.
The `make_*_images.R` scripts are used to make plots of the simulation results.
4. Run `realdata/submitter.sh` to analyze the mESC data. 
Individual R scripts in `realdata` are also responsible for generating plots -- `process_real.R` to generate barplots, `plottop_ESpresso.R` for the top DE genes, and `rank_ESpresso.R` to get the top DE pluripotency factors.

The `manuscript` directory contains all the LaTeX source code for the manuscript.
This can be compiled with `make`.

The other directories contain older analysis code that was not used in the final manuscript.
Nonetheless, some things may be interesting, e.g., `failsim/extra/dc_diagnosis.R` for a study of `duplicateCorrelation` behaviour, `variability` to examine how the results are affected by different types of variances.
