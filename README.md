# DE analyses of scRNA-seq data with plate effects

To run the code:

1. Download the count table from http://www.ebi.ac.uk/teichmann-srv/espresso, and store it in `reference/ESpresso`.
2. Run `reference/submitter.sh` to construct the simulation function based on real data.
3. Run `simulations/submitter.sh` to perform the simulations. This assumes you have an LSF system, otherwise remove the `bsub` preamble in front of each job.
4. Run `realdata/submitter.sh` to analyze the mESC data. Individual R scripts in `realdata` are also responsible for generating plots -- `process_real.R` to generate barplots, `plottop_ESpresso.R` for the top DE genes, and `rank_ESpresso.R` to get the top DE pluripotency factors.
