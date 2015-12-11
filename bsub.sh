bsub -R "rusage[mem=5000]" -n 1 -e "log.err" -o "log.out" R CMD BATCH --no-save quantile.R
bsub -R "rusage[mem=5000]" -n 1 -e "log.err" -o "log.out" R CMD BATCH --no-save power.R
bsub -R "rusage[mem=5000]" -n 6 -e "log.err" -o "log.out" R CMD BATCH --no-save failsim.R
bsub -R "rusage[mem=5000]" -n 1 -e "log.err" -o "log.out" R CMD BATCH --no-save alternatives.R
bsub -R "rusage[mem=5000]" -n 1 -e "log_de.err" -o "log_de.out" R CMD BATCH --no-save de_analyzer.R

