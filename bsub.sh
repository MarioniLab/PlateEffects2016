cd quantile
bsub -R "rusage[mem=5000]" -n 1 -e "log.err" -o "log.out" R CMD BATCH --no-save quantile.R
cd -

cd quantile/power
bsub -R "rusage[mem=5000]" -n 1 -e "log.err" -o "log.out" R CMD BATCH --no-save test.R
cd -

cd variability/power
bsub -R "rusage[mem=5000]" -n 1 -e "log.err" -o "log.out" R CMD BATCH --no-save power.R
cd -

cd variability/error
bsub -R "rusage[mem=5000]" -n 1 -e "log.err" -o "log.out" R CMD BATCH --no-save error.R
cd -

cd power
bsub -R "rusage[mem=5000]" -n 1 -e "log.err" -o "log.out" R CMD BATCH --no-save power.R
cd -

cd failsim
bsub -R "rusage[mem=5000]" -n 6 -e "log.err" -o "log.out" R CMD BATCH --no-save failsim.R
cd -

cd alternatives
bsub -R "rusage[mem=5000]" -n 1 -e "log.err" -o "log.out" R CMD BATCH --no-save alternatives.R
cd - 
 
cd realdata
bsub -R "rusage[mem=5000]" -n 1 -e "log_de.err" -o "log_de.out" R CMD BATCH --no-save de_analyzer.R
cd - 
