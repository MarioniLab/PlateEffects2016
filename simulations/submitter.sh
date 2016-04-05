mkdir ESpresso

for logfile in $(find ESpresso/ | grep "log_fail")
do
    rm $logfile
done

failargs="seed=10000 indir='../reference/results_ESpresso' outdir='ESpresso/results_failsim'"
bsub -R "rusage[mem=5000]" -n 6 -e "ESpresso/log_fail1.err" -o "ESpresso/log_fail1.out" R CMD BATCH --no-save "--args scenario=1 $failargs" failsim.R ESpresso/failsim_1.Rout
bsub -R "rusage[mem=5000]" -n 6 -e "ESpresso/log_fail2.err" -o "ESpresso/log_fail2.out" R CMD BATCH --no-save "--args scenario=2 $failargs" failsim.R ESpresso/failsim_2.Rout
bsub -R "rusage[mem=5000]" -n 6 -e "ESpresso/log_fail3.err" -o "ESpresso/log_fail3.out" R CMD BATCH --no-save "--args scenario=3 $failargs" failsim.R ESpresso/failsim_3.Rout
bsub -R "rusage[mem=5000]" -n 6 -e "ESpresso/log_fail4.err" -o "ESpresso/log_fail4.out" R CMD BATCH --no-save "--args scenario=4 $failargs" failsim.R ESpresso/failsim_4.Rout
bsub -R "rusage[mem=5000]" -n 6 -e "ESpresso/log_fail5.err" -o "ESpresso/log_fail5.out" R CMD BATCH --no-save "--args scenario=5 $failargs" failsim.R ESpresso/failsim_5.Rout
bsub -R "rusage[mem=5000]" -n 6 -e "ESpresso/log_fail6.err" -o "ESpresso/log_fail6.out" R CMD BATCH --no-save "--args scenario=6 $failargs" failsim.R ESpresso/failsim_6.Rout

for logfile in $(find ESpresso/ | grep "log_pow")
do
    rm $logfile
done

powargs="seed=10000 indir='../reference/results_ESpresso' outdir='ESpresso/results_power'"
bsub -R "rusage[mem=5000]" -n 1 -e "ESpresso/log_pow1.err" -o "ESpresso/log_pow1.out" R CMD BATCH --no-save "--args scenario=1 $powargs" power.R ESpresso/power_1.Rout
bsub -R "rusage[mem=5000]" -n 1 -e "ESpresso/log_pow2.err" -o "ESpresso/log_pow2.out" R CMD BATCH --no-save "--args scenario=2 $powargs" power.R ESpresso/power_2.Rout
bsub -R "rusage[mem=5000]" -n 1 -e "ESpresso/log_pow3.err" -o "ESpresso/log_pow3.out" R CMD BATCH --no-save "--args scenario=3 $powargs" power.R ESpresso/power_3.Rout
bsub -R "rusage[mem=5000]" -n 1 -e "ESpresso/log_pow4.err" -o "ESpresso/log_pow4.out" R CMD BATCH --no-save "--args scenario=4 $powargs" power.R ESpresso/power_4.Rout
