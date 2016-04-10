mkdir ESpresso

###########################
# Error simulations

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

# This takes intolerably long on a single run, so I'm going to split it up across many cores.
tmpdir=ESpresso/tmp7
mkdir $tmpdir
for i in {1..10}
do
    if [ $i -eq 1 ]
    then
        ncores=6
        longrun=TRUE
    else
        ncores=1
        longrun=FALSE
    fi
    mkdir $tmpdir/$i
    bsub -R "rusage[mem=10000]" -n $ncores -e "${tmpdir}/log_fail7_${i}.err" -o "${tmpdir}/log_fail7_${i}.out" \
        R CMD BATCH --no-save "--args seed=${i}0000 scenario=7 niter=1 do.long=${longrun} indir='../reference/results_ESpresso' outdir='ESpresso/${tmpdir}/${i}'"\
        failsim.R $tmpdir/failsim_7_${i}.Rout
done

############################
# Power simulations.

for logfile in $(find ESpresso/ | grep "log_pow")
do
    rm $logfile
done

powargs="seed=10000 indir='../reference/results_ESpresso' outdir='ESpresso/results_power'"
bsub -R "rusage[mem=5000]" -n 1 -e "ESpresso/log_pow1.err" -o "ESpresso/log_pow1.out" R CMD BATCH --no-save "--args scenario=1 $powargs" power.R ESpresso/power_1.Rout
bsub -R "rusage[mem=5000]" -n 1 -e "ESpresso/log_pow2.err" -o "ESpresso/log_pow2.out" R CMD BATCH --no-save "--args scenario=2 $powargs" power.R ESpresso/power_2.Rout
bsub -R "rusage[mem=5000]" -n 1 -e "ESpresso/log_pow3.err" -o "ESpresso/log_pow3.out" R CMD BATCH --no-save "--args scenario=3 $powargs" power.R ESpresso/power_3.Rout
bsub -R "rusage[mem=5000]" -n 1 -e "ESpresso/log_pow4.err" -o "ESpresso/log_pow4.out" R CMD BATCH --no-save "--args scenario=4 $powargs" power.R ESpresso/power_4.Rout

# Results aren't under version control, so they need to be pulled down manually:
# rsync -azv cruk:lustre/PlateEffects/simulations/ESpresso/ ESpresso/
