for i in {1..864}
do
  bsub -q long_normal -n 1 -R "span[hosts=1]" -o job_%J.out -e job_%J.err -W 24:00 -M 2G "Rscript Simulation_SMR.R $i"
done
