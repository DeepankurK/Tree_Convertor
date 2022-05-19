#!/bin/bash
#SBATCH --job-name=d7_50
#SBATCH --partition=commons
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --threads-per-core=1
#SBATCH --time=02:00:00
#SBATCH --mail-user=hz22@rice.edu

export SLURM_SUBMIT_DIR=/home/hz22/SiCloneFit/simData/noDoublet/100_cells/50_sites/10_clones/dataset7/SCITE_rslts
~/software/SCITE/scite -i ${SLURM_SUBMIT_DIR}/SCITE_noisy_genotype_dataset7.csv -m 100 -n 50 -fd 0.05 -r 4 -ad 0.2 0.2 -l 1000000 -e 0.25 -a 
python ~/SiCloneFit/helperFiles/processSCITEgv.py SCITE_noisy_genotype_dataset7_ml0.gv 100 50 
