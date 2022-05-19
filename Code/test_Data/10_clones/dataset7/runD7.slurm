#!/bin/bash
#SBATCH --job-name=d7_50
#SBATCH --partition=commons
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --threads-per-core=1
#SBATCH --time=08:00:00
#SBATCH --mail-user=hz22@rice.edu

module load Java/1.8.0_45
export SLURM_SUBMIT_DIR=/home/hz22/SiCloneFit/simData/noDoublet/100_cells/50_sites/10_clones/dataset7
java -jar ~/SiFit.jar -m 100 -n 50 -fp 0.05 -df 1 -ipMat ${SLURM_SUBMIT_DIR}/noisy_genotype_dataset7.txt -trueTree ${SLURM_SUBMIT_DIR}/Orig_tree_dataset7.txt -r 4 -fn 0.2 -iter 200000 -printIter 20000 > SiFit_rslt_noisy_dataset7.txt 
