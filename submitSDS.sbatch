#!/bin/bash
#SBATCH --job-name=mainSDS
#SBATCH --output=logs/mainSDS.out
#SBATCH --error=logs/mainSDS.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=2000
#SBATCH --partition=tier1q

module load gcc/12.1.0
module load python/3.10.5
module load plink/2.0
module load gcta/1.94.1
module load R
source ../myVilma/bin/activate 

echo "SLURM_JOBID="$SLURM_JOBID
#rm -r .snakemake/
cat snakefile_main_SDS
snakemake --unlock -s snakefile_main_SDS
snakemake --profile cluster-setup-SDS/ -s snakefile_main_SDS --rerun-incomplete
