#!/bin/bash
#SBATCH --account=BORODAVKA-SL3-CPU	# (-A)
#SBATCH --partition=cclake	# (-p)
#SBATCH --nodes=1 # (-N)
#SBATCH --ntasks=4  # (-n)
#SBATCH --cpus-per-task=8	# (-c)
#SBATCH --time=12:00:00	# (-t)
#SBATCH --mem=64GB

if [ $# -eq 0 ]; then
  echo "Error: Not enough arguments specified"
  echo "Usage: slurm_submit_tosca_SPLASH.sh <output_dir> <arg1> ..."
  exit 1
fi

OUTDIR=$1

shift

source ~/.bashrc

conda activate nextflow_env

nextflow run main.nf -ansi-log false --outdir ${OUTDIR} ${@}

mv slurm-${SLURM_JOB_ID}.out ${OUTDIR}