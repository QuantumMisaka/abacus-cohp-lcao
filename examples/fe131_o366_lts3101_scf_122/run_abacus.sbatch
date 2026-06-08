#!/bin/bash
#SBATCH --job-name=fe-cohp-hs
#SBATCH --partition=4V100
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --gpus-per-node=4
#SBATCH --qos=rush-gpu
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

set -euo pipefail

source /opt/sai_config/mps_mapping.d/${SLURM_JOB_PARTITION}.bash
module purge
module load abacus/LTSv3.10.1-sm70-auto

export OMP_NUM_THREADS=8
export OMP_PROC_BIND=spread
export OMP_PLACES=cores

nvidia-smi dmon -s pucvmte -o T > nvdmon_job-${SLURM_JOB_ID}.log &
mpirun -np $SLURM_NTASKS -map-by $MAP_OPT -mca coll_hcoll_enable 0 abacus > abacus.stdout 2> abacus.stderr
exit
