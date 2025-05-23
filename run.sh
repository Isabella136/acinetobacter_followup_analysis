#!/bin/bash

#SBATCH --job-name=acin_followup
#SBATCH --output=log/acin_followup.out.%j
#SBATCH --error=log/acin_followup.err.%j
#SBATCH --time=24:00:00
#SBATCH --qos=high
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --mem=16gb
#SBATCH --partition=cbcb
#SBATCH --account=cbcb

eval "$(conda shell.bash hook)"
conda activate snakemake

snakemake --directory workdir --unlock
snakemake --directory workdir --touch
snakemake --profile ./profile/default