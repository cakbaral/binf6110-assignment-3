#!/usr/bin/env bash
#SBATCH --account=def-cottenie
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=08:00:00
#SBATCH --job-name=kraken2_human_gut_class
#SBATCH --output=kraken2_%A_%a.out
#SBATCH --error=kraken2_%A_%a.err
#SBATCH --array=0-5

module load StdEnv/2023 kraken2/2.1.6

DB=/scratch/cakbaral/binf6110_assignment3/kraken_db
INDIR=/scratch/cakbaral/binf6110_assignment3/fastq_clean
OUTDIR=/scratch/cakbaral/binf6110_assignment3/taxonomic_classification

mkdir -p "$OUTDIR"

HUMAN=(SRR8146963 SRR8146952 SRR8146955 SRR8146975 SRR8146956 SRR6367588)
INPUT=${HUMAN[$SLURM_ARRAY_TASK_ID]}
R1="$INDIR/${INPUT}_1_clean.fastq.gz"
R2="$INDIR/${INPUT}_2_clean.fastq.gz"

kraken2 \
  --db "$DB" \
  --confidence 0.15 \
  --gzip-compressed \
  --paired "$R1" "$R2" \
  --report "$OUTDIR/${INPUT}_kraken2.report" \
  --output "$OUTDIR/${INPUT}_kraken2.kraken" \
  --threads $SLURM_CPUS_PER_TASK
