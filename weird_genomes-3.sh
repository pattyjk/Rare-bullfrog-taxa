#!/bin/bash
#SBATCH --job-name=jliv_anno
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=00-10:00:00
#SBATCH --mem=64gb
#SBATCH --error=%x-%j.err
#SBATCH --output=%x-%j.out
#SBATCH --partition=Intel6240,Intel6248
##SBATCH --partition=DGXA100
#SBATCH --mail-type=ALL
#SBATCH --mail-user=patrick.kearns@umb.edu

cd /hpcstor6/scratch01/p/patrick.kearns/genomes


# ---- CheckM ----
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate checkm
export CHECKM_DATA_PATH=/home/patrick.kearns/checkm
checkm lineage_wf aggregated_genomes checkM_results_$SLURM_JOB_ID -x fasta -t 48 --pplacer_threads 4 --tab_table -f checkM_results_$SLURM_JOB_ID.tsv
conda deactivate

#annotate with prokka
for i in *.fasta; do
    name="${i%.fasta}"
    prokka "$i" --outdir "$name" --prefix "$name" --force
done