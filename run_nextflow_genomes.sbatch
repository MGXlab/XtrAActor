#!/bin/bash
#SBATCH --job-name=nextflow
#SBATCH --output=tmp/nf_main_%j.out
#SBATCH --error=tmp/nf_main_%j.err
#SBATCH --cpus-per-task=20
#SBATCH --partition=standard,long,fat,short,interactive,gpu-veo,gpu-test,gpu
#SBATCH --mem=40G
#SBATCH --time=12:00:00

# Head job only: it submits and monitors the SLURM tasks defined in nextflow.config.
# It must run longer than the whole pipeline, otherwise all running tasks are orphaned.


source /vast/ri65fin/miniconda3/etc/profile.d/conda.sh
conda activate nf-env
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ from genomes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

nextflow run ../main.nf -resume \
  -c ../nextflow.config \
  --entry_point genome \
  --genome_dir /vast/ri65fin/nextflow_AAcidx_v2/genomes_test_run_3/data/fastas \
  --metadata_file /vast/ri65fin/nextflow_AAcidx_v2/genomes_test_run_3/data/metadata/freshwater_alkaline_metadata.tsv \
  --genome_type phage \
  --gene_predictor prodigal-gv \
  --pharokka_db "/vast/ri65fin/pharokka_db" \
  --id 0.5 --cov 0.5 \
  --biome1  "Environmental;Aquatic;Freshwater;Lake" \
  --biome2  "Environmental;Aquatic;Non-marine Saline and Alkaline;Alkaline" \
  --colour1 "#91ff00" \
  --colour2 "#ff7f0e" \
  --out_dir genomes_test_run_3_out
