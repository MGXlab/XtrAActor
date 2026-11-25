#!/bin/bash
#SBATCH --job-name=nextflow
#SBATCH --output=nf_main%j.out
#SBATCH --error=nf_main%j.err
#SBATCH --cpus-per-task=1
#SBATCH --partition=short,interactive,fat,long,standard
#SBATCH --mem=15G
#SBATCH --time 3:00:00

source /path/to/miniconda3/.../conda.sh

conda activate XtrAActor_env

# nextflow run pipeline_v2.nf -params-file user_inputs.yaml

nextflow run -resume pipeline_merged.nf \
  --fasta_dir  /genes_test/fastas\
  --metadata_dir genes_test/metadata \
  --biome1 "Environmental;Aquatic;Freshwater;Lake" \
  --biome2 "Environmental;Aquatic;Thermal springs;Hot (42-90C)" \
  --colour1 "#1f77b4" \
  --colour2 "#ff7f0e"\
  --out_dir genes_test_output


# BIOMES
# biome1 is always "Environmental;Aquatic;Freshwater;Lake"
# biome2:
#         "Environmental;Aquatic;Marine;Oceanic"
#         "Environmental;Aquatic;Thermal springs;Hot (42-90C)"
#         "Environmental;Aquatic;Non-marine Saline and Alkaline;Hypersaline"
#         "Environmental;Aquatic;Non-marine Saline and Alkaline;Alkaline"