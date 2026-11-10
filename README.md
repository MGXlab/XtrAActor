# XrAActor
This nextflow pipeline will help you identify which proteins and amino acid residues are important for the biome prediction of proteins. That is to say - we can investigate how genomes and proteins are adapted to different environments or phenotypes through their amino acid adaptations. 

You only need 2 "environments" + genomes or protein clusters/families  
e.g.  
biome1 = freshwater  
biome2 = marine  

the "environment" can also be a phenotype instead, for example, "disease" and "healthy"

**This pipeline is designed to be run with an HPC cluster like slurm** - change the `nextflow.config` file accordingly if using a different HPC or running locally (just ask your favourite Agent to do this)

## Overall workflow

![images/pipeline_1.png](https://github.com/MGXlab/XtrAActor/blob/main/images/pipeline_1.png)

## Installation (WIP)

### External programs used - stay tuned - I'm working on an easy way to get you all of these

You will require the following tools for running this pipeline and please make sure to cite them!

Pharokka:  
Pharokka database:  
Phold:  
Bakta:  
mmseqs2:  
foldseek with gpu:  
foldmason:   
mafft:  
phyloglm:   
Fasttree2-MP:   
Diamond Blastp:   

### conda environments required

There are two yml files included you can use to install all the python dependencies for both nextflow and for the scripts

```bash
conda env create -f XtAActor_env.yml
conda env create -f working-env-nf.yml
```

## Running the pipeline
There are two entry points to the pipeline

1. from genomes 
2. from protein clusters

### 1. Running with genomes

Create a new folder for your run, go into it, copy and edit the sbatch and run from there

```bash
mkdir my_run
cd my_run
cp ../run_xtraactor_genomes.sbatch ./

# create your input folders
mkdir -p data/fastas
mkdir -p data/fastas
```

Your input should be formatted like this:
```bash
genomes_test_run_3/
├── data
│   ├── fastas
│   │   ├── marine.fna
│   │   └── freshwater.fna
│   └── metadata
│       └── freshwater_marine_metadata.tsv 
```
See the examples folder for what these files should look like.  
[data/fastas](https://github.com/MGXlab/XtrAActor/tree/main/Examples/input_examples/from_genomes/data/fastas)  
[data/metadata](https://github.com/MGXlab/XtrAActor/tree/main/Examples/input_examples/from_genomes/data/metadata)  

Change the parameters in the sbatch file accordingly. You will need to change the paths for the inputs etc

Required params:

```bash
  --entry_point 'genome'
  --genome_dir path/to/data/fastas
  --metadata_file path/to/data/metadata/metadata.tsv
  --biome1 "Environmental;Aquatic;Freshwater"      # please make sure this matches whatever is in your metadata.tsv
  --biome2 "Environmental;Aquatic;Marine"
  --colour1 "#91ff00"
  --colour2 "#ff7f0e"               # colour1 and 2 are for biome1 and 2, for all the plots
  --out_dir results_out             # anyway it is the default
```

### 2. Running with protein clusters

Create a new folder for your run, go into it, copy and edit the sbatch and run from there

```bash
mkdir my_run
cd my_run
cp ../run_xtraactor_genomes.sbatch ./

# create your input folders (in this case, you'll also create the main output folder)
mkdir -p results_out/cluster_biomes
mkdir -p results_out/cluster_fastas
```

Your input should be formatted like this directly in the `results_out` directory:

```bash
├── results_out
│   ├── cluster_biomes
│   │   ├── cluster_1022537_biome.tsv
│   │   ├── cluster_1331451_biome.tsv
│   ├── cluster_fastas
│   │   ├── cluster_1022537.fasta
│   │   ├── cluster_1331451.fasta
```
See the examples folder for what these files should look like.

[cluster_fastas](https://github.com/MGXlab/XtrAActor/tree/main/Examples/input_examples/from_protein_clusters/results_out/cluster_fastas)  
[cluster_biomes](https://github.com/MGXlab/XtrAActor/tree/main/Examples/input_examples/from_protein_clusters/results_out/cluster_biomes)

Change the parameters in the sbatch file accordingly. You will need to change the paths for the inputs etc

Required params:

```bash
  --entry_point 'clusters'
  --metadata_dir results_out/cluster_biomes
  --fasta_dir results_out/cluster_fastas
  --biome1 "Environmental;Aquatic;Freshwater"      # please make sure this matches whatever is in your metadata.tsv
  --biome2 "Environmental;Aquatic;Marine"
  --colour1 "#91ff00"
  --colour2 "#ff7f0e"                         # colour1 and 2 are for biome1 and 2, for all the plots
  --out_dir results_out                       # anyway it is the default
```

## Welcome to parameter heaven

```bash
 N E X T F L O W   ~  version 25.04.7

Launching `main.nf` [desperate_volhard] DSL2 - revision: b01f1f5cd4


    ─────────────·◆·─◆─·◆·─────────────────────────────────────────
       Xtr·AA·ctor  ·  extracting amino-acid features across biomes
    ─────────────────────────────────────────·◆·─◆─·◆·─────────────

    Usage:
        # phage genomes (one multi-fasta per biome in --genome_dir)
        nextflow run main.nf --entry_point genome --genome_dir genomes \
            --metadata_file metadata.tsv --genome_type phage \
            --biome1 "..." --biome2 "..." --colour1 "#0bd600" --colour2 "#8e00fa"

        # one biome from genomes, the other from already dereplicated proteins
        nextflow run main.nf --entry_point genome --genome_dir genomes_new_biome \
            --precomputed_proteins "saved/freshwater_derep_rep_seq.fasta" \
            --metadata_file metadata.tsv --genome_type phage \
            --biome1 "..." --biome2 "..." --colour1 "#0bd600" --colour2 "#8e00fa"

        # non-phage genomes
        nextflow run main.nf --entry_point genome --genome_dir genomes \
            --metadata_file metadata.tsv --genome_type non-phage \
            --biome1 "..." --biome2 "..." --colour1 "#0bd600" --colour2 "#8e00fa"

        # pre-built clusters
        nextflow run main.nf --entry_point clusters --fasta_dir cluster_fastas \
            --metadata_dir metadata --biome1 "..." --biome2 "..." --colour1 "#0bd600" --colour2 "#8e00fa"

    General:
        --entry_point          'genome' or 'clusters' (required)                  [null]
        --out_dir              Output directory                                   [output_dir_1]
        --help                 Print this help and exit

    Entry point 'genome':
        --genome_dir           Directory with one multi-fasta per biome (*.fna, *.fa, *.fasta);
                               the file name is the biome set name                [null]
        --metadata_file        Combined metadata table: genome ID <whitespace> biome, all biomes
                                                                                  [null]
        --precomputed_proteins Already dereplicated protein fasta(s) of other biomes, added directly
                               to the clustering (glob allowed, put it in quotes) [null]
        --genome_type          'phage' (Pharokka) or 'non-phage' (Bakta) (required)
                                                                                  [null]
        --gene_predictor       Phage only: 'phanotate', 'prodigal' or 'prodigal-gv'
                                                                                  [prodigal-gv]
        --pharokka_db          Pharokka database directory                        [path/to/pharokka_db]
        --pharokka_args        Extra Pharokka arguments        [--meta --skip_extra_annotations --skip_mash]
        --bakta_db             Bakta database directory                         [path/to/bakta/v20240516/db]
        --bakta_args           Extra Bakta arguments                              [--skip-plot]
        --dereplicate          Dereplicate proteins within each biome first       [true]
        --derep_id             Dereplication min. sequence identity (fraction)    [0.99]
        --derep_cov            Dereplication coverage (fraction)                  [0.9]
        --id                   MMseqs2 clustering min. sequence identity (fraction)
                                                                                  [0.5]
        --cov                  MMseqs2 clustering coverage (fraction)             [0.5]

    Entry point 'clusters':
        --fasta_dir            Per-cluster fasta files (<cluster>.faa or .fasta)  [null]
        --metadata_dir         Per-cluster tables (<cluster>_biome.tsv or <cluster>_metadata.tsv)
                                                                                  [metadata]

    Downstream analysis:
        --biome1 / --biome2    Full biome labels, exactly as in the metadata (required)
                                                                                  [null / null]
        --colour1 / --colour2  Plot colours for biome1 / biome2 (required)       [null / null]
        --min_seqs             Skip clusters with fewer sequences                 [4]
        --do_trimal            Run trimAl on the MAFFT alignments                 [false]
        --n_chunks             Max phyloglm chunks per cluster                    [15]

    Tools and environments (usually left at their defaults):
        --mafft_bin            []
        --trimal_bin           []
        --fasttree_bin         []
        --conda_sh_veo         []
        --conda_sh_lydia       []
        --conda_sh_ri          []
        --pharokka_env         [pharokka]
        --bakta_env            [bakta_v1.9.3]
        --shared_bin           []

    Values in [brackets] are the current defaults (null = not set).
    Resources (cpus, memory, time, partitions) are set per process in nextflow.config.
```