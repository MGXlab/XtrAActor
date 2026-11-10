# XrAActor
This pipeline will help you identify which amino acid residues are important for the biome prediction of proteins

## Overall workflow

![images/pipeline_1.png](https://github.com/MGXlab/XtrAActor/blob/main/images/pipeline_1.png)

## Installation

### External programs used

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

## Running the pipeline
There are two entry points to the pipeline

1. from genomes 
2. from protein clusters

### 1. Running with genomes
Your input should be formatted like this:
```
genomes_test_run_3/
├── data
│   ├── fastas
│   │   ├── marine.fna
│   │   └── freshwater.fna
│   └── metadata
│       └── freshwater_marine_metadata.tsv 
```
### 2. Running with protein clusters
Your input should be formatted like this directly in the `results_out` directory:

```
├── genomes_test_run_3_out
│   ├── RF_results
│   │   └── raw
│   ├── annotation
│   │   └── pharokka
│   ├── cluster_biomes
│   │   ├── cluster_1022537_biome.tsv
│   │   ├── cluster_1331451_biome.tsv
│   ├── cluster_fastas
│   │   ├── cluster_1022537.fasta
│   │   ├── cluster_1331451.fasta
```
See the examples folder for what these files should look like.

[cluster_fastas]()  
[cluster_biomes]()