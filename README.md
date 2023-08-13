# AlphaFold-pairs

### This software tries to guess which proteins interact with a protein of interest using AlphaFold multimer

## Install

Install requirements:

* [python version 3.7.4 or newer](https://www.python.org)
* [Git](https://git-scm.com)

Install AlphaFold-pairs.

```shell
pip install git+https://github.com/benoitcoulombelab/https://github.com/benoitcoulombelab/alphafold-pairs.git
```

## Running AlphaFold-pairs

### Prepare FASTA files

```shell
mkdir fasta_pairs
fasta-pairs --baits baits.fasta --targets targets.fasta --output fasta_pairs -u -i
```

### Running AlphaFold using Nextflow

```shell
nextflow run alphafold-pairs.nf \
    --fasta 'fasta_pairs/*.fasta' \
    --outdir "$PWD/output"
```

### Scoring interactions

```shell
multi-interaction-score -w \
    -o interaction-scores.txt \
    output/alphafold/**/ranked_o.pdb
```
