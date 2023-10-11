# AlphaFold-pairs

### This software scores protein-protein interactions using AlphaFold-Multimer's output

[See AlphaFold-pairs paper](https://www.biorxiv.org/content/10.1101/2023.08.29.555151v1)

#### Contents

* [Install AlphaFold-pairs](#install)
* [Running AlphaFold-pairs](#running-alphafold-pairs)
    * [Prepare FASTA files](#prepare-fasta-files)
        * [Proteins with an ambigous amino acid](proteins-with-an-ambigous-amino-acid)
    * [Running AlphaFold-Multimer using Nextflow](#running-alphafold-multimer-using-nextflow)
    * [Scoring protein-protein interactions](#scoring-protein-protein-interactions)
* [Configuring AlphaFold-Multimer](#configuring-alphafold-multimer)
    * [Patching AlphaFold-Multimer](#patching-alphafold-multimer)
        * [Running only prepare or alphafold step](#running-only-prepare-or-alphafold-step)
    * [Alternative to patching AlphaFold-Multimer](#alternative-to-patching-alphafold-multimer)
* [Protein complexes](#protein-complexes)
    * [Interaction scores matrix](#interaction-scores-matrix)
    * [Zscore matrix](#zscore-matrix)
* [Extract interfaces of a protein-protein interaction](#extract-interfaces-of-a-protein-protein-interaction)
    * [Consensus interface](#consensus-interface)
* [Utilities](#utilities)
    * [Split FASTA file](#split-fasta-file)

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

If you have multiple baits or targets, you can use the `fasta-pairs` script.

```shell
mkdir fasta_pairs
fasta-pairs --baits baits.fasta --targets targets.fasta --output fasta_pairs -u -i
```

#### Proteins with an ambigous amino acid

AlphaFold does not work with proteins that have an ambigus amino acid.

You can delete FASTA files containing an ambigous amino acid using the following command:

```shell
delete-fasta -s fasta_pairs/*.fasta
```

To move the FASTA files instead of deleting them:

```shell
delete-fasta -s -b backup_folder fasta_pairs/*.fasta
```

### Running AlphaFold-Multimer using Nextflow

```shell
nextflow run alphafold-pairs.nf \
    --fasta 'fasta_pairs/*.fasta' \
    --outdir "$PWD/output"
```

### Scoring protein-protein interactions

```shell
multi-interaction-score -w \
    -o interaction-scores.txt \
    output/alphafold/**/ranked_o.pdb
```

For a single interaction:

```shell
interaction-score -w ranked_o.pdb
```

## Configuring AlphaFold-Multimer

The Nextflow pipeline relies on a script called `alphafold.sh` in the current directory to run AlphaFold.

An example of an implementation of this script using [Slurm](https://slurm.schedmd.com) is provided in
the [alphafold.sh](alphafold.sh) file.

### Patching AlphaFold-Multimer

The default Nextflow pipline expects to run AlphaFold in two steps.

1. Prepare
    * Run external programs to produce AlphaFold's input
    * Uses 8 CPUs and no GPU
2. AlphaFold
    * Run AlphaFold-mutlimer's AI program based on "prepare's" output
    * Uses 1 CPU and 1 GPU

To make it possible to run the "prepare" step separately, you must patch AlphaFold using the patch file available here:
[Patches for AlphaFold](https://github.com/benoitcoulombelab/modules/tree/main/alphafold)

#### Running only prepare or alphafold step

```shell
nextflow run alphafold-pairs.nf \
    --step prepare \
    --fasta 'fasta_pairs/*.fasta' \
    --outdir "$PWD/output"
```

```shell
nextflow run alphafold-pairs.nf \
    --step alphafold \
    --fasta 'fasta_pairs/*.fasta' \
    --outdir "$PWD/output"
```

This will allow you to optimize both steps independently to optimize HPC usage.

### Alternative to patching AlphaFold-Multimer

If you don't care about optimizing CPU/GPU usage, you can run the alternative pipeline that will run AlphaFold
normally (without patching).

```shell
nextflow run alphafold-pairs-noprepare.nf \
    --fasta 'fasta_pairs/*.fasta' \
    --outdir "$PWD/output"
```

## Protein complexes

To find all binary protein-protein interactions of a complex, start by running the normal AlphaFold-pairs pipeline.
[Running AlphaFold-pairs](#running-alphafold-pairs)

Once AlphaFold-Multimer completes for all FASTA files, use `multi-interaction-score` script.

```shell
multi-interaction-score -w \
    -o interaction-scores.txt \
    output/alphafold/**/ranked_o.pdb
```

### Interaction scores matrix

You can generate a matrix of protein-protein interaction scores using:

```shell
score-matrix -o score-matrix.txt -u \
    interaction-scores.txt
```

### Zscore matrix

You can also compute Zscores when generating the protein-protein interaction scores matrix:

```shell
score-matrix -o score-matrix.txt -u -z \
    interaction-scores.txt
```

## Extract interfaces of a protein-protein interaction

You can extract the interface(s) of a protein-protein interaction using `interaction-score`.

```shell
interaction-score -r 6 -R residue-pairs.txt -A atom-pairs.txt ranked_0.pdb
```

* The `residue-pairs.txt` file contains the residue-pairs that have an atom at a distance below the radius (`-s`
  parameter).
* The `atom-pairs.txt` file contains the atoms-pairs that are at a distance below the radius (`-s` parameter).

### Consensus interface

If you analyse similar bait-target pairs, you may want to extract a consensus interface.

```shell
consensus-interface -r residue-pairs/*.txt -o consensus-residue-pairs.txt \
    -c 0.5 -n "([\w-]+)__([\w-]+)" \
    -b aligned_baits.txt -t aligned_targets.txt -f clustal
```

Since baits and targets may slightly differ in their sequence, you can provide an alignment file to adjust the index of
the residues in different baits/targets.

## Utilities

### Split FASTA file

To split a FASTA file into a single file per sequence:

```shell
split-fasta -o fasta_per_sequence sequence.fasta
```
