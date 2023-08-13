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

## Configuring AlphaFold

The Nextflow pipeline relies on a script called `alphafold.sh` in the current directory to run AlphaFold.

An example of an implementation of this script using [Slurm](https://slurm.schedmd.com) is provided in
the [alphafold.sh](alphafold.sh) file.

### Patching AlphaFold

The default Nextflow pipline expects to run AlphaFold in two steps.

1. Prepare
    * Run external programs to produce AlphaFold's input
    * Uses 8 CPUs and no GPU
2. AlphaFold
    * Run AlphaFold-mutlimer's AI program based on "prepare's" output
    * Uses 1 CPU and 1 GPU

To make it possible to run the "prepare" step separately, you must patch AlphaFold using the patch file available here:
[Patches for AlphaFold](https://github.com/benoitcoulombelab/modules/tree/main/alphafold)

### Alternative to patching AlphaFold

If you don't care about optimizing CPU/GPU usage, you can run the alternative pipeline that will run AlphaFold
normally (without patching).

```shell
nextflow run alphafold-pairs-noprepare.nf \
    --fasta 'fasta_pairs/*.fasta' \
    --outdir "$PWD/output"
```
