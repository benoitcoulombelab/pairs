#!/bin/bash
#SBATCH --account=def-coulomb
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --mail-type=NONE
#SBATCH --output=alphafold-pairs-%A.out

# exit when any command fails
set -e

# load required modules
if [[ -n "$CC_CLUSTER" ]]
then
  module purge
  module load nextflow/22.10.6
  module load alphafold-pairs/1.0
fi

export NXF_OPTS="-Xms500M -Xmx8000M"

nextflow run alphafold-pairs.nf \
    --fasta "fasta_pairs/*.fasta" \
    --outdir "$PWD/output" \
    -c alliancecan.config \
    "$@"
