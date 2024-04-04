#!/bin/bash
#SBATCH --account=def-coulomb
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --mail-type=NONE
#SBATCH --output=pairs-%A.out

# exit when any command fails
set -e

# load required modules
if [[ -n "$CC_CLUSTER" ]]
then
  module purge
  module load StdEnv/2020
  module load nextflow/22.10.6
fi
if [[ "beluga" == "$CC_CLUSTER" ]]
then
  ulimit -v 40000000
fi

export NXF_OPTS="-Xms500M -Xmx8000M"

nextflow run pairs.nf \
    --fasta "fasta_pairs/*.fasta" \
    --outdir "$PWD/output" \
    -c alliancecan.config \
    "$@"
