#!/bin/bash
#SBATCH --account=def-coulomb
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --mail-type=NONE
#SBATCH --output=deepproteinconnector-%A.out

# exit when any command fails
set -e

# load required modules
if [[ -n "$CC_CLUSTER" ]]
then
  module purge
  module load nextflow/22.10.6
  module load deepproteinconnector/1.0
fi

export NXF_OPTS="-Xms500M -Xmx8000M"

nextflow run deepproteinconnector.nf \
    --baits baits.fasta --targets targets.fasta \
    --outdir "$PWD/output" \
    -c alliancecan.config \
    "$@"
