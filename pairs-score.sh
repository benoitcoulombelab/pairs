#!/bin/bash
#SBATCH --account=def-coulomb
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --mail-type=NONE
#SBATCH --output=pairs-score-%A.out

# exit when any command fails
set -e

# load required modules
if [ -n "$CC_CLUSTER" ]
then
  module -q purge
  module -q load pairs/1.0
fi

ranked_0_files_raw=$(find . -type f -name "ranked_0.pdb")
ranked_0_files=()
for ranked_0 in ${ranked_0_files_raw}
do
  ranked_0_files+=($ranked_0)
done

multi-interaction-score "$@" "${ranked_0_files[@]}"
