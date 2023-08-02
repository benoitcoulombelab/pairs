#!/bin/bash

# exit when any command fails
set -e

fasta=$1
output=$2
step=${3:-all}
max_template_date=$(date +"%Y-%m-%d")

echo -e "\n\nRun AlphaFold with step ${step} on fasta ${fasta}\n\n"

# load required modules
if [[ -n "$CC_CLUSTER" ]]
then
  module purge
  module load StdEnv/2020
  module load gcc/9.3.0 cuda/11.7
  module load apptainer/1.1
  module load alphafold/2.3.2
fi
if [[ -n "$SLURM_TMPDIR" ]]
then
  # Set TMPDIR to fast local storage.
  export TMPDIR="$SLURM_TMPDIR"
fi

### Check values of some environment variables
echo SLURM_JOB_ID="$SLURM_JOB_ID"
echo SLURM_JOB_GPUS="$SLURM_JOB_GPUS"
echo ALPHAFOLD_DIR="$ALPHAFOLD_DIR"
echo ALPHAFOLD_DATADIR="$ALPHAFOLD_DATADIR"

### Check values of variables
echo "output=${output}"
echo "max_template_date=${max_template_date}"


gpu_parameters=("--use_gpu" "--gpu_devices=${SLURM_JOB_GPUS}")
if [ "${step}" == "prepare" ]
then
  gpu_parameters=("--use_gpu=False")
fi

mkdir -p "${output}"
"${ALPHAFOLD}/venv/bin/python" "${ALPHAFOLD_DIR}/singularity/run_singularity.py" \
    --step="$step" \
    "${gpu_parameters[@]}" \
    --fasta_paths="$fasta" \
    --max_template_date="$max_template_date" \
    --data_dir="$ALPHAFOLD_DATADIR" \
    --model_preset=multimer \
    --db_preset=reduced_dbs \
    --output_dir="${output}"
alphafold_return=$?

echo "INFO: AlphaFold return code is ${alphafold_return}"
