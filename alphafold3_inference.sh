#!/bin/bash
#SBATCH --account=def-coulomb
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --gpus=a100:1
#SBATCH --mail-type=NONE
#SBATCH --output=alphafold3_inference-%A.out

# Exit when any command fails
set -e

json=$1
output=$2
database=$3
model_dir=${4:-${HOME/models}}

# Check values of some variables
echo SLURM_JOB_ID="$SLURM_JOB_ID"
echo SLURM_JOB_GPUS="$SLURM_JOB_GPUS"
echo "output=${output}"
echo "database=${database}"
echo "model_dir=${model_dir}"
echo
echo


echo -e "\n\nRun AlphaFold 3 inference pipeline on json ${json}\n\n"


if [[ -n "$CC_CLUSTER" ]]
then
  module purge
  module load StdEnv/2023 hmmer/3.4 rdkit/2024.03.5 python/3.12 cuda/12.2 cudnn/9.2

  echo "Create AlphaFold virtual environment in ${SLURM_TMPDIR}"
  venv="${SLURM_TMPDIR}/venv"
  virtualenv --no-download "$venv"
  source "${venv}/bin/activate"

  # Install alphafold and its dependencies
  pip install --no-index --upgrade pip
  pip install --no-index --requirement ~/alphafold3-requirements.txt

  # build data in $VIRTUAL_ENV
  build_data

  # https://github.com/google-deepmind/alphafold3/blob/main/docs/performance.md#compilation-time-workaround-with-xla-flags
  export XLA_FLAGS="--xla_gpu_enable_triton_gemm=false"

  # https://github.com/google-deepmind/alphafold3/blob/main/docs/performance.md#gpu-memory
  # For 3,000 tokens or less on A100 40GB.
  export XLA_PYTHON_CLIENT_PREALLOCATE=true
  export XLA_CLIENT_MEM_FRACTION=0.95
  # For more than 3000 tokens on A100 40GB.
  #export XLA_PYTHON_CLIENT_PREALLOCATE=false
  #export TF_FORCE_UNIFIED_MEMORY=true
  #export XLA_CLIENT_MEM_FRACTION=3.2
fi

mkdir -p "${output}"

echo "Start AlphaFold inference pipeline"
python run_alphafold.py \
    --db_dir="$database" \
    --model_dir="$model_dir" \
    --json_path="$json" \
    --output_dir="$output" \
    --jax_compilation_cache_dir="$HOME/.cache" \
    --norun_data_pipeline  # Run inference stage
