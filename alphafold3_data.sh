#!/bin/bash
#SBATCH --account=def-coulomb
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --mail-type=NONE
#SBATCH --output=alphafold3_data-%A.out

# Exit when any command fails
set -e

json=$1
output=$2
database=$3
max_template_date=${4:-${current_date}}
threads=${SLURM_CPUS_PER_TASK:-8}

# Check values of some variables
echo SLURM_JOB_ID="$SLURM_JOB_ID"
echo "output=${output}"
echo "database=${database}"
echo "threads=${threads}"
echo
echo


echo -e "\n\nRun AlphaFold 3 data pipeline on json ${json}\n\n"


if [[ -n "$CC_CLUSTER" ]]
then
  module purge
  module load StdEnv/2023 hmmer/3.4 rdkit/2024.03.5 python/3.12

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
fi

mkdir -p "${output}"

echo "Start AlphaFold data pipeline"
python run_alphafold.py \
    --db_dir="$database" \
    --json_path="$json" \
    --output_dir="$output" \
    --jax_compilation_cache_dir="$HOME/.cache" \
    --nhmmer_n_cpu="$threads" \
    --jackhmmer_n_cpu="$threads" \
    --max_template_date="$max_template_date" \
    --norun_inference  # Run data stage
