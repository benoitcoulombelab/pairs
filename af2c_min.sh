#!/bin/bash

# Exit when any command fails
set -e


targets=$1
structures=$2
output=$3
model=${4:-model_}
# Allow unified memory usage and pooling of memory for multiple GPUs.
export TF_FORCE_UNIFIED_MEMORY=${TF_FORCE_UNIFIED_MEMORY:-1}
export XLA_PYTHON_CLIENT_MEM_FRACTION=${XLA_PYTHON_CLIENT_MEM_FRACTION:-4.0}


echo -e "\n\nRun AF2Complex relaxation on target list ${targets}\n\n"


# Load required modules
if [[ -n "$CC_CLUSTER" ]]
then
  module purge
  module load StdEnv/2020
  module load gcc/9.3.0 openmpi/4.0.3 cuda/11.4 cudnn/8.2.0 kalign/2.03 hmmer/3.2.1 openmm-alphafold/7.5.1 \
              hh-suite/3.3.0 python/3.8
  module load alphafold/2.3.2
  module load af2complex/1.4.1-ca50922
fi
data_dir="$ALPHAFOLD_DATADIR"


# Check values of some variables
echo SLURM_JOB_ID="$SLURM_JOB_ID"
echo SLURM_JOB_GPUS="$SLURM_JOB_GPUS"
echo "structures=${structures}"
echo "output=${output}"
echo "model=${model}"
echo "data_dir=${data_dir}"

if [[ -n "$CC_CLUSTER" ]]
then
  echo "Create AlphaFold virtual environment in ${SLURM_TMPDIR}"
  venv="${SLURM_TMPDIR}/venv"
  virtualenv --no-download "$venv"
  source "${venv}/bin/activate"

  # Install alphafold and its dependencies
  pip install --no-index --upgrade pip
  pip install --no-index --requirement "${AF2COMPLEX}/af2complex-requirements.txt"
  pip install --no-index --requirement "${ALPHAFOLD}/alphafold-requirements.txt"
  pushd "${venv}/bin"
  git apply "${ALPHAFOLD}/alphafold-${ALPHAFOLD_VERSION}.patch"
  popd
  rsync -a --exclude=venv "${AF2COMPLEX}"/* "${SLURM_TMPDIR}"
fi

mkdir -p "${output}"

echo "Start AF2Complex relaxation using run_af2c_min.py"
python -u "${SLURM_TMPDIR}/src/run_af2c_min.py" \
  --target_lst_path="$targets" \
  --input_dir="$structures" \
  --output_dir="$output" \
  --model_str="$model" \
  --use_gpu_relax=True
