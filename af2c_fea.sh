#!/bin/bash

# Exit when any command fails
set -e

fasta=$1
output=$2
db_preset=${3:-uniprot}
feature_mode=${4:-monomer+organism+fullpdb}
max_template_date=${5:-$(date +"%Y-%m-%d")}

echo -e "\n\nRun AF2Complex feature generation on fasta ${fasta}\n\n"


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
pdb_mmcif_dir="$ALPHAFOLD_PDB_MMCIF"
pdb_seqres_dir="$ALPHAFOLD_PDB_SEQRES"


# Check values of some variables
echo SLURM_JOB_ID="$SLURM_JOB_ID"
echo SLURM_JOB_GPUS="$SLURM_JOB_GPUS"
echo "output=${output}"
echo "db_preset=${db_preset}"
echo "feature_mode=${feature_mode}"
echo "max_template_date=${max_template_date}"
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

echo "Start AF2Complex feature generation using run_af2c_fea.py"
python "${SLURM_TMPDIR}/src/run_af2c_fea.py" \
    --fasta_paths="$fasta" \
    --db_preset="$db_preset" \
    --feature_mode="$feature_mode" \
    --data_dir="$data_dir" \
    --output_dir="$output" \
    --max_template_date="$max_template_date" \
    --small_bfd_database_path="${data_dir}/small_bfd/bfd-first_non_consensus_sequences.fasta" \
    --mgnify_database_path="${data_dir}/mgnify/mgy_clusters_2022_05.fa" \
    --template_mmcif_dir="${pdb_mmcif_dir}/mmcif_files" \
    --obsolete_pdbs_path="${pdb_mmcif_dir}/obsolete.dat" \
    --pdb_seqres_database_path="${pdb_seqres_dir}/pdb_seqres.txt" \
    --uniprot_database_path="${data_dir}/uniprot/uniprot.fasta" \
    --uniref90_database_path="${data_dir}/uniref90/uniref90.fasta" \
    --hhblits_binary_path="${EBROOTHHMINSUITE}/bin/hhblits" \
    --hhsearch_binary_path="${EBROOTHHMINSUITE}/bin/hhsearch" \
    --jackhmmer_binary_path="${EBROOTHMMER}/bin/jackhmmer" \
    --hmmsearch_binary_path="${EBROOTHMMER}/bin/hmmsearch" \
    --hmmbuild_binary_path="${EBROOTHMMER}/bin/hmmbuild" \
    --kalign_binary_path="${EBROOTKALIGN}/bin/kalign" \
    --use_precomputed_msas=True
