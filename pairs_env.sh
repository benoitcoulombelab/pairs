#!/bin/bash

# load required modules
if [[ -n "$CC_CLUSTER" ]]
then
  module purge
  module load StdEnv/2023
  module load python/3.13.2
fi

python3 -m venv ./pairs-env
source ./pairs-env/bin/activate
pip install git+https://github.com/benoitcoulombelab/pairs.git@main

echo "Activate pairs virtual environment:"
echo "source pairs-env/bin/activate"
