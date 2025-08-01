#!/bin/bash

# load required modules
if [[ -n "$CC_CLUSTER" ]]
then
  module purge
  module load StdEnv/2023
  module load nextflow/25.04.6
fi
if [[ "beluga" == "$CC_CLUSTER" ]]
then
  ulimit -v 40000000
fi

export NXF_OPTS="-Xms500M -Xmx8000M"
