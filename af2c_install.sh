#!/bin/bash

# Exit when any command fails
set -e

rm -rf af2complex
git clone https://github.com/FreshAirTonight/af2complex.git
cd af2complex
git checkout 4c558034f0f1f784ac5ddd1b7234a3bc10df54b4

patch=.af2complex_path
{
  echo "diff --git a/src/alphafold/common/confidence.py b/src/alphafold/common/confidence.py"
  echo "index 68dca87..b2b96ec 100644"
  echo "--- a/src/alphafold/common/confidence.py"
  echo "+++ b/src/alphafold/common/confidence.py"
  echo "@@ -23,6 +23,7 @@"
  echo " from typing import Dict, Optional, Tuple, Union, List"
  echo " import numpy as np"
  echo " import scipy.special"
  echo "+import scipy.spatial"
  echo " import networkx as nx"
  echo ""
  echo ""
} >> "$patch"
git apply "$patch"
rm "$patch"
