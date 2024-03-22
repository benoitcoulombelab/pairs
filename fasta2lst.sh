#!/bin/bash

# Exit when any command fails
set -e

fasta=$1
fasta_base=$(basename "$fasta")
fasta_name="${fasta_base%.*}"
parent_dir="$(dirname -- "$(realpath -- "$fasta")")"
output="${parent_dir}/${fasta_name}.lst"

target=$(awk -F '__' '{print $1"/"$2}' <<< "${fasta_name}")
size=$(grep -v '>' "$fasta" | wc -m)

rm -f "${output}"
echo "${target} ${size} ${fasta_name}" > "${output}"
