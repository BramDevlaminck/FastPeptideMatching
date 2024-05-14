#!/bin/bash

#########################################################################################################
### This script is designed to run on the Ghent university HPC                                        ###
###                                                                                                   ###
### how to use:                                                                                       ###
### 1) Swap to the high-memory gallade cluster by executing `module swap cluster/gallade`             ###
### 2) Navigate the to root of the project                                                            ###
### 3) Submit the job to the queue with `qsub suffixarray/build_suffix_array.pbs`                     ###
#########################################################################################################

# go to cluster with high memory
module swap cluster/gallade

# define requested memory, cpu resources and email notifications
#PBS -m abe
#PBS -l walltime=10:00:00
#PBS -l mem=750gb
# ask for 1 node, 1 cpu (not more needed since we don't have parallelism)
#PBS -l nodes=1:ppn=1
#PBS -N suffix_array_construction_uniprot

# define output and error files
#PBS -o stdout.$PBS_JOBID
#PBS -e stderr.$PBS_JOBID

prefix="$VSC_DATA_VO/bram/"

# load Rust
module load Rust/1.75.0-GCCcore-12.3.0
module load Clang/16.0.6-GCCcore-12.3.0 # needed to build the bindings from Rust to C
module load CMake/3.26.3-GCCcore-12.3.0

# go to current working dir and execute
cd $PBS_O_WORKDIR

# compile
cargo build --release

# execute
./target/release/suffixarray_builder -d "$prefix"uniprot_protein_database_minimal.tsv -t "$prefix"taxons.tsv --sparseness-factor 3 --construction-algorithm lib-div-suf-sort -o "$prefix"uniprot_suffix_array_sparse3.bin
