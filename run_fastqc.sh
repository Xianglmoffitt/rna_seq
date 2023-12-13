#!/bin/bash

#PBS -N fastqc 
#PBS -l walltime=4:00:00
#PBS -l nodes=1:ppn=6,mem=24gb
#PBS -m bea
#PBS -M Xiang.Liu@moffitt.org


# run job parallelly on HPC using GNU parallel tool

export PROJECT_DIR="/share/lab_teng/xiangliu/lixin_rna"
export PARALLEL_WLIST="$PROJECT_DIR/fastqc_wlist"

module load fastqc/0.11.7

cd $PROJECT_DIR

parallel -j 6 -k < $PARALLEL_WLIST
