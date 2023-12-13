#!/bin/bash

#PBS -N salmon_gencode_filter 
#PBS -l walltime=8:00:00
#PBS -l nodes=1:ppn=6,mem=48gb
#PBS -m bea
#PBS -M Xiang.Liu@moffitt.org


# run job parallelly on HPC using GNU parallel tool

export PROJECT_DIR="/share/lab_teng/xiangliu/lixin_rna"
export PARALLEL_WLIST="$PROJECT_DIR/salmon_wlist_2"

#module load star/2.6.0c
module load samtools/1.9

cd $PROJECT_DIR

parallel -j 3 -k < $PARALLEL_WLIST
