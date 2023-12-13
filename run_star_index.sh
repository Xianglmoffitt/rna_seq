#!/bin/bash

#PBS -N star index 
#PBS -l walltime=4:00:00
#PBS -l nodes=1:ppn=6,mem=24gb
#PBS -m bea
#PBS -M Xiang.Liu@moffitt.org


# run job parallelly on HPC using GNU parallel tool

export PROJECT_DIR="/share/lab_teng/xiangliu/lixin_rna"
export PARALLEL_WLIST="$PROJECT_DIR/fastqc_wlist"

module load star/2.6.0c

cd $PROJECT_DIR

# make STAR index
STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir /share/lab_teng/xiangliu/hg38_star_index \
--genomeFastaFiles /share/lab_teng/xiangliu/ref_genome/GRCh38.p13.genome.fa \
--sjdbGTFfile /share/lab_teng/xiangliu/share_data/lixin_braf_rnaseq/gencode_filtered.gtf \
--sjdbOverhang 149
