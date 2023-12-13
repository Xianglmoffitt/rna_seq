#!/bin/bash
#PBS -N salmongc
#PBS -l walltime=48:00:00,mem=30gb
#PBS -l nodes=1:ppn=8
#PBS -o salmongc.out
#PBS -e salmongc.err

cd ~/workspace/cowork/derek/CICPT_1979_RNAseq_July2019
folder=salmongc
if [ ! -d $folder ]; then
    mkdir $folder
fi

fqf=fastq
fastqs=($(ls $fqf/*R1.fastq.gz))
for ((i=0; i<=${#fastqs[@]}-1; i++))
do
    name=$(echo "${fastqs[i]/_R1.fastq.gz/}")
    name=$(echo "${name/$fqf\//}")
    echo "salmon " $fqf/$name\_R1.fastq.gz  $fqf/$name\_R2.fastq.gz
    salmon quant -i ~/genomes/GRCh38/salmon-0.11.3/GRCh38.tx.index \
	-l ISF \
	--gcBias \
	-1 $fqf/$name\_R1.fastq.gz \
	-2 $fqf/$name\_R2.fastq.gz \
	-o $folder/$name \
	-g ~/genomes/GRCh38/gencode.v28.primary_assembly.annotation.gtf \
	-p 8
done
