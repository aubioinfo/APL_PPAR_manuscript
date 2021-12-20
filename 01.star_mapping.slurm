#!/bin/bash
#SBATCH -J STAR_mapping
#SBATCH -p small
#SBATCH --mail-type=end
#SBATCH --mail-user=xxxxxxxx@qq.com
#SBATCH -o /lustre/home/acct-medljm/medljm-jp/Logs/%j.out
#SBATCH -e /lustre/home/acct-medljm/medljm-jp/Logs/%j.err
#SBATCH -n 24
##SBATCH --exclusive

echo "Mission starts at `date`" &&\

### 1. ALIGNMENT with STAR-2PASS ###
####################################

inputdir=./rawdata
outdir=./star2pass.out
genomeindex=./index_star_2.7.8a_GRCh38_gencode.v22
annotation=./gencode.v22.annotation.gtf

thread=24
readFileComm=cat
R1=1.fastq.gz
R2=2.fastq.gz

cat ${1} |while read sample
do
STAR --genomeDir ${genomeindex} --readFilesCommand ${readFileComm} \
--readFilesIn ${inputdir}/${sample}_${R1} ${inputdir}/${sample}_${R2} \
--runThreadN ${thread} --alignIntronMax 500000 \
--outFilterMultimapScoreRange 1 \
--outFilterMultimapNmax 20 \
--outFilterMismatchNmax 10 \
--alignMatesGapMax 1000000 \
--sjdbScore 2 --alignSJDBoverhangMin 1 \
--genomeLoad NoSharedMemory --limitBAMsortRAM 700000000000 \
--outFilterMatchNminOverLread 0.33 \
--outFilterScoreMinOverLread 0.33 --sjdbOverhang 100 \
--outSAMstrandField intronMotif --outSAMattributes NH HI NM MD AS XS \
--sjdbGTFfile ${annotation} \
--limitSjdbInsertNsj 2000000 \
--outSAMmultNmax 1 \
--outSAMunmapped None \
--outFileNamePrefix ${outdir}/${sample}. \
--outBAMsortingThreadN ${thread} \
--outSAMtype BAM SortedByCoordinate \
--twopassMode Basic 
done