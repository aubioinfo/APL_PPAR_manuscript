#!/bin/bash
#SBATCH -J htseq
#SBATCH -p small
#SBATCH --mail-type=end
#SBATCH --mail-user=xxxxxx@qq.com
#SBATCH -o /lustre/home/acct-medljm/medljm-jp/Logs/%j.out
#SBATCH -e /lustre/home/acct-medljm/medljm-jp/Logs/%j.err
#SBATCH -n 1
##SBATCH --exclusive

echo "Mission starts at `date`" &&\

# conda activate base
annotation=./gencode.v22.annotation.gtf
inputdir=./star2pass.out
outdir=./htseq.out

sample=$1

htseq-count -f bam -r name -s no -a 10 -t exon -i gene_id -m intersection-nonempty \
${inputdir}/${sample}.Aligned.sortedByCoord.out.bam ${annotation} > ${outdir}/${sample}.count