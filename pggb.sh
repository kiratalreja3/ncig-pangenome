#!/bin/bash
#PBS -N pggb
#PBS -P te53
#PBS -q express
#PBS -l walltime=1:00:00
#PBS -l jobfs=400GB
#PBS -l mem=192GB
#PBS -l ncpus=48
#PBS -l wd
#PBS -l storage=gdata/if89+gdata/xl04+gdata/te53
#PBS -j oe


set -ex

module use -a /g/data/if89/shpcroot/modules
module load singularity
module load pangenome/pggb
unset -f which

cp ${fasta} ${fai} ${PBS_JOBFS}

fastabase=$(basename "${PBS_JOBFS}/${fasta}" .fasta.gz)
haps=$(cut -f1 "${PBS_JOBFS}/${fai}" | sed 's/\(#.*#\).*/\1/' | sort | uniq | wc -l)

pggb -i ${PBS_JOBFS}/${fasta} -o "${PBS_JOBFS}/${fastabase}.pggb.out" -n ${haps} -p 98 -s 100000 -k 331 -O 0.03 -m -A -S -V chm13,grch38 -t ${PBS_NCPUS} -T ${PBS_NCPUS}

mv ${PBS_JOBFS}/${fastabase}.pggb.out ${outdir}/
