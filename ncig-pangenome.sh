#!/bin/bash
#PBS -N ncig-pangenome
#PBS -P te53
#PBS -q normal
#PBS -l storage=gdata/xl04+gdata/if89+gdata/te53
#PBS -l walltime=48:00:00
#PBS -l mem=4GB
#PBS -l ncpus=1
#PBS -l wd
#PBS -j oe

module load nextflow

cd ${workdir}

nextflow run /g/data/te53/ka6418/refgen/workflows/ncig-pangenome/ncig-pangenome.nf -work-dir ${workdir} -config /g/data/te53/ka6418/refgen/workflows/ncig-pangenome/ncig-pangenome.config -profile NCI -resume -with-report "${workdir}/${wfname}_report.html" -with-trace "${workdir}/${wfname}_trace.txt" -with-timeline "${workdir}/${wfname}_timeline.html" -with-dag "${workdir}/${wfname}_graph.dot" 
