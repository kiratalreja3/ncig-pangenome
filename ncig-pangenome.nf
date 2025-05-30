inputs = Channel
    .fromPath("$params.metadata")
    .splitCsv(header: false)
    .map { row ->
        def (sample,assembly) = row
        tuple(sample,assembly)
    }


process map {

    input:
    tuple val (sample), val (assembly) 

    output:
    tuple val (sample), val(assembly), path("${sample}.primary.paf"), path("${sample}.unmapped.txt")

    script:

    """
    module load singularity samtools 

    samtools faidx ${assembly}
    samtools faidx $params.ref

    wfmash -t \${PBS_NCPUS} -m -N -s 50000 -p 90 $params.ref ${assembly} > ${sample}.primary.paf
    comm -23 <(cut -f 1 "${assembly}.fai" | sort) <(cut -f 1 ${sample}.primary.paf | sort) > ${sample}.unmapped.txt
    """

}

process unmapped_remap {

    input:
    tuple val (sample), val(assembly), path("${sample}.unmapped.txt")

    output:
    tuple val (sample), val(assembly), path("${sample}.rescue.paf")

    script:

    """

    samtools faidx ${assembly} -r "${sample}.unmapped.txt" > "${sample}.unmapped.fasta"
    samtools faidx "${sample}.unmapped.fasta"
    wfmash -t \${PBS_NCPUS} -m -s 50000 -p 90 $params.ref ${sample}.unmapped.fasta > ${sample}.unmapped.paf
    bash rescue.sh ${sample}.unmapped.paf > "${sample}.rescue.paf"

    """

}

process combinepaf {

    input:
    val (paf)

    output:
    path ("combinedfiltered.paf")


    script:
    

    """

    if [ -e "combined.paf" ]; then
        rm "combined.paf"
    fi

    for file in ${paf.join(' ')}; do
        cat \$file >> "combined.paf"
    done

    python3 filterpaf.py -i combined.paf -o combinedfiltered.paf

    

    """
}

process community_detection {

    publishDir "$params.outdir" , mode: 'copy', pattern:"*contigs*"

    input:
    val (paf)

    output:
    path ("*.contigs")

    script:
    """
    community_detection.sh ${paf}

    for chr in {1..22} X Y; do
        contigs_file="chr\${chr}.contigs"
        grep -w "chr\${chr}" "${params.ref}.fai" | cut -f1 | sed 's/^>//' >> \$contigs_file
    done

    cat chrX.contigs chrY.contigs > chrXY.contigs


    """

}

process processcontigs {

    publishDir "$params.outdir" , mode: 'copy', pattern:"*fasta*"

    input:
    val (contigs)

    output:
    tuple path ("*.fasta.gz"), path ("*.fai")

    script:

    """

    contigs=$contigs
    export contigs

    outdir=\${PWD}
    export outdir

    metadata=$params.metadataref
    export metadata

    bash processcontig_parallel.sh

    """

}

process pggb {

    publishDir "$params.outdir" , mode: 'copy', pattern:"*pggb.out*"

    input:
    tuple path (fasta), path (fai)

    output:
    path ("*pggb.out*")

    script:

    """

    fasta=${fasta}
    export fasta

    fai=${fai}
    export fai

    outdir=\${PWD}
    export outdir

    bash pggb.sh
    
    """


}


workflow {

    primary_and_unaligned = map(inputs)

    unmappedChannel = primary_and_unaligned.map { sample, assembly, primaryPaf, unmappedTxt ->
    tuple(sample, assembly, unmappedTxt)}
    
    unmapped_samples = unmappedChannel.filter { sample, assembly, unmappedTxt ->
        file(unmappedTxt).size() > 0
    }

    unmapped_remapChannel = unmapped_remap(unmapped_samples)

    rescuePafChannel = unmapped_remapChannel.filter { sample, assembly, rescuePaf ->
        file(rescuePaf).size() > 0
    }

    rescuePafChannel = rescuePafChannel.map { sample, assembly, rescuePaf ->
        rescuePaf
    }

    pafChannel = primary_and_unaligned.map { sample, assembly, primaryPaf, unmappedTxt ->
     primaryPaf}

    combinedPafChannel = rescuePafChannel.mix(pafChannel)
    combinedPafChannel = combinedPafChannel.collect()

    finalPafChannel = combinepaf(combinedPafChannel)
    
    contigfiles = community_detection(finalPafChannel)

    contigs = processcontigs(contigfiles.flatten())

    graph = pggb(contigs)

}