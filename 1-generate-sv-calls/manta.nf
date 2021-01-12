#!/usr/bin/env nextflow

params.reference = 'ref.fa'

reference_file = file(params.reference)

Channel
    .fromSRA(file('../library_list.txt').readLines())
    .set{reads}

process bwa_index {
    publishDir 'bwa_index'
    module 'bwa/bwa-0.7.17'

    input:
    file reference from reference_file

    output:
    file "${reference}.*" into reference_index

    """ bwa index ${reference} """
}

process samtools_faidx {
    module 'samtools/samtools-1.9'

    input:
    file reference from reference_file

    output:
    file "${reference}.fai" into faidx

    """ samtools faidx ${reference} """
}

process align {
    cpus 16
    module 'bwa/bwa-0.7.17:samtools/samtools-1.9'
    publishDir 'alignments'

    input:
    file ref from reference_file
    file index from reference_index
    set accession, file(both_ends) from reads

    output:
    file "${accession}.bam" into aligned
    file "${accession}.bam.bai" into aligned_index

    """
    bwa mem -t ${task.cpus} ${ref} ${both_ends} | samtools view -bh - | \
        samtools fixmate -m - - | samtools sort - | \
        samtools markdup -r - ${accession}.bam
    samtools index ${accession}.bam
    """
}

process manta {
    cpus 16
    module 'biocompute/biocompute-modules'
    module 'manta/manta-1.6.0'
    publishDir 'results'

    input:
    file 'ref.fa' from reference_file
    file 'ref.fa.fai' from faidx
    file bams from aligned.collect()
    file bais from aligned_index.collect()

    output:
    file "manta*" into results

    """
    bams=""; for bam in ${bams}; do bams+="--bam \$bam "; done
    configManta.py \$bams --referenceFasta ref.fa --runDir manta
    manta/runWorkflow.py -j ${task.cpus}
    """
}
