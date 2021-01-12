#!/usr/bin/env nextflow

params.reference = 'ref.fa'
params.scratch = '/local/scratch/esrbhb'

Channel
    .fromSRA(file('../library_list.txt').readLines())
    .set{reads}

reference_file = file(params.reference)

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
    set accession, "${accession}.bam*" into aligned

    """
    bwa mem -R "@RG\\tID:${accession}\\tSM:${accession}\\tPL:ILLUMINA" \
        -t ${task.cpus} ${ref} ${both_ends} | samtools view -bh - | \
        samtools fixmate -m - - | samtools sort - | \
        samtools markdup -r - ${accession}.bam
    samtools index ${accession}.bam
    """
}

aligned.into { aligned_for_smoove_call; aligned_for_smoove_genotype }
faidx.into { faidx_for_smoove_call; faidx_for_smoove_merge;
             faidx_for_smoove_genotype }

process smoove_call {
    container 'brentp/smoove:v0.2.3'
    publishDir 'unmerged'
    cpus 8

    input:
    file 'ref.fa' from reference_file
    file 'ref.fa.fai' from faidx_for_smoove_call
    set accession, file(bam) from aligned_for_smoove_call

    output:
    file "${accession}-smoove.genotyped.vcf.gz" into unmerged

    """
    smoove call --name ${accession} --fasta ref.fa -p ${task.cpus} \
        --genotype ${accession}.bam
    """
}

process smoove_merge {
    container 'brentp/smoove:v0.2.3'

    input:
    file 'ref.fa' from reference_file
    file 'ref.fa.fai' from faidx_for_smoove_merge
    file all_unmerged from unmerged.collect()

    output:
    file "merged.sites.vcf.gz" into merged

    """
    smoove merge --name merged -f ref.fa ${all_unmerged}
    """
}

merged_vcf_faidx_and_bams = merged.combine(faidx_for_smoove_genotype)
    .combine(aligned_for_smoove_genotype)

process smoove_genotype {
    container 'brentp/smoove:v0.2.3'

    input:
    file 'ref.fa' from reference_file
    set 'merged.sites.vcf.gz', 'ref.fa.fai', accession,
        file(bam) from merged_vcf_faidx_and_bams

    output:
    file "${accession}-joint-smoove.genotyped.vcf.gz" into joint_genotyped
    file "${accession}-joint-smoove.genotyped.vcf.gz.csi" into j_g_index

    """
    export TMPDIR=\$PWD
    smoove genotype -d -x --name ${accession}-joint --fasta ref.fa \
        --vcf merged.sites.vcf.gz ${accession}.bam
    echo "done!"
    """
}

process smoove_paste {
    container 'brentp/smoove:v0.2.3'
    publishDir 'output'

    input:
    file all_vcfs from joint_genotyped.collect()
    file all_indexes from j_g_index.collect()

    output:
    file "pasted.smoove.square.vcf.gz" into pasted

    """
    smoove paste --name pasted ${all_vcfs}
    """
}


