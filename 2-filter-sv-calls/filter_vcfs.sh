#!/bin/bash

# filter manta calls
cut -f2,4 ../sample_keys.tsv > sample2id.tsv
bcftools filter -Oz \
    -i 'INFO/SVTYPE="DEL" && INFO/SVLEN < -500 && INFO/SVLEN > -100000' \
    ../1-generate-sv-calls/manta/results/manta/results/variants/diploidSV.vcf.gz \
    | bcftools reheader -s sample2id.tsv - > manta_deletions_only.vcf.gz
rm sample2id.tsv

# filter lumpy calls
cut -f1,4 ../sample_keys.tsv > srr2id.tsv
bcftools filter -Oz \
    -i 'INFO/SVTYPE="DEL" && INFO/SVLEN < -500 && INFO/SVLEN > -100000' \
    ../1-generate-sv-calls/lumpy/output/pasted.smoove.square.vcf.gz | \
    bcftools reheader -s srr2id.tsv - > lumpy_deletions_only.vcf.gz
rm srr2id.tsv
