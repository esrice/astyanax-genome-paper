#!/bin/bash

# merge lumpy and manta deletions
./merge_deletions.py \
    ../vcf_filtering/lumpy_deletions_only.vcf.gz \
    ../vcf_filtering/manta_deletions_only.vcf.gz | bgzip \
    > merged_deletions.vcf.gz

# count number of deletions per sample
./count_variants_per_sample.py \
    ../2-filter-sv-calls/lumpy_deletions_only.vcf.gz > lumpy_deletion_counts.tsv
./count_variants_per_sample.py \
    ../2-filter-sv-calls/manta_deletions_only.vcf.gz > manta_deletion_counts.tsv
./count_variants_per_sample.py \
    merged_deletions.vcf.gz > merged_deletion_counts.tsv

# annotate the vcf
./annotate_deletions.py merged_deletions.vcf.gz ../AstMex.db \
    | bgzip > merged_annotated.vcf.gz

# filter the vcf to count only deletions of intronic sequence
bcftools filter -i intronic=1 -Oz merged_annotated.vcf.gz \
    > intronic_deletions.vcf.gz
./count_variants_per_sample.py \
    intronic_deletions.vcf.gz > intron_deletion_counts.tsv

# fitler the vcf to count only deletions of regulatory sequence
bcftools filter -i regulatory=1 -Oz merged_annotated.vcf.gz \
    > regulatory_deletions.vcf.gz
./count_variants_per_sample.py \
    regulatory_deletions.vcf.gz > regulatory_deletion_counts.tsv

# filter the vcf to count only deletions of coding sequence
bcftools filter -i coding=1 -Oz merged_annotated.vcf.gz \
    > coding_deletions.vcf.gz
./count_variants_per_sample.py \
    coding_deletions.vcf.gz > coding_deletion_counts.tsv
