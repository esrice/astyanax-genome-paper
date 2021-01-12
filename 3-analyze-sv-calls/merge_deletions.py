#!/usr/bin/env python3
"""
Given two vcf files containing deletions, output a new vcf file
containing only the deletions that appear in both inputs.
"""

import argparse
from itertools import starmap
import sys

import vcf


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-r",
        "--reciprocal-overlap",
        type=float,
        default=0.5,
        help="minimum reciprocal overlap to consider a "
        "deletion in both files to be the same deletion [0.5]",
    )
    parser.add_argument(
        "vcf1",
        help="the first vcf file to merge",
        type=lambda f: vcf.Reader(filename=f),
    )
    parser.add_argument(
        "vcf2",
        help="the second vcf file to merge",
        type=lambda f: vcf.Reader(filename=f),
    )
    return parser.parse_args()


def compare_records(record1, record2, contig_index):
    """
    Compares two vcf._Record instances by location to determine which
    would come earlier in a vcf file.

    Args:
        record1, record2 (vcf._Record): two vcf records to compare
        contig_index (dict): mapping CHROM names to their index in the
            VCF header

    Returns:
        compare_value (int): >0 if record1 > record2, <0 if record1 <
            record2; 0 if record1 and record2 overlap
    """
    if record1.CHROM != record2.CHROM:
        return contig_index[record1.CHROM] - contig_index[record2.CHROM]
    else:
        if (record1.POS >= record2.POS and record1.POS <= record2.sv_end) or (
            record1.sv_end >= record2.POS and record1.sv_end <= record2.sv_end
        ):
            return 0
        else:
            return record1.POS - record2.POS


def reciprocal_overlap(record1, record2, min_reciprocal_overlap=0.5):
    """
    Given two vcf._Record instances representing deletions, determines
    whether there is an overlap between them based on a minimum
    reciprocal overlap. If there is no overlap, returns None; if there
    is an overlap, returns the start and end points of that overlap.

    Args:
        record1, record2 (vcf._Record): the records to be intersected
        min_reciprocal_overlap(float): the minimum reciprocal overlap
            between the two records in order for them to be considered
            overlapping. Reciprocal overlap is defined as
                size_of_overlap / max(size_of_r1, size_of_r2)

    Returns:
        None if there is not an overlap between the two records, or
            there is an overlap but it is smaller than
            min_reciprocal_overlap
        (overlap_start, overlap_end) if there is an overlap of
            sufficient size between the two records
    """
    if record1.CHROM != record2.CHROM:
        return 0

    r1_deletion_size = record1.sv_end - record1.POS
    r2_deletion_size = record2.sv_end - record2.POS

    overlap_start = max(record1.POS, record2.POS)
    overlap_end = min(record1.sv_end, record2.sv_end)
    overlap_size = overlap_end - overlap_start
    reciprocal_overlap = min(
        overlap_size / r1_deletion_size, overlap_size / r2_deletion_size
    )

    return reciprocal_overlap


def call_type_key(call):
    if call.gt_type is None:
        return -1
    else:
        return call.gt_type


CallData = vcf.model.make_calldata_tuple(["GT"])


def merge_calls(call1, call2):
    """
    Given two calls on the same sample, outputs a new call with the
    lesser of the two genotypes; i.e., 0/0 < 0/1 < 1/1. Also simplifies
    the call data, leaving it with only a GT field.

    N.B. This function is not commutative, as the output call will have
    site = call1.site.
    """
    return vcf.model._Call(
        call1.site,
        call1.sample,
        CallData(min(call1, call2, key=call_type_key).data.GT),
    )


def make_sample_index_key(sample_index):
    def sample_index_key(call):
        return sample_index[call.sample]

    return sample_index_key


def merge_records(record1, record2, sample_index_key):
    """
    Takes two VCF records containing deletions and merges them into a
    single record. Rather than trying to split hairs over where the
    boundaries of the merged record should be when the called boundaries
    are imprecise anyway, this function just uses the metadata
    associated with the bigger deletion.

    Args:
        record1, record2 (vcf.model._Record): deletion records to
            merge. This function does not check to make sure it's
            reasonable to merge the two; it just takes the consensus
            for each of the individual sample calls.
        sample_index_key (function): function that maps _Call objects
            to the index of their sample. `make_sample_index_key()` can
            create one of these.
    """
    # figure out which record is bigger, so that we can use it as the
    # template for the new record
    small_record, big_record = sorted(
        [record1, record2], key=lambda r: r.sv_end - r.POS
    )

    # sort the calls by sample index, and then merge each pair of calls
    # for the same sample into a single call
    big_record.samples = list(
        starmap(
            merge_calls,
            zip(
                sorted(big_record.samples, key=sample_index_key),
                sorted(small_record.samples, key=sample_index_key),
            ),
        )
    )
    big_record.FORMAT = "GT"

    return big_record


def merge_all_deletions(reader1, reader2, min_reciprocal_overlap=0.5):
    # god help us if the VCFs have headers in different orders
    contig_index = {k: i for i, k in enumerate(reversed(reader1.contigs.keys()))}

    # keep the samples in the order in which they appear in the reader1
    # header
    sample_index = {k: i for i, k in enumerate(reader1.samples)}
    sample_index_key = make_sample_index_key(sample_index)

    try:
        record1, record2 = next(reader1), next(reader2)
        while True:  # TODO please don't do this
            while compare_records(record1, record2, contig_index) < 0:
                record1 = next(reader1)
            while compare_records(record1, record2, contig_index) > 0:
                record2 = next(reader2)
            overlap = reciprocal_overlap(record1, record2)
            if overlap >= min_reciprocal_overlap:
                yield merge_records(record1, record2, sample_index_key)
            record1 = next(reader1)
            record2 = next(reader2)
    except StopIteration:
        return


def main():
    args = parse_args()

    writer = vcf.Writer(sys.stdout, args.vcf1)
    # TODO fix output header
    for record in merge_all_deletions(args.vcf1, args.vcf2, args.reciprocal_overlap):
        writer.write_record(record)
        writer.flush()
    writer.close()


if __name__ == "__main__":
    main()
