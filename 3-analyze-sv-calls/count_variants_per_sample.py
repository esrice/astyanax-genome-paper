#!/usr/bin/env python3
"""
Count the number of called variants per sample in a VCF file.
"""

import argparse
import collections

import vcf


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "vcf", help="the vcf file to analyze", type=lambda f: vcf.Reader(filename=f)
    )
    return parser.parse_args()


def main():
    args = parse_args()

    call_counts = collections.Counter()
    hom_alt_counts = collections.Counter()
    het_counts = collections.Counter()
    for record in filter(lambda r: not r.is_filtered, args.vcf):
        for call in filter(lambda s: not s.is_filtered, record.samples):
            call_counts[call.sample] += 1
            if call.is_variant:
                if call.is_het:
                    het_counts[call.sample] += 1
                else:
                    hom_alt_counts[call.sample] += 1

    print("\t".join(["sample", "call_count", "hom_alt_count", "het_count"]))
    for sample in call_counts.keys():
        print(
            "\t".join(
                map(
                    str,
                    [
                        sample,
                        call_counts[sample],
                        hom_alt_counts[sample],
                        het_counts[sample],
                    ],
                )
            )
        )


if __name__ == "__main__":
    main()
