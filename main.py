#!/usr/bin/env python

"""
Decompose signatures using either:

    1. MAF with annotated Ref_Tri column giving nucleotide context (made with make_trinuc_maf)
    2. TSV mutational spectrum, i.e. base transition counts with column names in
    same order as Stratton_signatures.txt. Rows start with sample names and each
    value represents the number of times a certain base transition occurs in a
    sample.

"""

import argparse
from collections import OrderedDict
import signature
import sys


def parse_mutational_spectum_tsv(mutational_spectrum_file):
    """Return signature.make style mutational spectrum from tsv file. Assume
    subsitution order is same as stratton_file."""
    mut_spec = OrderedDict()

    with open(mutational_spectrum_file) as msf:
        for i, line in enumerate(msf):
            # skip header
            if i == 0:
                continue
            split_line = line.split()
            # sample name first column, counts other columns
            mut_spec[split_line[0]] = [int(v) for v in split_line[1:]]

    return mut_spec


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "stratton_file",
        type=str,
        help="Stratton Signatures file (included in repo)",
    )
    parser.add_argument(
        "in_file",
        type=str,
        help="Either MAF file or "
        "mutational spectrum tsv (set --spectrum-tsv "
        "see --help)",
    )
    parser.add_argument("out_file", type=str, help="Output file")
    parser.add_argument("--seed", type=int, default=None, help="Random seed")
    parser.add_argument(
        "--spectrum",
        action="store_true",
        default=False,
        help="input_file is " "mutational spectrum tsv instead of MAF",
    )
    parser.add_argument(
        "--spectrum_output",
        type=str,
        default=None,
        help="write out a " "mutational spectrum tsv",
    )
    args = parser.parse_args()

    print("Loading known signatures from %s" % args.stratton_file)
    stratton = signature.load_stratton_signatures(args.stratton_file)

    if args.spectrum:
        signatures = {}
        print("Parsing mutational spectrum tsv")
        signatures["signatures"] = parse_mutational_spectum_tsv(args.in_file)
    else:
        print("Making sample signatures from maf %s" % args.in_file)
        signatures = signature.make(
            args.in_file,
            substitution_order=stratton["substitution_order"],
            out_path=args.spectrum_output,
        )

        if signatures is None:
            print("Error during signature creation; quitting")
            sys.exit(0)

    print("Decomposing signatures and writing to %s" % args.out_file)
    signature.decompose_to_file(
        signatures["signatures"],
        stratton["signatures"],
        stratton["names"],
        args.out_file,
        random_seed=args.seed,
    )
