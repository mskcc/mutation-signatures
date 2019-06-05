#!/usr/bin/env python
from __future__ import print_function
import subprocess
import sys
from string import strip, join


def normalizeTrinuc(trinuc):
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    if not trinuc[1] in ["C", "T"]:
        return join([comp[x] for x in trinuc][::-1], "")
    else:
        return trinuc


def clean_up():
    print("Cleaning up")
    subprocess.call("rm -f ___tmp.maf".split(" "))
    subprocess.call("rm -f ___tmp.bed".split(" "))
    subprocess.call("rm -f ___tmp.tsv".split(" "))


def fetch_regions(ref_fasta):
    """
    Fetch regions from reference genome
    Original command:
    bedtools getfasta -tab -fi /ifs/depot/assemblies/H.sapiens/GRCh37/gr37.fasta -bed ___tmp.bed -fo ___tmp.tsv
    """

    print("Getting regions")
    args = [
        "bedtools",
        "getfasta",
        "-tab",
        "-fi",
        ref_fasta,
        "-bed",
        "___tmp.bed",
        "-fo",
        "___tmp.tsv",
    ]

    subprocess.call(args)


def make_bedfile():
    print("Making bed file")
    grep_ps = subprocess.Popen(
        ["grep", "-Ev", "^#|^Hugo", "___tmp.maf"], stdout=subprocess.PIPE
    )
    cut_ps = subprocess.Popen(
        "cut -f5-7".split(" "), stdin=grep_ps.stdout, stdout=subprocess.PIPE
    )
    awk_ps = subprocess.Popen(
        [
            "awk",
            'BEGIN{OFS="\\t"} {if ($1 == "M") {$1 = "MT"} ; print $1,$2-2,$3+1}',
        ],
        stdin=cut_ps.stdout,
        stdout=subprocess.PIPE,
    )

    # Write out to bed
    with open("___tmp.bed", "w") as bedf:
        bedf.write(awk_ps.communicate()[0])


def main(from_maf, to_maf, ref_fasta):
    """Run everything"""

    # First get rid of non-SNP mutations
    print("Ignoring non-SNP mutations")
    with open(from_maf, "r") as fromf:
        firstline = fromf.readline()
        while firstline.startswith("#"):
            firstline = fromf.readline()
        lines = list(map(strip, firstline.split("\t")))
        indels = []
        type_col = lines[0].index("Variant_Type")

        for line in fromf:
            split_line = list(map(strip, line.split("\t")))
            type = split_line[type_col]
            if type == "SNP":
                lines.append(split_line)
            if type != "SNP":
                indels.append(split_line)

    with open("___tmp.maf", "w") as tmpf:
        tmpf.write(join([join(x, "\t") for x in lines], "\n"))
        tmpf.flush()

    # Make bed file of regions
    make_bedfile()

    # Fetch regions subprocess
    fetch_regions(ref_fasta)

    # Add trinuc to lines
    print("Adding trinucs (normalized to start from C or T)")
    lines[0].append("Ref_Tri")

    with open("___tmp.tsv", "r") as trinucf:
        i = 0
        for line in trinucf:
            i += 1
            trinuc = line.split("\t")[1].strip()
            try:
                trinuc = normalizeTrinuc(trinuc)
            except:
                pass
            lines[i].append(trinuc)

    # What...?
    # Add back indels
    for line in indels:
        line.append("")

    lines.extend(indels)

    print("Writing to %s" % to_maf)
    with open(to_maf, "w") as tof:
        tof.write(join(map(lambda x: join(x, "\t"), lines), "\n"))

    # Clean up tmp files
    clean_up()


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print(
            "Usage: python make_trinuc_maf.py [path to ref_fasta] [source maf path] [target maf path]"
        )
        sys.exit(0)

    ref_fasta = sys.argv[1]
    from_maf = sys.argv[2]
    to_maf = sys.argv[3]

    main(from_maf, to_maf, ref_fasta)
