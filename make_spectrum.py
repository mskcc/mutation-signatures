import subprocess
import sys
import os
from string import strip, join

if len(sys.argv) < 3:
	print "Usage: python make_trinuc_maf.py [source maf path] [target maf path]"
	sys.exit(0)

from_maf = sys.argv[1]
to_maf = sys.argv[2]
spectrum_path = os.path.splitext(to_maf)[0] + ".spectrum.txt"

fromf = open(from_maf, 'r')
tof = open(to_maf, 'w')
tmpf = open("___tmp.maf",'w')

# First get rid of non-SNP mutations
print "Ignoring non-SNP mutations"
firstline = fromf.readline()
while firstline.startswith("#"):
	firstline = fromf.readline()

lines = [map(strip, firstline.split('\t'))]
indels = []
type_col = lines[0].index("Variant_Type")
	
for line in fromf:
	split_line = map(strip, line.split('\t'))
	type = split_line[type_col]
	if type == "SNP":
		lines.append(split_line)
	if type != "SNP": 
		indels.append(split_line)

tmpf.write(join(map(lambda x: join(x, '\t'), lines), '\n'))
tmpf.flush()

# Make bed file of regions
print "Making bed file"
grep_ps = subprocess.Popen(["grep", "-Ev", "^#|^Hugo", "___tmp.maf"], stdout=subprocess.PIPE)
cut_ps = subprocess.Popen("cut -f5-7".split(" "), stdin=grep_ps.stdout, stdout=subprocess.PIPE)
# awk_ps = subprocess.Popen(["awk", "{OFS=\"\\t\"; print $1,$2-2,$3+1}"], stdin=cut_ps.stdout, stdout=subprocess.PIPE)
awk_ps = subprocess.Popen(["awk", "BEGIN{OFS=\"\\t\"} {if ($1 == \"M\") {$1 = \"MT\"} ; print $1,$2-2,$3+1}"], stdin=cut_ps.stdout, stdout=subprocess.PIPE)
bedf = open("___tmp.bed","w")
bedf.write(awk_ps.communicate()[0])
bedf.close()

# Fetch regions
print "Getting regions"
subprocess.call("bedtools getfasta -tab -fi /ifs/depot/assemblies/H.sapiens/GRCh37/gr37.fasta -bed ___tmp.bed -fo ___tmp.tsv".split(" "))

# Add trinuc to lines
print "Adding trinucs (normalized to start from C or T)"
lines[0].append("Ref_Tri")
trinucf = open("___tmp.tsv","r")
i = 0

def normalizeTrinuc(trinuc):
	comp = {'A':'T',\
		'T':'A',\
		'C':'G',\
		'G':'C'}
	if not trinuc[1] in ['C','T']:
		return join(map(lambda x: comp[x], trinuc)[::-1], '')
	else:
		return trinuc

for line in trinucf:
	i += 1
	trinuc = line.split('\t')[1].strip()
	try:
		trinuc = normalizeTrinuc(trinuc)
	except:
		pass
	lines[i].append(trinuc)

# Add back indels
for line in indels:
	line.append('')

lines.extend(indels)

print "Writing to %s"%to_maf
tof.write(join(map(lambda x: join(x, '\t'), lines), '\n'))
tof.close()

nucleotides = ['A','C','G','T']
nucleotide_complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
transitions = ['CA','CG','CT','TA','TC','TG']
if substitution_order == None:
    substitution_order = []
    for change in transitions:
        for left in nucleotides:
            for right in nucleotides:
                substitution_order.append(left+change+right)

signatures = {}
line_number = 0
num_snps_counted = 0
num_snps_skipped = 0
num_nonsnps = 0
for line in fromf:
    line_number += 1
    sline = map(string.strip, line.split('\t'))
    if len(sline) != len(column_labels):
        print "Warning: Line length doesnt match number of columns." +\
                "skipping line %d"%line_number
        continue
    variant_type = sline[variant_type_index]
    if not variant_type.startswith("SNP"):
        # Skip
        num_nonsnps += 1
        continue
    num_snps_skipped += 1
    samp_id = sline[id_index]
    ref_allele = sline[ref_allele_index]
    tum_allele = sline[tum_allele_index]
    snp = ref_allele + tum_allele
    if not snp in transitions:
        snp = nucleotide_complement[ref_allele] + nucleotide_complement[tum_allele]
    ref_trinuc = sline[ref_tri_index]
    if ref_trinuc == "NA":
        print "Warning: Reference allele not available on "+\
                "line %d; skipping line"%line_number
        continue
    if not ref_trinuc[1] == snp[0]:
        print "Warning: Reference allele does not match reference "+\
                "trinucleotide; skipping line %d"%line_number
        continue
    snp_with_ctx = ref_trinuc[0] + snp + ref_trinuc[2]
    if not samp_id in signatures:
        signatures[samp_id] = [0 for i in substitution_order]
    if snp_with_ctx not in substitution_order:
        print "Warning: substitution on line " + \
                "%d is %s, not "%(line_number,snp_with_ctx) + \
                "found among possible substitutions. Skipping line."
        continue
    signatures[samp_id][substitution_order.index(snp_with_ctx)] += 1
    num_snps_counted += 1
    num_snps_skipped -= 1
fromf.close()
print "%d lines read; %d SNPS counted, %d SNPs skipped, %d non-SNPs skipped"%(line_number, num_snps_counted, num_snps_skipped, num_nonsnps)

out_f = open(spectrum_path, 'w')
out_f.write(string.join(["Tumor_Sample_Barcode"] + substitution_order,'\t'))
out_f.write('\n')
for samp_id in signatures:
    out_f.write(samp_id)
    out_f.write('\t')
    out_f.write(string.join(map(str, signatures[samp_id]), '\t'))
    out_f.write('\n')
out_f.close()


print "Cleaning up"
subprocess.call("rm -f ___tmp.maf".split(" "))
subprocess.call("rm -f ___tmp.bed".split(" "))
subprocess.call("rm -f ___tmp.tsv".split(" "))
