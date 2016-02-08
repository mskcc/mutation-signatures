import subprocess
import sys
from string import strip, join

if len(sys.argv) < 3:
	print "Usage: python make_trinuc_maf.py [source maf path] [target maf path]"
	sys.exit(0)

from_maf = sys.argv[1]
to_maf = sys.argv[2]

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
subprocess.call("bedtools getfasta -tab -fi /opt/common/CentOS_6-dev/vep/v81/homo_sapiens/81_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa -bed ___tmp.bed -fo ___tmp.tsv".split(" "))

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

print "Cleaning up"
subprocess.call("rm -f ___tmp.maf".split(" "))
subprocess.call("rm -f ___tmp.bed".split(" "))
subprocess.call("rm -f ___tmp.tsv".split(" "))
