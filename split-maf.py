#!/opt/common/CentOS_6-dev/bin/current/python
descr = "Chunkify MAF"
import os
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description = descr, formatter_class = argparse.RawTextHelpFormatter)
parser.add_argument('-m', '--maffile', help = 'Input MAF', required = True)
parser.add_argument('-o', '--outprefix', help = 'Prefix for output', required = False)
parser.add_argument('-c', '--chunks', help = 'Number of chunks', required = True, type = int)
args = parser.parse_args()

maf = args.maffile
outprefix = args.outprefix
chunks = args.chunks
if outprefix is None:
	outprefix = os.path.splitext(os.path.basename(maf))[0]

### Read MAF
maf = pd.read_table(maf,usecols = ['Tumor_Sample_Barcode', 'Reference_Allele', 'Tumor_Seq_Allele2', 'Variant_Type', 'Ref_Tri'])

### Determine chunk size
samples = maf['Tumor_Sample_Barcode'].unique().tolist()
sample_no = len(samples)
num = float(sample_no)/chunks
l = [samples[i:i + int(num)] for i in range(0, (chunks-1)*int(num), int(num))]
l.append(samples[(chunks-1)*int(num):])

for i in range(len(l)):
	mafChunk = maf[(maf['Tumor_Sample_Barcode'].isin(l[i])) & (maf['Variant_Type'] == 'SNP')]
	outname = outprefix+'_chunk'+str(i)+'.txt'
	mafChunk.to_csv(outname, sep = '\t', index = False)
