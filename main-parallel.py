#!/opt/common/CentOS_6-dev/bin/current/python
descr = "Parallelized signature decomposition"
import argparse
import signature
import sys
import os
import string

parser = argparse.ArgumentParser(description = descr, formatter_class = argparse.RawTextHelpFormatter)
parser.add_argument('-m', '--maffile', help = 'Input MAF', required = True)
parser.add_argument('-s', '--signatures', help = 'Predefined signatures', required = False,
	default = os.path.dirname(os.path.realpath(__file__))+'/Stratton_signatures30.txt')
parser.add_argument('-o', '--outfile', help = 'Name of output file', required = False)
parser.add_argument('-n', '--threads', help = 'Number of threads', required = False, type = int)
args = parser.parse_args()

maf = args.maffile
stratton_file = args.signatures
outfile = args.outfile
threads = args.threads
if outfile is None:
	outfile = os.path.splitext(os.path.basename(maf))[0]+'_sig_decomp.txt'
if threads is None:
	threads = 1

print "Loading known signatures from %s"%stratton_file
stratton = signature.load_stratton_signatures(stratton_file)

print "Making sample signatures from maf %s"%maf
signatures = signature.make(maf, substitution_order=stratton['substitution_order'])
if signatures == None:
    print "Error during signature creation; quitting"
    sys.exit(0)

print "Decomposing signatures and writing to %s"%outfile
signature.decompose_to_file_parallel(signatures['signatures'],
                            					stratton['signatures'],
                            					stratton['names'],
												outfile,
                            					threads)