import signature
import sys

if len(sys.argv) < 4:
    print "Usage: python main.py [stratton signature path] [maf file path] [decomposed output file path]"
    sys.exit(0)

stratton_file = sys.argv[1]
maf_file = sys.argv[2]
out_file = sys.argv[3]

print "Loading known signatures from %s"%stratton_file
stratton = signature.load_stratton_signatures(stratton_file)

print "Making sample signatures from maf %s"%maf_file
signatures = signature.make(maf_file, substitution_order=stratton['substitution_order'])
if signatures == None:
    print "Error during signature creation; quitting"
    sys.exit(0)

print "Decomposing signatures and writing to %s"%out_file
signature.decompose_to_file(signatures['signatures'], \
                            stratton['signatures'], \
                            stratton['names'], \
                            out_file)
