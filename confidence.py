import signature
import numpy as np
import sys

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print "Usage: python confidence.py [stratton signature path]"
        sys.exit(0)
    else:
        stratton = signature.load_stratton_signatures(sys.argv[1])	
	generate_decomposition_bank(stratton['signatures'        
        

def generate_decomposition_bank(stratton_signatures):
    def random_signature():
        # First draw random mixture
    raise "Not yet implemented - will randomly generate 10000 signatures "+\
            "and decompose them so that confidence can be quickly queried."

def confidence(num_mutations, decomposition, threshold):
    raise "Not yet implemented - will return a vector of p-values using an "+\
            "existing decomposition bank."
