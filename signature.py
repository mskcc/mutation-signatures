import numpy as np
import matplotlib.pyplot as plt
import string
from scipy.optimize import basinhopping
import pdb

def make(maf_path, out_path=None, substitution_order=None):
    id_col = "Tumor_Sample_Barcode"
    ref_allele_col = "Reference_Allele"
    variant_type_col = "Variant_Type"
    tum_allele_col = "Tumor_Seq_Allele2"
    ref_tri_col = "Ref_Tri"
    required_columns = [id_col, ref_allele_col, variant_type_col, \
                        tum_allele_col, ref_tri_col]

    maf_f = open(maf_path, 'r')
    column_labels = maf_f.readline().strip().split('\t')
    column_index = {}
    for i,x in enumerate(column_labels):
        column_index[x] = i
    for column_label in required_columns:
        if not column_label in column_labels:
            print "Fatal: Required column '%s' missing"%column_label
            return None

    id_index = column_index[id_col]
    ref_allele_index = column_index[ref_allele_col]
    variant_type_index = column_index[variant_type_col]
    tum_allele_index = column_index[tum_allele_col]
    ref_tri_index = column_index[ref_tri_col]

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
    for line in maf_f:
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
    maf_f.close()
    print "%d lines read; %d SNPS counted, %d SNPs skipped, %d non-SNPs skipped"%(line_number, num_snps_counted, num_snps_skipped, num_nonsnps)

    if out_path != None:
        out_f = open(out_path, 'w')
        out_f.write(string.join(substitution_order,'\t'))
        out_f.write('\n')
        for samp_id in signatures:
            f.write(samp_id)
            f.write('\t')
            f.write(string.join(map(str, signatures[samp_id]), '\t'))
            f.write('\n')
        f.close()

    return {'substitution_order':substitution_order, 'signatures':signatures}

def plot(array, ylim=0.2, ax=None):
    color = ((0.196,0.714,0.863),)*16 + ((0.102,0.098,0.098),)*16 + ((0.816,0.180,0.192),)*16 + ((0.777,0.773,0.757),)*16 + ((0.604,0.777,0.408),)*16 + ((0.902,0.765,0.737),)*16
    color = list(color)

    width = max(array.shape)
    x = np.arange(width)
    if ax == None:
        f,ax = plt.subplots(1)
    bars = ax.bar(x, array)
    for h in range(len(x)):
        bars[h].set_color(color[h])
    plt.ylim(0, ylim)
    plt.xlim(0, width)

def alogb(a,b):
    a = max(a, 0.0001)
    b = max(b, 0.0001)
    if a == 0 and b == 0:
        return 0
    elif b == 0:
        return -1*float("-inf")
    else:
        return a*np.log(b)

def div(to_array, from_array):
    to_l = to_array.tolist()
    from_l = from_array.tolist()
    return np.sum(map(alogb, to_l, to_l)) - np.sum(map(alogb, to_l, from_l))

def sym_div(a,b):
    return div(a,b) + div(b,a)

def load_stratton_signatures(path):
    f = open(path, 'r')
    substitution_order = map(string.strip, f.readline().split("\t"))
    sigs = []
    names = []
    for line in f:
        sline = line.split('\t')
        names.append(sline[0])
        sigs.append(map(float, sline[1:]))
    return {'signatures': np.array(sigs).T, \
            'names':names, \
            'substitution_order': substitution_order}

def load(path):
    f = open(path, 'r')
    f.readline() # skip substitution order
    ret = {}
    for line in f:
        sline = line.split('\t')
        ret[sline[0]] = np.array(map(float, sline[1:]))
    return ret


def decompose(target, sigs):
    num_muts = np.sum(target)
    if num_muts < 1:
        print "Warning: sample doesn't have any mutations; cancelling decomposition"
        return None

    num_sigs = sigs.shape[1]
    target = np.array(target)/float(np.sum(target))

    seed = np.random.multinomial(num_muts, np.ones((num_sigs,))/float(num_sigs), 1)[0].astype(float)
    seed = seed/np.linalg.norm(seed)

    class Step(object):
        def __init__(self, stepsize=5):
            self.stepsize = stepsize
        def __call__(self, x):
            x += np.random.uniform(-self.stepsize, self.stepsize, np.shape(x))
            x = x/np.linalg.norm(x)
            return x

    np.seterr(all="raise")
    def error(x):
        coeff = x*x
        approx = sigs.dot(coeff)
        try:
            return np.log(-target.dot(np.log(np.maximum(approx, 0.0001)))) + 100*(np.sum(coeff) - 1)**2 # take that extra log to prioritize the correct sum
        except:
            pdb.set_trace()

    result = basinhopping(error, seed, niter=3, disp=False, T=5, take_step=Step()).x
    return result*result

def decompose_to_file(targets, sigs, sigs_names, to_file):
    to_file = open(to_file, 'w')
    to_file.write(string.join(["Sample Name", "Number of Mutations"] + sigs_names, '\t'))
    to_file.write('\n')
    num_targets = len(targets.keys())
    num_targets_decomposed = 0
    for target_name in targets:
        if num_targets_decomposed%50 == 0:
            print "%d/%d decomposed"%(num_targets_decomposed, num_targets)
        decomposition = decompose(targets[target_name], sigs)
        num_targets_decomposed += 1
        if decomposition is None:
            print "Sample %s not decomposed"%target_name
            continue

        to_file.write(target_name)
        to_file.write('\t')
        to_file.write(str(np.sum(targets[target_name])))
        to_file.write('\t')
        to_file.write(string.join(map(str, decomposition), '\t'))
        to_file.write('\n')
    to_file.close()
