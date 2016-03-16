# mutation-signatures
Create mutation signatures from MAF's, and decompose them into Stratton signatures

### Required packages ###

You will need the python packages `numpy` and `scipy` for anything related to decomposition, and if you want to plot (using `plot.py`) you will also need `matplotlib`.

### Usage ###

Navigate to the project directory.

The SNPs in your MAF need to be sorted ([1-22, X, Y, M]) and annotated with trinucleotide contexts in a column called ```Ref_Tri```. If this is not already the case, you can use the following command:
```
python make_trinuc_maf.py <source-maf-path> <target-maf-path>
```

Next, use the following command, which will (1) Create SNP signatures for samples in the MAF (not saved to disk) (2) Decompose them and write the results to the given path.
```
python main.py Stratton_signatures.txt <maf-file-path> <decomposed-output-file-path>
```

In order to use this, you need the following columns in your MAF:  
```"Tumor_Sample_Barcode", "Reference_Allele", "Variant_Type", "Tumor_Seq_Allele2", "Ref_Tri"```

### Confidence Intervals ###

1) Resample each Tumor_Sample_Barcode in a maf file (with replacement) 1000 times, thusly:

    ./sigsig.R input.maf 1000 input.resamp.maf

Alternatively you can opt to split the output maf file by Tumor_Sample_Barcode with a 4th argument `split`:

    ./sigsig.R input.maf 1000 input_resamp_maf/ split

Note that the original maf is included at the top of the resampled file. Resamples get `:1`, `:2`, etc. at the end of the Tumor_Sample Barcode.

2) Run decomposition as usual:

    python main.py Stratton_signatures30.txt input.resamp.maf input.resamp.sig.txt
    
3) Calculate the (1 s.d.) confidence intervals and a quasi-pvalue for each signature:

    ./sigsig_conf_int.R input.resamp.sig.txt input.resamp.sig.conf_int.txt

Decompsed signatures without `:1`, `:2`, etc. at the end will appear as `observed_val`.

Example output:

Tumor_Sample_Barcode | Signature | observed_val | lower_val | median_val | upper_val | quasi_pvalue
--- | --- | --- | --- | --- | --- | ---
TCGA-05-4249 | 1 | 0.10231 | 0.013611 | 0.08426 | 0.13959 | 0.12098
TCGA-05-4249 | 2 | 2.953e-08 | 1.015e-11 | 2.709e-07 | 0.018698 | 0.50322
TCGA-05-4249 | 3 | 1.912e-10 | 1.131e-11 | 3.392e-10 | 5.492e-07 | 0.84942
TCGA-05-4249 | 4 | 0.38532 | 0.27567 | 0.37571 | 0.46896 | 0
TCGA-05-4249 | 5 | 3.432e-9 | 9.335e-12 | 3.725e-10 | 7.171e-09 | 0.95238

### References ###
Alexandrov, L. B., Nik-Zainal, S., Wedge, D. C., Campbell, P. J., & Stratton, M. R. (2013). Deciphering Signatures of Mutational Processes Operative in Human Cancer. Cell Reports, 3(1), 246–259. doi:10.1016/j.celrep.2012.12.008
