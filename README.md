# mutation-signatures
Create mutation signatures from MAF's, and decompose them into Stratton signatures

### Required packages ###

You will need the python packages `numpy` and `scipy` for anything related to decomposition, and if you want to plot (using `plot.py`) you will also need `matplotlib`.

### Usage ###

Navigate to the project directory.

The SNPs in your MAF need to be annotated with trinucleotide contexts in a column called ```Ref_Tri```. If this is not already the case, you can use the following command:  
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

In order to calculate the confidence intervals on mutational signatures, resample a maf file (with replacement) 1000 times thusly:

    ./sigsig.R input.maf 1000 input.resamp.maf

Run decomposition as usual:

    python main.py Stratton_signatures30.txt input.resamp.maf input.resamp.sig.txt
    
Then calculate the confidence intervals and a quasi-pvalue for each signature:

    ./sigsig_conf_int.R input.resamp.sig.txt


### References ###
Alexandrov, L. B., Nik-Zainal, S., Wedge, D. C., Campbell, P. J., & Stratton, M. R. (2013). Deciphering Signatures of Mutational Processes Operative in Human Cancer. Cell Reports, 3(1), 246–259. doi:10.1016/j.celrep.2012.12.008
