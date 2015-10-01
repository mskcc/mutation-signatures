# mutation-signatures
Create mutation signatures from MAF's, and decompose them into Stratton signatures

Usage: Navigate to the project directory.

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
