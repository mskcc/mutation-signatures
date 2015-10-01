# mutation-signatures
Create mutation signatures from MAF's, and decompose them into Stratton signatures

Usage: Navigate to the project directory and type the following command
```
python main.py Stratton_signatures.txt [maf file path] [decomposed output file path]
```

This will (1) Create SNP signatures for samples in the MAF (2) Decompose them and write the results to the given path.

In order to use this, you need the following columns in your MAF:  
```"Tumor_Sample_Barcode", "Reference_Allele", "Variant_Type", "Tumor_Seq_Allele2", "Ref_Tri"```
