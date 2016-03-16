#!/usr/bin/env Rscript
library(data.table)

write.maf <- function (...){
  write.table(..., quote = F, col.names = T, row.names = F,
              sep = "\t")
}

bootstrap_indices <- function(
  n_mutations = 10,
  n_draws = 10
){
  sapply(paste(1:n_draws),
         function(i) sample(n_mutations,
                            n_mutations,
                            replace = TRUE),
         USE.NAMES = T,
         simplify = F)
}

generate_bootstrap_maf <- function(maf, n_draws = 1000){
  if(! "Ref_Tri" %in% names(maf)){
    if("TriNuc" %in% names(maf)){
      maf[, Ref_Tri := TriNuc]
    } else {
      stop("must have either Ref_Tri or TriNuc column")
    }
  }

  maf <- maf[Variant_Type == "SNP"]
  maf <- maf[, list(Tumor_Sample_Barcode,
                    Reference_Allele,
                    Tumor_Seq_Allele2,
                    Ref_Tri,
                    Variant_Type)]

  ### for each Tumor_Sample_Barcode, create many (n_draws) bootstrap mafs
  ### and concatenate with the ID replicate
  set.seed(825)
  mafl <- maf[, rbindlist(idcol = "replicate",
                          lapply(
                            bootstrap_indices(n_mutations = .N,
                                              n_draws = n_draws),
                            function(i) .SD[i]
                          )),
              by = Tumor_Sample_Barcode]
  mafl[, Tumor_Sample_Barcode := paste0(Tumor_Sample_Barcode, ":", replicate)]
  mafl[, replicate := NULL]
  mafl
}

if(!interactive()){
  args <- commandArgs(TRUE)
  maf_filename <- args[1]; args <- args[-1]
  n_draws <- args[1]; args <- args[-1]
  output_maf_filename <- args[1]; args <- args[-1]
  split_files <- FALSE
  if(length(args)){
    if(args[1] == "split"){
      split_files <- TRUE
    }
  }


  maf <- suppressWarnings(fread(maf_filename))
  mafl <- generate_bootstrap_maf(maf, n_draws = n_draws)
  original_maf <- maf[, list(Tumor_Sample_Barcode,
                             Reference_Allele,
                             Tumor_Seq_Allele2,
                             Ref_Tri,
                             Variant_Type)]
  mafl <- rbind(original_maf,
                mafl)
  if(split_files){
    dir.create(output_maf_filename, showWarnings = FALSE)
    mafl[, sample := gsub(":.*$", "", Tumor_Sample_Barcode)]
    mafl[, write.maf(.SD, paste0(output_maf_filename, "/", sample, ".txt")), sample]
  } else {
    write.maf(mafl, output_maf_filename)
  }
}
