#!/usr/bin/env Rscript
library(data.table)
args <- commandArgs(TRUE)
sig_filename <- args[1]; args <- args[-1]
sig_conf_int_filename <- args[1]; args <- args[-1]

write.maf <- function (...){
  write.table(..., quote = F, col.names = T, row.names = F, 
              sep = "\t")
}

conf_int <- fread(sig_filename)
setnames(conf_int, names(conf_int), gsub("^Signature.", "", names(conf_int)))
conf_int[, c("Tumor_Sample_Barcode", "replicate") := tstrsplit(`Sample Name`, ":"), with = F]
conf_int_m <- melt.data.table(conf_int, 
                              id.vars = c("Sample Name", 
                                          "Number of Mutations", 
                                          "Tumor_Sample_Barcode",
                                          "replicate"),
                              value.name = "Fraction",
                              variable.name = "Signature")
conf_int_m[, Fraction := as.numeric(Fraction)]

conf_int_q <- conf_int_m[, {
  q <- quantile(Fraction, probs =  c(pnorm(-1), 0.5, pnorm(1)))
  names(q) <- c("lower_val", "median_val", "upper_val")
  q <- c(q, pval = mean(Fraction < 1e-5))
  lapply(q, function(x) as.numeric(formatC(x, digits = 3)))
}, 
by = list(Tumor_Sample_Barcode, Signature)]

write.maf(conf_int_q, sig_conf_int_filename)