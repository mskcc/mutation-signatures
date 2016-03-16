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
conf_int[, c("Tumor_Sample_Barcode", "replicate") := tstrsplit(`Sample Name`, ":", fill = NA), with = F]
conf_int_m <- melt.data.table(conf_int,
                              id.vars = c("Sample Name",
                                          "Number of Mutations",
                                          "Tumor_Sample_Barcode",
                                          "replicate"),
                              value.name = "Fraction",
                              variable.name = "Signature")
conf_int_m[, Fraction := as.numeric(Fraction)]
setkey(conf_int_m, Tumor_Sample_Barcode, Signature)

### for each combination of Tumor_Sample_Barcode and Signature
### calculate the median and range corresponding to 1 s.d. (68.3%)
conf_int_q <- conf_int_m[!is.na(replicate), {
  q <- quantile(Fraction, probs =  c(pnorm(-1), 0.5, pnorm(1)))
  names(q) <- c("lower_val", "median_val", "upper_val")
  q <- c(q, quasi_pval = mean(Fraction < 1e-5))
  lapply(q, function(x) as.numeric(formatC(x, digits = 3)))
},
by = list(Tumor_Sample_Barcode, Signature)]

conf_int_qm <- merge(conf_int_m[is.na(replicate),
                                list(Tumor_Sample_Barcode,
                                     Signature,
                                     observed_val = as.numeric(formatC(x, digits = 3)))],
                     conf_int_q,
                     by = c("Tumor_Sample_Barcode", "Signature"),
                     all = T)
write.maf(conf_int_q, sig_conf_int_filename)
