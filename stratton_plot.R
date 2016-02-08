library(data.table)
library(ggplot2)
old <- theme_set(theme_minimal(20))

summarize_maf_trinucleotide <- function(maf){
  adam_trinuc_colnames <-
    c("Sample Name", "Number of Mutations", "ACAA", "ACAC", "ACAG", "ACAT", "CCAA", "CCAC", "CCAG", "CCAT",
      "GCAA", "GCAC", "GCAG", "GCAT", "TCAA", "TCAC", "TCAG", "TCAT",
      "ACGA", "ACGC", "ACGG", "ACGT", "CCGA", "CCGC", "CCGG", "CCGT",
      "GCGA", "GCGC", "GCGG", "GCGT", "TCGA", "TCGC", "TCGG", "TCGT",
      "ACTA", "ACTC", "ACTG", "ACTT", "CCTA", "CCTC", "CCTG", "CCTT",
      "GCTA", "GCTC", "GCTG", "GCTT", "TCTA", "TCTC", "TCTG", "TCTT",
      "ATAA", "ATAC", "ATAG", "ATAT", "CTAA", "CTAC", "CTAG", "CTAT",
      "GTAA", "GTAC", "GTAG", "GTAT", "TTAA", "TTAC", "TTAG", "TTAT",
      "ATCA", "ATCC", "ATCG", "ATCT", "CTCA", "CTCC", "CTCG", "CTCT",
      "GTCA", "GTCC", "GTCG", "GTCT", "TTCA", "TTCC", "TTCG", "TTCT",
      "ATGA", "ATGC", "ATGG", "ATGT", "CTGA", "CTGC", "CTGG", "CTGT",
      "GTGA", "GTGC", "GTGG", "GTGT", "TTGA", "TTGC", "TTGG", "TTGT"
    )
  if(! "TriNuc" %in% names(maf)){
    if("Ref_Tri" %in% names(maf)){
      maf[, TriNuc := Ref_Tri]
    } else {
      stop("must have either Ref_Tri or TriNuc column")
    }
  }

  maf[Variant_Type == "SNP",
      Mut_Tri := paste0(substr(TriNuc, 1, 2),
                        Tumor_Seq_Allele2,
                        substr(TriNuc, 3, 3))]
  trinuc_counts <- dcast.data.table(maf[!is.na(Mut_Tri),
                                        .N,
                                        by=list(Tumor_Sample_Barcode,
                                                Mut_Tri)],
                                    Tumor_Sample_Barcode ~ Mut_Tri,
                                    value.var = "N",
                                    fill=0)
  trinuc_counts_norm_to_one <- with(trinuc_counts,
                                    data.table("Sample Name"=Tumor_Sample_Barcode,
                                               "Number of Mutations" = rowSums(trinuc_counts[,-1,with=FALSE]),
                                               trinuc_counts[,-1,with=FALSE] / rowSums(trinuc_counts[,-1,with=FALSE])))
  setnames(trinuc_counts_norm_to_one, gsub("[^[:alnum:] ]", "", names(trinuc_counts_norm_to_one)))
  new_names <- setdiff(adam_trinuc_colnames, names(trinuc_counts_norm_to_one))
  trinuc_counts_norm_to_one[,new_names := 0, with = F]
  bad_names <- setdiff(names(trinuc_counts_norm_to_one), adam_trinuc_colnames)
  if(length(bad_names > 0)){
    cat("Bad trinucleotide context:", bad_names, "\n")
    trinuc_counts_norm_to_one[, bad_names := NULL, with = F]
  }
  setcolorder(trinuc_counts_norm_to_one, adam_trinuc_colnames)
}


stratton_plot <- function(maf){
  count_summary <- summarize_maf_trinucleotide(maf)
  count_summary_m <- melt.data.table(count_summary,
                                     id.vars = c("Sample Name", "Number of Mutations"),
                                     variable.name = "Mut_Tri",
                                     value.name = "Fraction")
  count_summary_m[, subtype := factor(paste0(substr(Mut_Tri, 2, 2),
                                             " > ",
                                             substr(Mut_Tri, 3, 3)
  ))]
  ggplot(count_summary_m,
         aes(Mut_Tri, Fraction, fill = subtype)) +
    geom_bar(stat="identity") +
    facet_grid(`Sample Name`~.) +
    scale_fill_manual("", values = c("#1EBFF0", "#050708", "#E62725", "#CBCACB", "#A1CF64", "#EDC8C5")) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y=element_text(vjust=1),
          panel.margin = unit(0, "cm")) +
    xlab("") +
    ylab("Fraction of Mutations") +
    theme(legend.text=element_text(size=16, family="Courier"))
  # ggsave(file=paste0(gsub(".txt$", "_stratton_plot.pdf", summ_file)),
  #                    width = 15, height = 3)
}
