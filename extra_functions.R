suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(grid))
suppressMessages(library(cowplot))
suppressMessages(library(stringr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(Cairo))
'%nin%' = Negate('%in%')

####################################################
### ggplot theme
####################################################
theme_bwmin = function(base_size = 12, base_family = 'Helvetica')
{
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(fill = NA, colour = 'black', size = 0.5),
          strip.background = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_line(size = 0.5),
          legend.key = element_rect(color = 'white'),
          legend.key.size = unit(.25, 'cm')
          )
}
####################################################

####################################################
### Parse files
####################################################
parse_maf = function(filepath) {
  maf = fread(filepath) %>%
    tbl_df()
  if ('Ref_Tri' %nin% colnames(maf)) { stop('No Ref_Tri column in MAF.')
  } else { maf }
}

parse_signatures = function(filepath) {
  fread(filepath) %>%
    tbl_df() %>%
    rename(Tumor_Sample_Barcode = `Sample Name`) %>%
    select(-`Number of Mutations`) %>%
    select(Tumor_Sample_Barcode, everything())
}
####################################################

####################################################
### Stratton plot matrix
####################################################
trinuc_context = c(
  "ACAA", "ACAC", "ACAG", "ACAT", "CCAA", "CCAC", "CCAG", "CCAT",
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

nt_comp = list('G'='C', 'A'='T', 'C'='G', 'T'='A')

stratton_plot = function(maf, output_prefix, save_output = TRUE) {

  ### Count fractions of transitions in given trinucleotide context
  trinuc_maf = filter(maf, Variant_Type == 'SNP') %>%
      mutate(TriNuc_Context = ifelse(
          Reference_Allele %in% c('C','T'),
          paste0(str_sub(Ref_Tri,1,2),Tumor_Seq_Allele2,str_sub(Ref_Tri,3,3)),
          paste0(str_sub(Ref_Tri,1,1),nt_comp[Reference_Allele],nt_comp[Tumor_Seq_Allele2],str_sub(Ref_Tri,3,3)))) %>%
      filter(str_sub(TriNuc_Context,2,2) != str_sub(TriNuc_Context,3,3)) %>%
      group_by(Tumor_Sample_Barcode, TriNuc_Context) %>%
      summarize(TriNuc_Count = n()) %>%
      mutate(TriNuc_Frac = TriNuc_Count/sum(TriNuc_Count))

  ### Plot
  trinuc_maf = mutate(trinuc_maf, TriNuc_Context = factor(TriNuc_Context, levels = trinuc_context, ordered = T)) %>%
      mutate(Transition = paste0(str_sub(TriNuc_Context,2,2),'>',str_sub(TriNuc_Context,3,3)))

  outplot = ggplot(trinuc_maf, aes(TriNuc_Context, TriNuc_Frac, fill = Transition)) +
    geom_bar(stat = 'identity') +
    labs(x = '', y = 'Fraction') +
    scale_x_discrete(labels = paste0(str_sub(levels(trinuc_maf$TriNuc_Context),1,2),str_sub(levels(trinuc_maf$TriNuc_Context),4,4))) +
    scale_fill_manual(values = c("#1EBFF0", "#050708", "#E62725", "#CBCACB", "#A1CF64", "#EDC8C5")) +
    scale_y_continuous(labels = scales::percent, breaks = seq(0, max(as.numeric(trinuc_maf$TriNuc_Frac)), .05), expand = c(0,0)) +
    geom_hline(yintercept = seq(0, max(as.numeric(trinuc_maf$TriNuc_Frac)), .05), col = 'darkgrey', size = .5) +
    facet_wrap(~Tumor_Sample_Barcode) +
    theme_bwmin() +
    theme(text = element_text(size = 8), axis.ticks.x = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = .5, hjust = 0),
      panel.background = element_blank(), legend.title = element_blank(),
      strip.text = element_text(size = 10))

  CairoPDF(outplot, file = paste0(output_prefix, '_stratton_plots.pdf'), w = 8, h = 8)

  ### Save output
  # if (save_output == T) saveRDS(trinuc_maf, file = '')
}
####################################################

####################################################
### Composite plot
####################################################

indels = c(
  'Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins','frameshift_deletion',
  'frameshift_insertion', 'nonframeshift_deletion', 'nonframeshift_insertion'
  )

### Annotate known signatures
sign_map = c(
  '1'='Age','2'='APOBEC.1','3'='BRCA1/2','4'='Smoking','6'='MMR.1','7'='UV','9'='IGHV_hypermutation',
  '10'='POLE','11'='TMZ','13'='APOBEC.2','15'='MMR.2','20'='MMR.3','22'='Aristolochic_acid','26'='MMR.4','29'='Tobacco'
  )

### Signature colors // grey for unknowns
grey_pal = colorRampPalette(brewer.pal(9, 'Greys'))(16) # Greys except for hand-picked signatures
sign_cols = c(
  "#e41a1c","#377eb8","#4daf4a","#984ea3",grey_pal[2],"#ff7f00","#ffff33",grey_pal[3],"#a65628",
  "#f781bf","#b2df8a",grey_pal[4],"#377eb8",grey_pal[5],"#ff7f00",grey_pal[6:9],"#ff7f00",
  grey_pal[10],"#cab2d6",grey_pal[11:13],"#ff7f00",grey_pal[14:15],"#984ea3",grey_pal[16]
  )


composite_plot = function(maf, sign, print_genes = NULL, output_prefix, save_output = TRUE) {

  ### Count somatic mutations
  mut_count = group_by(maf, Tumor_Sample_Barcode) %>%
    dplyr::summarize(
      Total_Mutation_Count = n(),
      Indel_Count = sum(Variant_Classification %in% indels),
      Synonymous_SNV_Count = sum(Variant_Classification == 'Silent'),
      Nonsynonymous_SNV_Count = sum(Variant_Classification %nin% c(indels, 'Silent'))
      )

  ### Cluster by signatures
  sample_clust = hclust(dist(sign[,-c(1:2)]))
  sample_order = sign$Tumor_Sample_Barcode[sample_clust$order]

  ### Plot indel counts
  mut_count = melt(mut_count, id.vars = 'Tumor_Sample_Barcode')
  mut_count$Tumor_Sample_Barcode = factor(mut_count$Tumor_Sample_Barcode, levels = sample_order, ordered = T)
  mut_indels = filter(mut_count, variable == 'Indel_Count') %>% droplevels()
  mut_snvs = filter(mut_count, variable %in% c('Synonymous_SNV_Count', 'Nonsynonymous_SNV_Count'))
  barplot_indels = ggplot(mut_indels, aes(x = Tumor_Sample_Barcode, y = value, fill = '#636363')) +
      geom_bar(stat = 'identity') +
      labs(y ='#Indels', x = '') +
      scale_y_continuous(expand = c(0,0)) + geom_hline(yintercept = 0) +
      guides(fill = guide_legend()) +
      scale_fill_manual(values = '#9C824A') +
      theme_bwmin() +
      theme(axis.text.y = element_text(size = 6),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_blank(),
        plot.margin = unit(c(0,.25,0,.25), 'in'),
        legend.position = 'none')

  ### Plot SNV counts
  barplot_snvs = ggplot(mut_snvs, aes(x = Tumor_Sample_Barcode, y = value, fill = variable)) +
      geom_bar(stat = 'identity') +
      scale_fill_manual(values = c('#023474', '#EF0107')) +
      labs(y = '#SNVs', x = '') +
      guides(fill = guide_legend(keywidth = .5, keyheight = .5, title = '')) +
      scale_y_continuous(expand = c(0,0)) +
      geom_hline(yintercept = 0) +
      theme_bwmin() +
      theme(axis.text.y = element_text(size = 6),
        legend.key = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_blank(),
        plot.margin = unit(c(0,.25,0,.25), 'in'),
        legend.position = 'none')

  ### Plot signatures
  sign_melt = melt(sign, id.vars = c('Tumor_Sample_Barcode', 'Patient'), variable.name = 'Signature', value.name = 'Fraction') %>%
    mutate(Signature = str_replace(Signature, 'Signature.', ''))
  sign_melt$Signature = factor(sign_melt$Signature, levels = seq(1,30), ordered = T)
  sign_melt$Signature = revalue(sign_melt$Signature, sign_map)
  sign_melt$Tumor_Sample_Barcode = factor(sign_melt$Tumor_Sample_Barcode, levels = sample_order, ordered = T)

  sample_sig = ggplot(sign_melt, aes(x = Tumor_Sample_Barcode, fill = Fraction)) +
      geom_bar(aes(weight = Fraction, fill = Signature)) +
      guides(fill = guide_legend(keywidth = .5, keyheight = .5, ncol = 2, title = 'Signature')) +
      labs(y = 'Fraction of mutations', x = '') +
      scale_fill_manual(values = sign_cols) +
      scale_y_continuous(expand = c(0,0)) +
      theme_bwmin() +
      theme(legend.key = element_blank(),
        axis.text.y = element_text(size = 6),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(0,.25,0,.25), 'in'),
        legend.position = 'none')

  if (!is.null(print_genes)) {

    ### Subset on selected genes
    mutations = filter(maf, Hugo_Symbol %in% print_genes)

    ### Make oncoprint
    mutations$Tumor_Sample_Barcode = factor(mutations$Tumor_Sample_Barcode, levels = sample_order, ordered = T)
    mutation_heat = ggplot(mutations, aes(x = Tumor_Sample_Barcode, y = Hugo_Symbol)) +
        geom_tile(data = filter(mutations, Mutation_Status == 'Somatic', Variant_Classification %in% c('Frame_Shift_Del', 'Frame_Shift_Ins')), fill = 'darkorange') +
        geom_tile(data = filter(mutations, Mutation_Status == 'Somatic', Variant_Classification %in% c('In_Frame_Del', 'In_Frame_Ins', 'Missense_Mutation', 'Splice_Site')), fill = 'darkgreen') +
        geom_tile(data = filter(mutations, Mutation_Status == 'Somatic', Variant_Classification %in% 'Nonsense_Mutation'), fill = 'black') +
        geom_tile(data = filter(mutations, Mutation_Status == 'Germline'), fill = 'red', col = 'red') +
        labs(y = '', x = '')
        theme_bwmin() +
        theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 8),
          plot.margin = unit(c(0,.25,0,.25), 'in'))

    outplot = plot_grid(barplot_snvs, barplot_indels, sample_sig, mutation_heat, ncol = 1, align = 'v', rel_heights = c(1,1,4,2))
    CairoPDF(outplot, file = paste0(output_prefix, '_composite_plot.pdf'), w = 8, h = 8)

  } else {
    outplot = plot_grid(barplot_snvs, barplot_indels, sample_sig, ncol = 1, align = 'v', rel_heights = c(1,1,4))
    print(outplot)
    CairoPDF(outplot, file = paste0(output_prefix, '_composite_plot.pdf'), w = 8, h = 8)
  }

  if (save_output) print('saving') ### Save output in tabular format
}
####################################################