
## Title: This Script defines the functions for the omics analyses
## Date: 2023
## Author: Linda Dieckmann


library(RColorBrewer)
library(colorspace)
library(plotrix)
library(stringr)
library(ggpubr)
library(ggplot2)
library(tidyverse) 
library(scales)
library(Polychrome)
library(extrafont)

# color-blind friendly palette 
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

orange_nice <- "#CC6677"
yellow_nice <- "#DDCC77"
green_nice <- "#44AA99"

three_colors_colorblind <- c(orange_nice, yellow_nice, green_nice)

brown_red_cb <- "#661100"

# scales::show_col(cbPalette)
# scales::show_col(safe_colorblind_palette)

# distinct colors
n <- 60
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',] 
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
cols_30 <- c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", 
             "#1B9E77", "#D95F02", "#E6AB02", "#B15928", "#FBB4AE", 
             "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", 
             "#E5D8BD", "#FDDAEC", "#F2F2F2", "#B3E2CD", "#FDCDAC", 
             "#CBD5E8", "#F4CAE4", "#E6F5C9", "#D9D9D9", "#BC80BD", 
             "#FB8072", "#FFED6F", "#999999", "#A6761D", "#8DD3C7")


set.seed(18)
colors_35 <- createPalette(35,  c("#010101", "#ff0000"), M=1000, range = c(30, 85))
colors_35 <- lighten(colors_35, 0.4)

greysPalette <- c("#BDBDBD", "#969696", "#737373", "#525252")
colors_placenta_cvs_greys <- c("gray49", "gray84")

# define some basic colors for cvs and placenta 
color_cvs <- "skyblue"
color_placenta <- "steelblue4"
colors_placenta_cvs_colored <- c(color_placenta, color_cvs)

colors_cvs_3levels <- c(lighten(color_cvs, 0.4), lighten(color_cvs, 0.1), darken(color_cvs, 0.2))
colors_placenta_3levels <- c(lighten(color_placenta, 0.4), lighten(color_placenta, 0.1), darken(color_placenta, 0.2))

# colors for qtls
# rna + genotype -> orange + green -> eQTLs -> brown
eqtl_color <- "#8B4513"
eqtl_color_2levels <- c(lighten(eqtl_color, 0.4), lighten(eqtl_color, 0.2))

# methylation + rna -> blue + orange -> cyan -> eQTMs
eqtm_color <- 	"darkcyan"
eqtm_color_2levels <- c(lighten(eqtm_color, 0.4), lighten(eqtm_color, 0.2))

# methylation + geno -> blue + green -> brown/khaki -> choose gold
meqtl_color <- 	"gold4"
meqtl_color_2levels <- c(lighten(meqtl_color, 0.4), lighten(meqtl_color, 0.2))

# need colors for other plot good to distinguish:
qtl_color <- c("coral1", "cyan3", "gold2")

# define some figure parameters
standard_textsize <- 8
standard_dpi_fig <- 600
standard_width_mm_single <- 85
standard_width_mm_double <- 170
standard_height_mm <- 190

# written for Preparation Scripts -----------------------------------------

lmp <- function(modelobject) {
  # extracts the overall ANOVA p-value out of a linear model object
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1], f[2], f[3], lower.tail = F)
  attributes(p) <- NULL
  return(p)
}

meth_BtoM <- function(x) {
  # transforms methylation Beta to M values
  # input: methylation Beta matrix (CpGs in rows and samples in columns)
  m <- log2(x/(1-x))
  return(m)
}

lm_p_data_confounders <- function(y, preds, dataframe) {
  # takes as input a list of predictors and a vector of dependent variables, runs lm model and extracts p-values with lmp function
  p_data_confounders <- data.frame(matrix(NA, nrow = length(y), ncol = length(preds)))
  names(p_data_confounders) <- preds
  for (i in preds) {
    p_model_PC_confounders <- unlist(lapply(y, function(x) {
      lmp(lm(
        formula = paste0("`", x, "` ~", i),
        data = dataframe, na.action = na.omit
      ))
    }))
    p_data_confounders[, i] <- p_model_PC_confounders
  }
  return(p_data_confounders)
}

Get_Annotation_RNA <- function(gene_names_in_rna_vec, ensembl){
  # takes as input a vector of ensembl gene names and a biomaRt ensembl object and gives out genes with available location information
  # print number of genes in the beginning
  print(paste0("start length of genes: ", length(gene_names_in_rna_vec)))
  # get position of genes using biomaRt
  annotation_genes_rna <- getBM(
    attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "hgnc_symbol"), #"transcription_start_site"
    filters = "ensembl_gene_id",
    values = gene_names_in_rna_vec,
    mart = ensembl) %>% unique()
  
  # print genes not found and number of rows
  na_genes <- setdiff(gene_names_in_rna_vec,annotation_genes_rna$ensembl_gene_id)
  print(paste0("missing: ",na_genes))
  print(paste0("number of genes remaining: ",nrow(annotation_genes_rna)))
  
  # filter for chromosomes = include only 1-22
  chromosome_list <- as.character(1:22)
  annotation_genes_rna <- annotation_genes_rna[annotation_genes_rna$chromosome_name %in% chromosome_list, ]
  print(paste0("number of genes remaining with chromosomes 1-22: ",nrow(annotation_genes_rna)))
  
  return(annotation_genes_rna)
}

pvalue_sig_plot <- function(datamatrix, dataset) {
  # plots color-coded p values from matrix of p-values (melted to long format)
  plotcolors <- brewer.pal(4, "Greys") # "Blues"
  ggplot(datamatrix, aes(as.factor(Var2), as.factor(Var1))) +
    geom_tile(aes(fill = value), colour = "transparent") +
    theme_bw() +
    scale_fill_manual(values = c(plotcolors[4], plotcolors[3], plotcolors[2], plotcolors[1]), name = "significance level", guide = guide_legend(reverse = TRUE)) +
    theme(
      axis.line = element_line(colour = "grey"),
      panel.background = element_blank(),
      text = element_text(family = "Helvetica", size = standard_textsize),
      axis.title = element_text(family = "Helvetica", size = standard_textsize),
      axis.text = element_text(family = "Helvetica", size = standard_textsize),
      legend.text = element_text(family = "Helvetica", size = standard_textsize),
      legend.position="none"
    ) +
    labs(x = "Principal Component", y = dataset)
}


# outlier identification for boxplots
is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}

## Basic Plots
basic_gg_point_plot <- function(data, x, y, title) {
  ggplot(data, aes_string(x, y)) +
  geom_point()+
  theme_bw() +
  theme(
    axis.line = element_line(colour = "grey"),
    panel.background = element_blank(),
    text = element_text(family = "Helvetica", size = standard_textsize),
    axis.title = element_text(family = "Helvetica", size = standard_textsize),
    axis.text = element_text(family = "Helvetica", size = standard_textsize),
    legend.text = element_text(family = "Helvetica", size = standard_textsize),
  ) +
  labs(x = x, y = y, title = title)
}

basic_gg_box_plot <- function(data, x, y, title) {
  ggplot(data, aes_string(x, y)) +
  geom_boxplot()+
  theme_bw() +
  theme(
    axis.line = element_line(colour = "grey"),
    panel.background = element_blank(),
    text = element_text(family = "Helvetica", size = standard_textsize),
    axis.title = element_text(family = "Helvetica", size = standard_textsize),
    axis.text = element_text(family = "Helvetica", size = standard_textsize),
    legend.text = element_text(family = "Helvetica", size = standard_textsize),
  ) +
  # geom_text(aes(label=ID), na.rm=TRUE, size=2, nudge_x=0.15) + # to plot with ID name, not for public (data security)
  labs(x = "", y = "", title = title)
}

## function for desaturating colors
desatcolor <- function(cols, sat=0.5) {
  X <- diag(c(1, sat, 1)) %*% rgb2hsv(col2rgb(cols))
  hsv(X[1,], X[2,], X[3,])
}


# written for Characteristics ---------------------------------------------


embold_vector <- function(src, boulder) {
  if (!is.factor(src)) src <- factor(src)                   
  src_levels <- levels(src)                                
  brave <- boulder %in% src_levels                         
  if (all(brave)) {                                       
    b_pos <- purrr::map_int(boulder, ~which(.==src_levels))
    b_vec <- rep("plain", length(src_levels))              
    b_vec[b_pos] <- "bold"                                
    b_vec                                                   
  } else {
    stop("All elements of 'boulder' must be in src")
  }
}

# to display QTL results on different levels
barplot_counts_function <- function(hits_eqtl_dataframe, tissue_name, bold_level, colorpalette, ylimits) {
  # takes as input a prepared data frame with numbers of hits and a character saying for which tissue and which level should be bold
  # plots the counts as barplot
  hits_eqtl_dataframe %>% ggplot(aes(y = cnt, x = method, fill = as.factor(as.numeric(method)))) + 
    geom_col() + 
    geom_text(aes(label = comma(cnt, accuracy = 1L), y = cnt), stat = "identity", vjust=-0.2, size = 5*0.36) + 
    labs(title = tissue_name, y = "Count", x = " ") + facet_wrap(~number, ncol = 3) + 
    theme(legend.position = "none", 
          legend.title = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), 
          plot.title = element_text(family = "Helvetica", size = standard_textsize *1.3), axis.title.y = element_blank(), axis.title.x = element_blank(),
          axis.text.x= element_text(family = "Helvetica", size=standard_textsize, angle=45, hjust=1, face=embold_vector(hits_eqtl_dataframe$method, c(bold_level))), 
          axis.text.y = element_blank(), axis.ticks.y=element_blank(), strip.text = element_text(family = "Helvetica", size=standard_textsize),
          plot.margin = margin(t = 0.7, r = 0, b = 0, l = 1.8, unit = "cm")) +
    scale_y_continuous(labels = scientific, expand = expansion(mult = c(0, 0.3)), limits=ylimits) +
    scale_fill_manual(values = colorpalette) 
}

volcano_plot_matrix <- function(eqtl_statistic_dataframe, plotname) {
  # takes as input a data frame with eqtl statistics and a character saying for which tissue and plots effect vs. p-value
  ggplot(eqtl_statistic_dataframe, aes(x = beta, y = -log10(pvalue),  colour = FDR < 0.05)) + 
    geom_point(alpha = 1.5, size = 1.2) + 
    geom_vline(xintercept = 0, colour = "#990000", linetype = "dashed") + 
    theme_bw()+
    theme(plot.title = element_text(family = "Helvetica", size = standard_textsize), 
          axis.title = element_text(family = "Helvetica", size = standard_textsize),
          legend.text = element_text(family = "Helvetica", size = standard_textsize)) + 
    scale_colour_manual(name = 'FDR < 0.05', values = setNames(c('darkgreen','darkred'),c(T, F))) +
    labs(x = "beta", y="-log10 nominal p-value", title = plotname)
} 

volcano_plot_boot <- function(eqtl_statistic_dataframe, plotname) {
  # takes as input a data frame with eqtl statistics and a character saying for which tissue and plots effect vs. p-value
  ggplot(eqtl_statistic_dataframe, aes(x = corrected_beta, y = -log10(nominal_pval))) + 
    geom_point(alpha = 1.5, size = 1.2) + 
    geom_vline(xintercept = 0, colour = "#990000", linetype = "dashed") + 
    theme_bw()+
    theme(plot.title = element_text(family = "Helvetica", size = standard_textsize), 
          axis.title = element_text(family = "Helvetica", size = standard_textsize),
          legend.text = element_text(family = "Helvetica", size = standard_textsize)) + 
    labs(x = "corrected beta", y="-log10 nominal p-value", title = plotname)
} 

# Histogram
hist_plot <- function(data, columnname, xtitle, ytitle, title) {
  ggplot(data, aes_string(x=columnname)) + 
    geom_histogram(fill="darkgrey", binwidth = 20000) + 
    #geom_density() +
    theme_classic() +
    theme(
      plot.title = element_text(family = "Helvetica", size = standard_textsize), 
      axis.title = element_text(family = "Helvetica", size = standard_textsize), 
      axis.text.x = element_text(family = "Helvetica",size = standard_textsize, angle = 45, hjust=1),
      axis.text.y = element_text(family = "Helvetica", size = standard_textsize),
      legend.text = element_text(family = "Helvetica", size = standard_textsize)) + 
    labs(x = xtitle, y=ytitle, title = title)
} 

### change values / groupings of gene regions
group_as_chipseeker <- function(vector_regions) {
  vector_regions[grep("exon 1 of", vector_regions)]   <- "1st Exon"
  vector_regions[grep("Exon \\(", vector_regions)]    <- "Other Exon"
  vector_regions[grep("intron 1 of", vector_regions)] <- "1st Intron"
  vector_regions[grep("Intron \\(", vector_regions)]  <- "Other Intron"
  vector_regions[grep("Downstream", vector_regions)]  <- "Downstream (<=300)"
  vector_regions[grep("^Distal", vector_regions)]     <- "Distal Intergenic"
  return(vector_regions)
} 
group_exons_introns <- function(vector_regions) {
  vector_regions[grep("1st Exon", vector_regions)]   <- "Exon"
  vector_regions[grep("Other Exon", vector_regions)]    <- "Exon"
  vector_regions[grep("1st Intron", vector_regions)] <- "Intron"
  vector_regions[grep("Other Intron", vector_regions)]  <- "Intron"
  return(vector_regions)
} 
group_transcribed_regions <- function(vector_regions) {
  vector_regions[grep("5' UTR", vector_regions)]   <- "Transcribed"
  vector_regions[grep("Exon", vector_regions)]    <- "Transcribed"
  vector_regions[grep("Intron", vector_regions)] <- "Transcribed"
  vector_regions[grep("3' UTR", vector_regions)]  <- "Transcribed"
  vector_regions[grep("Promoter", vector_regions)]  <- "Promoter"
  vector_regions[grep("Downstream (<=300)", vector_regions)]  <- "Intergenic"
  vector_regions[grep("Distal Intergenic", vector_regions)]  <- "Intergenic"
  return(vector_regions)
} 

group_body <- function(vector_regions) {
  vector_regions[grep("5' UTR", vector_regions)]   <- "UTR"
  vector_regions[grep("Exon", vector_regions)]    <- "Gene Body"
  vector_regions[grep("Intron", vector_regions)] <- "Gene Body"
  vector_regions[grep("3' UTR", vector_regions)]  <- "UTR"
  vector_regions[grep("Promoter", vector_regions)]  <- "Promoter"
  vector_regions[grep("Downstream (<=300)", vector_regions)]  <- "Intergenic"
  vector_regions[grep("Distal Intergenic", vector_regions)]  <- "Intergenic"
  return(vector_regions)
} 

# Barplot to show percentage of SNPs in regions; for ChipSeeker; transcribed region is colored
barplot_snps_regions <- function(data, group, yvalues, xvalues, plottitle) {
  ggplot(data, aes(fill={{group}}, y={{yvalues}}, x={{xvalues}})) + 
    geom_rect(fill = "lavender",xmin = 1.5,xmax = 5.5, ymin = -Inf,ymax = Inf, alpha = 0.02, color = "grey") +
    geom_rect(fill = "lavenderblush",xmin = 0,xmax = 1.5, ymin = -Inf,ymax = Inf, alpha = 0.02, color = "grey") +
    geom_bar(position="dodge", stat="identity") + 
    scale_fill_grey(start=0.8, end=0.4) +
    scale_y_continuous(breaks=seq(0,1, by = 0.1), labels = scales::percent_format(accuracy = 1)) +
    theme_classic() +
    theme(
      plot.title = element_text(family = "Helvetica", size = standard_textsize), 
      axis.title = element_text(family = "Helvetica", size = standard_textsize), 
      axis.text.x = element_text(family = "Helvetica", angle = 45, hjust=1, size = standard_textsize), 
      axis.text.y = element_text(family = "Helvetica", size = standard_textsize),
      legend.text = element_text(family = "Helvetica", size = standard_textsize),
      legend.title=element_blank()) + 
    labs(x = "", y="", title = plottitle)
} 

# Barplot to show percentage of CpGs in regions for Methylation annotation
barplot_cpgs_regions <- function(data, group, yvalues, xvalues, plottitle) {
  ggplot(data, aes(fill={{group}}, y={{yvalues}}, x={{xvalues}})) + 
    geom_bar(position="dodge", stat="identity") + 
    scale_fill_grey(start=0.8, end=0.4) +
    scale_y_continuous(breaks=seq(0,1, by = 0.1), labels = scales::percent_format(accuracy = 1)) +
    theme_classic() +
    theme(
      plot.title = element_text(family = "Helvetica", size = standard_textsize), 
      axis.title = element_text(family = "Helvetica", size = standard_textsize), 
      axis.text.x = element_text(family = "Helvetica", angle = 45, hjust=1, size = standard_textsize), 
      axis.text.y = element_text(family = "Helvetica", standard_textsize),
      legend.text = element_text(family = "Helvetica", size = standard_textsize),
      legend.title=element_blank()) + 
    labs(x = "", y="", title = plottitle)
} 

# test enrichment with permuatation
  # this function tests enrichment of 'feature' among vector 'own', considering 'background'
  # own and background should be distinct sets
enrichment_permed <- function(own, background, feature, nperm){
  nsample <- length(own)
  
  overlap      <- length(own[which(own==feature)])  
  non_overlap  <- length(own) - overlap
  
  overlap_bkgr      <- length(background[which(background==feature)])  
  non_overlap_bkgr  <- length(background) - overlap_bkgr
  
  conf_mtrx        <- matrix(c(overlap, overlap_bkgr, non_overlap, non_overlap_bkgr), 2, 2, byrow = TRUE)
  fisher_test_rslt <- fisher.test(conf_mtrx)
  p_value <- fisher_test_rslt$p.value
  or      <- fisher_test_rslt$estimate
  
  resampling <- lapply(1:nperm, function(x){
    sample_gr          <- sample(background, nsample)     
    sample_overlap     <- length(sample_gr[which(sample_gr==feature)])
    sample_non_overlap <- nsample - sample_overlap
    conf_mtrx          <- matrix(c(overlap, sample_overlap, non_overlap, sample_non_overlap), 2, 2, byrow = TRUE)
    fisher_test_rslt   <- fisher.test(conf_mtrx)
    c(p_value = fisher_test_rslt$p.value, fisher_test_rslt$estimate, sample_overlap = sample_overlap)
  }
  )
  
  p_value_emp <- (sum(unlist(resampling)[attr(unlist(resampling), "names") == "sample_overlap"] >= overlap)+1) / (nperm+1)
  p_value_permutat <- unlist(resampling)[attr(unlist(resampling), "names") == "p_value"] %>% mean() 
  or_permutat <- unlist(resampling)[attr(unlist(resampling), "names") == "odds ratio"] %>% mean() 
  
  return(c(number_feature = overlap, or, p_val = p_value, or_perm_mean = or_permutat, p_val_perm_mean = p_value_permutat, p_val_emp = p_value_emp)) 
}

# test enrichment without permuatation
# this function tests enrichment of 'feature' among vector 'own', considering 'background'
# own and background should be distinct sets
enrichment_basic <- function(own, background, feature){
  nsample <- length(own)
  
  overlap      <- length(own[which(own==feature)])  
  non_overlap  <- length(own) - overlap
  
  overlap_bkgr      <- length(background[which(background==feature)])  
  non_overlap_bkgr  <- length(background) - overlap_bkgr
  
  conf_mtrx        <- matrix(c(overlap, overlap_bkgr, non_overlap, non_overlap_bkgr), 2, 2, byrow = TRUE)
  fisher_test_rslt <- fisher.test(conf_mtrx)
  p_value <- fisher_test_rslt$p.value
  or      <- fisher_test_rslt$estimate
  
  return(c(number_feature = overlap, or, p_val = p_value)) 
}


# written for Scripts Placenta-CVS ----------------------------------------

# pie chart showing the % of same/different direction of effects in placenta and cvs
piechart_signs <- function(data, percent, name, plottitle, colors) {
  ggplot(data, aes(x="", y=percent, fill=name)) +
    geom_bar(stat="identity", width=1, color="white") + 
    coord_polar("y", start=0)+
    geom_text(aes(x = 1,label = paste0(round(percent), "%")), position = position_stack(vjust = 0.5), size = 4)+
    labs(x = NULL, y = NULL, fill = NULL, title = plottitle) +
    scale_fill_manual(values = colors) +
    theme_classic()+
    theme(
      axis.title.x = element_blank(),
      axis.text.x=element_blank(),
      axis.title.y = element_blank(),
      panel.border = element_blank(),
      axis.line = element_blank(),
      panel.grid=element_blank(),
      axis.ticks = element_blank(),
      plot.title=element_text(family = "Helvetica", size = standard_textsize, face="bold"),
      legend.text=element_text(family = "Helvetica", size=standard_textsize))
} 

# boxplot for betas
betasboxplot <- function(meltdata, color, xlabels, ylabtitle) {
  ggboxplot(meltdata, x = "variable", y = "value", fill = color) +
    scale_x_discrete(labels = xlabels) + 
    scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
    #geom_jitter(shape=16, position=position_jitter(0.2)) +
    xlab("") + ylab(ylabtitle) +
    theme_classic()+
    theme(
      axis.title.x = element_blank(),
      axis.text.x= element_text(family = "Helvetica", size=standard_textsize),
      axis.title.y = element_text(family = "Helvetica", size=standard_textsize),
      plot.title=element_text(family = "Helvetica", size=standard_textsize, face="bold"))
} 

betasboxplot_signed <- function(meltdata, color, xlabels, ylabtitle) {
  ggplot(meltdata, aes(x = variable, y = value, fill = sign_effect)) + 
    geom_boxplot()+
    scale_fill_manual(values = color) +
    scale_x_discrete(labels = xlabels) + 
    scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
    #geom_jitter(shape=16, position=position_jitter(0.2)) +
    xlab("") + ylab(ylabtitle) +
    theme_classic()+
    theme(
      axis.title.x = element_blank(),
      axis.text.x= element_text(family = "Helvetica", size=standard_textsize),
      axis.title.y = element_text(family = "Helvetica", size=standard_textsize),
      plot.title=element_text(family = "Helvetica", size=standard_textsize, face="bold"),
      legend.text = element_text(family = "Helvetica", size = standard_textsize),
      legend.title=element_blank())
} 

scaled_betas_bar <- function(data, color, plottitle) {
  ggplot(data, aes(x=variable, y=beta_mean)) +
    geom_bar(stat="identity", aes(fill = sign_effect), width = 0.3) +
    geom_errorbar(aes(ymin=beta_mean-se_beta_mean, ymax=beta_mean+se_beta_mean), width=.1) +
    scale_fill_manual(values = color) +
    xlab("") + ylab("mean +/- mean standard error") +
    labs(title = plottitle) +
    scale_y_continuous(breaks = seq(-0.4, 0.4, 0.2), limits = c(-0.6, 0.6)) +
    scale_x_discrete(labels = c(expression(paste("standardized betas ", bold("Placenta"))),expression(paste("standardized betas ", bold("CVS"))))) + 
    theme_classic()+
    theme(
      axis.title.x = element_blank(),
      axis.text.x= element_text(family = "Helvetica", size=standard_textsize),
      axis.title.y = element_text(family = "Helvetica", size=standard_textsize),
      plot.title=element_text(family = "Helvetica", size=standard_textsize, face="bold"),
      legend.text = element_text(family = "Helvetica", size = standard_textsize),
      legend.title=element_blank())
} 

snp_vs_expression_plot <- function(gene, snp, genesymbol, snpdataset_placenta, genedataset_placenta, snpdataset_cvs, genedataset_cvs) {
  # takes as input a snp and a gene and plots expression in cvs vs. placenta
  # note that data frames containing the prepared info about eqtls need to be loaded, including geno_info file
  top_snp1_data_placenta <- snpdataset_placenta[rownames(snpdataset_placenta) == snp, ]
  top_snp1_data_placenta_info <- geno_info[rownames(geno_info) == snp, ]
  top_snp1_data_placenta_alleles <- factor(str_replace_all(top_snp1_data_placenta, c("0" = paste(top_snp1_data_placenta_info$allele.1, top_snp1_data_placenta_info$allele.1, sep=""), "2" = paste(top_snp1_data_placenta_info$allele.2, top_snp1_data_placenta_info$allele.2, sep=""), "1" = paste(top_snp1_data_placenta_info$allele.1,  top_snp1_data_placenta_info$allele.2, sep=""))), levels = c(paste(top_snp1_data_placenta_info$allele.1, top_snp1_data_placenta_info$allele.1, sep=""), paste(top_snp1_data_placenta_info$allele.1,  top_snp1_data_placenta_info$allele.2, sep=""), paste(top_snp1_data_placenta_info$allele.2, top_snp1_data_placenta_info$allele.2, sep="")))
  top_gene1_data_placenta <- t(genedataset_placenta[rownames(genedataset_placenta) == gene, ])
  
  top_snp1_data_cvs = snpdataset_cvs[rownames(snpdataset_cvs) == snp, ]
  top_snp1_data_cvs_info <- geno_info[rownames(geno_info) == snp, ]
  top_snp1_data_cvs_alleles <- factor(str_replace_all(top_snp1_data_cvs, c("0" = paste(top_snp1_data_cvs_info$allele.1, top_snp1_data_cvs_info$allele.1, sep=""), "2" = paste(top_snp1_data_cvs_info$allele.2, top_snp1_data_cvs_info$allele.2, sep=""), "1" = paste(top_snp1_data_cvs_info$allele.1,  top_snp1_data_cvs_info$allele.2, sep=""))), levels = c(paste(top_snp1_data_cvs_info$allele.1, top_snp1_data_cvs_info$allele.1, sep=""), paste(top_snp1_data_cvs_info$allele.1,  top_snp1_data_cvs_info$allele.2, sep=""), paste(top_snp1_data_cvs_info$allele.2, top_snp1_data_cvs_info$allele.2, sep="")))
  top_gene1_data_cvs = t(genedataset_cvs[rownames(genedataset_cvs) == gene, ])
  
  plot_data_top1_placenta = data.frame(top_snp1_data_placenta_alleles, top_gene1_data_placenta)
  colnames(plot_data_top1_placenta) = c("snp", "gene_expr")
  plot_data_top1_placenta$tissue <- rep("Placenta", nrow(plot_data_top1_placenta))
  
  plot_data_top1_cvs = data.frame(top_snp1_data_cvs_alleles, top_gene1_data_cvs)
  colnames(plot_data_top1_cvs) = c("snp", "gene_expr")
  plot_data_top1_cvs$tissue <- rep("CVS", nrow(plot_data_top1_cvs))
  
  plot_data_top1 <- rbind(plot_data_top1_cvs, plot_data_top1_placenta)
  plot_data_top1$tissue <- as.factor(plot_data_top1$tissue)
  plot_data_top1$tissue <- relevel(plot_data_top1$tissue, "Placenta")
  
  filter(plot_data_top1, !is.na(snp)) %>% 
    ggplot(aes(x = snp, y = gene_expr, fill = tissue)) + #fill = tissue fill = eqtl_color_2levels[2]
    scale_y_continuous(limits = c(-4, 4), breaks = seq(-4, 4, by = 2)) +
    geom_boxplot(fill = orange_nice) +
    facet_grid(~tissue) +
    theme_classic()+
    theme(plot.title = element_text(family = "Helvetica", size = standard_textsize), 
          axis.title = element_text(family = "Helvetica", size = standard_textsize), 
          axis.text = element_text(family = "Helvetica", size = standard_textsize), 
          strip.text.x = element_text(family = "Helvetica", size = standard_textsize), 
          legend.position="none") + 
    labs(x = "Genotype", y="Gene Expression in TMM values \n(corrected for covariates)", title = paste(genesymbol, snp, sep="~")) 
}

snp_vs_methylation_plot <- function(cpg, snp, snpdataset_placenta, cpgdataset_placenta, snpdataset_cvs, cpgdataset_cvs) {
  # takes as input a snp and a cpg and plots methylation in cvs vs. placenta
  # note that data frames containing the prepared info about meqtls need to be loaded, including geno_info file
  top_snp1_data_placenta <- snpdataset_placenta[rownames(snpdataset_placenta) == snp, ]
  top_snp1_data_placenta_info <- geno_info[rownames(geno_info) == snp, ]
  top_snp1_data_placenta_alleles <- factor(str_replace_all(top_snp1_data_placenta, c("0" = paste(top_snp1_data_placenta_info$allele.1, top_snp1_data_placenta_info$allele.1, sep=""), "2" = paste(top_snp1_data_placenta_info$allele.2, top_snp1_data_placenta_info$allele.2, sep=""), "1" = paste(top_snp1_data_placenta_info$allele.1,  top_snp1_data_placenta_info$allele.2, sep=""))), levels = c(paste(top_snp1_data_placenta_info$allele.1, top_snp1_data_placenta_info$allele.1, sep=""), paste(top_snp1_data_placenta_info$allele.1,  top_snp1_data_placenta_info$allele.2, sep=""), paste(top_snp1_data_placenta_info$allele.2, top_snp1_data_placenta_info$allele.2, sep="")))
  top_cpg1_data_placenta <- t(cpgdataset_placenta[rownames(cpgdataset_placenta) == cpg, ])
  
  top_snp1_data_cvs = snpdataset_cvs[rownames(snpdataset_cvs) == snp, ]
  top_snp1_data_cvs_info <- geno_info[rownames(geno_info) == snp, ]
  top_snp1_data_cvs_alleles <- factor(str_replace_all(top_snp1_data_cvs, c("0" = paste(top_snp1_data_cvs_info$allele.1, top_snp1_data_cvs_info$allele.1, sep=""), "2" = paste(top_snp1_data_cvs_info$allele.2, top_snp1_data_cvs_info$allele.2, sep=""), "1" = paste(top_snp1_data_cvs_info$allele.1,  top_snp1_data_cvs_info$allele.2, sep=""))), levels = c(paste(top_snp1_data_cvs_info$allele.1, top_snp1_data_cvs_info$allele.1, sep=""), paste(top_snp1_data_cvs_info$allele.1,  top_snp1_data_cvs_info$allele.2, sep=""), paste(top_snp1_data_cvs_info$allele.2, top_snp1_data_cvs_info$allele.2, sep="")))
  top_cpg1_data_cvs = t(cpgdataset_cvs[rownames(cpgdataset_cvs) == cpg, ])
  
  plot_data_top1_placenta = data.frame(top_snp1_data_placenta_alleles, top_cpg1_data_placenta)
  colnames(plot_data_top1_placenta) = c("snp", "cpg_meth")
  plot_data_top1_placenta$tissue <- rep("Placenta", nrow(plot_data_top1_placenta))
  
  plot_data_top1_cvs = data.frame(top_snp1_data_cvs_alleles, top_cpg1_data_cvs)
  colnames(plot_data_top1_cvs) = c("snp", "cpg_meth")
  plot_data_top1_cvs$tissue <- rep("CVS", nrow(plot_data_top1_cvs))
  
  plot_data_top1 <- rbind(plot_data_top1_cvs, plot_data_top1_placenta)
  plot_data_top1$tissue <- as.factor(plot_data_top1$tissue)
  plot_data_top1$tissue <- relevel(plot_data_top1$tissue, "Placenta")
  
  filter(plot_data_top1, !is.na(snp)) %>% 
    ggplot(aes(x = snp, y = cpg_meth)) +
    geom_boxplot(fill = orange_nice) +
    facet_wrap(~tissue) +
    scale_y_continuous(limits = c(-4, 8), breaks = seq(-4, 8, by = 2)) +
    theme_classic()+
    theme(plot.title = element_text(family = "Helvetica", size = standard_textsize), 
          axis.title = element_text(family = "Helvetica", size = standard_textsize), 
          axis.text = element_text(family = "Helvetica", size = standard_textsize), 
          strip.text.x = element_text(family = "Helvetica", size = standard_textsize), 
          legend.position="none") + 
    labs(x = "Genotype", y="DNAm in M values \n(corrected for covariates)", title = paste(cpg, snp, sep="~")) 
}

cpg_vs_expression_plot <- function(gene, cpg, genesymbol, cpgdataset_placenta, genedataset_placenta, cpgdataset_cvs, genedataset_cvs) {
  # takes as input a cpg and a gene and plots expression in cvs vs. placenta
  # note that data frames containing the prepared info about eqtms need to be loaded
  top_cpg1_data_placenta <- cpgdataset_placenta[rownames(cpgdataset_placenta) == cpg, ]
  top_gene1_data_placenta <- t(genedataset_placenta[rownames(genedataset_placenta) == gene, ])
  
  top_cpg1_data_cvs = cpgdataset_cvs[rownames(cpgdataset_cvs) == cpg, ]
  top_gene1_data_cvs = t(genedataset_cvs[rownames(genedataset_cvs) == gene, ])
  
  plot_data_top1_placenta = data.frame(top_cpg1_data_placenta, top_gene1_data_placenta)
  colnames(plot_data_top1_placenta) = c("cpg", "gene_expr")
  plot_data_top1_placenta$tissue <- rep("Placenta", nrow(plot_data_top1_placenta))
  
  plot_data_top1_cvs = data.frame(top_cpg1_data_cvs, top_gene1_data_cvs)
  colnames(plot_data_top1_cvs) = c("cpg", "gene_expr")
  plot_data_top1_cvs$tissue <- rep("CVS", nrow(plot_data_top1_cvs))
  
  plot_data_top1 <- rbind(plot_data_top1_cvs, plot_data_top1_placenta)
  plot_data_top1$tissue <- as.factor(plot_data_top1$tissue)
  plot_data_top1$tissue <- relevel(plot_data_top1$tissue, "Placenta")

    plot_data_top1 %>%
    ggplot(aes(x = cpg, y = gene_expr, color= orange_nice)) +
    geom_point() +
    facet_wrap(~tissue) +
    scale_color_manual(values = orange_nice) +
    scale_y_continuous(limits = c(-4, 4), breaks = seq(-4, 4, by = 2)) +
    scale_x_continuous(limits = c(-3, 3), breaks = seq(-3, 3, by = 1)) +
    theme_classic()+
    theme(plot.title = element_text(family = "Helvetica", size = standard_textsize), 
          axis.title = element_text(family = "Helvetica", size = standard_textsize), 
          axis.text = element_text(family = "Helvetica", size = standard_textsize), 
          strip.text.x = element_text(family = "Helvetica", size = standard_textsize), 
          legend.position="none") + 
    labs(x = "DNAm", y="Gene Expression in TMM values \n(corrected for covariates)", title = paste(genesymbol, cpg, sep="~")) 
}

scatterplot_signs <- function(data, plottitle) {
  ggplot(data, aes(x = beta.cvs, y = beta.placenta, color = type_sig, shape = type_clumped)) +
    geom_point(size = 0.8, position=position_jitter(h=0.05,w=0.05), alpha = 0.5) +
    scale_color_manual(values = three_colors_colorblind) +
    scale_shape_manual(values = c(8, 7, 0, 4)) +
    xlab("(nominal) beta CVS") + ylab("(nominal) beta Placenta") +
    labs(title = plottitle, color = "significant in", shape = "in clumped results in") +
    theme_classic()+
    theme(
      axis.title.x = element_text(family = "Helvetica", size=standard_textsize),
      axis.text.x= element_text(family = "Helvetica", size=standard_textsize),
      axis.title.y = element_text(family = "Helvetica", size=standard_textsize),
      plot.title=element_text(family = "Helvetica", size=standard_textsize, face="bold"),
      legend.title=element_text(family = "Helvetica", size=standard_textsize, face = "italic"),
      legend.text = element_text(family = "Helvetica", size = standard_textsize),
      legend.position = "bottom") +
    geom_abline(intercept = 0, slope = 1, linetype=2, color="darkgrey", linewidth=0.5)
} 

scatterplot_signs_eqtms <- function(data, plottitle) {
  ggplot(data, aes(x = beta.cvs, y = beta.placenta, color = type_sig)) +
    geom_point(size = 0.8, position=position_jitter(h=0.05,w=0.05), alpha = 0.5) +
    scale_color_manual(values = three_colors_colorblind) +
    xlab("(nominal) beta CVS") + ylab("(nominal) beta Placenta") +
    labs(title = plottitle, color = "significant in", shape = "in clumped results in") +
    theme_classic()+
    theme(
      axis.title.x = element_text(family = "Helvetica", size=standard_textsize),
      axis.text.x= element_text(family = "Helvetica", size=standard_textsize),
      axis.title.y = element_text(family = "Helvetica", size=standard_textsize),
      plot.title=element_text(family = "Helvetica", size=standard_textsize, face="bold"),
      legend.title=element_text(family = "Helvetica", size=standard_textsize, face = "italic"),
      legend.text = element_text(family = "Helvetica", size = standard_textsize),
      legend.position = "bottom") +
    geom_abline(intercept = 0, slope = 1, linetype=2, color="darkgrey", linewidth=0.5)
} 

reformat_long_plot <- function(data) {
  df <- data %>%
  melt(value.name = "absolute.beta", variable.name="tissue") %>%
  mutate(x = ifelse(tissue == "CVS", 1, 2)) %>%
  mutate(xj = jitter(x, amount=.03))
  return(df)
} 

boxplot_paired_data <- function(data, textlabel, plottitle) {
  ggplot(data=data, aes(y=absolute.beta)) +
    geom_line(aes(x=xj, group=id_qtl), color = "grey", linewidth = 0.3) +
    geom_boxplot(aes(x=x, group=tissue), width=0.2, outlier.shape = NA, fill = orange_nice) +
    geom_point(aes(x=xj), color = "grey", size = 0.5) +
    xlab("Tissue Sample") + ylab("absolute standardized (nominal) beta value") +
    scale_x_continuous(breaks=c(1,2), labels=c("CVS", "Placenta"), limits=c(0.8, 2.3)) +
    scale_y_continuous(limits=c(0, 1.0), breaks= seq(0,1, by = 0.2)) +
    theme_classic()+
    annotate("text",x=0.8, y=1.0, label= textlabel, hjust = 0, vjust = 1, family = "Helvetica", size = standard_textsize/.pt) +
    labs(title = plottitle) +
    # labs(subtitle = get_test_label(t.test_common_eqtms, detailed= TRUE))+
    theme(
      axis.title.x = element_text(family = "Helvetica", size=standard_textsize),
      axis.text.x= element_text(family = "Helvetica", size=standard_textsize),
      axis.title.y = element_text(family = "Helvetica", size=standard_textsize),
      axis.text.y= element_text(family = "Helvetica", size=standard_textsize),
      plot.subtitle=element_text(family = "Helvetica", size=standard_textsize, face="bold"),
      legend.text = element_text(family = "Helvetica", size = standard_textsize),
      plot.title=element_text(family = "Helvetica", size=standard_textsize, face="bold"))
} 





