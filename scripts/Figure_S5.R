
library(dplyr)
library(GSVA)
library(doParallel)
library(meta)
library(ggplot2)
library(ggpubr)
library(gridExtra)

# Output folders
ssgsea_dir <- "./oxphos/data/ssgseas/FigS5/"
figures_dir <- "./oxphos/figures/"
stats_dir <- "./oxphos/results/"

# Variables to correlate
corr_a <- c("PDK1")
corr_b <- c("BACH1",
            "HIF1A", 
            "HIF1A-SPECIFIC TARGETS")

# Create folders (if they don't exist already)
dir.create(ssgsea_dir, recursive = TRUE)
dir.create(figures_dir, recursive = TRUE)
dir.create(stats_dir, recursive = TRUE)

# List of all TCGA files
tcga_dir <- "./oxphos/data/TCGA/FPKM-UQ_filtered/"
tcga_files <- list.files(path = tcga_dir, recursive = TRUE, full.names = TRUE)

# Read all gene lists 
geneset_paths <- list(
  KEGG_OXIDATIVE_PHOSPHORYLATION = "./oxphos/data/genesets/HIF1A_specific_targets.txt"
)
signatures_figS5 <- lapply(geneset_paths, scan, what = character())

# Update gene set names
names(signatures_figS5) <- c("HIF1A-SPECIFIC TARGETS")

# Function to compute ssGSEA score using the GSVA library
compute_ssgsea <- function(tcga_file, ssgsea_dir, signatures_figS5){
  tcga <- gsub(".*/(TCGA-[A-Z]+)_rnaSeq\\.rds", "\\1", tcga_file)
  rnaSeq <- readRDS(tcga_file) %>% as.data.frame()
  rownames(rnaSeq) <- rnaSeq$hgnc_symbol
  rnaSeq <- rnaSeq %>% dplyr::select(-ensembl_probe_id, -ensembl_gene_id, -hgnc_symbol)
  
  ssGSEAs <- gsva(as.matrix(rnaSeq), (signatures_figS5), 
                  method = "ssgsea", 
                  kcdf = "Gaussian", 
                  parallel.sz = detectCores()-2)
  ssGSEAs <- as.data.frame(t(ssGSEAs))
  # Save it
  saveRDS(ssGSEAs, paste0(ssgsea_dir, tcga, "_ssGSEA.rds"))
}

# Function to create an empty DF to fill later
empty_df <- function(){
  return(data.frame(tcga = character(), n = integer(), corr = double(), 
                    lower.ci = double(), upper.ci = double(), p_val = double(),
                    color = character()))
}

# Function that computes the Pearson correlation between two variables and 
# prepare values for a Forest plot
calc_corr_ci <- function(x, y) {
  corr <- cor.test(x, y, method="pearson")
  pval <- corr$p.value
  lower.ci <- corr$conf.int[[1]]
  upper.ci <- corr$conf.int[[2]]
  color <- ifelse(upper.ci < 0, "blue", ifelse(lower.ci > 0, "red", "grey"))
  return(c(corr = corr$estimate, lower.ci = lower.ci, upper.ci = upper.ci, pval = pval, color = color))
}

# Generate a row with the statistics needed for the forest plot
forest_df <- function(tcga_file, ssgsea_dir, rnaseq, corr_var_1, corr_var_2){

  tcga <- gsub(".*/(TCGA-[A-Z]+)_rnaSeq\\.rds", "\\1", tcga_file)
  ssGSEAs <- readRDS(paste0(ssgsea_dir, tcga, "_ssGSEA.rds"))
  ssGSEAs <- merge(ssGSEAs, rnaseq, by = 0)
  corr_res <- calc_corr_ci(ssGSEAs[,corr_var_1],ssGSEAs[,corr_var_2])
  
  return(data.frame(tcga = paste0(tcga, " (",nrow(ssGSEAs),")"),
                    n = nrow(ssGSEAs),
                    corr = as.double(corr_res[[1]]), 
                    lower.ci = as.double(corr_res[[2]]), 
                    upper.ci = as.double(corr_res[[3]]), 
                    pval = corr_res[[4]],
                    color = corr_res[[5]]))
  
}

add_to_df <- function(df1, df2){
  return(rbind(df1, df2))
}

bach1_forest <- empty_df()
hif1a_forest <- empty_df()
hif1a_st_forest <- empty_df()

# Generate the data frames for the forest plot for each TCGA
generate_df_forest <- function(tcga_file) {

  rnaSeq <- readRDS(tcga_file)
  rnaSeq <- rnaSeq %>% dplyr::filter(hgnc_symbol %in% c(corr_a, corr_b))
  rownames(rnaSeq) <- rnaSeq$hgnc_symbol
  rnaSeq <- rnaSeq %>% dplyr::select(-ensembl_probe_id, -ensembl_gene_id, -hgnc_symbol)
  rnaSeq <- as.data.frame(t(rnaSeq))
  
  tcga_short <- gsub(".*/TCGA-([A-Z]+)_rnaSeq\\.rds", "\\1", tcga_file)
  compute_ssgsea(tcga_file, ssgsea_dir, signatures_figS5)
  
  # Add the correlation values to dataframes to forest plot
  bach1_forest <<- add_to_df(bach1_forest, forest_df(tcga_file, ssgsea_dir, rnaSeq, corr_a, corr_b[1]))
  hif1a_forest <<- add_to_df(hif1a_forest, forest_df(tcga_file, ssgsea_dir, rnaSeq, corr_a, corr_b[2]))
  hif1a_st_forest <<- add_to_df(hif1a_st_forest, forest_df(tcga_file, ssgsea_dir, rnaSeq, corr_a, corr_b[3]))
  
}

lapply(tcga_files, generate_df_forest)

# Function to compute the statistics for the summary
compute_summary <- function(n, corr, lower.ci, upper.ci){
  summary_res <- metacor(cor = corr, n = n)
  return(c(n = sum(summary_res$n), 
           corr = tanh(summary_res$TE.random), 
           lower.ci = tanh(summary_res[["lower.random"]]), 
           upper.ci = tanh(summary_res[["upper.random"]]), 
           pval = summary_res[["pval.random"]]))
}

add_summary <- function(data.corr){
  
  summary.corr <- compute_summary(data.corr$n, 
                                  data.corr$corr)
  tcga_name <- paste0("Summary (", summary.corr[[1]] ,")")
  
  data.corr.sub <- data.frame(tcga = tcga_name,
                              n = summary.corr[[1]],
                              corr = summary.corr[[2]], 
                              lower.ci = summary.corr[[3]], 
                              upper.ci = summary.corr[[4]], 
                              pval = summary.corr[[5]],
                              color = "black")
  data.corr <- rbind(data.corr, data.corr.sub)
  data.corr$tcga <- factor(data.corr$tcga, 
                           levels = rev(data.corr$tcga))
  return(data.corr)
  
}

# Add summary
bach1_forest <- add_summary(bach1_forest)
hif1a_forest <- add_summary(hif1a_forest)
hif1a_st_forest <- add_summary(hif1a_st_forest) 

# Generate the forestplots

create_forest <- function(forest.df, subtitle){
  return(ggplot(forest.df, aes(x=corr, y=tcga, color=color)) +
           geom_point() +
           geom_errorbarh(aes(xmin=lower.ci, xmax=upper.ci), height=0.2) +
           xlim(-1, 1) +
           scale_color_identity() + 
           theme_minimal() +
           labs(x="Pearson Correlation Coefficient", y="TCGA Project",
                title="Forest Plot of Pearson Correlations across TCGA primary tumors",
                subtitle=subtitle) +
           geom_vline(xintercept=0, linetype="dashed", color="purple") + 
           theme(legend.position="none")) # Hide legend
  
}


pdf("./oxphos/figures/FigS5.pdf")
create_forest(bach1_forest,
              "PDK1 gene expression vs BACH1")
create_forest(hif1a_forest,
              "PDK1 gene expression vs HIF1A gene expression")
create_forest(hif1a_st_forest,
              "PDK1 gene expression vs HIF1A-specific targets")
dev.off()

# Save statistics as a TSV file
write.table(bach1_forest, "./oxphos/results/FigS5_PDK1-BACH1_statistics.txt",
            sep="\t", quote = F, row.names = F)
write.table(hif1a_forest, "./oxphos/results/FigS5_PDK1-HIF1A_statistics.txt",
            sep="\t", quote = F, row.names = F)
write.table(hif1a_st_forest, "./oxphos/results/FigS5_HIF1A-targets_statistics.txt",
            sep="\t", quote = F, row.names = F)



