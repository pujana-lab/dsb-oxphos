
library(dplyr)
library(GSVA)
library(doParallel)
library(meta)
library(ggplot2)
library(ggpubr)
library(gridExtra)

# Output folders
ssgsea_dir <- "./oxphos/data/ssgseas/Fig1A/"
figures_dir <- "./oxphos/figures/"
stats_dir <- "./oxphos/results/"

# Variables to correlate
corr_a <- c("AltEJ")
corr_b <- c("KEGG OXIDATIVE PHOSPHORYLATION",
            "MOOTHA VOXPHOS", 
            "RESPIRATORY CHAIN COMPLEX I")

# Create folders (if they don't exist already)
dir.create(ssgsea_dir, recursive = TRUE)
dir.create(figures_dir, recursive = TRUE)
dir.create(stats_dir, recursive = TRUE)

# List of all TCGA files
tcga_dir <- "./oxphos/data/TCGA/FPKM-UQ_filtered/"
tcga_files <- list.files(path = tcga_dir, recursive = TRUE, full.names = TRUE)

# Read all gene lists 
geneset_paths <- list(
  AltEJ = "./oxphos/data/genesets/AltEJ.txt",
  KEGG_OXIDATIVE_PHOSPHORYLATION = "./oxphos/data/genesets/KEGG_OXIDATIVE_PHOSPHORYLATION.txt",
  MOOTHA_VOXPHOS = "./oxphos/data/genesets/MOOTHA_VOXPHOS.txt",
  RESPIRATORY_CHAIN_COMPLEX_I = "./oxphos/data/genesets/RESPIRATORY_CHAIN_COMPLEX_I.txt"
)
signatures_fig1A <- lapply(geneset_paths, scan, what = character())

# Update gene names
signatures_fig1A$AltEJ[signatures_fig1A$AltEJ=="MRE11A"] <- "MRE11"
signatures_fig1A$AltEJ[signatures_fig1A$AltEJ=="C19orf40"] <- "FAAP24"

# Update gene set names
names(signatures_fig1A) <- c("AltEJ",
                             "KEGG OXIDATIVE PHOSPHORYLATION",
                             "MOOTHA VOXPHOS", 
                             "RESPIRATORY CHAIN COMPLEX I")

# Function to compute ssGSEA score using the GSVA library
compute_ssgsea <- function(tcga_file, ssgsea_dir, signatures_fig1A){
  tcga <- gsub(".*/(TCGA-[A-Z]+)_rnaSeq\\.rds", "\\1", tcga_file)
  rnaSeq <- readRDS(tcga_file) %>% as.data.frame()
  rownames(rnaSeq) <- rnaSeq$hgnc_symbol
  rnaSeq <- rnaSeq %>% dplyr::select(-ensembl_probe_id, -ensembl_gene_id, -hgnc_symbol)
  
  ssGSEAs <- gsva(as.matrix(rnaSeq), (signatures_fig1A), 
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
forest_df <- function(tcga_file, ssgsea_dir, corr_var_1, corr_var_2){

  tcga <- gsub(".*/(TCGA-[A-Z]+)_rnaSeq\\.rds", "\\1", tcga_file)
  ssGSEAs <- readRDS(paste0(ssgsea_dir, tcga, "_ssGSEA.rds"))
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

kegg_oxphos_forest <- empty_df()
mootha_oxphos_forest <- empty_df()
respiratory_chain_forest <- empty_df()

# Generate the data frames for the forest plot for each TCGA
generate_df_forest <- function(tcga_file) {

  tcga_short <- gsub(".*/TCGA-([A-Z]+)_rnaSeq\\.rds", "\\1", tcga_file)
  compute_ssgsea(tcga_file, ssgsea_dir, signatures_fig1A)
  
  # Add the correlation values to dataframes to forest plot
  kegg_oxphos_forest <<- add_to_df(kegg_oxphos_forest, forest_df(tcga_file, ssgsea_dir, corr_a, corr_b[1]))
  mootha_oxphos_forest <<- add_to_df(mootha_oxphos_forest, forest_df(tcga_file, ssgsea_dir, corr_a, corr_b[2]))
  respiratory_chain_forest <<- add_to_df(respiratory_chain_forest, forest_df(tcga_file, ssgsea_dir, corr_a, corr_b[3]))
  
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
kegg_oxphos_forest <- add_summary(kegg_oxphos_forest)
mootha_oxphos_forest <- add_summary(mootha_oxphos_forest)
respiratory_chain_forest <- add_summary(respiratory_chain_forest) 

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


pdf("./oxphos/figures/Fig1A.pdf")
create_forest(kegg_oxphos_forest,
              "AltEJ vs KEGG Oxidative Phosphorylation")
create_forest(mootha_oxphos_forest,
              "AltEJ vs Mootha Voxphos")
create_forest(respiratory_chain_forest,
              "AltEJ vs Respiratory Chain Complex I")
dev.off()

# Save statistics as a TSV file
write.table(kegg_oxphos_forest, "./oxphos/results/Fig1A_AltEJ-KEGGoxphos_statistics.txt",
            sep="\t", quote = F, row.names = F)
write.table(mootha_oxphos_forest, "./oxphos/results/Fig1A_AltEJ-MoothaOxphos_statistics.txt",
            sep="\t", quote = F, row.names = F)
write.table(respiratory_chain_forest, "./oxphos/results/Fig1A_AltEJ-RespChain_statistics.txt",
            sep="\t", quote = F, row.names = F)









