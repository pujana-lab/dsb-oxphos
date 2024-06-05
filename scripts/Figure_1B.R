
library(dplyr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(GSVA)
library(doParallel)

# Output folders
ssgsea_dir <- "./oxphos/data/ssgseas/Fig1B/"
figures_dir <- "./oxphos/figures/"
dir.create(ssgsea_dir, recursive = TRUE)
dir.create(figures_dir, recursive = TRUE)

## Datasets
#Gene expression an IDs downloaded from https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources/Home.html)
ic50_cells <- data.table::fread("./oxphos/data/datasets/Cell_line_RMA_proc_basalExp.txt")
ic50_ids <- data.table::fread("./oxphos/data/datasets/ic50_cells_cosmic_IDs.csv")
# IC50 of Rehydroroterone downloaded from # https://www.cancerrxgene.org/compound/Dihydrorotenone/1827/overview/ic50?screening_set=GDSC2
ic50_rehydro <- data.table::fread("./oxphos/data/datasets/IC50_dihydrorotenone.csv")

## Signatures
# Read all gene lists 
geneset_paths <- list(
  RESPIRATORY_CHAIN_COMPLEX_I = "./oxphos/data/genesets/RESPIRATORY_CHAIN_COMPLEX_I.txt",
  HALLMARK_GLYCOLYSIS = "./oxphos/data/genesets/HALLMARK_GLYCOLYSIS.txt",
  AltEJ = "./oxphos/data/genesets/AltEJ.txt"
)
signatures_fig1B <- lapply(geneset_paths, scan, what = character())

# Update gene set names
names(signatures_fig1B) <- c("RESPIRATORY CHAIN COMPLEX I",
                             "HALLMARK GLYCOLYSIS", 
                             "AltEJ")

# Rename colnames of gene expr. so they match the ones from ic50_ids file
colnames(ic50_cells)[3:ncol(ic50_cells)] <- substr(colnames(ic50_cells)[3:ncol(ic50_cells)], 
                                                   6, nchar(colnames(ic50_cells))[3:ncol(ic50_cells)])
# Remove genes without gene symbol name
ic50_cells <- ic50_cells %>% dplyr::filter(GENE_SYMBOLS != "") %>% as.data.frame()
# Leave only gene expr. for ssGSEA computing
rownames(ic50_cells) <- ic50_cells$GENE_SYMBOLS
ic50_cells <- ic50_cells %>% dplyr::select(-GENE_SYMBOLS, -GENE_title)

# Compute ssGSEAs
ssGSEAs <- gsva(as.matrix(ic50_cells), (signatures_fig1B), 
                method = "ssgsea", 
                kcdf = "Gaussian", 
                parallel.sz = detectCores()-2)

ssGSEAs <- as.data.frame(t(ssGSEAs))

# Leave only useful variables to IC50 of Dihydroroterone
ic50_rehydro <- ic50_rehydro %>% dplyr::select(Sample_Name = `Cell line`, IC50, 
                                               TCGA = `TCGA classification`)
# Merge with metadata file
ic50_rehydro <- merge(ic50_rehydro, ic50_ids, by = "Sample_Name")

# Merge with ssGSEAs by COSMIC ID
all_data <- merge(ssGSEAs, ic50_rehydro, by.x=0, by.y = "COSMIC_ID")

# Scale ssGSEA data for better visualization. Apply log10 to IC50
all_data[names(signatures_fig1B)] <- lapply(all_data[names(signatures_fig1B)], scale)
all_data["IC50"] <- log10(all_data$IC50)

write.table(all_data, paste0(ssgsea_dir, "IC50_ssGSEAs.txt"),
            sep="\t", quote = F, row.names = F)

# Function to generate the scatterplots
generate_plots <- function(df, signature_names) {
  
  plot_list <- list()
  
  # Loop through each column
  for(i in 1:length(signature_names)){
    col_name <- signature_names[i]
    
    # Create the base plot
    plot_list[[i]] <- ggplot(df, aes(x = IC50, y = .data[[col_name]])) +
      geom_point(colour="#A9A9A9") +
      
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour="black") +
      ggtitle(paste("Correlation between IC50 and", col_name)) + 
      theme_minimal() + 
      theme(plot.title = element_text(size = 10, face = "bold"))  + 
      xlab("Dihydrorotenone IC50 (\U03bcM)")
    
    # Manually calculate correlations
    all_corr_pearson <- cor.test(df[[col_name]], df$IC50, method = "pearson")
    all_corr_spearman <- cor.test(df[[col_name]], df$IC50, method = "spearman")

    # Add annotations for correlations
    labels_all <- paste0("All Pearson: ", round(all_corr_pearson$estimate, 4), 
                         "\nAll Spearman: ", round(all_corr_spearman$estimate, 4),
                         "\np-values: ", ifelse(all_corr_pearson$p.value >= 1e-3, 
                                                round(all_corr_pearson$p.value,4 ), 
                                                sprintf("%.2e", all_corr_pearson$p.value)), 
                         "; ", ifelse(all_corr_spearman$p.value >= 1e-3, 
                                      round(all_corr_spearman$p.value,4 ), 
                                      sprintf("%.2e", all_corr_spearman$p.value)))
    plot_list[[i]] <- plot_list[[i]] + annotate("text", x = Inf, y = Inf, 
                                                label = labels_all, hjust = 1, 
                                                vjust = 1, color = "black", size = 3)
    
  }
  
  return(plot_list)
}

# Generate and save the Figure
plot_corr_list <- generate_plots(all_data, names(signatures_fig1B))
pdf(paste0(figures_dir, "Fig1B.pdf"), width = 17, height = 7.5)
ggarrange(plotlist=plot_corr_list, ncol=3, nrow=1)
dev.off()
