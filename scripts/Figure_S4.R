
library(dplyr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(GSVA)
library(doParallel)

# Output folders
ssgsea_dir <- "./oxphos/data/ssgseas/FigS4/"
figures_dir <- "./oxphos/figures/"
dir.create(ssgsea_dir, recursive = TRUE)
dir.create(figures_dir, recursive = TRUE)

## Datasets
# Data downloaded from https://depmap.org/portal/data_page/?tab=allData
drugs_df <- data.table::fread("./oxphos/data/datasets/Repurposing_Public_23Q2_Extended_Primary_Data_Matrix.csv")
metadata <- data.table::fread("./oxphos/data/datasets//Repurposing_Public_23Q2_Cell_Line_Meta_Data.csv")
cell_lines <- data.table::fread("./oxphos/data/datasets/Expression_Public_23Q4.csv") %>% as.data.frame()

## Signatures
# Read all gene lists 
geneset_paths <- list(
  RESPIRATORY_CHAIN_COMPLEX_I = "./oxphos/data/genesets/RESPIRATORY_CHAIN_COMPLEX_I.txt",
  HALLMARK_GLYCOLYSIS = "./oxphos/data/genesets/HALLMARK_GLYCOLYSIS.txt",
  AltEJ = "./oxphos/data/genesets/AltEJ.txt"
)
signatures_figS4 <- lapply(geneset_paths, scan, what = character())

# Update gene set names
names(signatures_figS4) <- c("RESPIRATORY CHAIN COMPLEX I",
                             "HALLMARK GLYCOLYSIS", 
                             "AltEJ")

# Select which drugs to analyze
sel_drugs <- data.frame(drug_id = c("BRD:BRD-K00003105-003-01-9",
                                    "BRD:BRD-A81541225-001-06-8"),
                        drug_name = c("IACS-10759", "OLIGOMYCIN-A"))
drugs_df <- merge(sel_drugs, drugs_df, by.x = "drug_id", by.y = "V1")

# Rename drugs
drugs_df <- drugs_df %>%
  as.data.frame() %>%
  column_to_rownames(var = "drug_name") %>%
  select(-drug_id) %>%
  mutate(drug = rownames(.)) %>%
  t() %>% as.data.frame()
drugs_df$depmap_id <- rownames(drugs_df)

metadata <- metadata %>% dplyr::select(depmap_id)

# Reformat gene expression dataset for ssGSEA
cell_lines <- cell_lines %>%
  as.data.frame() %>%
  column_to_rownames(var = "V1")

# Compute ssGSEA scores
ssGSEAs <- gsva(t(cell_lines), (signatures_figS4), 
                method = "ssgsea", 
                kcdf = "Gaussian", 
                parallel.sz = detectCores()-2)
ssGSEAs <- as.data.frame(t(ssGSEAs))
ssGSEAs$depmap_id <- rownames(ssGSEAs)

# Merge all data
all_data <- merge(metadata, ssGSEAs, by = "depmap_id")
all_data <- merge(drugs_df, all_data, by = "depmap_id")

# Scale ssGSEA data for better visualization. 
all_data[names(signatures_figS4)] <- lapply(all_data[names(signatures_figS4)], scale)
all_data[c("OLIGOMYCIN-A", "IACS-10759")] <- lapply(all_data[c("OLIGOMYCIN-A", "IACS-10759")], as.double)

write.table(all_data, paste0(ssgsea_dir, "Depmap_ssGSEAs.txt"),
            sep="\t", quote = F, row.names = F)

# Function to generate the scatterplots
generate_plots <- function(df, drug_name, signature_names) {

  plot_list <- list()

  # Loop through each COL_i to create a plot
  for(i in 1:length(signature_names)) {
   
    col_name <- signature_names[i]
    
    # Create the base plot
    plot_list[[i]] <- ggplot(df, aes(x = .data[[drug_name]], y = .data[[col_name]])) +
      geom_point(colour="#A9A9A9") +
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour="black") +
      ggtitle(paste("Correlation between",drug_name,"\nand", col_name)) + 
      theme_minimal() + 
      theme(plot.title = element_text(size = 10, face = "bold"))  + 
      xlab(paste0(drug_name))
    
    # Manually calculate correlations for all 
    all_corr_pearson <- cor.test(df[[col_name]], df[[drug_name]], method = "pearson")
    all_corr_spearman <- cor.test(df[[col_name]], df[[drug_name]], method = "spearman")

    # Add annotations for correlations
    labels_all <- paste0("All Pearson: ", round(all_corr_pearson$estimate, 4), 
                         "\nAll Spearman: ", round(all_corr_spearman$estimate, 4),
                         "\np-values: ", ifelse(all_corr_pearson$p.value >= 1e-3, 
                                                round(all_corr_pearson$p.value,4 ), 
                                                sprintf("%.2e", all_corr_pearson$p.value)), 
                         "; ", ifelse(all_corr_spearman$p.value >= 1e-3, 
                                      round(all_corr_spearman$p.value,4 ), 
                                      sprintf("%.2e", all_corr_spearman$p.value)))

    # Add annotations
    plot_list[[i]] <- plot_list[[i]] + annotate("text", x = Inf, y = Inf, 
                                                label = labels_all, hjust = 1, 
                                                vjust = 1, color = "black", size = 3)
   
  }
  
  return(plot_list)
}

# Generate and save the Figure
plot_corr_list_1 <- generate_plots(all_data, "OLIGOMYCIN-A", names(signatures_figS4))
plot_corr_list_2 <- generate_plots(all_data, "IACS-10759", names(signatures_figS4))

pdf(paste0(figures_dir, "FigS4.pdf"), width = 17, height = 7.5)
ggarrange(plotlist=plot_corr_list_1, ncol=3, nrow=1)
ggarrange(plotlist=plot_corr_list_2, ncol=3, nrow=1)
dev.off()


