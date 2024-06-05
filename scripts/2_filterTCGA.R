
library(dplyr)
library(tibble)

# FPKM-UQ data can be downloaded from GDC and from Xena UCSC: 
# https://gdc.xenahubs.net/

# Read selected protein genes
sel_genes <- readRDS("./oxphos/data/selected_protein_genes.rds")

# Define directory paths
output_dir <- "./oxphos/data/TCGA/FPKM-UQ_filtered/"
tcga_directory_path  <- "./oxphos/data/tcga/FPKM-UQ/"
dir.create(output_dir, recursive = TRUE)

# List of all TCGA files
tcga_files <- list.files(path = tcga_directory_path, recursive = TRUE, full.names = TRUE)

# Function to process each TCGA file
process_tcga <- function(file_path) {
  
  rnaSeq <- data.table::fread(file_path)
  
  # Get TCGA:
  tcga <- gsub(".*/(TCGA-[A-Z]+)\\..*", "\\1", file_path)
  
  # Select protein coding genes
  rnaSeq <- rnaSeq %>% 
    mutate(Ensembl_ID = gsub("\\..*", "", Ensembl_ID)) %>%
    filter(Ensembl_ID %in% sel_genes$ensembl_gene_id) %>%
    column_to_rownames(var = "Ensembl_ID")
  
  # Remove non primary tumour samples
  primary_samples <- colnames(rnaSeq)[substr(colnames(rnaSeq), 14, 15) == "01"]
  rnaSeq <- rnaSeq %>%
    select(all_of(primary_samples))
  
  # Filter genes with more than 20% of zeroes
  zero_counts <- rowSums(rnaSeq == 0)
  rnaSeq <- rnaSeq[zero_counts <= floor(0.2 * ncol(rnaSeq)),]
  
  # Add gene info to the TCGA matrix
  all_data <- merge(sel_genes, rnaSeq, by.x = "ensembl_gene_id", by.y = "row.names")
  
  # Save filtered TCGA
  saveRDS(all_data, paste0("./oxphos/data/TCGA/FPKM-UQ_filtered/",tcga,"_rnaSeq.rds"))
  
}

# Process each TCGA
lapply(tcga_files, process_tcga)


