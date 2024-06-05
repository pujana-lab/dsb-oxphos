
library(dplyr)
library(biomaRt)

# Read the FPKM-UQ RNAseq data to get all genes to filter
rnaSeq <- data.table::fread("./oxphos/data/tcga/FPKM-UQ/TCGA-ACC.htseq_fpkm-uq.tsv")
rnaSeq$Ensembl_gene_ID <- gsub("\\..*", "", rnaSeq$Ensembl_ID)
rnaSeq <- rnaSeq %>% dplyr::select(ensembl_probe_id = Ensembl_ID, ensembl_gene_id = Ensembl_gene_ID)

# BiomaRt query
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
sel_genes <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
                    filters = "ensembl_gene_id",
                    values = unique(rnaSeq$ensembl_gene_id),
                    mart = mart)

# Select protein coding genes with a gene name
sel_genes <- sel_genes %>%
  filter(gene_biotype == "protein_coding" & hgnc_symbol != "")

# Remove ensembl gene ids from duplicated gene names
duplicated_genes <- sel_genes[duplicated(sel_genes$hgnc_symbol),]
sel_genes <- sel_genes %>% dplyr::filter(!ensembl_gene_id %in% 
                                             c("ENSG00000279195",
                                               "ENSG00000254093"))
# Remove unneeded columns
sel_genes <- merge(sel_genes, rnaSeq, by = "ensembl_gene_id") %>% 
  dplyr::select(-gene_biotype)

# Save the list of genes to filter
saveRDS(sel_genes, "./oxphos/data/selected_protein_genes.rds")

