library(mygene) # installed with bioconductor
library(AUCell) # installed with bioconductor
library(GSEABase)
library(tidyverse)

# reading in counts and gene list
setwd("~/glab/cell_line_selection")
gene_list <- read_delim('gene_list.txt',delim ='\t')
counts <- read_table2('CCLE_RNAseq_rsem_genes_tpm_20180929.zip')

################################################################################
# ----------- adding in ENS identifiers to gene list ---------------------------

gene_list <- as_tibble(queryMany(gene_list,scopes='symbol',fields='all',species='human')) # annotating with ensembl ids
gene_list <- gene_list %>% dplyr::select(c('query','entrezgene','ensembl.gene'))
write_delim(gene_list,'gene_listUpdated.txt',delim='\t')

################################################################################
# --------------------- AUCell use ---------------------------------------------
counts <- counts %>% dplyr::select(-transcript_ids)
counts$gene_id <- sub('\\.[0-9]*$','',counts$gene_id) # removing version information from gene ids

exprMatrix <- as.matrix(counts %>% select(-gene_id)) # saving only the tissue data in the matrix
rownames(exprMatrix) <- counts$gene_id # setting the rownames of the matrix to the gene IDs  

length(intersect(gene_list$ensembl.gene,rownames(exprMatrix))) # checking overlap of genes counted and gene_list

# formatting the gene_list tbl to a GeneSet using the Ensembl gene ids 
geneset <- GeneSet(gene_list$ensembl.gene, setName='geneSet1')

cells_rankings <- AUCell_buildRankings(exprMatrix,plotStats=TRUE) # Build gene expression rankings for each cell
cells_AUC <- AUCell_calcAUC(geneset, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05) #Calculates the 'AUC' for each gene-set in each cell.
cells_assignment <- AUCell_exploreThresholds(cells_AUC,plotHist=TRUE,nCores=1,assign=TRUE)

length(cells_assignment[["geneSet1"]][["assignment"]]) # number of cells that are assigned to the gene list

################################################################################
# -------------------- production of ranked table ------------------------------
# saving and filtering rankings to show only genes in the gene list 
ranked <- as.data.frame(cells_rankings@assays@data@listData[["ranking"]]) 
ranked <- ranked %>% tibble::rownames_to_column("genes")
length(intersect(ranked$genes,gene_list$ensembl.gene)) # checking to see if our genes are present 
ranked <- ranked %>% filter(genes %in% gene_list$ensembl.gene) # filtering for our genes
ranked <- ranked %>% left_join(gene_list,by=c('genes'='ensembl.gene'))
ranked <- ranked %>% dplyr::select(-c('entrezgene','genes')) 
ranked <- ranked %>% rename('genes'='query')

# data reformatting so that the column names are the rank # and the 
    # values are the name of the cell line
pivoted_rank <- ranked %>% pivot_longer(!genes,names_to='cell_line',values_to='rank')
pivoted_rank1 <- pivoted_rank %>% pivot_wider(names_from='genes',values_from='rank')
write_delim(pivoted_rank1,'ranks_unformatted.txt',delim='\t')

# filtering for AUC == 0 
AUC <-  pivot_longer(as.data.frame(cells_AUC@assays@data@listData[["AUC"]]),cols=everything(),names_to='cell_line',values_to = 'AUC') %>% 
  filter(AUC == 0)
ranked_filtered <- pivoted_rank %>% filter(cell_line %in% AUC$cell_line) #filtering our ranked list with the cell lines created above
ranked_filtered_pivoted <- ranked_filtered %>% pivot_wider(names_from='genes',values_from='rank') 
write_delim(ranked_filtered,'ranks_filtered.txt',delim='\t')

################################################################################
# -------------- creation of filtered TPM table --------------------------------
# replacing rank with tpm value in filtered table
  # goal: filter oritinal file for our specific gene names and for the cell liens discovered by AUC filtering
ranked_filtered_tpm <- counts %>% filter(gene_id %in% gene_list$ensembl.gene)
ranked_filtered_tpm_graphable <- ranked_filtered_tpm %>% pivot_longer(!gene_id,values_to='TPM',names_to='cell_line')
#chaning ensembl ids to gene names 
ranked_filtered_tpm <- ranked_filtered_tpm_graphable %>% left_join(gene_list,by=c('gene_id'='ensembl.gene'))
ranked_filtered_tpm <- ranked_filtered_tpm %>% dplyr::select(-c('entrezgene','gene_id')) 
ranked_filtered_tpm <- ranked_filtered_tpm %>% rename('genes'='query')
# fixing table structure
ranked_filtered_tpm <- ranked_filtered_tpm %>% pivot_wider(names_from='genes',values_from='TPM')

write_delim(ranked_filtered_tpm,'tpm_filtered.txt',delim='\t')

################################################################################
# ----------------------------- extra visuals ----------------------------------
# graphing tpm values from the AUC filtered dataset
ggplot(ranked_filtered_tpm_graphable, aes(x=TPM)) + geom_histogram(binwidth=1) +xlim(0,300)+ylim(0,2000) +
  ggtitle('hist of TPM for AUC filtered cell lines')
ggplot(ranked_filtered_tpm_graphable, aes(y=TPM)) + geom_boxplot() +ggtitle('boxplot of TPM for AUC filtered cell lines')

