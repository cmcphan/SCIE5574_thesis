library(reticulate)
use_condaenv('r-reticulate') # Default environment created through conda_create()
library(anndata)
library(SingleCellExperiment)
library(biomaRt)
library(hash)


#########
# Exploration/Testing

hum_sce_final = readRDS("sepp_data/original_data/hum_sce_final.rds")
unique(hum_sce_final$batch) # 38 batches
dim(assay(hum_sce_final)) # Access actual matrix using assay() function
# This data has counts stored as a sparse matrix with dimensions 27715x180956
unique(hum_sce_final$Tissue) # All tissues are cerebellar
unique(hum_sce_final$`Capture System`) # Either Chromium or Chromium v3
# Differences between these that may impact data handling?
head(names(rowRanges(hum_sce_final))) # Get rowNames for matrix (i.e. gene names)
# Gene names for this data given as ENSG num
# or...
head(rownames(rowData(hum_sce_final))) # rowData returns DataFrame
assays(hum_sce_final) # Only one matrix available, called 'umi'
colData(hum_sce_final) # `obs` equivalent - first col is index, other metadata follows after
colnames(colData(hum_sce_final)) # `obsnames`
assay(hum_sce_final)
test_row = assay(hum_sce_final)['ENSG00000188290',] # Grab an example row to 
# test against to make sure transposition works properly
transposed_matrix = t(assay(hum_sce_final))
unique(transposed_matrix[,'ENSG00000188290'] == test_row) # Seems to have worked

# Look for gene datasets
ensembl = useEnsembl(biomart='genes', dataset='hsapiens_gene_ensembl')
searchDatasets(mart=ensembl, pattern='hsapiens') # Human, Homo sapiens, version GRCh38.p14
searchDatasets(mart=ensembl, pattern='mmus') # Mouse, Mus musculus, version GRCm39
searchDatasets(mart=ensembl, pattern='mdom') # Opossum, Monodelphis domestica, version ASM229v1
# Databases accessed 15/11/24

# Verify everything works as intended
unique(ad$obs == colData(hum_sce_final))
# All TRUE except for some NA entries as described below
ad$obs['cell_type'] == colData(hum_sce_final)['cell_type']
colData(hum_sce_final)['SN308_HUM_newborn_CCGGACAAGATACTGA-1', 'cell_type']
colData(hum_sce_final)['SN308_HUM_newborn_CCGGACAAGATACTGA-1', 'subtype']
ad['SN308_HUM_newborn_CCGGACAAGATACTGA-1']$obs$cell_type == colData(hum_sce_final)['SN308_HUM_newborn_CCGGACAAGATACTGA-1', 'cell_type']
# Matching NAs return NA
unique(ad$var_names == rownames(rowData(hum_sce_final)))
# Data seems to match as intended based on col/row checks
sum(data.frame(ad$X[999,]) == data.frame(assay(hum_sce_final)[,999]))
sum(data.frame(ad$X[,999]) == data.frame(assay(hum_sce_final)[999,]))

######
# Convert to Anndata

# Unabbreviated species names with ensembl database names and gene symbol filter names
species_info = hash('hum'=c('human', 'hsapiens_gene_ensembl', 'hgnc_symbol'), 
                    'mou'=c('mouse', 'mmusculus_gene_ensembl', 'mgi_symbol'),
                    'opo'=c('opossum', 'mdomestica_gene_ensembl', 'hgnc_symbol'))
# Mouse gene symbols from the Mouse Genome Informatics (MGI database)
# Opossum gene symbols are *projected* HGNC symbols

for(dataset in list.files(path='sepp_data/original_data/')){
  data = readRDS(paste("sepp_data/original_data/",dataset, sep=''))
  species = strsplit(dataset, '_')[[1]][1]
  # Get gene names from ensembl
  ensembl = useEnsembl(biomart='genes', dataset=species_info[[species]][2])
  gene_names = getBM(attributes=c('ensembl_gene_id', species_info[[species]][3]), filters='ensembl_gene_id',
                     values=rownames(rowData(data)), mart=ensembl, useCache=FALSE)
  names(gene_names) = c('ensembl_gene_id', 'gene_symbol') # Make naming consistent across species
  # Convert to dictionary
  gene_names_dict = hash()
  for(i in 1:length(gene_names$ensembl_gene_id)){ # Build dictionary of gene names
    row = gene_names[i,]
    id = row$ensembl_gene_id
    name = row$gene_symbol
    if(is.null(gene_names_dict[[id]])){ # If there's already an entry for this ID, skip it
      gene_names_dict[[id]] = name
    }
  }
  # Build anndata object
  transposed_matrix = t(assay(data))
  cols = data.frame(colData(data))
  rows = data.frame(row.names='ensembl_gene_id')
  for(gene in rownames(rowData(data))){
    rows = rbind(rows, c(gene, gene_names_dict[[gene]]))
  }
  colnames(rows) = c('ensembl_gene_id', 'gene_symbol')
  row.names(rows) = rows$ensembl_gene_id
  ad = AnnData(
    X = transposed_matrix,
    obs = cols,
    var = rows,
  )
  # Verify conversion worked
  if(!unique(ad$var_names == rownames(rowData(data)))){
    print(paste('var_names do not match for',species_info[[species]][1]))
    return(1)
  }
  if(!unique(ad$obs_names == rownames(colData(data)))){
    print(paste('obs_names do not match for',species_info[[species]][1]))
    return(1)
  }
  # Choose 10 random indices from anywhere in the dataset and check the counts match between the 
  # new anndata object and the original dataset in both directions
  for(i in sample(1:length(rownames(colData(data))), 10)){
    if(!unique(data.frame(ad$X[i,]) == data.frame(assay(data)[,i]))){
      print(paste('Matrices do not match at ad row index',i,'for',species_info[[species]][1]))
      return(1)
    }
  }
  for(i in sample(1:length(rownames(rowData(data))), 10)){
    if(!unique(data.frame(ad$X[,i]) == data.frame(assay(data)[i,]))){
      print(paste('Matrices do not match at ad column index',i,'for',species_info[[species]][1]))
      return(1)
    }
  }
  write_h5ad(ad, paste('sepp_data/sepp_', species_info[[species]][1], '.h5ad', sep=''))
}
