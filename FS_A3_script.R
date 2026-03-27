install.packages("vcfR")
install.packages("HardyWeinberg")
install.packages("LDlinkR")
install.packages("randomForest")

library(vcfR)
library(readr)
library(HardyWeinberg)
library(LDlinkR)
library(vegan)
library(ggplot2)
library(ggpubr)
library(randomForest)
library(glue)

# Load Data ----
# Load VCF files
snps_DPB1 <- read.vcfR("HLADPB1.vcf.gz")
snps_DRB1 <- read.vcfR("HLADQA1.vcf.gz")
snps_DQA1 <- read.vcfR("HLADRB1.vcf.gz")

# Load population metadata file
pop_info <- read_tsv("igsr_samples.tsv")
df_pop_info <- as.data.frame(pop_info)

# Define populations of interest: 
pops_of_interest <- c("YRI", "JPT", "FIN", "CHB")

# Define Helper Functions ----
# EDA Stats
run_EDA <- function(x){
  str(x)
  cat("Dimensions:", dim(x), "\n")
  print(summary(x))
}

# Remove invariant columns
invar_omit <- function(x){
  cat("Dataframe of dim",dim(x), "processed...\n")
  sds <- apply(x, 2, sd, na.rm = TRUE)
  i_var0 <- which(sds == 0)
  if(length(i_var0) > 0){
    x <- x[, -i_var0]
  }
  return(x)                      
}

# Extract percent variation explained for PCs
PCA_variation <- function(pca_summary, PCs = 2){
  var_explained <- pca_summary$importance[2,1:PCs]*100
  var_explained <- round(var_explained, 1)
  return(var_explained)
}

# Run PCA for each gene: HLA-DPB1, HLA-DRB1, HLA-DQA1 ----
# Function to run analysis ----
run_pca <- function(snp_data, gene_name, df_pop_info, pops_of_interest) {
  # 1. Extract the genotype data from the VCF file
  snps_num <- vcfR::extract.gt(snp_data, 
                   element = "GT", 
                   IDtoRowNames = F, 
                   as.numeric = T, 
                   convertNA = T, 
                   return.alleles = F)
 
  # 2. Transpose → rows = individuals and check data
  snps_num_t <- as.data.frame(t(snps_num))
  run_EDA(snps_num_t) # use EDA function to look at structure and summary of data

  # 3. Look at NA values 
  total_NAs <- sum(is.na(snps_num_t))
  print(glue('Number of NAs in {gene_name}: {total_NAs}'))
  
  # 4. Filter for populations of interest
  # Define population names to keep
  pop_lookup <- setNames(df_pop_info$`Population code`, df_pop_info$`Sample name`)
  keep <- rownames(snps_num_t)[pop_lookup[rownames(snps_num_t)] %in% pops_of_interest]
  snps_num_t <- snps_num_t[keep, ]

  # 5. Remove invariants
  snps_no_invar <- invar_omit(snps_num_t)
  
  # 6. Scale SNPs
  snps_scaled <- scale(snps_no_invar)
  
  # 7. Run PCA
  pca_scaled <- prcomp(snps_scaled)
  pca_summary <- summary(pca_scaled)
  
  # 8. Extract variation explained and scores for each PC
  var_explained <- PCA_variation(pca_summary, PCs = 10)
  pca_scores <- as.data.frame(vegan::scores(pca_scaled))
  pca_scores$pop_id <- pop_lookup[rownames(pca_scores)]
  pca_scores$gene <- gene_name

  return(list(
    pca = pca_scaled, 
    scores = pca_scores, 
    var = var_explained, 
    gene = gene_name, 
    n_snps = ncol(snps_scaled)
  ))
}

## HLA-DPB1 ----
results_DPB1 <- run_pca(snps_DPB1, "HLA-DPB1", df_pop_info, pops_of_interest)

### Plot results ----
# Scree plot for scaled PCA with relative importance
screeplot(results_DPB1$pca, 
          ylab = "Relative Importance", 
          main = "Scree plot - A. DPB1 (Relative Importance)")

# Bar plot for PCA with percent variance explained
barplot(results_DPB1$var,
        main = "Percent variation Scree plot - HLA-DPB1",
        ylab = "Percent variation explained")
abline(h = 1/results_DPB1$n_snps * 100, col = 2, lwd = 2) # Eigenvalue above 1 (this line) captures meaningful data

# PCA plot 
ggpubr::ggscatter(data = results_DPB1$scores,
                  y = "PC2",
                  x = "PC1",
                  color = "pop_id",
                  shape = "pop_id",
                  xlab = "PC1",
                  ylab = "PC2",
                  main = "SNPs PCA Scatterplot")


## HLA-DRB1 ----
results_DRB1 <- run_pca(snps_DRB1, "HLA-DRB1", df_pop_info, pops_of_interest)

### Plot results ----
# Scree plot for scaled PCA with relative importance
screeplot(results_DRB1$pca, 
          ylab = "Relative Importance", 
          main = "Scree plot - B. DRB1 (Relative Importance)")

# Bar plot for PCA with percent variance explained
barplot(results_DRB1$var,
        main = "Percent variation Scree plot - HLA-DRB1",
        ylab = "Percent variation explained")
abline(h = 1/results_DRB1$n_snps * 100, col = 2, lwd = 2) # Eigenvalue above 1 (this line) captures meaningful data

# PCA plot 
ggpubr::ggscatter(data = results_DRB1$scores,
                  y = "PC2",
                  x = "PC1",
                  color = "pop_id",
                  shape = "pop_id",
                  xlab = "PC1",
                  ylab = "PC2",
                  main = "SNPs PCA Scatterplot")

## HLA-DQA1 ----
results_DQA1 <- run_pca(snps_DQA1, "HLA-DQA1", df_pop_info, pops_of_interest)

### Plot results ----
# Scree plot for scaled PCA with relative importance
screeplot(results_DQA1$pca, 
          ylab = "Relative Importance", 
          main = "Scree plot - C. DQA1 (Relative Importance)")

# Bar plot for PCA with percent variance explained
barplot(results_DQA1$var,
        main = "Percent variation Scree plot - HLA-DQA1",
        ylab = "Percent variation explained")
abline(h = 1/results_DRB1$n_snps * 100, col = 2, lwd = 2) # Eigenvalue above 1 (this line) captures meaningful data

# PCA plot 
ggpubr::ggscatter(data = results_DQA1$scores,
                  y = "PC2",
                  x = "PC1",
                  color = "pop_id",
                  shape = "pop_id",
                  xlab = "PC1",
                  ylab = "PC2",
                  main = "SNPs PCA Scatterplot")
