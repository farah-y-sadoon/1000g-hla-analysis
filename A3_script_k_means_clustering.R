##Assignment 03

## April 02, 2026

library(vcfR)
library(readr)
library(SNPRelate)
library(HardyWeinberg)
library(LDlinkR)
library(vegan)
library(ggplot2)
library(ggpubr)
library(randomForest)
library(glue)
library(factoextra)
library(dartR)
library(NbClust)
library(cluster)
library(gridExtra)

# Load Data ----

# Load population metadata file
pop_info <- read_tsv("igsr_samples.tsv")
df_pop_info <- as.data.frame(pop_info)

metadata <- df_pop_info %>%
  rename(id = `Sample name`, pop = `Population code`) %>%
  select(id, pop)

# metadata write to csv
write.csv(metadata, "metadata.csv")

# Define populations of interest: 
pops_of_interest <- c("YRI", "JPT", "FIN", "CHB")

# read in vcfs
gl_HLADPB1 <- gl.read.vcf(vcffile = "HLADPB1.vcf.gz", ind.metafile = "metadata.csv")
gl_HLADRB1 <- gl.read.vcf(vcffile = "HLADRB1.vcf.gz", ind.metafile = "metadata.csv")
gl_HLADQA1 <- gl.read.vcf(vcffile = "HLADQA1.vcf.gz", ind.metafile = "metadata.csv")

#HWE filtering
#HLA-DPB1
HLADPB1_filtered_HWE <- gl.filter.hwe(gl_HLADPB1, subset = "each", n.pop.threshold = 1, method_sig = "Exact", multi_comp = FALSE, multi_comp_method = "BH", alpha_val = 0.05, pvalue_type = "midp", cc_val = 0.5, min_sample_size = 5, verbose = 5)

# HLA-DRB1
HLADRB1_filtered_HWE <- gl.filter.hwe(gl_HLADRB1, subset = "each", n.pop.threshold = 1, method_sig = "Exact", multi_comp = FALSE, multi_comp_method = "BH", alpha_val = 0.05, pvalue_type = "midp", cc_val = 0.5, min_sample_size = 5, verbose = 5)

# HLA-DQA1
HLADQA1_filtered_HWE <- gl.filter.hwe(gl_HLADQA1, subset = "each", n.pop.threshold = 1, method_sig = "Exact", multi_comp = FALSE, multi_comp_method = "BH", alpha_val = 0.05, pvalue_type = "midp", cc_val = 0.5, min_sample_size = 5, verbose = 5)

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
run_pca <- function(snp_data, gene_name, metadata, pops_of_interest) {
  # 1. Extract the genotype data from the VCF file
  snps_num <- as.matrix(snp_data)
  
  # 2. Transpose → rows = individuals and check data
  snps_num_t <- as.data.frame(snps_num)
  
  # 3. Look at NA values 
  total_NAs <- sum(is.na(snps_num_t))
  print(glue('Number of NAs in {gene_name}: {total_NAs}'))
  
  # 4. Filter for populations of interest
  # Define population names to keep
  pop_lookup <- setNames(metadata$pop, metadata$id)
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
    scaled_snps = snps_scaled,
    scores = pca_scores, 
    var = var_explained, 
    gene = gene_name, 
    n_snps = ncol(snps_scaled)
  ))
}

## HLA-DPB1 ----
results_DPB1 <- run_pca(HLADPB1_filtered_HWE, "HLA-DPB1", metadata, pops_of_interest)

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
                  xlab = paste0("PC1", " [", results_DPB1$var[1], "%", "]"),
                  ylab = paste0("PC2", " [", results_DPB1$var[2], "%", "]"),
                  main = "SNPs PCA Scatterplot for HLA-DPB1")

#One outlier? need to figure out which one it is
outlier_id <- rownames(results_DPB1$scores)[which.min(results_DPB1$scores$PC1)]

#generate clean df_pop_info with NA18915 sample removed
metadata_clean <- metadata %>%
  filter(id != outlier_id)

sum(metadata_clean$id== outlier_id) #is 0 which is good

# Remove NA18915 from post-HW filter
gl_clean_HLADPB1 <- HLADPB1_filtered_HWE[!indNames(HLADPB1_filtered_HWE) %in% outlier_id]

#Rerun function again
results_DPB1_2 <- run_pca(gl_clean_HLADPB1, "HLA-DPB1", metadata, pops_of_interest)

#Screenplot without outlier
screeplot(results_DPB1_2$pca, 
          ylab = "Relative Importance", 
          main = "Scree plot - A. DPB1 (Relative Importance)")

# Bar plot for PCA with percent variance explained
barplot(results_DPB1_2$var,
        main = "Percent variation Scree plot - HLA-DPB1",
        ylab = "Percent variation explained")
abline(h = 1/results_DPB1_2$n_snps * 100, col = 2, lwd = 2) # Eigenvalue above 1 (this line) captures meaningful data

#plot with new dataset without outlier
ggpubr::ggscatter(data = results_DPB1_2$scores,
                  y = "PC2",
                  x = "PC1",
                  color = "pop_id",
                  shape = "pop_id",
                  xlab = paste0("PC1", " [", results_DPB1_2$var[1], "%", "]"),
                  ylab = paste0("PC2", " [", results_DPB1_2$var[2], "%", "]"),
                  main = "SNPs PCA Scatterplot for HLA-DPB1")

## HLA-DRB1 ----
results_DRB1 <- run_pca(HLADRB1_filtered_HWE, "HLA-DRB1", metadata, pops_of_interest)

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
                  xlab = paste0("PC1", " [", results_DRB1$var[1], "%", "]"),
                  ylab = paste0("PC2", " [", results_DRB1$var[2], "%", "]"),
                  main = "SNPs PCA Scatterplot for HLADRB1")

## HLA-DQA1 ----
results_DQA1 <- run_pca(HLADQA1_filtered_HWE, "HLA-DQA1", metadata, pops_of_interest)

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
                  xlab = paste0("PC1", " [", results_DQA1$var[1], "%", "]"),
                  ylab = paste0("PC2", " [", results_DQA1$var[2], "%", "]"),
                  main = "SNPs PCA Scatterplot for HLADQA1")


#############################################
# Non - Hierarchical clustering using k means
#############################################

## Find optimal k means with several methods

# For HLA-DPB1
# 1. Using NbClust visualization method on scaled values
set.seed(321)
fviz_nbclust(results_DPB1_2$scaled_snps, kmeans, method = "wss") +
  geom_vline(xintercept = 3, linetype = 2)+
  labs(subtitle = "Elbow method HLA-DPB1") #elbow plot seems to dip at 3 

# 2. Using the within cluster sums of square (wss) method
set.seed(321)
wss_DPB1 <- (nrow(results_DPB1_2$scaled_snps)-1)*sum(apply(results_DPB1_2$scaled_snps,2,var))

wss_DPB1[2:15] <- unlist(
  lapply(2:15,function(i){sum(kmeans(results_DPB1_2$scaled_snps,centers=i)$withinss)})
)

par(mfrow=c(1,1))
plot(1:15, wss_DPB1, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")
# 3-4 clusters, dip at 3-4

# 3. Using gap statistics to determine the best k number for clustering
set.seed(321)

#Using scaled SNPs
nclust_gapstat_DPB1 <- clusGap(results_DPB1_2$scaled_snps, 
                               kmeans, nstart = 25, d.power=2,
                          K.max = 10, B = 100)

#plot gap statistics plot
fviz_gap_stat(nclust_gapstat_DPB1)
#10 clusters is the best?

#Using the first 10 PCAs?
#nclust_gapstat_DPB1_pca <- clusGap(results_DPB1_2$scores[,1:10], 
                         #      kmeans, nstart = 25, d.power=2,
                       #        K.max = 10, B = 100)

#plot gap statistics plot
#fviz_gap_stat(nclust_gapstat_DPB1_pca)
#1 cluster is the best now?

# Use NBClust from NbClust library method use on first 10 PCAs instead of scaled values, cuz scaled snps keeps giving an error
nbout_DPB1 <- NbClust(results_DPB1_2$scores[, 1:10], method="kmeans")

barplot(table(nbout_DPB1$Best.nc[1,]), ylab = "Indexes", 
        xlab = "Number of clusters (k)",
        main = "Optimal Number of Clusters for HLA-DPB1")
#optimal number of cluster is 3

#plot to see the best k cluster number
#Compare and contrast number of k means clustering between 2-4 and 10 from the methods determined above and plot
p_DPB1_k2 <- fviz_cluster(kmeans(results_DPB1_2$scaled_snps,2,nstart=10),data=results_DPB1_2$scaled_snps) + theme_minimal()

p_DPB1_k3 <- fviz_cluster(kmeans(results_DPB1_2$scaled_snps,3,nstart=10),data=results_DPB1_2$scaled_snps) + theme_minimal()

p_DPB1_k4 <- fviz_cluster(kmeans(results_DPB1_2$scaled_snps,4,nstart=10),data=results_DPB1_2$scaled_snps) + theme_minimal()

p_DPB1_k10 <- fviz_cluster(kmeans(results_DPB1_2$scaled_snps,10,nstart=10),data=results_DPB1_2$scaled_snps) + theme_minimal()

#compare and contrast 
plots_DPB1 <- grid.arrange(p_DPB1_k2, p_DPB1_k3, p_DPB1_k4, p_DPB1_k10, ncol = 2, nrow = 2)
#From the plots above, it seems that k means of 3 may be ideal since it separates the groups nicely
