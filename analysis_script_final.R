# BINF 6970 Assignment 3
# Rebekah Hest, Eva Innocente, Farah Sadoon, Liona Vu
# April 03, 2023

# Load Libraries ----
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
library(mlbench)
library(caret)
library(class)
library(e1071)

# Load Data ----

# Load population metadata file
pop_info <- read_tsv("igsr_samples.tsv")
df_pop_info <- as.data.frame(pop_info)

metadata <- df_pop_info %>%
  rename(id = `Sample name`, pop = `Population code`) %>%
  dplyr::select(id, pop)

# Write metadata to CSV
# write.csv(metadata, "metadata.csv")

# Define populations of interest: 
pops_of_interest <- c("YRI", "JPT", "FIN", "CHB")

# Read in VCF files
gl_HLADPB1 <- gl.read.vcf(vcffile = "HLADPB1.vcf.gz", ind.metafile = "metadata.csv")
gl_HLADRB1 <- gl.read.vcf(vcffile = "HLADRB1.vcf.gz", ind.metafile = "metadata.csv")
gl_HLADQA1 <- gl.read.vcf(vcffile = "HLADQA1.vcf.gz", ind.metafile = "metadata.csv")

# Filter SNPs with HWE ----
# Function to filter with HWE 
HWE_filter <- function(gl_gene) {
  gene_filtered_HWE <- gl.filter.hwe(gl_gene, subset = "each", n.pop.threshold = 1, 
                                     method_sig = "Exact", multi_comp = FALSE, 
                                     multi_comp_method = "BH", alpha_val = 0.05, 
                                     pvalue_type = "midp", cc_val = 0.5, 
                                     min_sample_size = 5, verbose = 5)
  return(gene_filtered_HWE)
}

# Filter HLA-DPB1
HLADPB1_filtered_HWE <- HWE_filter(gl_HLADPB1)

# Filter HLA-DRB1
HLADRB1_filtered_HWE <- HWE_filter(gl_HLADRB1)

# Filter HLA-DQA1
HLADQA1_filtered_HWE <- HWE_filter(gl_HLADQA1)

# Define Helper Functions for PCA ----
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

# PRINCIPAL COMPONENTS ANALYSIS ----

# Function to run analysis
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

# Plot results
# Scree plot for scaled PCA with relative importance
screeplot(results_DPB1$pca, 
          npcs = 50,
          ylab = "Relative Importance", 
          main = "Scree plot - DPB1 (Relative Importance)")

# Bar plot for PCA with percent variance explained
barplot(results_DPB1$var,
        main = "Percent variation Scree plot - HLA-DPB1",
        ylab = "Percent variation explained", 
        cex.names = 0.96)
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

# One outlier - need to figure out which one it is
outlier_id <- rownames(results_DPB1$scores)[which.min(results_DPB1$scores$PC1)]

# Generate clean df_pop_info with NA18915 sample removed
metadata_clean <- metadata %>%
  filter(id != outlier_id)

sum(metadata_clean$id== outlier_id) #is 0 which is good

# Remove NA18915 from post-HW filter
gl_clean_HLADPB1 <- HLADPB1_filtered_HWE[!indNames(HLADPB1_filtered_HWE) %in% outlier_id]

# Re-run function again
results_DPB1_2 <- run_pca(gl_clean_HLADPB1, "HLA-DPB1", metadata, pops_of_interest)

# Scree plot without the outlier
screeplot(results_DPB1_2$pca,
          npcs = 50,
          ylab = "Relative Importance", 
          main = "Scree plot - A. DPB1 (Relative Importance)")

# Bar plot for PCA with percent variance explained
barplot(results_DPB1_2$var,
        main = "A. HLA-DPB1",
        ylab = "Percent variation explained",
        cex.names = 0.96)
abline(h = 1/results_DPB1_2$n_snps * 100, col = 2, lwd = 2)

# Plot with new data set without the outlier
ggpubr::ggscatter(data = results_DPB1_2$scores,
                  y = "PC2",
                  x = "PC1",
                  color = "pop_id",
                  shape = "pop_id",
                  xlab = paste0("PC1", " [", results_DPB1_2$var[1], "%", "]"),
                  ylab = paste0("PC2", " [", results_DPB1_2$var[2], "%", "]"),
                  main = "A. HLA-DPB1",
                  legend.title = "Population")

## HLA-DRB1 ----
results_DRB1 <- run_pca(HLADRB1_filtered_HWE, "HLA-DRB1", metadata, pops_of_interest)

# Plot results
# Scree plot for scaled PCA with relative importance
screeplot(results_DRB1$pca,
          ylab = "Relative Importance", 
          main = "Scree plot - B. DRB1 (Relative Importance)")

# Bar plot for PCA with percent variance explained
barplot(results_DRB1$var,
        main = "B. HLA-DRB1",
        ylab = "Percent variation explained", 
        cex.names = 0.96)
abline(h = 1/results_DRB1$n_snps * 100, col = 2, lwd = 2)

# PCA plot 
ggpubr::ggscatter(data = results_DRB1$scores,
                  y = "PC2",
                  x = "PC1",
                  color = "pop_id",
                  shape = "pop_id",
                  xlab = paste0("PC1", " [", results_DRB1$var[1], "%", "]"),
                  ylab = paste0("PC2", " [", results_DRB1$var[2], "%", "]"),
                  main = "B. HLA-DRB1",
                  legend.title = "Population")

## HLA-DQA1 ----
results_DQA1 <- run_pca(HLADQA1_filtered_HWE, "HLA-DQA1", metadata, pops_of_interest)

#Plot results
# Scree plot for scaled PCA with relative importance
screeplot(results_DQA1$pca, 
          ylab = "Relative Importance", 
          main = "Scree plot - DQA1 (Relative Importance)",
          cex.names = 0.96)

# Bar plot for PCA with percent variance explained
barplot(results_DQA1$var,
        main = "C. HLA-DQA1",
        ylab = "Percent variation explained", 
        cex.names = 0.96)
abline(h = 1/results_DQA1$n_snps * 100, col = 2, lwd = 2)

# PCA plot 
ggpubr::ggscatter(data = results_DQA1$scores,
                  y = "PC2",
                  x = "PC1",
                  color = "pop_id",
                  shape = "pop_id",
                  xlab = paste0("PC1", " [", results_DQA1$var[1], "%", "]"),
                  ylab = paste0("PC2", " [", results_DQA1$var[2], "%", "]"),
                  main = "C. HLA-DQA1",
                  legend.title = "Population")
# Look at cumulative variation for first 10 PCs of each gene: 
cumsum(results_DPB1_2$var)[10]
cumsum(results_DRB1$var)[10]
cumsum(results_DQA1$var)[10]
 

# NON-HIERARCHICAL CLUSTERING: K-MEANS ----

# Function to find the optimal K values
find_k <- function(gene_pca_results, gene_name) {
  set.seed(321)
  #1. Use Within-Cluster Sums of Squares (WSS) Visualization Method
  k_elbow_plot <- fviz_nbclust(gene_pca_results$scores[, 1:10], kmeans, method = "wss") +
    labs(subtitle = glue("{gene_name}"))
  print(k_elbow_plot)
  
  #2. Use Nbclust method
  nbout_gene <- NbClust(gene_pca_results$scores[, 1:10], method = "kmeans")
  k_barplot <- barplot(table(nbout_gene$Best.nc[1,]), ylab = "Indexes", 
                       xlab = "Number of clusters (k)",
                       main = glue("B. {gene_name}"))
  
  #3. Use Gap Statistics
  gapstat_gene <- clusGap(gene_pca_results$scores[, 1:10], 
                          FUN = kmeans, nstart = 25, d.power = 2, K.max = 10, B = 100)
  k_gapstat_plot <- fviz_gap_stat(gapstat_gene)
  print(k_gapstat_plot)
  
  # Extract Best K from NbClust
  best_k <- as.integer(names(which.max(table(nbout_gene$Best.nc[1,]))))
  
  return(data.frame(gene = gene_name, best_k_nbclust = best_k))
}

# Best k for HLA-DPB1
best_k_DPB1 <- find_k(results_DPB1_2, "HLA-DPB1")

# Best k for HLA-DRB1
best_k_DRB1 <- find_k(results_DRB1, "HLA-DRB1")

# Best k for HLA-DQA1
best_k_DQA1 <- find_k(results_DQA1, "HLA-DQA1")

# Combine best k values into a table 
best_k_table <- rbind(best_k_DPB1, best_k_DRB1, best_k_DQA1)

# Perform Clustering Analysis
## First 10 PCs ----
### HLA-DPB1 ----
kmeans_results_DPB1_10PC <- kmeans(results_DPB1_2$scores[, 1:10], best_k_DPB1$best_k_nbclust, nstart = 10)

# Plot results
fviz_cluster(kmeans_results_DPB1_10PC, data = results_DPB1_2$scores[, 1:10])

### HLA-DRB1 ----
kmeans_results_DRB1_10_PC <- kmeans(results_DRB1$scores[, 1:10], best_k_DRB1$best_k_nbclust, nstart = 10)

# Plot results
fviz_cluster(kmeans_results_DRB1_10_PC, data = results_DRB1$scores[, 1:10])

### HLA-DQA1 ----
kmeans_results_DQA1_10PC <- kmeans(results_DQA1$scores[, 1:10], best_k_DQA1$best_k_nbclust, nstart = 10)

# Plot results
fviz_cluster(kmeans_results_DQA1_10PC, data = results_DQA1$scores[, 1:10])

## All Scaled SNPs ----
### HLA-DPB1 ----
kmeans_results_DPB1 <- kmeans(results_DPB1_2$scaled_snps, best_k_DPB1$best_k_nbclust, nstart = 10)

# Plot results
clust_plot_DPB1 <- fviz_cluster(kmeans_results_DPB1, data = results_DPB1_2$scaled_snps, geom = "point", main = "A. HLA-DPB1")

### HLA-DRB1 ----
kmeans_results_DRB1 <- kmeans(results_DRB1$scaled_snps, best_k_DRB1$best_k_nbclust, nstart = 10)

# Plot results
clust_plot_DRB1 <- fviz_cluster(kmeans_results_DRB1, data = results_DRB1$scaled_snps, geom = "point", main = "A. HLA-DRB1")

### HLA-DQA1 ----
kmeans_results_DQA1 <- kmeans(results_DQA1$scaled_snps, best_k_DQA1$best_k_nbclust, nstart = 10)

# Plot results
clust_plot_DQA1 <- fviz_cluster(kmeans_results_DQA1, data = results_DQA1$scaled_snps, geom = "point", main = "A. HLA-DQA1")



# EDIT!!! MAKE FIGURES LOOK PRETTY ####
# EDIT!!! WE COULD ALSO TRY TO PLOT THE POPULATION GROUPS ON TO THESE CLUSTERING FIGURES BUT IDK HOW ####

# CLASSIFICATION: K NEAREST NEIGHBOURS ----
set.seed(1234)

# Function to convert genlights into data frames
genlight_to_df <- function(hwe_filtered_df) {
  gene_df <- data.frame(population = as.factor(pop(hwe_filtered_df)), 
                        as.matrix(hwe_filtered_df), check.names = FALSE)
  return(gene_df)
}

# Apply function
HLADRB1_df <- genlight_to_df(HLADRB1_filtered_HWE)
HLADPB1_df <- genlight_to_df(HLADPB1_filtered_HWE) # EDIT!!! Do we want to use the one without the outlier?
HLADQA1_df <- genlight_to_df(HLADQA1_filtered_HWE)

# Remove rows with abnormal population code
HLADRB1_df <- HLADRB1_df %>% 
  filter(!grepl(",", population)) %>%
  filter(population %in% pops_of_interest) %>%
  droplevels()

HLADPB1_df <-  HLADPB1_df %>% 
  filter(!grepl(",", population)) %>%
  filter(population %in% pops_of_interest) %>%
  droplevels()

HLADQA1_df <- HLADQA1_df %>% 
  filter(!grepl(",", population)) %>%
  filter(population %in% pops_of_interest) %>%
  droplevels()

# Split into test and training data
data_split <- function(gene_df){
  gene_train_index <- createDataPartition(gene_df$population, p = 0.75, list = FALSE)
  train_data <- gene_df[gene_train_index, ]
  test_data  <- gene_df[-gene_train_index, ]
  return(list(index = gene_train_index, training = train_data, test = test_data))
}

# Apply function
HLADRB1_list <- data_split(HLADRB1_df)
HLADRB1_train <- HLADRB1_list$training
HLADRB1_test  <- HLADRB1_list$test

HLADPB1_list <- data_split(HLADPB1_df)
HLADPB1_train <- HLADPB1_list$training
HLADPB1_test  <- HLADPB1_list$test

HLADQA1_list <- data_split(HLADQA1_df)
HLADQA1_train <- HLADQA1_list$training
HLADQA1_test  <- HLADQA1_list$test

# CV to choose best k

set.seed(123)

# Function to select best k for each data

cv_choosek <- function(train_data){
  control <- trainControl(method = "cv", number = 10)
  knn_model <- train(population ~ ., data = train_data, method = "knn", trControl = control, tuneLength = 20, preProcess = c("center", "scale"))
  return(knn_model)
}

# Apply function
HLADRB1_bestk <- cv_choosek(HLADRB1_train)
plot(HLADRB1_bestk)
HLADRB1_bestk_results <- HLADRB1_bestk$bestTune$k

HLADPB1_bestk <- cv_choosek(HLADPB1_train)
plot(HLADPB1_bestk)
HLADPB1_bestk_results <- HLADPB1_bestk$bestTune$k

HLADQA1_bestk <- cv_choosek(HLADQA1_train)
plot(HLADQA1_bestk)
HLADQA1_bestk_results <- HLADQA1_bestk$bestTune$k 

# knn function

knn_function <- function(train_data, test_data, bestk){
  model_train <- train_data[,-1]
  model_test <- test_data[,-1]
  k_neighbours <- class::knn(train = model_train, test = model_test, cl = train_data$population, k = bestk)
  return(k_neighbours)
}

# Apply function
HLADRB1_knnmod <- knn_function(HLADRB1_train, HLADRB1_test, HLADRB1_bestk_results)
HLADPB1_knnmod <- knn_function(HLADPB1_train, HLADPB1_test, HLADPB1_bestk_results)
HLADQA1_knnmod <- knn_function(HLADQA1_train, HLADQA1_test, HLADQA1_bestk_results)

# Function to calculate classification performance metrics using confusion matrix
classification_performance <- function(pred, actual) {
  
  # Confusion matrix
  c_matrix <- confusionMatrix(pred, actual)
  
  # Simple table
  con_matrix <- table(Predicted = pred, Actual = actual)
  print(con_matrix)
  
  # Extract metrics
  accuracy <- round(c_matrix$overall['Accuracy'],3)
  sensitivity <- round(c_matrix$byClass[, 'Sensitivity'],3)
  specificity <- round(c_matrix$byClass[, 'Specificity'],3)

  print(accuracy)
  print(sensitivity)
  print(specificity)
  
  # Return as a list (useful later if you want a table)
  return(list(
    confusion_matrix = con_matrix,
    accuracy = accuracy,
    sensitivity = sensitivity,
    specificity = specificity
  ))
}

# Apply function
metrics_DRB1 <- classification_performance(HLADRB1_knnmod, HLADRB1_test$population)
metrics_DPB1 <- classification_performance(HLADPB1_knnmod, HLADPB1_test$population)
metrics_DQA1 <- classification_performance(HLADQA1_knnmod, HLADQA1_test$population)
