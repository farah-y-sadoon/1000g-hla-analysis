# Install the packages

#BiocManager::install("dartR")
library(dartR)
library(tidyr)
library(readr)

# convert metadata to .csv, must contain sample names as column "id" and population as "pop"
tsv <- read_tsv("igsr_samples.tsv")
names(tsv)
metadata <- tsv %>%
  rename(id = `Sample name`, pop = `Population code`) %>%
  select(id, pop)

# metadata write to csv
write.csv(metadata, "metadata.csv")

# read in vcfs

gl_HLADPB1 <- gl.read.vcf(vcffile = "HLADPB1.vcf.gz", ind.metafile = "metadata.csv")

## repeat for other genes
gl_HLADRB1 <- gl.read.vcf(vcffile = "HLADRB1.vcf.gz", ind.metafile = "metadata.csv")
gl_HLADQA1 <- gl.read.vcf(vcffile = "HLADQA1.vcf.gz", ind.metafile = "metadata.csv")

#HWE filtering
library(HardyWeinberg)

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
    scores = pca_scores, 
    var = var_explained, 
    gene = gene_name, 
    n_snps = ncol(snps_scaled)
  ))
}

# pops of interest
pops_of_interest <- c("YRI", "JPT", "FIN", "CHB")

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
                  xlab = "PC1",
                  ylab = "PC2",
                  main = "SNPs PCA Scatterplot")

### Remove outlier ----
outlier <- results_DPB1$scores[which.min(results_DPB1$scores$PC1), ]
outlier

outlier_id <- rownames(results_DPB1$scores)[which.min(results_DPB1$scores$PC1)]
outlier_id

# Remove NA18915 from post-HW filter
gl_clean_HLADPB1 <- HLADPB1_filtered_HWE[!indNames(HLADPB1_filtered_HWE) %in% outlier_id]

## HLA-DPB1 no outlier ----
results_DPB1_clean <- run_pca(gl_clean_HLADPB1, "HLA-DPB1", metadata, pops_of_interest)

### Plot results ----
# Scree plot for scaled PCA with relative importance
screeplot(results_DPB1_clean$pca, 
          ylab = "Relative Importance", 
          main = "Scree plot - A. DPB1 (Relative Importance)")

# Bar plot for PCA with percent variance explained
barplot(results_DPB1_clean$var,
        main = "Percent variation Scree plot - HLA-DPB1",
        ylab = "Percent variation explained")
abline(h = 1/results_DPB1_clean$n_snps * 100, col = 2, lwd = 2) # Eigenvalue above 1 (this line) captures meaningful data

# PCA plot 
ggpubr::ggscatter(data = results_DPB1_clean$scores,
                  y = "PC2",
                  x = "PC1",
                  color = "pop_id",
                  shape = "pop_id",
                  xlab = "PC1",
                  ylab = "PC2",
                  main = "SNPs PCA Scatterplot - HLA-DPB1")

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
                  xlab = "PC1",
                  ylab = "PC2",
                  main = "SNPs PCA Scatterplot - HLA-DRB1")

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
                  xlab = "PC1",
                  ylab = "PC2",
                  main = "SNPs PCA Scatterplot - HLA-DQA1")




### Classification

# Load libraries
library(mlbench)
library(caret)
library(class)
library(e1071)

set.seed(1234)

# Function to convert genlights into data frames
genlight_to_df <- function(hwe_filtered_df) {
  gene_df <- data.frame(population = as.factor(pop(hwe_filtered_df)), 
          as.matrix(hwe_filtered_df), check.names = FALSE)
  return(gene_df)
}

# Apply function
HLADRB1_df <- genlight_to_df(HLADRB1_filtered_HWE)
HLADPB1_df <- genlight_to_df(HLADPB1_filtered_HWE)
HLADQA1_df <- genlight_to_df(HLADQA1_filtered_HWE)
  
  
# Remove rows with abnormal population code
  HLADRB1_df <- HLADRB1_df %>% 
    filter(!grepl(",", population)) %>%
    droplevels()
  
  HLADPB1_df <-  HLADPB1_df %>% 
    filter(!grepl(",", population)) %>%
    droplevels()
  
  HLADQA1_df <- HLADQA1_df %>% 
    filter(!grepl(",", population)) %>%
    droplevels()
  

# split inot test and rtraining daat
  
data_split <- function(gene_df){
  gene_train_index <- createDataPartition(gene_df$population, p = 0.75, list = FALSE)
  train_data <- gene_df[gene_train_index, ]
  test_data  <- gene_df[-gene_train_index, ]
  return(list(index = gene_train_index, training = train_data, test = test_data))
}

# apply function
HLADRB1_list <- data_split(HLADRB1_df)
HLADRB1_train <- HLADRB1_list$training
HLADRB1_test  <- HLADRB1_list$test

# Apply to other 2 genes...

# CV to choose best k

set.seed(123)

# function to select best k for each data

cv_choosek <- function(train_data){
  control <- trainControl(method = "cv", number = 10)
  knn_model <- train(population ~ ., data = train_data, method = "knn", trControl = control, tuneLength = 20, preProcess = c("center", "scale"))
  return(knn_model)
}

# apply
HLADRB1_bestk <- cv_choosek(HLADRB1_train)
plot(HLADRB1_bestk)
HLADRB1_bestk_results <- HLADRB1_bestk$bestTune$k 

# knn function

knn_function <- function(train_data, test_data, bestk){
  model_train <- train_data[,-1]
  model_test <- test_data[,-1]
  k_neighbours <- class::knn(train = model_train, test = model_test, cl = train_data$population, k = bestk)
  return(k_neighbours)
}

# apply
HLADRB1_knnmod <- knn_function(HLADRB1_train, HLADRB1_test, HLADRB1_bestk_results)

#worked

# evaluate model
# convert to factor if needed
#test_data$population <- as.factor(test_data$population)
#k_neighbours <- as.factor(k_neighbours)

#levels(k_neighbours)
#levels(test_data$population)

#class(k_neighbours)
#class(test_data$population)

# had to force labels to be the same so the confusiob matrix would work
#k_neighbours <- factor(k_neighbours, levels = levels(test_data$population))


# Confusion matrix
c_matrix <- confusionMatrix(HLADRB1_knnmod, HLADRB1_test$population)

con_matrix <- table(Predicted = HLADRB1_knnmod, Actual = HLADRB1_test$population)
print(con_matrix)


c_matrix$overall['Accuracy']
c_matrix$byClass['Sensitivity'] # NA why
c_matrix$byClass['Specificity'] # NA why?

# repeat for other genes, write function