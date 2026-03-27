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

# Load VCF
snps <- read.vcfR("HLADPB1.vcf.gz")
snps_HLADRB1 <- read.vcfR("HLADQA1.vcf.gz")
snps_HLADQA1 <- read.vcfR("HLADRB1.vcf.gz")
# Load population metadata
pop_info <- read_tsv("igsr_samples.tsv")
df <- as.data.frame(pop_info)

# Extract genotype strings (0/0, 0/1, etc.)
snps_num <- vcfR::extract.gt(snps, 
                             element = "GT",
                             IDtoRowNames  = F,
                             as.numeric = T,
                             convertNA = T,
                             return.alleles = F)

# Transpose → rows = individuals
snps_num_t <- as.data.frame(t(snps_num))


# Remove NAs
find_NAs <- function(x){
  NAs_TF <- is.na(x)
  i_NA <- which(NAs_TF == TRUE)
  N_NA <- length(i_NA)
  
  cat("Results:",N_NA, "NAs present\n.")
  return(i_NA)
}

# N_rows
# number of rows (individuals)
N_rows <- nrow(snps_num_t) #2504

# N_NA
# vector to hold output (number of NAs)
N_NA   <- rep(x = 0, times = N_rows)

# N_SNPs
# total number of columns (SNPs)
N_SNPs <- ncol(snps_num_t) # 582

# the for() loop
for(i in 1:N_rows){
  
  # for each row, find the location of
  ## NAs with snps_num_t()
  i_NA <- find_NAs(snps_num_t[i,]) 
  
  # then determine how many NAs
  ## with length()
  N_NA_i <- length(i_NA)
  
  # then save the output to 
  ## our storage vector
  N_NA[i] <- N_NA_i
}

# O NA 

# Remove invariant columns
invar_omit <- function(x){
  cat("Dataframe of dim",dim(x), "processed...\n")
  sds <- apply(x, 2, sd, na.rm = TRUE)
  i_var0 <- which(sds == 0)
  
  
  cat(length(i_var0),"columns removed\n")
  
  if(length(i_var0) > 0){
    x <- x[, -i_var0]
  }
  
  ## add return()  with x in it
  return(x)                      
}

snps_no_invar <- invar_omit(snps_num_t)

# Scale snps
snps_scaled <- scale(snps_no_invar)

# Run PCA
pca_scaled <- prcomp(snps_scaled)

# Generate Scree Plot
screeplot(pca_scaled, 
          ylab  = "Relative importance",
          main = "Scree plot of PC for SNP data")

summary <- summary(pca_scaled)

PCA_variation <- function(pca_summary, PCs = 2){
  var_explained <- pca_summary$importance[2,1:PCs]*100
  var_explained <- round(var_explained,1)
  return(var_explained)
}

var_out <- PCA_variation(summary,PCs = 10)
var_out

N_columns <- ncol(snps_scaled)

barplot(var_out,
        main = "Percent variation Scree plot",
        ylab = "Percent variation explained")
abline(h = 1/N_columns*100, col = 2, lwd = 2)

biplot(pca_scaled)

pca_scores <- vegan::scores(pca_scaled)

# Attach ID columns to merge by ID 
snps_no_invar$ID <- rownames(snps_no_invar)
df$ID <- df$`Sample name`

# Merge genotype data with population metadata using Sample ID
merged <- merge(snps_no_invar, df, by = "ID")

pop_id <- merged$`Population code`

pca_scores2 <- data.frame(pop_id,
                          pca_scores)

ggpubr::ggscatter(data = pca_scores2,
                  y = "PC2",
                  x = "PC1",
                  color = "pop_id",
                  shape = "pop_id",
                  xlab = "PC1 (14.4% variation)",
                  ylab = "PC2 (13.0% variation)",
                  main = "PCA results from SNPs with real population identity")
