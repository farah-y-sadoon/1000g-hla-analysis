# Install the packages

#BiocManager::install("dartR")
library(dartR)

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