library(TwoSampleMR)
library(R.utils)
library(dplyr)
library(jsonlite)
library(httr)
library(MRPRESSO)
library(RadialMR)
library(biomaRt)
exposure_dat <- read_exposure_data(
  filename = "GCST010144.h.tsv",
  sep = "\t",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "standard_error",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "effect_allele_frequency",
  pval_col = "p_value"
)

header <- read.table("pcos222.tsv", sep = "\t", nrows = 1, header = TRUE)
colnames(header)
head(header)
readLines("pcos222.tsv", n = 15)

outcome_dat <- read_outcome_data(
  filename = "pcos222.tsv",
  sep = "\t",
  snp_col = "rs_id",
  beta_col = "beta",
  se_col = "standard_error",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "effect_allele_frequency",
  pval_col = "p_value"
)

harmonised_dat <- harmonise_data(
  exposure_dat = exposure_dat,
  outcome_dat = outcome_dat,
)
results = mr(harmonised_dat) # added afterwards to see if this code works
exposure_dat <- exposure_dat[exposure_dat$pval.exposure < 5e-8, ]
radial_data <- format_radial(
  BXG = harmonised_dat$beta.exposure,
  BYG = harmonised_dat$beta.outcome,
  seBXG = harmonised_dat$se.exposure,
  seBYG = harmonised_dat$se.outcome,
)

results_radial <- ivw_radial(radial_data, alpha = 0.05)
str(results_radial$outliers)
outlier_indices <- results_radial$outliers$SNP
outlier_snps <- harmonised_dat$SNP[outlier_indices]
cleaned_data <- harmonised_dat[!harmonised_dat$SNP %in% outlier_snps, ]

results = mr(cleaned_data)
mr_scatter_plot(results, cleaned_data)
print(results)
print(mr_pleiotropy_test(cleaned_data))
print(mr_heterogeneity(cleaned_data))
write.table(results, "PCOSMR_resultsfinal.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
nrow(read.table("pcos222.tsv", header=TRUE, sep="\t"))

table(cleaned_data$mr_keep)
