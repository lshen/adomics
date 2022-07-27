args = commandArgs(TRUE)

ind = as.numeric(args)

# simulations = read.csv("simulation.csv", as.is=TRUE)
# simulations = read.csv("simulation_beta_dist.csv", as.is=TRUE)  # varying SNP effect sizes
simulations = read.csv("simulation_real_genotypes.csv", as.is=TRUE)  # varying SNP effect sizes
simulation_setup = simulations[ind, ]

outcome = simulation_setup$outcome
sample_size = 500
IV_strength = simulation_setup$IV_strength
pleiotropy = simulation_setup$pleiotropy
effect = simulation_setup$effect
overlap = simulation_setup$overlap
r2 = simulation_setup$r2  # not used
xQTL_from_one_sample = simulation_setup$xQTL_from_one_sample
number_of_SNPs = simulation_setup$number_of_SNPs   # not used

stopifnot(sample_size == 500)
stopifnot(is.na(r2))
stopifnot(xQTL_from_one_sample)
stopifnot(is.na(number_of_SNPs))

ncores = 7

# 1000 Genome genotype data downloaded from:
# http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/integrated_call_sets/
# ALL.chr19.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz
#
# Extract APOE gene and pad the genomic region with 1m bp on each side.
# https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000130203;r=19:45409011-45412650
# 
################################################################################
# Exposure and outcome simulated from a dosage matrix
# with 1092 individuals (cols) and 6157 SNPs (rows).
################################################################################

library(ieugwasr)
dosage_matrix_path =  "../data/APOE_region_dosage_matrix.rds"
stopifnot(file.exists(dosage_matrix_path))
dosage_matrix = readRDS(dosage_matrix_path)

EUR_path = "../data/EUR"
if (!file.exists(paste0(EUR_path, ".bed"))) {
  system("cd ../data/ && wget http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz &&  xzvf 1kg.v3.tgz EUR.bed EUR.bim EUR.fam")
}

plink_path = "/appl/plink-1.9-20210416/plink"
is_plink_available = file.exists(plink_path)
if (!is_plink_available) {
  plink_path = EUR_path = NULL
}

library(harmonicmeanp)
source("SMR_functions.R")  # get_p_values_from_summary()

library(matrixStats)

# Table format: see table IV in:
# "Bias due to participant overlap in
# two-sample Mendelian randomization" by Burgess et al. 2016.

# Columns:
# alpha (non-zero corresponding to IV strength), casual effect of SNPs on exposure
#       signal to noise ratio measured in:
#       (mean F (and mean R^2); 3 choices around 5, 10, 15)
#       Mean F   for each exposure
#       Mean R^2 for each exposure
# Power
# (the followings are not applicable -- Median estimate (Median standard errors) [Power])
#
# Rows:
# 
# horizontal vs. vertical pleiotropy (2 choices)
# positive vs. null casual effect (\beta_X)
# sample overlap (percent overlap 0%, 20%, 40%, 60%, 80%)
# linkage disequilibrium: r2 = 0, 0.01, 0.1
# number of SNPs used as instrumental variables (5, 10)
#
# "Bias due to participant overlap in
# two-sample Mendelian randomization" by Burgess et al. 2016:
# > The ordinary least squares (OLS, also known as standard least
# > squares regression) estimate is obtained by regressing the outcome
# > on the risk factor. This “observational” estimate is typically
# > biased due to confounding. The relative bias of the 2SLS
# > estimate —the bias of the 2SLS estimate divided by the bias
# > of the OLS estimate—is approximately and asymptotically
# > equal to 1∕E(F) (Staiger & Stock, 1997), where E(F) is the
# > “F parameter.” An F parameter of 10 therefore corresponds
# > to a 1∕10 = 10% relative bias of the 2SLS estimate compared
# > to the OLS estimate.

library(doMC)
library(doRNG)
date()
# 2. Generate x1, x2, x3, y, z1, z2 100 times with different ways to generate y
times = 10000
set.seed(1234)
stopifnot(overlap == 0 | xQTL_from_one_sample)
method_names = c(
  # multiple SNPs
  "Multivariable", "Multivariable_tdist",
  "MultivariableLD", "MultivariableLD_tdist",
  "IVW_Cauchy_cauchy", "IVW_tdist_Cauchy_cauchy",
  "IVW_HMP", "IVW_tdist_HMP", 
  "IVW_Modality1", "IVW_tdist_Modality1", 
  "IVW_Modality2", "IVW_tdist_Modality2", 
  "IVW_Modality3", "IVW_tdist_Modality3",
  "IVW_MinP", "IVW_tdist_MinP", 
  "IVW_Fisher_chisq", "IVW_tdist_Fisher_chisq", 
  "WGLR_Cauchy_cauchy","WGLR_tdist_Cauchy_cauchy",
  "WGLR_HMP", "WGLR_tdist_HMP",
  "WGLR_Modality1", "WGLR_tdist_Modality1",
  "WGLR_Modality2", "WGLR_tdist_Modality2",
  "WGLR_Modality3", "WGLR_tdist_Modality3",
  "WGLR_MinP", "WGLR_tdist_MinP",
  "WGLR_Fisher_chisq", "WGLR_tdist_Fisher_chisq",
  "WGLR_Fisher_gamma", "WGLR_tdist_Fisher_gamma",
  "GSMR_Cauchy_cauchy", "GSMR_HMP", 
  "GSMR_Modality1", "GSMR_Modality2", "GSMR_Modality3",
  "GSMR_MinP", "GSMR_Fisher_chisq", 
  "SMR_allSNPs_Cauchy_cauchy", "SMR_allSNPs_HMP", 
  "SMR_allSNPs_MinP", "SMR_allSNPs_Fisher_chisq", 
  # single SNP
  "SMR_singleSNP_Cauchy_cauchy", "SMR_singleSNP_HMP", 
  "SMR_Modality1", "SMR_Modality2", "SMR_Modality3",
  "SMR_singleSNP_Fisher_chisq",
  "Multivariable_beta1",        
  "Multivariable_beta2",         "Multivariable_beta3",        
  "MultivariableLD_beta1",       "MultivariableLD_beta2",      
  "MultivariableLD_beta3",       "IVW_beta1",                  
  "IVW_beta2",                   "IVW_beta3",                  
  "WGLR_beta1",                  "WGLR_beta2",                 
  "WGLR_beta3"                
)
statistics_names = c("OLS_beta1", "OLS_beta2", "OLS_beta3",
                     "F_1", "R2_1", "F_2", "R2_2", "F_3", "R2_3", "nsnps")
# total number of sample (case + control)
# we will simulate 2*n samples for a two-sample MR
# For a case-control study we need to simulate a larger number
# and then take a subset.
n = sample_size

if (outcome %in% c("binary_use_control_in_xQTL", "binary_use_all_in_xQTL")) {
  stop("We are simulating using real genotypes with a limited sample size. Binary outcome is not supported.") 
}
registerDoMC(ncores)
RNGkind("L'Ecuyer-CMRG")
set.seed(1234)

# combined_array = matrix(ncol=length(c(method_names, statistics_names)), nrow=times)
# colnames(combined_array) = c(method_names, statistics_names)
  
combined_array = foreach(i = 1:times, .combine = rbind,
                         .packages = c("ieugwasr", "harmonicmeanp", "matrixStats")) %dorng% {
# for (i in 1:times) {
  if (i %% 100 == 0) cat(sprintf("i=%d, %s\n", i, date()))
  # sample 1,000 subjects from 1,092 subjects
  z_matrix = t(dosage_matrix[, sample(ncol(dosage_matrix), 2*n, replace=FALSE)])
  # z_matrix[] = sample(as.numeric(z_matrix)); 
  # z_matrix = apply(z_matrix, 2, scale)  # scaling does not work
  u = rnorm(2*n)

  # determine the number of the SNPs
  p = ncol(z_matrix)
  
  # SNP-exposure effect
  generate_SNP_effects = function(zmat) {
    number_of_SNPs_in_laplace_dist = 100
    zmat[] = t(t(zmat) * sample(c(
      rep(0, ncol(zmat) - number_of_SNPs_in_laplace_dist),
      rexp(number_of_SNPs_in_laplace_dist, rate = 1) *
        sample(c(-1,1), number_of_SNPs_in_laplace_dist, replace=TRUE)
      )))
    rowSums(zmat)
  }

  z = generate_SNP_effects(z_matrix)
  if (pleiotropy == "horizontal") {
    z_first_half = generate_SNP_effects(z_matrix[,1:(floor(p/2))])
    z_second_half = generate_SNP_effects(z_matrix[,(floor(p/2)+1):p])
    x1 = 1/10 * IV_strength * z + u + rnorm(2*n) 
    x2 = - 1/10 * IV_strength * z_first_half - u + rnorm(2*n)
    x3 = 1/10 * IV_strength * z_second_half - u + rnorm(2*n)
    v = effect * (x1 + x2 + x3) + u
    odds = exp(-2 + v)
  } else if (pleiotropy == "vertical") {
    x1 = 1/10 * IV_strength * z + u + rnorm(2*n)
    x2 = 2 * x1 - u + rnorm(2*n)
    x3 = - x2 + u + rnorm(2*n)
    v = effect * x3 - u + rnorm(2*n)
    odds = exp(-2 + v)
  } else stop("Unrecognized pleiotropy")
  

  if (xQTL_from_one_sample) {
    if (outcome == "continuous") {
      y = v + rnorm(2*n)
      samples_in_GWAS = 1:n
      samples_in_xQTL = c(
        seq_len(round(n*overlap)),
        n + seq_len(n-round(n*overlap)))
    }
    samples_in_xQTL1 = samples_in_xQTL2 = samples_in_xQTL3 = samples_in_xQTL
  } else {
    if (outcome == "continuous") {
      y = v + rnorm(2*n)
      samples_in_GWAS = 1:n
      samples_in_xQTL1 = n + seq_len(n)
      samples_in_xQTL2 = n*2 + seq_len(n)
      samples_in_xQTL3 = n*3 + seq_len(n)
    }
  }
  
  x1 = x1[samples_in_xQTL1]
  x2 = x2[samples_in_xQTL2]
  x3 = x3[samples_in_xQTL3]
  y0 = y[samples_in_GWAS]
  z1 = z_matrix[samples_in_xQTL1, ]
  z2 = z_matrix[samples_in_xQTL2, ]
  z3 = z_matrix[samples_in_xQTL3, ]
  z0 = z_matrix[samples_in_GWAS, ]
  
  zx_est = zx_se = matrix(NA, nrow=3, ncol=p)
  zy_est = zy_se = rep(NA, p)
  for (k in 1:p) {
    x1_lm_summary = summary(lm(x1~z1[,k]))
    x2_lm_summary = summary(lm(x2~z2[,k]))
    x3_lm_summary = summary(lm(x3~z3[,k]))
    
    x1_lm_coefficients = x1_lm_summary$coefficients
    x2_lm_coefficients = x2_lm_summary$coefficients
    x3_lm_coefficients = x3_lm_summary$coefficients
    
    zx_est[1, k] = x1_lm_coefficients[2, "Estimate"]
    zx_se[1, k] = x1_lm_coefficients[2, "Std. Error"]
    zx_est[2, k] = x2_lm_coefficients[2, "Estimate"]
    zx_se[2, k] = x2_lm_coefficients[2, "Std. Error"]
    zx_est[3, k] = x3_lm_coefficients[2, "Estimate"]
    zx_se[3, k] = x3_lm_coefficients[2, "Std. Error"]
    
    if (outcome == "continuous") {
      y0_lm_coefficients = summary(lm(y0~z0[,k]))$coefficients
      
      zy_est[k] = y0_lm_coefficients[2, "Estimate"]
      zy_se[k] = y0_lm_coefficients[2, "Std. Error"]
    }
  }
  
  # select instrumental variables using ieugwasr
  # we have plink_path
  # get r2 matrix using ieugwasr::ld_matrix() and local references
  min_p = colMins(pnorm(zx_est/zx_se))
  
  # local version
  # http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz
  ld_clump_result = ieugwasr::ld_clump(dat = data.frame(rsid = colnames(z_matrix), pval = min_p),
                             clump_kb = 10000,
                             clump_r2 = 1e-4,
                             clump_p = max(5e-3, min(min_p)),
                             bfile = EUR_path,
                             plink_bin = plink_path)
  ld_matrix_signed = ieugwasr::ld_matrix(ld_clump_result$rsid,
                               with_allele = FALSE,
                               pop = "EUR",
                               bfile = EUR_path,
                               plink_bin = plink_path)
  selected_SNPs = ld_clump_result$rsid
  which_selected_SNPs = match(ld_clump_result$rsid, colnames(z_matrix))
  
  # assess instrument strength
  x1_lm_multi_summary = summary(lm(x1~z1[,which_selected_SNPs]))
  x2_lm_multi_summary = summary(lm(x2~z2[,which_selected_SNPs]))
  x3_lm_multi_summary = summary(lm(x3~z3[,which_selected_SNPs]))
  
  # OLS
  x1_OLS_lm = lm(y[samples_in_xQTL1]~x1)
  x2_OLS_lm = lm(y[samples_in_xQTL2]~x2)
  x3_OLS_lm = lm(y[samples_in_xQTL3]~x3)
  
  # F and R2 statistics
  statistics = c(
    coef(x1_OLS_lm)[2],
    coef(x2_OLS_lm)[2],
    coef(x3_OLS_lm)[2],
    x1_lm_multi_summary$fstatistic[1],
    x1_lm_multi_summary$r.squared,
    x2_lm_multi_summary$fstatistic[1],
    x2_lm_multi_summary$r.squared,
    x3_lm_multi_summary$fstatistic[1],
    x3_lm_multi_summary$r.squared,
    nrow(ld_clump_result)
    )
  names(statistics) = c("OLS_beta1", "OLS_beta2", "OLS_beta3",
                        "F_1", "R2_1", "F_2", "R2_2", "F_3", "R2_3", "nsnps")
  
  pval_array_i = get_p_values_from_summary(zx_est[,which_selected_SNPs, drop=FALSE], 
                                           zx_se[,which_selected_SNPs, drop=FALSE],
                                           zy_est[which_selected_SNPs], 
                                           zy_se[which_selected_SNPs],
                                           r2 = ld_matrix_signed,
                                           include_coef_estimates = TRUE,
                                           include_number_of_SNPs = FALSE,
                                           need_ld_clump = FALSE)
  
  c(pval_array_i, statistics)
  # combined_array[i, ] = c(pval_array_i, statistics)
}

pval_array = combined_array[,!grepl("_beta|^F_|^R2_|nsnps", colnames(combined_array))]
statistics_array = combined_array[,grepl("_beta|^F_|^R2_|nsnps", colnames(combined_array))]

powers = apply(pval_array, 2, function(x) sum(x < 0.05) / times)

means = apply(statistics_array, 2, mean)
names(means) = paste0(colnames(statistics_array), "_mean")

sds = apply(statistics_array, 2, sd)
names(sds) = paste0(colnames(statistics_array), "_sd")

simulation_details = c(simulation_setup, powers, means, sds)

# cat real_genotypes_0*.csv > result_real_genotypes.csv
if (ind == 1) {
  write.table(simulation_details, file = sprintf("real_genotypes_%03d.csv", ind),  
              row.names=FALSE, col.names=TRUE, sep=",")
} else {
  write.table(simulation_details, file = sprintf("real_genotypes_%03d.csv", ind),  
              row.names=FALSE, col.names=FALSE, sep=",")
}

write.table(combined_array, file = sprintf("real_genotypes_combined_array_%03d.csv", ind),
            row.names=FALSE, col.names=TRUE, sep=",")