args = commandArgs(TRUE)

ind = as.numeric(args)

simulations = read.csv("simulation.csv", as.is=TRUE)
# simulations = read.csv("simulation_beta_dist.csv", as.is=TRUE)  # varying SNP effect sizes
simulation_setup = simulations[ind, ]

outcome = simulation_setup$outcome
sample_size = 500
IV_strength = simulation_setup$IV_strength
pleiotropy = simulation_setup$pleiotropy
effect = simulation_setup$effect
overlap = simulation_setup$overlap
r2 = simulation_setup$r2
xQTL_from_one_sample = simulation_setup$xQTL_from_one_sample
number_of_SNPs = simulation_setup$number_of_SNPs

ncores = 1

library(harmonicmeanp)
source("SMR_functions.R")  # get_p_values_from_summary()

# set desired r2 by setting the covariance matrix in haplosim()
if (FALSE) {
# Investigate the relation of r2 vs. covariance matrix in haplosim() empirically.
# We would like to be able to simulate genotype data with a
# prespecified number of SNPs, MAF, and r2.
  model_filename = "sqrt_r2_model.rds"
  library(snpStats)
  library(hapsim)
  set.seed(1234)
  p = 5
  MAF = 0.3
  n_large = 100000
  correlation_in_normal_dist_vec = seq(0, 0.9, by = 0.1)
  r2_vec = rep(NA, length(correlation_in_normal_dist_vec))
  for (i in seq_along(correlation_in_normal_dist_vec)) {
    haplo_cov = matrix(correlation_in_normal_dist_vec[i], nrow = p, ncol = p)
    diag(haplo_cov) = rep(1, p)
    haplo_freqs = rep(1 - MAF, p)
    names(haplo_freqs) = paste0("SNP", 1:p)
    haplo_data = list(freqs = haplo_freqs, cov = haplo_cov)
    repeat {
      simdata_hap1 = haplosim(n_large*2, haplo_data, summary=FALSE)
      if (sd(simdata_hap1$freqs) > 0) break
    }
    repeat {
      simdata_hap2 = haplosim(n_large*2, haplo_data, summary=FALSE)
      if (sd(simdata_hap2$freqs) > 0) break
    }
    ldmat = ld(new("SnpMatrix", as.matrix(simdata_hap1$data + simdata_hap2$data)), depth=p-1, stats=c("R.squared"))
    r2_vec[i] = median(ldmat[upper.tri(ldmat)])
  }
  sqrt_r2_model = smooth.spline(correlation_in_normal_dist_vec ~ sqrt(r2_vec))
  predict(sqrt_r2_model, sqrt(c(0.01, 0.05, 0.1, 0.2, 0.5)))
  saveRDS(sqrt_r2_model, model_filename)
  #$x
  #[1] 0.1000000 0.2236068 0.3162278 0.4472136 0.7071068
  #$y
  #[1] 0.29631021  0.52706164  0.66044081  0.79026279  0.97147630
}
  
if (abs(r2 - 0) < 1e-10) {
  correlation_in_normal_dist = 0
} else if (abs(r2 - 0.01) < 1e-10) {
  correlation_in_normal_dist = 0.29631021
} else if (abs(r2 - 0.1) < 1e-10) {
  correlation_in_normal_dist = 0.66044081
} else if (abs(r2 - 0.2) < 1e-10) {
  correlation_in_normal_dist = 0.79026279
} else stop("Unrecognized r2")

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
library(hapsim)
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
  "IVW_Fisher_gamma", "IVW_tdist_Fisher_gamma",
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
  "SMR_singleSNP_Fisher_chisq"
)

# total number of sample (case + control)
# we will simulate 2*n samples for a two-sample MR
# For a case-control study we need to simulate a larger number
# and then take a subset.
n = sample_size

if (outcome %in% c("binary_use_control_in_xQTL", "binary_use_all_in_xQTL")) {
  n_max = 10*n  # we need more samples for equal-number of cases and controls
} else if (outcome == "continuous") {
  n_max = n
}
if (!xQTL_from_one_sample) {
  n_max = n_max * 2
}
# simulate a single SNP with MAF = 0.3
MAF = 0.3
# p: number of SNPs as instrumental variables
p = number_of_SNPs
registerDoMC(ncores)
RNGkind("L'Ecuyer-CMRG")
set.seed(1234)
combined_array = foreach(i = 1:times, .combine = rbind) %dorng% {
  if (i %% 100 == 0) cat(sprintf("i=%d, %s\n", i, date()))
  haplo_cov = matrix(correlation_in_normal_dist, nrow = p, ncol = p)
  diag(haplo_cov) = rep(1, p)
  haplo_freqs = rep(1 - MAF, p)
  names(haplo_freqs) = paste0("SNP", 1:p)
  haplo_data = list(freqs = haplo_freqs, cov = haplo_cov)
  repeat {
    simdata_hap1 = haplosim(n_max*2, haplo_data, summary=FALSE)
    if (sd(simdata_hap1$freqs) > 0) break
  }
  repeat {
    simdata_hap2 = haplosim(n_max*2, haplo_data, summary=FALSE)
    if (sd(simdata_hap2$freqs) > 0) break
  }
  
  z_matrix = simdata_hap1$data + simdata_hap2$data
  u = rnorm(2*n_max)

  generate_SNP_effects = function(zmat) {
    # zmat[] = zmat * (0.5 + rbeta(length(zmat), shape1 = 5, shape2 = 5))
    zmat[] = zmat * (2 * rbeta(length(zmat), shape1 = 5, shape2 = 5))
    rowSums(zmat)
  }

  z = generate_SNP_effects(z_matrix)
  if (pleiotropy == "horizontal") {
    z_first_half = generate_SNP_effects(z_matrix[,1:(floor(p/2))])
    z_second_half = generate_SNP_effects(z_matrix[,(floor(p/2)+1):p])
    x1 = 1 / p * IV_strength * z + u + rnorm(2*n_max) 
    x2 = - 1 / p * IV_strength * z_first_half - u + rnorm(2*n_max)
    x3 = 1 / p * IV_strength * z_second_half - u + rnorm(2*n_max)
    v = effect * (x1 + x2 + x3) + u
    odds = exp(-2 + v)
  # } else if (pleiotropy == "horizontal6") {
  #   z2_first_half = generate_SNP_effects(z_matrix[,1:(floor(p/2))])
  #   z3_second_half = generate_SNP_effects(z_matrix[,(floor(p/2)+1):p])
  #   z5_first_half = generate_SNP_effects(z_matrix[,1:(floor(p/2))])
  #   z6_second_half = generate_SNP_effects(z_matrix[,(floor(p/2)+1):p])
  #   x1 = 1 / p * IV_strength * z + u + rnorm(2*n_max) 
  #   x2 = - 1 / p * IV_strength * z2_first_half - u + rnorm(2*n_max)
  #   x3 = 1 / p * IV_strength * z3_second_half - u + rnorm(2*n_max)
  #   x4 = 1 / p * IV_strength * z + u + rnorm(2*n_max) 
  #   x5 = - 1 / p * IV_strength * z5_first_half - u + rnorm(2*n_max)
  #   x6 = 1 / p * IV_strength * z6_second_half - u + rnorm(2*n_max)
  #   v = effect * (x1 + x2 + x3) + u
  #   odds = exp(-2 + v)
  } else if (pleiotropy == "vertical") {
    x1 = 1 / p * IV_strength * z + rnorm(2*n_max) + u
    x2 = 2 * x1 - u + rnorm(2*n_max)
    x3 = - x2 + u + rnorm(2*n_max)
    v = effect * x3 - u + rnorm(2*n_max)
    odds = exp(-2 + v)
  } else stop("Unrecognized pleiotropy")
  

  if (xQTL_from_one_sample) {
    if (outcome %in% c("binary_use_control_in_xQTL", "binary_use_all_in_xQTL")) {
      y = rbinom(2*n_max, size=1, prob=odds / (1+odds))
      controls = which(y == 0)[1:round(n*0.5)]
      cases = which(y == 1)[1:round(n*0.5)]
      controls_in_xQTL = which(y == 0)[c(
        seq_len(round(n*0.5*overlap)),
        round(n*0.5) + seq_len(round(n*0.5)-round(n*0.5*overlap))
      )]
      cases_in_xQTL = which(y == 1)[c(
        seq_len(round(n*0.5*overlap)),
        round(n*0.5) + seq_len(round(n*0.5)-round(n*0.5*overlap))
      )]
      samples_in_GWAS = c(cases, controls)
      
      if (outcome == "binary_use_control_in_xQTL") {
        samples_in_xQTL = controls_in_xQTL
      } else if (outcome == "binary_use_all_in_xQTL") {
        samples_in_xQTL = c(controls_in_xQTL, cases_in_xQTL)
      }
      
    } else if (outcome == "continuous") {
      y = v + rnorm(2*n_max)
      samples_in_GWAS = 1:n
      samples_in_xQTL = c(
        seq_len(round(n*overlap)),
        n + seq_len(n-round(n*overlap)))
    }
    samples_in_xQTL1 = samples_in_xQTL2 = samples_in_xQTL3 = samples_in_xQTL
  } else {
    if (outcome %in% c("binary_use_control_in_xQTL", "binary_use_all_in_xQTL")) {
      y = rbinom(2*n_max, size=1, prob=odds / (1+odds))
      controls = which(y == 0)[1:round(n*0.5)]
      cases = which(y == 1)[1:round(n*0.5)]
      controls_in_xQTL1 = which(y == 0)[round(n*0.5) + seq_len(round(n*0.5))]
      controls_in_xQTL2 = which(y == 0)[round(n*0.5)*2 + seq_len(round(n*0.5))]
      controls_in_xQTL3 = which(y == 0)[round(n*0.5)*3 + seq_len(round(n*0.5))]
      cases_in_xQTL1 = which(y == 1)[round(n*0.5) + seq_len(round(n*0.5))]
      cases_in_xQTL2 = which(y == 1)[round(n*0.5)*2 + seq_len(round(n*0.5))]
      cases_in_xQTL3 = which(y == 1)[round(n*0.5)*3 + seq_len(round(n*0.5))]
      samples_in_GWAS = c(cases, controls)
      
      if (outcome == "binary_use_control_in_xQTL") {
        samples_in_xQTL1 = controls_in_xQTL1
        samples_in_xQTL2 = controls_in_xQTL2
        samples_in_xQTL3 = controls_in_xQTL3
      } else if (outcome == "binary_use_all_in_xQTL") {
        samples_in_xQTL1 = c(controls_in_xQTL1, cases_in_xQTL1)
        samples_in_xQTL2 = c(controls_in_xQTL2, cases_in_xQTL2)
        samples_in_xQTL3 = c(controls_in_xQTL3, cases_in_xQTL3)
      }
      
    } else if (outcome == "continuous") {
      y = v + rnorm(2*n_max)
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
    
    if (outcome %in% c("binary_use_control_in_xQTL", "binary_use_all_in_xQTL")) {
      y0_glm_coefficients = summary(glm(y0~z0[,k], family="binomial"))$coefficients
      
      zy_est[k] = y0_glm_coefficients[2, "Estimate"]
      zy_se[k] = y0_glm_coefficients[2, "Std. Error"]
    } else if (outcome == "continuous") {
      y0_lm_coefficients = summary(lm(y0~z0[,k]))$coefficients
      
      zy_est[k] = y0_lm_coefficients[2, "Estimate"]
      zy_se[k] = y0_lm_coefficients[2, "Std. Error"]
    }
  }
  
  # assess instrument strength
  x1_lm_multi_summary = summary(lm(x1~z1))
  x2_lm_multi_summary = summary(lm(x2~z2))
  x3_lm_multi_summary = summary(lm(x3~z3))

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
    x3_lm_multi_summary$r.squared
    )
  names(statistics) = c("OLS_beta1", "OLS_beta2", "OLS_beta3",
                        "F_1", "R2_1", "F_2", "R2_2", "F_3", "R2_3")

  pval_array_i = get_p_values_from_summary(zx_est, zx_se, zy_est, zy_se, r2 = r2,
                                           include_coef_estimates = TRUE,
                                           include_number_of_SNPs = FALSE)
  
  c(pval_array_i, statistics)
}

pval_array = combined_array[,!grepl("_beta|^F_|^R2_", colnames(combined_array))]
statistics_array = combined_array[,grepl("_beta|^F_|^R2_", colnames(combined_array))]

powers = apply(pval_array, 2, function(x) sum(x < 0.05) / times)

means = apply(statistics_array, 2, mean)
names(means) = paste0(colnames(statistics_array), "_mean")

sds = apply(statistics_array, 2, sd)
names(sds) = paste0(colnames(statistics_array), "_sd")

simulation_details = c(simulation_setup, powers, means, sds)

if (ind == 1) {
  write.table(simulation_details, file = sprintf("simulation_%03d.csv", ind),  
              row.names=FALSE, col.names=TRUE, sep=",")
} else {
  write.table(simulation_details, file = sprintf("simulation_%03d.csv", ind),  
              row.names=FALSE, col.names=FALSE, sep=",")
}
