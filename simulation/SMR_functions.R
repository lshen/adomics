#' compute p-values with summary data as input
#'
#' @param zx_est estimate of regression coefficient of x1, x2, and x3 on instrumental variable(s) z.
#'        The matrix dimension is 3 (number of modalities) x p (number of instrumental variables).
#' @param zx_se  SE estimate of regression coefficient of x1, x2, and x3 on instrumental variable(s) z.
#'        The matrix dimension is 3 (number of modalities) x p (number of instrumental variables).
#' @param zy_est estimate of regression coefficient of z on instrumental variable(s) z
#'        A vector of length p (number of instrumental variables).
#' @param zy_se  SE estimate of regression coefficient of z on instrumental variable(s) z
#'        A vector of length p (number of instrumental variables).
#' @param r2     r^2 that is assumed to be known and the same between all instrumental variables
get_p_values_from_summary = function(zx_est, zx_se, zy_est, zy_se, r2 = 0, bp = NULL, max_r2 = 0.2, min_dist_bp = 1e5,
                                     include_coef_estimates = FALSE, include_number_of_SNPs = FALSE, need_ld_clump = TRUE) {
  m = nrow(zx_est)  # only m=3 has been implemented
  p = ncol(zx_est)
  
  # function to perform LD clumping
  # Here we do LD clumping separately for each modality
  # returns a vector of SNP indices, with the first one being the one with the
  # most significant statistics.
  LD_clumping = function(est, se, r2, bp, max_r2 = 0.2, min_dist_bp = 1e5) {
    # est and se may be a m by p matrix or a vector of length p.
    statistics = abs(est / se)
    if (is.matrix(statistics)) {
      statistics = apply(statistics, 2, max)
    }
    p = length(statistics)
    
    is_SNP_selected = rep(NA, p) # is a SNP selected?
    
    NA_SNPs = which(is.na(is_SNP_selected))  # undecided SNPs
    while (length(NA_SNPs) > 0) {
      # select the probe with min GWAS p-value
      selected_SNP_ind = which.max(statistics[NA_SNPs])[1]
      selected_SNP = NA_SNPs[selected_SNP_ind]
      is_SNP_selected[selected_SNP] = TRUE
      
      # 2nd round: label SNPs that need to be excluded
      NA_SNPs = NA_SNPs[NA_SNPs != selected_SNP]
      for (NA_SNP in which(is.na(is_SNP_selected))) {
        if ((abs(bp[selected_SNP] - bp[NA_SNP]) < min_dist_bp) || 
            (r2[selected_SNP, NA_SNP] > max_r2))
          is_SNP_selected[NA_SNP] = FALSE
      }
      NA_SNPs = which(is.na(is_SNP_selected))  # undecided SNPs
    }
    SNPs_selected = which(is_SNP_selected)
    SNPs_selected[order(statistics[SNPs_selected], decreasing = TRUE)]
  }
  
  SNPs = list()
  if (is.null(bp) || !need_ld_clump) {
    SNPs$MVMR = SNPs[[3]] = SNPs[[2]] = SNPs[[1]] = seq_len(p)
  } else {
    # LD clumping separately for each modality
    SNPs[[1]] = LD_clumping(zx_est[1, ], zx_se[1, ], r2 = r2, bp = bp, max_r2 = max_r2, min_dist_bp = min_dist_bp)
    SNPs[[2]] = LD_clumping(zx_est[2, ], zx_se[2, ], r2 = r2, bp = bp, max_r2 = max_r2, min_dist_bp = min_dist_bp)
    SNPs[[3]] = LD_clumping(zx_est[3, ], zx_se[3, ], r2 = r2, bp = bp, max_r2 = max_r2, min_dist_bp = min_dist_bp)
    # LD clumping together for Multivariable MR methods
    SNPs$MVMR = LD_clumping(zx_est, zx_se, r2 = r2, bp = bp, max_r2 = max_r2, min_dist_bp = min_dist_bp)
  }
   
  # LD matrix, unsigned, not squared, see https://mrcieu.github.io/TwoSampleMR/reference/ld_matrix.html
  if (need_ld_clump) {
    if (is.data.frame(r2)) {
      r2 = as.matrix(r2)
      LDmat = sqrt(r2)
    } else if (is.matrix(r2)) {
      LDmat = sqrt(r2)
    } else {
      LDmat = matrix(sqrt(r2), nrow=p, ncol=p)
      diag(LDmat) = 1
    }
  } else {
    LDmat = r2
  }
  
  # Under a fixed model and a weighted regression analysis, the residual standard error should be 1.
  # See the Supplementary Materials of:
  # Bowden, Jack, George Davey Smith, and Stephen Burgess. “Mendelian Randomization with Invalid Instruments: 
  #  Effect Estimation and Bias Detection through Egger Regression.” International Journal of Epidemiology 
  #  44, no. 2 (April 1, 2015): 512–25. https://doi.org/10.1093/ije/dyv080.
  # Burgess, Stephen, Frank Dudbridge, and Simon G. Thompson. “Re: ‘Multivariable Mendelian Randomization:
  #  The Use of Pleiotropic Genetic Variants to Estimate Causal Effects.’” American Journal of Epidemiology 181,
  #  no. 4 (February 15, 2015): 290–91. https://doi.org/10.1093/aje/kwv017.
  #
  # "Standard errors are calculated by only constraining the residual standard error to be 1
  # when it is less than 1 (which would imply under-dispersion), and otherwise allowing the 
  # residual standard error to take its estimated value."
  
  # According to the R package MendelianRandomization:
  # "the model used for estimation can be either random-effects or fixed-effect.
  #  The default option is to use a fixed-effect model when there are three or
  #  fewer genetic variants, and a random-effects model when there are four or more.
  #  The (multiplicative) random-effects model allows for heterogeneity between 
  #  the causal estimates targeted by the genetic variants by allowing
  #  over-dispersion in the regression model. Under-dispersion is not permitted
  #  (in case of under-dispersion, the residual standard error is set to 1, as
  #  in a fixed-effect analysis)."
  
  allow_overdispersion = function(p) ifelse(p >= 4, TRUE, FALSE)
  
  # Another option for distribution is "t-dist". The default option used in 
  # the R package MendelianRandomization is "normal"
  # distribution = "normal" 
  # distribution = "t-dist" 
  
  # Multivariable. Make sure you have at least as many genetic instruments as exposures.
  coef_multivariable = coef_multivariable_LD = rep(NA, m)
  if (length(SNPs$MVMR) <= m) {
    Multivariable = MultivariableLD = Multivariable_tdist = MultivariableLD_tdist = 1
  } else {
    # Most of the times the variants need to be uncorrelated when using summary data in multivariable MR
    # In the future we need to check the implementation in: https://github.com/qingyuanzhao/mr.raps/tree/multivariate
    summary_multivariable = summary(lm(zy_est[SNPs$MVMR]~0+t(zx_est)[SNPs$MVMR, ], weights = 1/as.numeric(zy_se[SNPs$MVMR])^2))
    if (allow_overdispersion(length(SNPs$MVMR))) {
      stat_multivariable = summary_multivariable$fstatistic[1] * min(1, summary_multivariable$sigma^2)
    } else {
      stat_multivariable = summary_multivariable$fstatistic[1] * summary_multivariable$sigma^2
    }
    coef_multivariable = coef(summary_multivariable)[,"Estimate"]
  
    # MultivariableLD
    stat_multivariable_LD = tryCatch({
      chol_LDmat_MVMR = chol(LDmat[SNPs$MVMR,SNPs$MVMR])
      zy_est_LD = backsolve(chol_LDmat_MVMR, zy_est[SNPs$MVMR] / zy_se[SNPs$MVMR], transpose = TRUE)
      zx_est_LD = t(backsolve(chol_LDmat_MVMR, t(zx_est)[SNPs$MVMR, ] / zy_se[SNPs$MVMR], transpose = TRUE))
      summary_multivariable_LD = summary(lm(zy_est_LD~0+t(zx_est_LD)))
      if (allow_overdispersion(length(SNPs$MVMR))) {
        summary_multivariable_LD$fstatistic[1] * min(1, summary_multivariable_LD$sigma^2)
      } else {
        summary_multivariable_LD$fstatistic[1] * summary_multivariable_LD$sigma^2
      }
    }, error = function(e) {NA})
    if (!is.na(stat_multivariable_LD)) {
      coef_multivariable_LD = coef(summary_multivariable_LD)[,"Estimate"]
    }
    
    # if (distribution == "normal") {
    Multivariable = pchisq(stat_multivariable * m, m, lower.tail = FALSE)
    MultivariableLD = pchisq(stat_multivariable_LD * m, m, lower.tail = FALSE)
    # } else {
    Multivariable_tdist = pf(stat_multivariable, m, length(SNPs$MVMR)-m, lower.tail = FALSE)
    MultivariableLD_tdist = pf(stat_multivariable_LD, m, length(SNPs$MVMR)-m, lower.tail = FALSE)
  }

  # IVW
  summary_x = beta_x = se_x = list()
  pval_x = pval_tdist_x = coef_IVW = rep(NA, m)
  for (k in seq_len(m)) {
    if (length(SNPs[[k]]) < 2) {
      pval_x[k] = pval_tdist_x[k] = 1
      coef_IVW[k] = NA
      next
    }
    summary_x[[k]] = summary(lm(zy_est[SNPs[[k]]]~0+zx_est[k,SNPs[[k]]], weights = 1/as.numeric(zy_se[SNPs[[k]]]^2)))
    beta_x[[k]] = summary_x[[k]]$coef[1,1]
    if (allow_overdispersion(length(SNPs[[k]]))) {
      se_x[[k]]   = summary_x[[k]]$coef[1,2] / min(1, summary_x[[k]]$sigma)
    } else {
      se_x[[k]]   = summary_x[[k]]$coef[1,2] / summary_x[[k]]$sigma
    }
    # if (distribution == "normal") {
    pval_x[k] = 2*pnorm(abs(beta_x[[k]]/se_x[[k]]), lower.tail=FALSE)
    # } else {
    pval_tdist_x[k] = 2*pt(abs(beta_x[[k]]/se_x[[k]]), length(SNPs[[k]])-1, lower.tail=FALSE)
    coef_IVW[k] = coef(summary_x[[k]])[,"Estimate"]
  }
  stat_Fisher = -2 * sum(log(pval_x))
  stat_Cauchy = mean(ifelse(pval_x < 1e-15, 1/pval_x/pi, tan((0.5 - pval_x) * pi)))
  stat_MinP = min(pval_x)
  stat_tdist_Fisher = -2 * sum(log(pval_tdist_x))
  stat_tdist_Cauchy = mean(ifelse(pval_tdist_x < 1e-15, 1/pval_tdist_x/pi, tan((0.5 - pval_tdist_x) * pi)))
  stat_tdist_MinP = min(pval_tdist_x)
  
  # IVW
  IVW_MinP = stat_MinP
  IVW_Fisher_chisq = pchisq(stat_Fisher, df = 2*m, lower.tail = FALSE)
  # # Fisher combination function: considers the correlation between modalities
  # mu = 2*m
  # rho = rep(NA, m)
  # # cosine similarity matrix; see https://stats.stackexchange.com/q/367216
  # beta_x_matrix = matrix(t(zx_est) / zy_se, ncol = m)
  # for (k in seq_len(m)) {
  #   beta_x_matrix[-SNPs[[k]], k] = 0
  # }
  # beta_x_matrix_scaled = beta_x_matrix / sqrt(rowSums(beta_x_matrix^2))
  # beta_x_matrix_scaled[!is.finite(beta_x_matrix_scaled)] = 0
  # # TODO: the similarity matrix below has variance term != 1 because of LD
  # sim = t(beta_x_matrix_scaled) %*% LDmat %*% beta_x_matrix_scaled
  # rho = sim[upper.tri(sim)]
  # r = rho #* (1 + (1-rho^2)/(2*(n-3)))
  # delta = 3.9081*r^2 + 0.0313*r^4 + 0.1022*r^6 - 0.1378*r^8 + 0.0941*r^10 #- 3.9081/n*(1-r^2)^2
  # sigma2 = 4*m + 2*sum(delta)
  # IVW_Fisher_gamma = pgamma(stat_Fisher, shape=mu^2/sigma2, scale=sigma2/mu, lower.tail=FALSE)
  IVW_Modality1 = pval_x[1]
  IVW_Modality2 = pval_x[2]
  IVW_Modality3 = pval_x[3]
  IVW_HMP =  harmonicmeanp::p.hmp(pval_x + .Machine$double.xmin, L = m)
  IVW_Cauchy_cauchy = pcauchy(stat_Cauchy, lower.tail = FALSE) # 0.5 - atan(stat_Cauchy)/pi
  
  IVW_tdist_MinP = stat_tdist_MinP
  IVW_tdist_Fisher_chisq = pchisq(stat_tdist_Fisher, df = 2*m, lower.tail = FALSE)
  # IVW_tdist_Fisher_gamma = pgamma(stat_tdist_Fisher, shape=mu^2/sigma2, scale=sigma2/mu, lower.tail=FALSE)
  IVW_tdist_Modality1 = pval_tdist_x[1]
  IVW_tdist_Modality2 = pval_tdist_x[2]
  IVW_tdist_Modality3 = pval_tdist_x[3]
  IVW_tdist_HMP =  ifelse(!is.finite(stat_tdist_Fisher[1]), NA, harmonicmeanp::p.hmp(pval_tdist_x + .Machine$double.xmin, L = m))
  IVW_tdist_Cauchy_cauchy = pcauchy(stat_tdist_Cauchy, lower.tail = FALSE) # 0.5 - atan(stat_tdist_Cauchy)/pi
  
  # WGLR
  # when a LD matrix is provided, account for correlation between SNPs:
  summary_x_LD = beta_x_LD = se_x_LD = list()
  pval_x_LD = pval_tdist_x_LD = coef_WGLR = rep(NA, m)
  chol_LDmat_list = list()
  pval_LDs = tryCatch({
    for (k in seq_len(m)) {
      if (length(SNPs[[k]]) < 2) {
        pval_x_LD[k] = pval_tdist_x_LD[k] = 1
        chol_LDmat_list[[k]] = matrix(1, nrow=1, ncol=1)
        coef_WGLR[k] = NA
        next
      }
      chol_LDmat_list[[k]] = chol(LDmat[SNPs[[k]],SNPs[[k]]])
      zy_est_LD_k = backsolve(chol_LDmat_list[[k]], zy_est[SNPs[[k]]] / zy_se[SNPs[[k]]], transpose = TRUE)
      zx_est_LD_k = backsolve(chol_LDmat_list[[k]], t(zx_est)[SNPs[[k]], k] / zy_se[SNPs[[k]]], transpose = TRUE)
      summary_x_LD[[k]] = summary(lm(zy_est_LD_k~0+zx_est_LD_k))
      beta_x_LD[[k]] = summary_x_LD[[k]]$coef[1,1]
      if (allow_overdispersion(length(SNPs[[k]]))) {
        se_x_LD[[k]]   = summary_x_LD[[k]]$coef[1,2] / min(1, summary_x_LD[[k]]$sigma)
      } else {
        se_x_LD[[k]]   = summary_x_LD[[k]]$coef[1,2] / summary_x_LD[[k]]$sigma
      }
      coef_WGLR[k] = coef(summary_x_LD[[k]])[,"Estimate"]
      
      # if (distribution == "normal") {
      pval_x_LD[k] = 2*pnorm(abs(beta_x_LD[[k]]/se_x_LD[[k]]), lower.tail=FALSE)
      # } else {
      pval_tdist_x_LD[k] = 2*pt(abs(beta_x_LD[[k]]/se_x_LD[[k]]), length(SNPs[[k]])-1, lower.tail=FALSE)
    }
    c(pval_x_LD, pval_tdist_x_LD)
  }, error = function(e) {NA})
  if (is.na(pval_LDs[1])) {
    pval_x_LD = pval_tdist_x_LD = NA
  } else {
    pval_x_LD = pval_LDs[seq_len(m)]
    pval_tdist_x_LD = pval_LDs[m + seq_len(m)]
  }
  stat_Fisher_LD = -2 * sum(log(pval_x_LD))
  stat_Cauchy_LD = mean(ifelse(pval_x_LD < 1e-15, 1/pval_x_LD/pi, tan((0.5 - pval_x_LD) * pi)))
  stat_MinP_LD = min(pval_x_LD)
  stat_tdist_Fisher_LD = -2 * sum(log(pval_tdist_x_LD))
  stat_tdist_Cauchy_LD = mean(ifelse(pval_x_LD < 1e-15, 1/pval_x_LD/pi, tan((0.5 - pval_tdist_x_LD) * pi)))
  stat_tdist_MinP_LD = min(pval_tdist_x_LD)
  
  # WGLR
  if (is.na(pval_LDs[1])) {
    WGLR_Cauchy_cauchy = WGLR_tdist_Cauchy_cauchy =
      WGLR_HMP = WGLR_tdist_HMP =
      WGLR_Modality1 = WGLR_tdist_Modality1 =
      WGLR_Modality2 = WGLR_tdist_Modality2 =
      WGLR_Modality3 = WGLR_tdist_Modality3 =
      WGLR_MinP = WGLR_tdist_MinP =
      WGLR_Fisher_chisq = WGLR_tdist_Fisher_chisq =
      WGLR_Fisher_gamma = WGLR_tdist_Fisher_gamma =
      GSMR_Cauchy_cauchy = GSMR_HMP = 
      GSMR_Modality1 = GSMR_Modality3 = GSMR_Modality3 = NA
  } else {
    WGLR_Cauchy_cauchy = pcauchy(stat_Cauchy_LD, lower.tail = FALSE) # 0.5 - atan(stat_Cauchy_LD)/pi
    WGLR_MinP = stat_MinP_LD
    WGLR_Fisher_chisq = pchisq(stat_Fisher_LD, df = 2*m, lower.tail = FALSE)
    # Fisher combination function: considers the correlation between modalities
    mu = 2*m
    rho = rep(NA, m)
    # cosine similarity matrix; see https://stats.stackexchange.com/q/367216
    beta_x_matrix = matrix(t(zx_est) / zy_se, ncol = m)
    beta_x_list = beta_1_list = list()
    for (k in seq_len(m)) {
      beta_x_list[[k]] = backsolve(chol_LDmat_list[[k]], beta_x_matrix[SNPs[[k]], k], transpose = TRUE)
      beta_1_list[[k]] = backsolve(chol_LDmat_list[[k]], diag(1, length(SNPs[[k]])), transpose = TRUE)
    }
    
    sim = matrix(nrow = m, ncol = m)
    for (j in seq_len(m)) {
      for (i in seq_len(j-1)) {
        sim[i, j] =  ((t(beta_x_list[[i]]) %*% beta_1_list[[i]]) / as.numeric(sqrt(t(beta_x_list[[i]]) %*% beta_x_list[[i]]))) %*%
                    LDmat[SNPs[[i]], SNPs[[j]]] %*% 
                    t((t(beta_x_list[[j]]) %*% beta_1_list[[j]]) / as.numeric(sqrt(t(beta_x_list[[j]]) %*% beta_x_list[[j]])))
      }
    }
    rho = sim[upper.tri(sim)]
    r = rho #* (1 + (1-rho^2)/(2*(n-3)))
    delta = 3.9081*r^2 + 0.0313*r^4 + 0.1022*r^6 - 0.1378*r^8 + 0.0941*r^10 #- 3.9081/n*(1-r^2)^2
    sigma2 = 4*m + 2*sum(delta)
    WGLR_Fisher_gamma = pgamma(stat_Fisher_LD, shape=mu^2/sigma2, scale=sigma2/mu, lower.tail=FALSE)
    WGLR_Modality1 = pval_x_LD[1]
    WGLR_Modality2 = pval_x_LD[2]
    WGLR_Modality3 = pval_x_LD[3]
    WGLR_HMP = ifelse(!is.finite(stat_Fisher_LD[1]), NA, harmonicmeanp::p.hmp(pval_x_LD + .Machine$double.xmin, L = m))
    
    WGLR_tdist_Cauchy_cauchy = pcauchy(stat_tdist_Cauchy_LD, lower.tail = FALSE) # 0.5 - atan(stat_tdist_Cauchy_LD)/pi
    WGLR_tdist_MinP = stat_tdist_MinP_LD
    WGLR_tdist_Fisher_chisq = pchisq(stat_tdist_Fisher_LD, df = 2*m, lower.tail = FALSE)
    WGLR_tdist_Fisher_gamma = pgamma(stat_tdist_Fisher_LD, shape=mu^2/sigma2, scale=sigma2/mu, lower.tail=FALSE)
    WGLR_tdist_Modality1 = pval_tdist_x_LD[1]
    WGLR_tdist_Modality2 = pval_tdist_x_LD[2]
    WGLR_tdist_Modality3 = pval_tdist_x_LD[3]
    WGLR_tdist_HMP = ifelse(!is.finite(stat_tdist_Fisher_LD[1]), NA, harmonicmeanp::p.hmp(pval_tdist_x_LD + .Machine$double.xmin, L = m))
  }
  
  # GSMR
  # The GSMR method involve inverting a covariance matrix,
  # which tends to give numerical problems when the IV is not strong enough.
  # The following code block follows gsmr:::bxy_gsmr()
  # https://cnsgenomics.com/software/gsmr/

  # In the supplementary materials of the SMR paper (2016 Zhu et al.),
  # The authors stated that z is selected to be strongly associated with x.
  # In the gsmr package, they specified that in terms of 1 df chisq distribution,
  # the association statistic needs to be larger than 10 (approx. pval < 1e-3)
  eps = 1e-6
  stat_GSMR = rep(NA, m)
  
  for (k in seq_len(m)) {
    if (length(SNPs[[k]]) < 2) {
      stat_GSMR[k] = 0
      next
    }
    xy_est_k = zy_est[SNPs[[k]]] / zx_est[k, SNPs[[k]]]
    
    V = diag(length(SNPs[[k]]))
    if (length(SNPs[[k]]) > 1) {
      V = LDmat[SNPs[[k]],SNPs[[k]]] * (tcrossprod(zy_se[SNPs[[k]]]) / tcrossprod(zx_est[k, SNPs[[k]]]) +
                   median(xy_est_k)^2 / tcrossprod(zx_est[k, SNPs[[k]]] / zx_se[k, SNPs[[k]]]))
    }
    diag(V) = diag(V) + eps
    V_eigen = eigen(V, symmetric = TRUE)
    V_eigenvalues = as.numeric(V_eigen$values)
    V_eigenvectors = V_eigen$vectors
    if (min(V_eigenvalues) < eps) {
      message("The covariance matrix is not invertible in GSMR. The p-value is set to NA.")
      stat_GSMR[k] = NA
    }
    V_inv = V_eigenvectors %*% diag(1/V_eigenvalues) %*% t(V_eigenvectors)
    vec_1 = rep(1, length(V_eigenvalues))
    var_GSMR_k = as.numeric(solve(t(vec_1) %*% V_inv %*% vec_1))
    est_GSMR_k = as.numeric(var_GSMR_k * as.numeric(t(vec_1) %*% V_inv) %*% xy_est_k)
    stat_GSMR[k] = est_GSMR_k^2 / var_GSMR_k
  }
  
  pval_GSMR = pchisq(stat_GSMR, 1, lower.tail = FALSE)
  stat_GSMR_Cauchy = mean(ifelse(pval_GSMR < 1e-15, 1/pval_GSMR/pi, tan((0.5 - pval_GSMR) * pi)))
  GSMR_Cauchy_cauchy = pcauchy(stat_GSMR_Cauchy, lower.tail = FALSE) # 0.5 - atan(stat_GSMR_Cauchy)/pi
  
  GSMR_HMP = ifelse(!is.finite(stat_GSMR[1]), NA, harmonicmeanp::p.hmp(pval_GSMR + .Machine$double.xmin, L = m))
  
  GSMR_Modality1 = pval_GSMR[1]
  GSMR_Modality2 = pval_GSMR[2]
  GSMR_Modality3 = pval_GSMR[3]
  
  GSMR_MinP = min(pval_GSMR)
  
  GSMR_Fisher_chisq = pchisq(-2*sum(log(pval_GSMR)), df = 2*m, lower.tail = FALSE)
  
  # directly combine the single SNP and single modality statistics (SMR)
  # SMR_statistics = (t(zx_est) / t(zx_se))^2 * (zy_est / zy_se)^2 / ((t(zx_est) / t(zx_se))^2 + (zy_est / zy_se)^2)
  SMR_statistics = c((zx_est[1, SNPs[[1]]] / zx_se[1, SNPs[[1]]])^2 * (zy_est[SNPs[[1]]] / zy_se[SNPs[[1]]])^2 /
                       ((zx_est[1, SNPs[[1]]] / zx_se[1, SNPs[[1]]])^2 + (zy_est[SNPs[[1]]] / zy_se[SNPs[[1]]])^2),
                     (zx_est[2, SNPs[[2]]] / zx_se[2, SNPs[[2]]])^2 * (zy_est[SNPs[[2]]] / zy_se[SNPs[[2]]])^2 /
                       ((zx_est[2, SNPs[[2]]] / zx_se[2, SNPs[[2]]])^2 + (zy_est[SNPs[[2]]] / zy_se[SNPs[[2]]])^2),
                     (zx_est[3, SNPs[[3]]] / zx_se[3, SNPs[[3]]])^2 * (zy_est[SNPs[[3]]] / zy_se[SNPs[[3]]])^2 /
                       ((zx_est[3, SNPs[[3]]] / zx_se[3, SNPs[[3]]])^2 + (zy_est[SNPs[[3]]] / zy_se[SNPs[[3]]])^2)
                     )
  stat_SMR_Fisher = -2 * sum(log(pchisq(SMR_statistics, 1, lower.tail = FALSE)))
  SMR_allSNPs_Fisher_chisq = pchisq(stat_SMR_Fisher, df = 2*length(SMR_statistics), lower.tail = FALSE)
  pval_SMR = pchisq(SMR_statistics, 1, lower.tail = FALSE)
  stat_SMR_Cauchy = mean(ifelse(pval_SMR < 1e-15, 1/pval_SMR/pi, tan((0.5 - pval_SMR) * pi)))
  SMR_allSNPs_Cauchy_cauchy = pcauchy(stat_SMR_Cauchy, lower.tail = FALSE) # 0.5 - atan(stat_SMR_Cauchy)/pi
  stat_SMR_MinP = min(pval_SMR)
  SMR_allSNPs_MinP = stat_SMR_MinP
  SMR_allSNPs_HMP = ifelse(!is.finite(SMR_statistics[1]), NA, harmonicmeanp::p.hmp(pval_SMR + .Machine$double.xmin, L = m * p))
  
  # single SNP
  # GWAS_multipleSNPs = pnorm(zy_est / zy_se) * 2
  # topSNP = which.min(GWAS_multipleSNPs)
  # GWAS_singleSNP = GWAS_multipleSNPs[topSNP]
  # topSNP -- select the top one according to LD clumping using the strongest xQTL for each modality
  statistics_SMR_topSNP = c((zx_est[1, SNPs[[1]][1]] / zx_se[1, SNPs[[1]]][1])^2 * (zy_est[SNPs[[1]]][1] / zy_se[SNPs[[1]]][1])^2 /
                              ((zx_est[1, SNPs[[1]][1]] / zx_se[1, SNPs[[1]]][1])^2 + (zy_est[SNPs[[1]]][1] / zy_se[SNPs[[1]]][1])^2),
                            (zx_est[2, SNPs[[2]]][1] / zx_se[2, SNPs[[2]]][1])^2 * (zy_est[SNPs[[2]]][1] / zy_se[SNPs[[2]]][1])^2 /
                              ((zx_est[2, SNPs[[2]]][1] / zx_se[2, SNPs[[2]]][1])^2 + (zy_est[SNPs[[2]]][1] / zy_se[SNPs[[2]]][1])^2),
                            (zx_est[3, SNPs[[3]]][1] / zx_se[3, SNPs[[3]]][1])^2 * (zy_est[SNPs[[3]]][1] / zy_se[SNPs[[3]]][1])^2 /
                              ((zx_est[3, SNPs[[3]]][1] / zx_se[3, SNPs[[3]]][1])^2 + (zy_est[SNPs[[3]]][1] / zy_se[SNPs[[3]]][1])^2)
                           )
  stat_topSNP_Fisher = -2 * sum(log(pchisq(statistics_SMR_topSNP, 1, lower.tail = FALSE)))
  SMR_singleSNP_Fisher_chisq = pchisq(stat_topSNP_Fisher, df = 2*length(statistics_SMR_topSNP), lower.tail = FALSE)
  # # Fisher combination function using gamma distribution
  # # This does not work as intended because we don't know what is the
  # # distribution of the GWAS Z statistics under the null hypothesis.
  # # considers the correlation between modalities
  # mu = 2*m
  # delta = rep(0.194, m)  # The number 0.194 is based on simulation when the z values from different xQTLs are independent
  # # suppose we don't have individual level measurement of modalities
  # # and we use the association between GWAS and xQTL summary data instead
  # # TODO: we need to properly calculate the measurement of correlation
  # sigma2 = 4*m + 2*sum(delta)
  # SMR_singleSNP_Fisher_gamma = pgamma(stat_topSNP_Fisher, shape=mu^2/sigma2, scale=sigma2/mu, lower.tail=FALSE)
  
  pval_topSNP = pchisq(statistics_SMR_topSNP, 1, lower.tail = FALSE)
  stat_topSNP_Cauchy = mean(ifelse(pval_topSNP < 1e-15, 1/pval_topSNP/pi, tan((0.5 - pval_topSNP) * pi)))
  SMR_singleSNP_Cauchy_cauchy =  pcauchy(stat_topSNP_Cauchy, lower.tail = FALSE) # 0.5 - atan(stat_topSNP_Cauchy)/pi
  
  SMR_singleSNP_HMP = ifelse(!is.finite(statistics_SMR_topSNP[1]), NA, harmonicmeanp::p.hmp(pval_topSNP + .Machine$double.xmin, L = m))
  
  SMR_Modality1 = pval_topSNP[1]
  SMR_Modality2 = pval_topSNP[2]
  SMR_Modality3 = pval_topSNP[3]
  
  vec = c(
    # multiple SNPs
    Multivariable, Multivariable_tdist,
    MultivariableLD, MultivariableLD_tdist,
    IVW_Cauchy_cauchy, IVW_tdist_Cauchy_cauchy,
    IVW_HMP, IVW_tdist_HMP, 
    IVW_Modality1, IVW_tdist_Modality1, 
    IVW_Modality2, IVW_tdist_Modality2, 
    IVW_Modality3, IVW_tdist_Modality3,
    IVW_MinP, IVW_tdist_MinP,
    IVW_Fisher_chisq, IVW_tdist_Fisher_chisq, 
    # IVW_Fisher_gamma, IVW_tdist_Fisher_gamma, 
    WGLR_Cauchy_cauchy, WGLR_tdist_Cauchy_cauchy,
    WGLR_HMP, WGLR_tdist_HMP,
    WGLR_Modality1, WGLR_tdist_Modality1,
    WGLR_Modality2, WGLR_tdist_Modality2,
    WGLR_Modality3, WGLR_tdist_Modality3,
    WGLR_MinP, WGLR_tdist_MinP,
    WGLR_Fisher_chisq, WGLR_tdist_Fisher_chisq,
    WGLR_Fisher_gamma, WGLR_tdist_Fisher_gamma,
    GSMR_Cauchy_cauchy, GSMR_HMP, 
    GSMR_Modality1, GSMR_Modality3, GSMR_Modality3,
    GSMR_MinP, GSMR_Fisher_chisq, 
    SMR_allSNPs_Cauchy_cauchy, SMR_allSNPs_HMP, 
    SMR_allSNPs_MinP, SMR_allSNPs_Fisher_chisq, 
    # single SNP
    SMR_singleSNP_Cauchy_cauchy, SMR_singleSNP_HMP, 
    SMR_Modality1, SMR_Modality2, SMR_Modality3,
    SMR_singleSNP_Fisher_chisq)
  
  names(vec) = c(
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
    # "IVW_Fisher_gamma", "IVW_tdist_Fisher_gamma",
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
  
  if (include_coef_estimates) {
    coef_estimates = c(coef_multivariable,
                       coef_multivariable_LD,
                       coef_IVW,
                       coef_WGLR)
    if (length(coef_estimates) != 12) {
      write.csv(data.frame(coefs=coef_estimates), "coef_estimate_error.out")
      coef_estimates = rep(NA, 12)
    }
    names(coef_estimates) = c("Multivariable_beta1",
                              "Multivariable_beta2",
                              "Multivariable_beta3",
                              "MultivariableLD_beta1",
                              "MultivariableLD_beta2",
                              "MultivariableLD_beta3",
                              "IVW_beta1", "IVW_beta2", "IVW_beta3",
                              "WGLR_beta1", "WGLR_beta2", "WGLR_beta3"
                              )
    vec = c(vec, coef_estimates)
  }
  
  if (include_number_of_SNPs) {
    number_of_SNPs = sapply(SNPs, length)
    names(number_of_SNPs) = c(paste0("number_of_SNPs_", seq_len(m)), "number_of_SNPs_MVMR")
    vec = c(vec, number_of_SNPs)
  }
  
  vec
}
