library(tidyverse)

# Generate simulation settings
simulations = expand.grid(
  outcome = c("continuous", "binary_use_control_in_xQTL", "binary_use_all_in_xQTL"),
  pleiotropy = c("horizontal", "vertical"),
  #sample_size = c(500, 1000, 2000),
  IV_strength = c(0.5, 1, 2),
  effect = c(0, 0.1),
  overlap = c(0, 0.5, 1),
  r2 = c(0, 0.01, 0.2),
  number_of_SNPs = c(5, 20),
  xQTL_from_one_sample = c(TRUE, FALSE)
)

r2_values = c(0, 0.01, 0.2)

data_all = read.csv("result.csv")

data_all$overlap = sprintf("%s; xQTLs from same sample", data_all$overlap)
data_all$overlap[!data_all$xQTL_from_one_sample] = "0; separate xQTL samples"

library(tidyr)
library(wesanderson)
library(ggplot2)

for (number_of_SNPs in c(5, 20)) {
  
  data = data_all[data_all$number_of_SNPs == number_of_SNPs,]
  
  data_long = gather(data, "method", "power", Multivariable:SMR_singleSNP_Fisher_chisq, factor_key=TRUE)
  #data_long = gather(data, "method", "power", second_order_Fisher_chisq:first_order_MinP, factor_key=TRUE)
  #data_long = gather(data, "method", "power", p_Multivariable:p_SMR_topSNP_MinP, factor_key=TRUE)
  # data_long
  # data_long = data_long[!grepl("GSMR|Multivariable", data_long$method), ]
  data_long$effect = factor(data_long$effect)
  data_long$positive_rate = data_long$power
  data_long$positive_rate[data_long$effect == 0] = - data_long$positive_rate[data_long$effect == 0]
  
  data_long_power = data_long[data_long$effect == 0.1, ]
  data_long_type_I_error = data_long[data_long$effect == 0, ]
  
  all.equal(data_long_power$outcome, data_long_type_I_error$outcome)
  all.equal(data_long_power$pleiotropy, data_long_type_I_error$pleiotropy)
  all.equal(data_long_power$IV_strength, data_long_type_I_error$IV_strength)
  all.equal(data_long_power$overlap, data_long_type_I_error$overlap)
  all.equal(data_long_power$r2, data_long_type_I_error$r2)
  
  data_long_power$type_I_error = data_long_type_I_error$power
  
  # order the barplot by mean type I error
  aggr = aggregate(power ~ method, data=data_long_type_I_error, mean)
  aggr = aggr[order(aggr$power), ]
  data_long$method = factor(data_long$method, levels = rev(levels(aggr$method)))
  method_color = wes_palette(length(levels(data_long$method)),
                                       name = "Darjeeling1", type = "continuous")
  names(method_color) = levels(data_long$method)
  
  # select a certain outcome, pleiotropy, and r2
  outcome = "continuous"
  pleiotropy = "horizontal"
  r2 = 0.2
  
  # TODO: more simulations -- 10,000 would be ideal
  
  # TODO: some LD methods are not working as intended and we need to remove them right now.
  # data_long_power = data_long_power[!data_long_power$method %in% c("p_MultivariableLD", "p_WGLR_Cauchy_cauchy", "p_WGLR_oneExposure_Cauchy_cauchy", "p_WGLR_twoExposures_Cauchy_cauchy"),]
  
  # scale_factor = 3
  # scale_function = function(x) {0.5*(scale_factor+1)*x-0.5*(scale_factor-1)*abs(x)}
  # data_long$positive_rate = scale_function(pmax(-0.2, data_long$positive_rate))
  
  pdf(paste0("SNP_",number_of_SNPs,"_simulation_barplot.pdf"), width=8.5, height=16.5)
  for (rtwo in r2_values) {
    for (o in c("continuous", "binary_use_control_in_xQTL", "binary_use_all_in_xQTL")) {
      for (p in c("horizontal", "vertical")) {
        print(ggplot(subset(data_long, outcome == o & pleiotropy == p & r2 == rtwo)
                     , aes(x = method, y = positive_rate, fill = method)) +
                geom_bar(position="stack", stat="identity") +
                facet_grid(overlap ~ IV_strength) +
                geom_vline(xintercept=0.05, linetype="solid", color = "black", size = 0.2) +
                geom_rect(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=0, fill="grey93", color=NA, alpha=0.05) +
                geom_hline(yintercept=0, linetype="solid", color = "black") +
                geom_hline(yintercept=-0.05, linetype="dotted", color = "red") +
                #geom_hline(yintercept=-0.05*scale_factor, linetype="dotted", color = "red") +
                #scale_y_continuous(labels = scale_function, limits = scale_function(c(-0.2, 1))) +
                scale_fill_manual(values=method_color) +
                theme_bw() +
                coord_flip() +
                theme(axis.line = element_line(colour = "black"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      plot.margin = margin(t = .35, r = .1, b = 0, l = .1, unit = "in"),
                      legend.position = "none",
                      legend.text=element_text(size=8),
                      axis.text.x = element_text(size=8, angle=90),
                      axis.text.y = element_text(size=8)) +
                guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
                ggtitle(sprintf("outcome type = %s, pleiotropy = %s, r2 = %.2f
                          IV strength: {0.5, 1, 2}, overlap: {0, 0.5, 1}", o, p, rtwo)))  
      }
    }
  }
  dev.off()
  
  pdf(paste0("SNP_",number_of_SNPs,"_simulation_barplot_subset.pdf"), width=8.5, height=8.5)
  for (rtwo in r2_values) {
    for (o in c("continuous", "binary_use_control_in_xQTL", "binary_use_all_in_xQTL")) {
      for (p in c("horizontal", "vertical")) {
        print(ggplot(subset(data_long, outcome == o & pleiotropy == p & r2 == rtwo &
                              !grepl("^IVW|Fisher|Modality|MinP|^Multivariable$|^Multivariable_tdist$", method))
                     , aes(x = method, y = positive_rate, fill = method)) +
                geom_bar(position="stack", stat="identity") +
                facet_grid(overlap ~ IV_strength) +
                geom_vline(xintercept=0.05, linetype="solid", color = "black", size = 0.2) +
                geom_rect(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=0, fill="grey93", color=NA, alpha=0.05) +
                geom_hline(yintercept=0, linetype="solid", color = "black") +
                geom_hline(yintercept=-0.05, linetype="dotted", color = "red") +
                #geom_hline(yintercept=-0.05*scale_factor, linetype="dotted", color = "red") +
                #scale_y_continuous(labels = scale_function, limits = scale_function(c(-0.2, 1))) +
                scale_fill_manual(values=method_color) +
                theme_bw() +
                coord_flip() +
                theme(axis.line = element_line(colour = "black"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      plot.margin = margin(t = .35, r = .1, b = 0, l = .1, unit = "in"),
                      legend.position = "none",
                      legend.text=element_text(size=8),
                      axis.text.x = element_text(size=8, angle=90),
                      axis.text.y = element_text(size=8)) +
                guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
                ggtitle(sprintf("outcome type = %s, pleiotropy = %s, r2 = %.2f
                          IV strength: {0.5, 1, 2}, overlap: {0, 0.5, 1}", o, p, rtwo)))  
      }
    }
  }
  dev.off()
  
  
  pdf(paste0("SNP_",number_of_SNPs,"_simulation_barplot_outcome.pdf"), width=8.5, height=12.5)
  for (rtwo in r2_values) {
    for (s in c(0.5, 1, 2)) {
      for (p in c("horizontal", "vertical")) {
        print(ggplot(subset(data_long, IV_strength == s & pleiotropy == p & r2 == rtwo)
                     , aes(x = method, y = positive_rate, fill = method)) +
                geom_bar(position="stack", stat="identity") +
                facet_grid(overlap ~ outcome) +
                geom_vline(xintercept=0.05, linetype="solid", color = "black", size = 0.2) +
                geom_rect(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=0, fill="grey93", color=NA, alpha=0.05) +
                geom_hline(yintercept=0, linetype="solid", color = "black") +
                geom_hline(yintercept=-0.05, linetype="dotted", color = "red") +
                #geom_hline(yintercept=-0.05*scale_factor, linetype="dotted", color = "red") +
                #scale_y_continuous(labels = scale_function, limits = scale_function(c(-0.2, 1))) +
                scale_fill_manual(values=method_color) +
                theme_bw() +
                coord_flip() +
                theme(axis.line = element_line(colour = "black"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      plot.margin = margin(t = .35, r = .1, b = 0, l = .1, unit = "in"),
                      legend.position = "none",
                      legend.text=element_text(size=8),
                      axis.text.x = element_text(size=8, angle=90),
                      axis.text.y = element_text(size=8)) +
                guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
                ggtitle(sprintf("IV strength = %s, pleiotropy = %s, r2 = %.2f
                                outcome type: x axis, overlap: {0, 0.5, 1}", s, p, rtwo)))
      }
    }
  }
  dev.off()
  
  
  pdf(paste0("SNP_",number_of_SNPs,"_simulation_barplot_outcome_subset.pdf"), width=8.5, height=8.5)
  for (rtwo in r2_values) {
    for (s in c(0.5, 1, 2)) {
      for (p in c("horizontal", "vertical")) {
        print(ggplot(subset(data_long, IV_strength == s & pleiotropy == p & r2 == rtwo &
                              !grepl("^IVW|Fisher|OneModality|MinP|^Multivariable$", method))
                     , aes(x = method, y = positive_rate, fill = method)) +
                geom_bar(position="stack", stat="identity") +
                facet_grid(overlap ~ outcome) +
                geom_vline(xintercept=0.05, linetype="solid", color = "black", size = 0.2) +
                geom_rect(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=0, fill="grey93", color=NA, alpha=0.05) +
                geom_hline(yintercept=0, linetype="solid", color = "black") +
                geom_hline(yintercept=-0.05, linetype="dotted", color = "red") +
                #geom_hline(yintercept=-0.05*scale_factor, linetype="dotted", color = "red") +
                #scale_y_continuous(labels = scale_function, limits = scale_function(c(-0.2, 1))) +
                scale_fill_manual(values=method_color) +
                theme_bw() +
                coord_flip() +
                theme(axis.line = element_line(colour = "black"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      plot.margin = margin(t = .35, r = .1, b = 0, l = .1, unit = "in"),
                      legend.position = "none",
                      legend.text=element_text(size=8),
                      axis.text.x = element_text(size=8, angle=90),
                      axis.text.y = element_text(size=8)) +
                guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
                ggtitle(sprintf("IV strength = %s, pleiotropy = %s, r2 = %.2f
                                outcome type: x axis, overlap: {0, 0.5, 1}", s, p, rtwo)))
      }
    }
  }
  dev.off()
  
  ############################################################
  # multiple vs. single SNP
  ############################################################
  # Under horizontal pleiotropy, single SNP methods perform much worse,
  # while the difference is lower under vertical pleiotropy.
  # Note that there is inflated type I error due to pleiotropy for single SNP methods.
  method_col = method_color
  names(method_col)[names(method_col) =="IVW_Cauchy_cauchy"] = "multipleSNPs_IVW"
  names(method_col)[names(method_col) =="WGLR_Cauchy_cauchy"] = "multipleSNPs_WGLR"
  names(method_col)[names(method_col) =="SMR_allSNPs_Cauchy_cauchy"] = "multipleSNPs_SMR"
  names(method_col)[names(method_col) =="SMR_singleSNP_Cauchy_cauchy"] = "singleSNP_SMR"
  pdf(paste0("SNP_",number_of_SNPs,"_simulation_barplot_multiple_vs_single_SNP.pdf"), width=6, height=4)
  rtwo = 0.2
  o = "binary_use_all_in_xQTL"
  print(ggplot(subset(data_long, outcome == o & r2 == rtwo & overlap == "0; xQTLs from same sample" &
                                 grepl("cauchy$", method)) %>%
                 dplyr::mutate(method = case_when(
                   method == "IVW_Cauchy_cauchy" ~ "multipleSNPs_IVW",
                   method == "WGLR_Cauchy_cauchy" ~ "multipleSNPs_WGLR",
                   method == "SMR_allSNPs_Cauchy_cauchy" ~ "multipleSNPs_SMR",
                   method == "SMR_singleSNP_Cauchy_cauchy" ~ "singleSNP_SMR"
                   ))
               , aes(x = method, y = positive_rate, fill = method)) +
          geom_bar(position="stack", stat="identity") +
          facet_grid(pleiotropy ~ IV_strength) +
          geom_vline(xintercept=0.05, linetype="solid", color = "black", size = 0.2) +
          geom_rect(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=0, fill="grey80", color=NA, alpha=0.05) +
          geom_hline(yintercept=0, linetype="solid", color = "black") +
          geom_hline(yintercept=-0.05, linetype="dotted", color = "red") +
          #geom_hline(yintercept=-0.05*scale_factor, linetype="dotted", color = "red") +
          #scale_y_continuous(labels = scale_function, limits = scale_function(c(-0.2, 1))) +
          #scale_fill_scico_d() +
          scale_fill_manual(values=method_col) +
          theme_bw() +
          coord_flip() +
          theme(axis.line = element_line(colour = "black"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                plot.margin = margin(t = .35, r = .1, b = 0, l = .1, unit = "in"),
                legend.position = "none",
                legend.text=element_text(size=8),
                axis.text.x = element_text(size=8, angle=90),
                axis.text.y = element_text(size=8)) +
          guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
          ggtitle(paste0(
            "Cauchy combination function\n",
            "outcome type = binary_use_all_in_xQTL, r2 = 0.1\n",
            "overlap: 0; xQTLs from same sample")))
  
  dev.off()
  
  # pdf(paste0("SNP_",number_of_SNPs,"simulation_linegraph_multiple_vs_single_SNP.pdf"), width=6, height=4)
  # rtwo = 0.2
  # o = "binary_use_all_in_xQTL"
  # print(ggplot(subset(data_long, outcome == o & r2 == rtwo & overlap == "0; xQTLs from same sample" &
  #                       grepl("cauchy$", method)) %>%
  #                dplyr::mutate(method = case_when(
  #                  method == "IVW_Cauchy_cauchy" ~ "multipleSNPs_IVW",
  #                  method == "WGLR_Cauchy_cauchy" ~ "multipleSNPs_WGLR",
  #                  method == "SMR_allSNPs_Cauchy_cauchy" ~ "multipleSNPs_SMR",
  #                  method == "SMR_singleSNP_Cauchy_cauchy" ~ "singleSNP_SMR"
  #                ))
  #              , aes(x = effect, y = power, group = method, color = method)) +
  #         geom_line() +
  #         facet_grid(pleiotropy ~ IV_strength) +
  #         geom_vline(xintercept=0.05, linetype="solid", color = "black", size = 0.2) +
  #         geom_rect(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=0, fill="grey80", color=NA, alpha=0.05) +
  #         geom_hline(yintercept=0, linetype="solid", color = "black") +
  #         geom_hline(yintercept=0.05, linetype="dotted", color = "red") +
  #         #geom_hline(yintercept=-0.05*scale_factor, linetype="dotted", color = "red") +
  #         #scale_y_continuous(labels = scale_function, limits = scale_function(c(-0.2, 1))) +
  #         #scale_color_scico_d() +
  #         scale_color_manual(values=method_col) +
  #         theme_bw() +
  #         #coord_flip() +
  #         theme(axis.line = element_line(colour = "black"),
  #               panel.grid.major = element_blank(),
  #               panel.grid.minor = element_blank(),
  #               panel.background = element_blank(),
  #               plot.margin = margin(t = .35, r = .1, b = 0, l = .1, unit = "in"),
  #               legend.position = "right",
  #               legend.text=element_text(size=8),
  #               axis.text.x = element_text(size=8, angle=90),
  #               axis.text.y = element_text(size=8)) +
  #         guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  #         ggtitle(paste0(
  #           "Cauchy combination function\n",
  #           "outcome type = binary_use_all_in_xQTL, r2 = 0.1\n",
  #           "overlap: 0; xQTLs from same sample")))
  # 
  # dev.off()
  
  ############################################################
  # p-value combinations (WGLR)
  ############################################################
  
  pdf(paste0("SNP_",number_of_SNPs,"_simulation_barplot_pval_combination_WGLR.pdf"), width=6, height=4)
  rtwo = 0.2
  o = "binary_use_all_in_xQTL"
  print(ggplot(subset(data_long, outcome == o & r2 == rtwo & overlap == "0; xQTLs from same sample" &
                        grepl("^WGLR", method))
               , aes(x = method, y = positive_rate, fill = method)) +
          geom_bar(position="stack", stat="identity") +
          facet_grid(pleiotropy ~ IV_strength) +
          geom_vline(xintercept=0.05, linetype="solid", color = "black", size = 0.2) +
          geom_rect(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=0, fill="grey80", color=NA, alpha=0.05) +
          geom_hline(yintercept=0, linetype="solid", color = "black") +
          geom_hline(yintercept=-0.05, linetype="dotted", color = "red") +
          #geom_hline(yintercept=-0.05*scale_factor, linetype="dotted", color = "red") +
          #scale_y_continuous(labels = scale_function, limits = scale_function(c(-0.2, 1))) +
          #scale_fill_scico_d() +
          scale_fill_manual(values=method_color) +
          theme_bw() +
          coord_flip() +
          theme(axis.line = element_line(colour = "black"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                plot.margin = margin(t = .35, r = .1, b = 0, l = .1, unit = "in"),
                legend.position = "none",
                legend.text=element_text(size=8),
                axis.text.x = element_text(size=8, angle=90),
                axis.text.y = element_text(size=8)) +
          guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
          ggtitle(paste0(
            "Different p-value combination methods for WGLR\n",
            "outcome type = binary_use_all_in_xQTL, r2 =", rtwo, "\n",
            "overlap: 0; xQTLs from same sample")))
  
  dev.off()
  
  
  
  ############################################################
  # p-value combinations (SMR)
  ############################################################
  
  pdf(paste0("SNP_",number_of_SNPs,"_simulation_barplot_pval_combination_SMR_singleSNP.pdf"), width=6, height=4)
  rtwo = 0.2
  o = "binary_use_all_in_xQTL"
  print(ggplot(subset(data_long, outcome == o & r2 == rtwo & overlap == "0; xQTLs from same sample" &
                        grepl("SMR_singleSNP_Cauchy_cauchy|SMR_singleSNP_Fisher_chisq|SMR_singleSNP_HMP|SMR_OneModality", method))
               , aes(x = method, y = positive_rate, fill = method)) +
          geom_bar(position="stack", stat="identity") +
          facet_grid(pleiotropy ~ IV_strength) +
          geom_vline(xintercept=0.05, linetype="solid", color = "black", size = 0.2) +
          geom_rect(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=0, fill="grey80", color=NA, alpha=0.05) +
          geom_hline(yintercept=0, linetype="solid", color = "black") +
          geom_hline(yintercept=-0.05, linetype="dotted", color = "red") +
          #geom_hline(yintercept=-0.05*scale_factor, linetype="dotted", color = "red") +
          #scale_y_continuous(labels = scale_function, limits = scale_function(c(-0.2, 1))) +
          #scale_fill_scico_d() +
          scale_fill_manual(values=method_color) +
          theme_bw() +
          coord_flip() +
          theme(axis.line = element_line(colour = "black"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                plot.margin = margin(t = .35, r = .1, b = 0, l = .1, unit = "in"),
                legend.position = "none",
                legend.text=element_text(size=8),
                axis.text.x = element_text(size=8, angle=90),
                axis.text.y = element_text(size=8)) +
          guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
          ggtitle(paste0(
            "Different p-value combination methods for SMR\n",
            "outcome type = binary_use_all_in_xQTL, r2 = ", rtwo, "\n",
            "overlap: 0; xQTLs from same sample")))
  
  dev.off()
  
  ############################################################
  # outcome and overlap
  ############################################################
  
  
  pdf(paste0("SNP_",number_of_SNPs,"_simulation_barplot_outcome_and_overlap.pdf"), width=8.5, height=8.5)
  s = 1
  rtwo = 0.2
  for (p in c("horizontal", "vertical")) {
    print(ggplot(subset(data_long, IV_strength == s & pleiotropy == p & r2 == rtwo & method %in%
                          c("WGLR_Cauchy_cauchy", "SMR_singleSNP_Cauchy_cauchy") )
                 , aes(x = method, y = positive_rate, fill = method)) +
            geom_bar(position="stack", stat="identity") +
            facet_grid(overlap ~ outcome) +
            geom_vline(xintercept=0.05, linetype="solid", color = "black", size = 0.2) +
            geom_rect(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=0, fill="grey93", color=NA, alpha=0.05) +
            geom_hline(yintercept=0, linetype="solid", color = "black") +
            geom_hline(yintercept=-0.05, linetype="dotted", color = "red") +
            #geom_hline(yintercept=-0.05*scale_factor, linetype="dotted", color = "red") +
            #scale_y_continuous(labels = scale_function, limits = scale_function(c(-0.2, 1))) +
            scale_fill_manual(values=method_color) +
            theme_bw() +
            coord_flip() +
            theme(axis.line = element_line(colour = "black"),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  plot.margin = margin(t = .35, r = .1, b = 0, l = .1, unit = "in"),
                  legend.position = "none",
                  legend.text=element_text(size=8),
                  axis.text.x = element_text(size=8, angle=90),
                  axis.text.y = element_text(size=8)) +
            guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
            ggtitle(paste0(
              "Different outcome and overlap situations\n",
              "outcome type = binary_use_all_in_xQTL, r2 = ", rtwo, "\n",
              "IV Strength = 1, pleiotropy = ", p)))
  }
  dev.off()

}

# pdf("simulation.pdf", width=7.5, height=5.5)
# for (rtwo in r2_values) {
#   for (o in c("continuous", "binary_use_control_in_xQTL", "binary_use_all_in_xQTL")) {
#     for (p in c("horizontal", "vertical")) {
#        print(ggplot(subset(data_long_power, outcome == o & pleiotropy == p & r2 == rtwo)
# , aes(x = type_I_error, y = power)) +
#         geom_point(aes(col = method, shape = method), stroke = 0.2) +
#         facet_grid(overlap ~ IV_strength) +
#         geom_vline(xintercept=0.05, linetype="solid", color = "black", size = 0.2) +
#         geom_rect(xmin=0.05, xmax=2, ymin=-1, ymax=2, fill="grey80", color=NA, alpha=0.05) +
#         scale_x_continuous(limits = c(0, max(data_long_power$type_I_error))) +
#         scale_shape_manual(values=1:nlevels(data_long_power$method)) +
#         scale_color_manual(values=wes_palette(name="Darjeeling1", n = length(levels(data_long$method)), type = "continuous")) +
#         theme_bw() +
#         theme(axis.line = element_line(colour = "black"),
#               panel.grid.major = element_blank(),
#               panel.grid.minor = element_blank(),
#               panel.background = element_blank(),
#               plot.margin = margin(t = .35, r = .1, b = 0, l = .1, unit = "in"),
#               # legend.position = "bottom",
#               legend.text=element_text(size=8),
#               axis.text.x = element_text(size=8, angle=90),
#               axis.text.y = element_text(size=8)) +
#         guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
#         ggtitle(sprintf("outcome type = %s, pleiotropy = %s, r2 = %.2f
#                         IV strength: {0.5, 1, 2}, overlap: {0, 0.5, 1}", o, p, rtwo)))    
#       }
#   }
# }
# dev.off()
# 
# pdf("simulation_type_I_error_max_0.1.pdf", width=7.5, height=5.5)
# for (rtwo in r2_values) {
#   for (o in c("continuous", "binary_use_control_in_xQTL", "binary_use_all_in_xQTL")) {
#     for (p in c("horizontal", "vertical")) {
#       print(ggplot(subset(data_long_power, outcome == o & pleiotropy == p & r2 == rtwo)
#                    , aes(x = type_I_error, y = power)) +
#               geom_point(aes(col = method, shape = method), stroke = 0.2) +
#               facet_grid(overlap ~ IV_strength) +
#               geom_vline(xintercept=0.05, linetype="solid", color = "black", size = 0.2) +
#               geom_rect(xmin=0.05, xmax=2, ymin=-1, ymax=2, fill="grey80", color=NA, alpha=0.05) +
#               scale_x_continuous(limits = c(0, 0.1)) +
#               scale_shape_manual(values=1:nlevels(data_long_power$method)) +
#               scale_color_manual(values=wes_palette(name="Darjeeling1", n = length(levels(data_long$method)), type = "continuous")) +
#               theme_bw() +
#               theme(axis.line = element_line(colour = "black"),
#                     panel.grid.major = element_blank(),
#                     panel.grid.minor = element_blank(),
#                     panel.background = element_blank(),
#                     plot.margin = margin(t = .35, r = .1, b = 0, l = .1, unit = "in"),
#                     # legend.position = "bottom",
#                     legend.text=element_text(size=8),
#                     axis.text.x = element_text(size=8, angle=90),
#                     axis.text.y = element_text(size=8)) +
#               guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
#               ggtitle(sprintf("outcome type = %s, pleiotropy = %s, r2 = %.2f
#                               IV strength: {0.5, 1, 2}, overlap: {0, 0.5, 1}", o, p, rtwo)))
#     }
#   }
# }
# dev.off()
# 
# pdf("simulation_outcome_type_I_error_max_0.1.pdf", width=7.5, height=5.5)
# for (rtwo in r2_values) {
#   for (s in c(0.5, 1, 2)) {
#     for (p in c("horizontal", "vertical")) {
#       print(ggplot(subset(data_long_power, IV_strength == s & pleiotropy == p & r2 == rtwo)
#                    , aes(x = type_I_error, y = power)) +
#               geom_point(aes(col = method, shape = method), stroke = 0.2) +
#               facet_grid(overlap ~ outcome) +
#               geom_vline(xintercept=0.05, linetype="solid", color = "black", size = 0.2) +
#               geom_rect(xmin=0.05, xmax=2, ymin=-1, ymax=2, fill="grey80", color=NA, alpha=0.05) +
#               scale_x_continuous(limits = c(0, 0.1)) +
#               scale_shape_manual(values=1:nlevels(data_long_power$method)) +
#               scale_color_manual(values=wes_palette(name="Darjeeling1", n = length(levels(data_long$method)), type = "continuous")) +
#               theme_bw() +
#               theme(axis.line = element_line(colour = "black"),
#                     panel.grid.major = element_blank(),
#                     panel.grid.minor = element_blank(),
#                     panel.background = element_blank(),
#                     plot.margin = margin(t = .35, r = .1, b = 0, l = .1, unit = "in"),
#                     # legend.position = "bottom",
#                     legend.text=element_text(size=8),
#                     axis.text.x = element_text(size=8, angle=90),
#                     axis.text.y = element_text(size=8)) +
#               guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
#               ggtitle(sprintf("IV strength = %s, pleiotropy = %s, r2 = %.2f
#                               outcome type: x axis, overlap: {0, 0.5, 1}", s, p, rtwo)))
#     }
#   }
# }
# dev.off()
# 
# 
# pdf("simulation_pleiotropy_type_I_error_max_0.1.pdf", width=7.5, height=5.5)
# for (rtwo in c(0, 0.01, 0.1)) {
#   for (s in c(0.5, 1, 2)) {
#     for (o in c("continuous", "binary_use_control_in_xQTL", "binary_use_all_in_xQTL")) {
#       print(ggplot(subset(data_long_power, IV_strength == s & outcome == o & r2 == rtwo)
#                    , aes(x = type_I_error, y = power)) +
#               geom_point(aes(col = method, shape = method), stroke = 0.2) +
#               facet_grid(overlap ~ pleiotropy) +
#               geom_vline(xintercept=0.05, linetype="solid", color = "black", size = 0.2) +
#               geom_rect(xmin=0.05, xmax=2, ymin=-1, ymax=2, fill="grey80", color=NA, alpha=0.05) +
#               scale_x_continuous(limits = c(0, 0.1)) +
#               scale_shape_manual(values=1:nlevels(data_long_power$method)) +
#               scale_color_manual(values=wes_palette(name="Darjeeling1", n = length(levels(data_long$method)), type = "continuous")) +
#               theme_bw() +
#               theme(axis.line = element_line(colour = "black"),
#                     panel.grid.major = element_blank(),
#                     panel.grid.minor = element_blank(),
#                     panel.background = element_blank(),
#                     plot.margin = margin(t = .35, r = .1, b = 0, l = .1, unit = "in"),
#                     # legend.position = "bottom",
#                     legend.text=element_text(size=8),
#                     axis.text.x = element_text(size=8, angle=90),
#                     axis.text.y = element_text(size=8)) +
#               guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
#               ggtitle(sprintf("IV strength = %s, outcome type = %s, r2 = %.2f
#                               pleiotropy: x axis, overlap: {0, 0.5, 1}", s, o, rtwo)))
#     }
#   }
# }
# dev.off()
