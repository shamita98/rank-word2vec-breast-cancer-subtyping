library(survminer)
library(survival)
library(tidyverse)

# load files
counts_res = read.csv("C:/Users/User/Documents/master_thesis_project_analysis/pam50_classification_counts/counts_model_evaluation/counts_prediction_results/tcga_brca_counts_rf_predictions.csv",
                  header = 1, row.names = 1)
rank_res = read.csv("C:/Users/User/Documents/master_thesis_project_analysis/pam50_classification_ranks/ranks_model_evaluation/ranks_prediction_results/tcga_brca_ranks_svm_predictions.csv",
                  header = 1, row.names = 1)
word2vec_res = read.csv("C:/Users/User/Documents/master_thesis_project_analysis/pam50_classification_word2vec/word2vec_model_evaluation/word2vec_prediction_results/tcga_brca_word2vec_svm_predictions.csv",
                  header = 1, row.names = 1)

# remove a 'Basal' sample with negative overall survival in all files
# row index: "TCGA-PL-A8LV-01A-21R-A41B-07"
counts_res = counts_res[rownames(counts_res) != "TCGA-PL-A8LV-01A-21R-A41B-07", ]
rank_res = rank_res[rownames(rank_res) != "TCGA-PL-A8LV-01A-21R-A41B-07", ]
word2vec_res = word2vec_res[rownames(word2vec_res) != "TCGA-PL-A8LV-01A-21R-A41B-07", ]

# keep only LumA__LumA, LumA__LumB, LumB__LumB in actual_pred column
counts_res = counts_res[counts_res$actual_pred %in% c("LumA__LumA", "LumA__LumB", "LumB__LumB"), ]
rank_res = rank_res[rank_res$actual_pred %in% c("LumA__LumA", "LumA__LumB", "LumB__LumB"), ]
word2vec_res = word2vec_res[word2vec_res$actual_pred %in% c("LumA__LumA", "LumA__LumB", "LumB__LumB"), ]

# remove patients with survival years greater than 12.5, because that is when the plateau starts
# 39 patients (29 alive and 10 dead were removed)
# distribution: 8 Basal, 2 Her2, 24 LumA, 5 LumB
counts_res_filtered = counts_res[counts_res$overall_survival_years<=10, ]
rank_res_filtered = rank_res[rank_res$overall_survival_years<=10, ]
word2vec_res_filtered = word2vec_res[word2vec_res$overall_survival_years<=10, ]

# function to create a survival plot
predicted_survivalPlot = function(df, survival_event, 
                               alive_indicator, plot_title)
{
  df$deceased = ifelse(df[[survival_event]] == alive_indicator, FALSE, TRUE)
  survival_fit = survfit(Surv(overall_survival_years, deceased) ~ actual_pred, 
                         data = df)
  print('Proportion of deceased and alive:')
  print(table(df$deceased))
  
  plot = ggsurvplot(survival_fit,
                    ggtheme = theme_bw(),
                    data = df,
                    break.time.by = 1,
                    pval = TRUE,
                    # risk.table = TRUE,
                    title = plot_title,
                    pval.method = TRUE,
                    pval.size = 7,
                    pval.coord = c(0,0.05),
                    pval.method.coord = c(0,0.12),
                    censor = T,
                    censor.size = 3,
                    legend.labs = c('LumA/LumA', 'LumA/LumB', 'LumB/LumB'),
                    xlab = "Survival time (years)",
                    ylab = "Overall survival probability")
  
  plot$plot <- plot$plot +
    labs(title = paste(plot_title), fill = NULL, color = "Actual/Predicted",
         x = "Follow-up time (years)") +
    scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, by=0.1), labels = scales::percent_format())+
    guides(color = guide_legend(override.aes = list(shape = NA,
                                                    linewidth=1,
                                                    size=5)))+
    theme(panel.border = element_rect(color="black", fill=NA,
                                      linewidth=0.6),
          plot.title = element_text(size=20, hjust = 0.5, face = 'bold'),
          legend.justification = c(1,1),
          legend.key.size = unit(5,"mm"),
          legend.spacing.y = unit(0.5,"mm"),
          legend.background = element_rect(color="black", fill=NA,
                                           linewidth = 0.6*0.5),
          legend.title = element_text(size=18, hjust=0.5),
          legend.text = element_text(size=18),
          legend.position = c(0.99,0.22),
          legend.box.margin = margin(t = 0.6*(.pt*72.27/96/4),
                                     r = 0.6*(.pt*72.27/96/4)),
          axis.ticks = element_line(color="black"),
          axis.text = element_text(color="black", size=18),
          axis.title = element_text(size = 19, color="black"),
          plot.background = element_rect(fill="white"),
          plot.tag.position = c(0,1))+
    geom_vline(xintercept = 5, col = "black", lwd = 0.6, lty = 2)

  print('Plot is generated.')
  
  return(plot)
  
}


# survival plot using counts-logreg labels
counts_rf_plot = predicted_survivalPlot(df = counts_res_filtered, 
                                        survival_event = "vital_status", 
                                        alive_indicator = "Alive", 
                                        plot_title = "Predicted Luminal Categories \n in TCGA-BRCA (Counts-RF)")
print(counts_rf_plot)

ggsave(counts_rf_plot$plot, filename = "tcga_brca_luminal_counts_rf_labels_survival.png", 
       dpi = 300, height = 7, width = 8)


# survival plot using rank-svm labels
rank_svm_plot = predicted_survivalPlot(df = rank_res_filtered, 
                                            survival_event = "vital_status", 
                                            alive_indicator = "Alive", 
                                            plot_title = "Predicted Luminal Categories \n in TCGA-BRCA (Rank-SVM)")
print(rank_svm_plot)

ggsave(rank_svm_plot$plot, filename = "tcga_brca_luminal_rank_svm_labels_survival.png", 
       dpi = 300, height = 7, width = 8)


# survival plot using word2vecpca-svm labels
word2vec_svm_plot = predicted_survivalPlot(df = word2vec_res_filtered, 
                                       survival_event = "vital_status", 
                                       alive_indicator = "Alive", 
                                       plot_title = "Predicted Luminal Categories \n in TCGA-BRCA (Word2vec-SVM)")
print(word2vec_svm_plot)

ggsave(word2vec_svm_plot$plot, filename = "tcga_brca_luminal_word2vec_svm_labels_survival.png", 
       dpi = 300, height = 7, width = 8)


