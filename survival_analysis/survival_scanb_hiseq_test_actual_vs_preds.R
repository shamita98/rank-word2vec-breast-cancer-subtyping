#--------------------------------------------------------------------
# Title: Survival Analysis of SCAN-B HiSeq Test Set (20% Held-Out Test Set)
# Comparison between Actual Subtype Labels and Model-Predicted Subtype labels
#--------------------------------------------------------------------

library(survminer)
library(survival)
library(tidyverse)

#--------------------------------------------------------------------
# load files
#--------------------------------------------------------------------
counts_res = read.csv("C:/Users/User/Documents/master_thesis_project_analysis/pam50_classification_counts/counts_model_evaluation/counts_prediction_results/scanb_hiseq_test_counts_rf_predictions.csv",
                  header = 1, row.names = 1)
rank_res = read.csv("C:/Users/User/Documents/master_thesis_project_analysis/pam50_classification_ranks/ranks_model_evaluation/ranks_prediction_results/scanb_hiseq_test_ranks_svm_predictions.csv",
                  header = 1, row.names = 1)
word2vec_res = read.csv("C:/Users/User/Documents/master_thesis_project_analysis/pam50_classification_word2vec/word2vec_model_evaluation/word2vec_prediction_results/scanb_hiseq_test_word2vec_svm_predictions.csv",
                  header = 1, row.names = 1)


#--------------------------------------------------------------------
# function to create a survival plot
#--------------------------------------------------------------------
#--------------------------------------------------------------------
# function for creating survival plots based on actual subtype labels
#--------------------------------------------------------------------
#--------------------------------------------------------------------
actual_survivalPlot = function(df, survival_event, 
                        alive_indicator, plot_title)
{
  df$deceased = ifelse(df[[survival_event]] == alive_indicator, FALSE, TRUE)
  survival_fit = survfit(Surv(overall_survival_years, deceased) ~ subtype, 
                         data = df)
  print('Proportion of deceased and alive:')
  print(table(df$deceased))
  
  plot = ggsurvplot(survival_fit,
                    ggtheme = theme_bw(),
                    data = df,
                    break.time.by = 1,
                    pval = TRUE,
                    title = plot_title,
                    pval.method = TRUE,
                    pval.size = 7,
                    pval.coord = c(0,0.05),
                    pval.method.coord = c(0,0.12),
                    censor = T,
                    censor.size = 3,
                    legend.labs = c('Basal', 'HER2', 'LumA', 'LumB'),
                    ylab = "Overall survival probability")
  
  plot$plot <- plot$plot +
    labs(title = paste(plot_title), fill = NULL, color = 'Subtypes',
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
          legend.position = c(0.99,0.27),
          legend.box.margin = margin(t = 0.6*(.pt*72.27/96/4),
                                     r = 0.6*(.pt*72.27/96/4)),
          axis.ticks = element_line(color="black"),
          axis.text = element_text(color="black", size=18),
          axis.title = element_text(size = 19, color="black"),
          plot.background = element_rect(fill="white"),
          plot.tag.position = c(0,1)) +
    geom_vline(xintercept = 5, col = "black", lwd = 0.6, lty = 2) 
  
  print('Plot is generated.')
  
  return(plot)
                    
}

#--------------------------------------------------------------------
#--------------------------------------------------------------------
# function for creating survival plots based on model-predicted subtype labels
#--------------------------------------------------------------------
#--------------------------------------------------------------------

predicted_survivalPlot = function(df, survival_event, 
                               alive_indicator, plot_title)
{
  df$deceased = ifelse(df[[survival_event]] == alive_indicator, FALSE, TRUE)
  survival_fit = survfit(Surv(overall_survival_years, deceased) ~ predicted_subtype, 
                         data = df)
  print('Proportion of deceased and alive:')
  print(table(df$deceased))
  
  plot = ggsurvplot(survival_fit,
                    ggtheme = theme_bw(),
                    data = df,
                    break.time.by = 1,
                    pval = TRUE,
                    title = plot_title,
                    pval.method = TRUE,
                    pval.size = 7,
                    pval.coord = c(0,0.05),
                    pval.method.coord = c(0,0.12),
                    censor = T,
                    censor.size = 3,
                    legend.labs = c('Basal', 'HER2', 'LumA', 'LumB'),
                    ylab = "Overall survival probability")
  
  plot$plot <- plot$plot +
    labs(title = paste(plot_title), fill = NULL, color = 'Subtypes',
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
          legend.position = c(0.99,0.27),
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

#--------------------------------------------------------------------
# Survival plots using actual labels
#--------------------------------------------------------------------
actual_label_plot = actual_survivalPlot(df = counts_res, 
                                 survival_event = "overall_survival_event", 
                                 alive_indicator = 0, 
                                 plot_title = "SCAN-B HiSeq Test Set \n with Actual Subtype Labels (n=551)")

print(actual_label_plot)
ggsave(actual_label_plot$plot, filename = "scanb_hiseq_test_actual_labels_survival.png", 
       dpi = 500, height = 7, width = 8)


#--------------------------------------------------------------------
# survival plot using counts-rf labels
#--------------------------------------------------------------------
counts_rf_plot = predicted_survivalPlot(df = counts_res, 
                                        survival_event = "overall_survival_event", 
                                        alive_indicator = 0, 
                                        plot_title = "SCAN-B HiSeq Test Set \n with Counts-RF Predicted Labels (n=551)")
print(counts_rf_plot)

ggsave(counts_rf_plot$plot, filename = "scanb_hiseq_test_counts_rf_labels_survival.png", 
       dpi = 300, height = 7, width = 8)

#--------------------------------------------------------------------
# survival plot using rank-svm labels
#--------------------------------------------------------------------
rank_svm_plot = predicted_survivalPlot(df = rank_res, 
                                            survival_event = "overall_survival_event", 
                                            alive_indicator = 0, 
                                            plot_title = "SCAN-B HiSeq Test Set \n with Rank-SVM Predicted Labels (n=551)")
print(rank_svm_plot)

ggsave(rank_svm_plot$plot, filename = "scanb_hiseq_test_rank_svm_labels_survival.png", 
       dpi = 300, height = 7, width = 8)

#--------------------------------------------------------------------
# survival plot using word2vec-svm labels
#--------------------------------------------------------------------
word2vec_svm_plot = predicted_survivalPlot(df = word2vec_res, 
                                       survival_event = "overall_survival_event", 
                                       alive_indicator = 0, 
                                       plot_title = "SCAN-B HiSeq Test Set \n with Word2vec-SVM Predicted Labels (n=551)")
print(word2vec_svm_plot)

ggsave(word2vec_svm_plot$plot, filename = "scanb_hiseq_test_word2vec_svm_labels_survival.png", 
       dpi = 300, height = 7, width = 8)


