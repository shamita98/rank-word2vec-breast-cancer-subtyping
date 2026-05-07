#--------------------------------------------------------------------
# Title: Survival Analysis of SCAN-B HiSeq Training Set using Actual Subtype Labels
#--------------------------------------------------------------------

library(survminer)
library(survival)
library(tidyverse)

#--------------------------------------------------------------------
# load files
#--------------------------------------------------------------------
scanb_train = read.csv("C:/Users/User/Documents/master_thesis_project_analysis/datasets/SCANB_GSE202203/scanb_hiseq_train_test_sets/train_test_80_20/scanb_hiseq_train80_survival_subtype.csv",
                  header = 1, row.names = 1)

#--------------------------------------------------------------------
# function to create a survival plot
#--------------------------------------------------------------------
#--------------------------------------------------------------------
# function to create a survival plot using actual subtype labels
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
          legend.position = c(0.99,0.25),
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
# survival plot using actual subtype labels
#--------------------------------------------------------------------
actual_label_plot = actual_survivalPlot(df = scanb_train, 
                                 survival_event = "overall_survival_event", 
                                 alive_indicator = 0, 
                                 plot_title = "SCAN-B HiSeq Training Set (n=2203)")

ggsave(actual_label_plot$plot, filename = "scanb_hiseq_train_actual_labels_survival.png", 
       dpi = 300, height = 7, width = 8)
print(actual_label_plot)




