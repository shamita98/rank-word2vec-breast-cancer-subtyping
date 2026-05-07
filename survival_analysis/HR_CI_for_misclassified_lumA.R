#---------------------------------------------------------------------------------------------------------------
# Title: hazard ratios (HR) and 95% confidence intervals (CI) computation
# Pairwise comparisons between misclassified luminal A cases and true positive luminal A and luminal B cases
#---------------------------------------------------------------------------------------------------------------

library(survminer)
library(survival)
library(tidyverse)

#---------------------------------------------------------------------------------------------------------------
# load files
#---------------------------------------------------------------------------------------------------------------
counts_res = read.csv("C:/Users/User/Documents/master_thesis_project_analysis/pam50_classification_counts/counts_model_evaluation/counts_prediction_results/tcga_brca_counts_rf_predictions.csv",
                  header = 1, row.names = 1)
rank_res = read.csv("C:/Users/User/Documents/master_thesis_project_analysis/pam50_classification_ranks/ranks_model_evaluation/ranks_prediction_results/tcga_brca_ranks_svm_predictions.csv",
                  header = 1, row.names = 1)
word2vec_res = read.csv("C:/Users/User/Documents/master_thesis_project_analysis/pam50_classification_word2vec/word2vec_model_evaluation/word2vec_prediction_results/tcga_brca_word2vec_svm_predictions.csv",
                  header = 1, row.names = 1)

#---------------------------------------------------------------------------------------------------------------
# remove a 'Basal' sample with negative overall survival in all files
# row index: "TCGA-PL-A8LV-01A-21R-A41B-07"
counts_res = counts_res[rownames(counts_res) != "TCGA-PL-A8LV-01A-21R-A41B-07", ]
rank_res = rank_res[rownames(rank_res) != "TCGA-PL-A8LV-01A-21R-A41B-07", ]
word2vec_res = word2vec_res[rownames(word2vec_res) != "TCGA-PL-A8LV-01A-21R-A41B-07", ]

# filter to relevant groups
filter_luminal <- function(df) {
  # keep only LumA__LumA, LumA__LumB, LumB__LumB in actual_pred column
  df_filtered = df[df$actual_pred %in% c("LumA__LumA", "LumA__LumB", "LumB__LumB"), ]
  
  # remove patients with survival years greater than 10
  # because that is when the plateau starts
  # 39 patients (29 alive and 10 dead were removed)
  # distribution: 8 Basal, 2 Her2, 24 LumA, 5 LumB
  df_filtered = df_filtered[df_filtered$overall_survival_years <= 10, ]
  
  # create survival event column
  df_filtered$event = ifelse(df_filtered$vital_status == "Alive", 0, 1)
  return(df_filtered)
}

counts_res_filtered = filter_luminal(counts_res)
rank_res_filtered = filter_luminal(rank_res)
word2vec_res_filtered = filter_luminal(word2vec_res)

#---------------------------------------------------------------------------------------------------------------
# Function to calculate HRs using Cox model with broom::tidy()
#---------------------------------------------------------------------------------------------------------------
calculate_HR = function(df, reference_group) {
  # Ensure factor levels: reference group first
  df$actual_pred = factor(df$actual_pred, 
                           levels = c(reference_group, 
                                      setdiff(unique(df$actual_pred), reference_group)))
  
  # Survival object
  surv_obj = Surv(time = df$overall_survival_years, event = df$event)
  
  # Fit Cox model
  cox_model = coxph(surv_obj ~ actual_pred, data = df)
  
  # Tidy output with HR, CI, p-value
  hr_table = broom::tidy(cox_model, exponentiate = TRUE, conf.int = TRUE) %>%
    select(term, estimate, conf.low, conf.high, p.value) %>%
    rename(Comparison = term,
           HR = estimate,
           CI_lower = conf.low,
           CI_upper = conf.high,
           p_value = p.value)
  
  return(hr_table)
}

#---------------------------------------------------------------------------------------------------------------
# Calculate HRs for all three approaches vs Luminal A
#---------------------------------------------------------------------------------------------------------------
hr_counts_vs_lumA = calculate_HR(counts_res_filtered, "LumA__LumA")
hr_rank_vs_lumA = calculate_HR(rank_res_filtered, "LumA__LumA")
hr_word2vec_vs_lumA = calculate_HR(word2vec_res_filtered, "LumA__LumA")

#---------------------------------------------------------------------------------------------------------------
# Calculate HRs for all three approaches vs Luminal B
#---------------------------------------------------------------------------------------------------------------
hr_counts_vs_lumB = calculate_HR(counts_res_filtered, "LumB__LumB")
hr_rank_vs_lumB = calculate_HR(rank_res_filtered, "LumB__LumB")
hr_word2vec_vs_lumB = calculate_HR(word2vec_res_filtered, "LumB__LumB")

#---------------------------------------------------------------------------------------------------------------
# Combine into one tidy table
#---------------------------------------------------------------------------------------------------------------
hr_summary = bind_rows(
  hr_counts_vs_lumA %>% mutate(Method = "Counts-RF", Reference = "Luminal A"),
  hr_rank_vs_lumA %>% mutate(Method = "Rank-SVM", Reference = "Luminal A"),
  hr_word2vec_vs_lumA %>% mutate(Method = "Word2Vec-SVM", Reference = "Luminal A"),
  
  hr_counts_vs_lumB %>% mutate(Method = "Counts-RF", Reference = "Luminal B"),
  hr_rank_vs_lumB %>% mutate(Method = "Rank-SVM", Reference = "Luminal B"),
  hr_word2vec_vs_lumB %>% mutate(Method = "Word2Vec-SVM", Reference = "Luminal B")
) %>%
  select(Method, Reference, Comparison, HR, CI_lower, CI_upper, p_value)

# Print final summary table
print(hr_summary)
