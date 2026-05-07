# import libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# gene clusters
gene_cluster_file = read.csv('hierarchical_clusters_pam50_word2vec_embedding.csv')

# get clusters
cluster_1 = gene_cluster_file[gene_cluster_file$hier_cluster==1,'new_gene_symbol']
cluster_2 = gene_cluster_file[gene_cluster_file$hier_cluster==2,'new_gene_symbol']
cluster_3 = gene_cluster_file[gene_cluster_file$hier_cluster==3,'new_gene_symbol']
cluster_4 = gene_cluster_file[gene_cluster_file$hier_cluster==4,'new_gene_symbol']
cluster_5 = gene_cluster_file[gene_cluster_file$hier_cluster==5,'new_gene_symbol']
cluster_6 = gene_cluster_file[gene_cluster_file$hier_cluster==6,'new_gene_symbol']

## GO Biological Process Enrichment Analysis

#------------------------------------------------------------------------------

# cluster 1
go_bp_term_cluster_1 = enrichGO(gene=cluster_1, OrgDb = "org.Hs.eg.db",
                        keyType = "SYMBOL",pvalueCutoff = 0.05, 
                        pAdjustMethod = "BH",
                        qvalueCutoff = 0.05,ont = "BP")

simplified_cluster_1 = simplify(go_bp_term_cluster_1, cutoff = 0.7, by = "p.adjust", select_fun  = min)  

go_bp_term_cluster_1_df = as.data.frame((go_bp_term_cluster_1))
simplified_cluster_1_df = as.data.frame((simplified_cluster_1))

#------------------------------------------------------------------------------

# cluster 2
go_bp_term_cluster_2 = enrichGO(gene = cluster_2, OrgDb = "org.Hs.eg.db",
                                keyType = "SYMBOL",pvalueCutoff = 0.05, 
                                pAdjustMethod = "BH",
                                qvalueCutoff = 0.05,ont = "BP")

go_bp_term_cluster_2_df = as.data.frame((go_bp_term_cluster_2))

#------------------------------------------------------------------------------

# cluster 3
go_bp_term_cluster_3 = enrichGO(gene = cluster_3,OrgDb = "org.Hs.eg.db",
                                keyType = "SYMBOL",pvalueCutoff = 0.05, 
                                pAdjustMethod = "BH",
                                qvalueCutoff = 0.05,ont = "BP")

go_bp_term_cluster_3_df = as.data.frame((go_bp_term_cluster_3))

#------------------------------------------------------------------------------

# cluster 4
go_bp_term_cluster_4 = enrichGO(gene = cluster_4,OrgDb = "org.Hs.eg.db",
                                keyType = "SYMBOL",pvalueCutoff = 0.05, 
                                pAdjustMethod = "BH",
                                qvalueCutoff = 0.05,ont = "BP")

go_bp_term_cluster_4_df = as.data.frame((go_bp_term_cluster_4))

#------------------------------------------------------------------------------

# cluster 5
go_bp_term_cluster_5 = enrichGO(gene = cluster_5,OrgDb = "org.Hs.eg.db",
                                keyType = "SYMBOL",pvalueCutoff = 0.05, 
                                pAdjustMethod = "BH",
                                qvalueCutoff = 0.05,ont = "BP")

go_bp_term_cluster_5_df = as.data.frame((go_bp_term_cluster_5))

simplified_cluster_5 = simplify(go_bp_term_cluster_5, cutoff = 0.5, by = "p.adjust", select_fun = min, 
                                measure = "Wang", semData = NULL)

simplified_cluster_5_df = as.data.frame((simplified_cluster_5))

#------------------------------------------------------------------------------

# cluster 6
go_bp_term_cluster_6 = enrichGO(gene = cluster_6,OrgDb = "org.Hs.eg.db",
                                keyType = "SYMBOL",pvalueCutoff = 0.05, 
                                pAdjustMethod = "BH",
                                qvalueCutoff = 0.05,ont = "BP")

go_bp_term_cluster_6_df = as.data.frame((go_bp_term_cluster_6))

#------------------------------------------------------------------------------

# # save files
# write.csv(go_bp_term_cluster_1_df, 'cluster_1_go_bp.csv')
# write.csv(go_bp_term_cluster_2_df, 'cluster_2_go_bp.csv')
# write.csv(go_bp_term_cluster_3_df, 'cluster_3_go_bp.csv')
# write.csv(go_bp_term_cluster_4_df, 'cluster_4_go_bp.csv')
# write.csv(go_bp_term_cluster_5_df, 'cluster_5_go_bp.csv')
# write.csv(go_bp_term_cluster_6_df, 'cluster_6_go_bp.csv')


