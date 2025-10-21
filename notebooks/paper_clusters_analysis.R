library(tidyverse)
library(cowplot)
library(ggsankey)

# To run this script, open the imoc_pdac.Rproj file in the main folder (imoc_pdac)


papers_cluster_df <- read_csv("data/TCGA/comparison_papers/clusters_papers.csv")
papers_cluster_df$`My Clusters`[papers_cluster_df$`My Clusters` == 1] <- "Cluster 1"
papers_cluster_df$`My Clusters`[papers_cluster_df$`My Clusters` == 2] <- "Cluster 2"

sankey_bailey <- papers_cluster_df %>% make_long(`My Clusters`, bailey)

p_bailey <- ggplot(sankey_bailey, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  geom_sankey(flow.alpha = 0.33) +
  scale_fill_manual(values = c("#F8766D", "#00BFC4", rep("gray50",4))) +
  geom_sankey_label(size = 4, color = "white", alpha = 0.75) +
  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5)) +
  ggtitle("Bailey subtypes")
p_bailey


sankey_collisson <- papers_cluster_df %>% make_long(`My Clusters`, collisson)

p_collisson <- ggplot(sankey_collisson, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  geom_sankey(flow.alpha = 0.33) +
  scale_fill_manual(values = c("#F8766D", "#00BFC4",rep("gray50",3))) +
  geom_sankey_label(size = 4, color = "white", alpha = 0.75) +
  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5)) +
  ggtitle("Collisson subtypes")
p_collisson

sankey_moffitt <- papers_cluster_df %>% make_long(`My Clusters`, moffitt)

p_moffitt <- ggplot(sankey_moffitt, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  geom_sankey(flow.alpha = 0.33) +
  scale_fill_manual(values = c("#F8766D", "#00BFC4", rep("gray50",2))) +
  geom_sankey_label(size = 4, color = "white", alpha = 0.75) +
  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5)) +
  ggtitle("Moffitt subtypes")
p_moffitt


plot_grid(p_bailey, p_collisson, p_moffitt, nrow = 1)
#ggsave("figures/pdac_subtypes.pdf", width = 12, height = 4)

bailey_tbl <- table(papers_cluster_df$`My Clusters`, papers_cluster_df$bailey)

rcompanion::pairwiseNominalIndependence(
  bailey_tbl,
  compare = "column",
  fisher = T,
  gtest = F,
  chisq = F,
  method = "fdr",
  cramer = TRUE
)

collisson_tbl <- table(papers_cluster_df$`My Clusters`, papers_cluster_df$collisson)

rcompanion::pairwiseNominalIndependence(
  collisson_tbl,
  compare = "column",
  fisher = T,
  gtest = F,
  chisq = F,
  method = "fdr",
  cramer = TRUE
)


