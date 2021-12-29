library(qiime2R)
library(ggplot2)
library(decontam)
library(phyloseq)
library(tidyverse)
library(microbiome)
library(ANCOMBC)
library(patchwork)
library(data.table)
library(readxl)
library(gghalves)
library(ggpubr)
library(stats)
library(rstatix)
##
physeqtocontam <- qza_to_phyloseq(
  features = "/home/irshad/Schilling_Project_075/processed/shi7_learning/fastqs_stitched/fastafiles/dada2-single-end-sepp-filtered-table.qza",
  tree="/home/irshad/Schilling_Project_075/processed/shi7_learning/fastqs_stitched/fastafiles/insertion-tree-silva-128.qza",
  "/home/irshad/Schilling_Project_075/processed/shi7_learning/fastqs_stitched/fastafiles/taxonomy-silva-138-515-806-nb.qza", 
  metadata = "/home/irshad/Schilling_Project_075/processed/shi7_learning/fastqs_stitched/fastafiles/metadata-decontam.tsv"
)

View(physeqtocontam)


### To identify contaminants using prevalence method
sample_data(physeqtocontam)$is.neg <- sample_data(physeqtocontam)$Sample_or_Control == "Control Sample"
contamdf.prev <- isContaminant(physeqtocontam, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev$contaminant)
View(contamdf.prev)
head(which(contamdf.prev$contaminant))

### To make phyloseq object of presence-absence in negative controls and true samples
physeqtocontam.prev <- transform_sample_counts(physeqtocontam, function(abund) 1*(abund>0))
physeqtocontam.prev.neg <- prune_samples(sample_data(physeqtocontam.prev)$Sample_or_Control == "Control Sample", physeqtocontam.prev)
physeqtocontam.prev.pos <- prune_samples(sample_data(physeqtocontam.prev)$Sample_or_Control == "True Sample", physeqtocontam.prev)

### To make data.frame of prevalence in positive and negative samples
df.physeqtocontam <- data.frame(physeqtocontam.prev.pos=taxa_sums(physeqtocontam.prev.pos), physeqtocontam.prev.neg=taxa_sums(physeqtocontam.prev.neg), contaminant=contamdf.prev$contaminant)
ggplot(df.physeqtocontam, aes(x=physeqtocontam.prev.neg, y=physeqtocontam.prev.pos, color=contaminant)) + geom_point()+
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
###

### Beta diversity unweighted unifrac PCoA figure_02

setwd("/home/irshad/Schilling_Project_075/processed/shi7_learning/fastqs_stitched/fastafiles/") # to set working directory

metadata <-read_q2metadata("metadata_new_fig_3.tsv")
uwunifrac<-read_qza("core-metrics-results-1189/unweighted_unifrac_pcoa_results.qza")
shannon<-read_qza("core-metrics-results-1189/shannon_vector.qza")$data %>% rownames_to_column("SampleID") 

unwe_2 <- uwunifrac$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>%
  left_join(shannon) %>%
  ggplot(aes(x=PC1, y=PC2, color=`Category`, size=shannon_entropy)) +
  geom_point(alpha=0.4) +
  scale_shape_manual(values=c(16,1)) + 
  scale_size_continuous(name="Shannon diversity") +
  scale_color_discrete(name="Category")+
  theme_q2r()+
  theme(
    axis.title.y = element_text(size = 8, face="bold"),
    axis.text = element_text(size=6, face = "bold"),
    axis.title.x = element_text(size = 8, face="bold"),
    legend.position = "none")

unwe_2

### Beta diversity weighted unifrac PCoA figure_02

metadata <-read_q2metadata("metadata_new_fig_3.tsv")
wunifrac<-read_qza("core-metrics-results-1189/weighted_unifrac_pcoa_results.qza")
shannon<-read_qza("core-metrics-results-1189/shannon_vector.qza")$data %>% rownames_to_column("SampleID") 

we_2 <- wunifrac$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>%
  left_join(shannon) %>%
  ggplot(aes(x=PC1, y=PC2, color=`Category`, size=shannon_entropy)) +
  geom_point(alpha=0.4) + 
  scale_shape_manual(values=c(16,1)) + 
  scale_size_continuous(name="Shannon diversity") +
  scale_color_discrete(name="Category")+
  theme_q2r()+
  theme(
    axis.title.y = element_text(size = 8, face="bold"),
    axis.text = element_text(size=6, face = "bold"),
    axis.title.x = element_text(size = 8, face="bold"),
    legend.text = element_text(face = "italic", size = 6),
    legend.title = element_text(size = 8),
    legend.key.height = unit(0.1, "lines"))
    
we_2

fig_02 <- unwe_2 + we_2
fig_02 + plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face="bold"))
##
ggsave("figures/weighted_unweighted_UniFrac_Fig_02.tiff", width = 8, height = 3.5, compression="lzw", dpi = 600)


### Alpha diversity visualization Supplementary Figure_01
setwd("/home/irshad/Schilling_Project_075/processed/shi7_learning/fastqs_stitched/fastafiles/figures/New_raw_data_alpha_diversity/")

qiimealpha <- read_xlsx("alpha_diversity.xlsx")
qiimealpha
adp <- qiimealpha

p1 <- ggplot(adp, aes(x=category, y=Shannon_entropy, fill=category))+
  stat_summary(geom="bar", fun.data=mean_se, color="black")+
  geom_jitter(shape=21, width=0.2, height=0) +
  theme_q2r()+
  labs(x = "", y = "Shannon entropy") +
  theme(
    legend.title = element_blank(),
    axis.text = element_text(size=10, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12, face="bold"),
    axis.text.x = element_blank(),
    legend.text = element_text(face="italic", size = 12),
    legend.position = "none"
  )
p1
p2 <- ggplot(adp, aes(x=category, y=Faith_pd, fill=category))+
  stat_summary(geom="bar", fun.data=mean_se, color="black")+
  geom_jitter(shape=21, width=0.2, height=0) +
  theme_q2r()+
  labs(x = "", y = "Faith's pd") +
  theme(
    legend.title = element_blank(),
    axis.text = element_text(size=10, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12, face="bold"),
    axis.text.x = element_blank(),
    legend.text = element_text(face="italic", size = 12),
    legend.position = "none"
  )
p2

p3 <- ggplot(adp, aes(x=category, y=Observed_features, fill=category))+
  stat_summary(geom="bar", fun.data=mean_se, color="black") + 
  geom_jitter(shape=21, width=0.2, height=0) +
  theme_q2r()+
  labs(x = "", y = "Observed features") +
  theme(
    legend.title = element_blank(),
    axis.text = element_text(size=10, face = "bold"),
    axis.title.y = element_text(size = 12, face="bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    legend.text = element_text(face="italic", size = 12),
    legend.position = "bottom",
    legend.key.size = unit(0.4, 'cm')
  )
p3
p4 <- ggplot(adp, aes(x=category, y=Pielou_evenness, fill=category))+
  stat_summary(geom="bar", fun.data=mean_se, color="black") + 
  geom_jitter(shape=21, width=0.2, height=0) +
  theme_q2r()+
  labs(x = "", y = "Pielou's evenness") +
  theme(
    legend.title = element_blank(),
    axis.text = element_text(size=10, face = "bold"),
    axis.title.y = element_text(size = 12, face="bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    legend.text = element_text(face="italic", size = 12),
    legend.position = "bottom",
    legend.key.size = unit(0.4, 'cm')
  )
p4
###
combind <- (p1 + p2) / (p3 + p4)
combind + plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face="bold"))
##
ggsave("Supplementary_figure_1_alpha_diversity.tiff", width = 8.2, height = 9, compression="lzw", dpi = 600)

### Analysis of Compositions of Microbiomes with Bias Correction (ANCOM-BC)
setwd("/home/irshad/Schilling_Project_075/processed/shi7_learning/fastqs_stitched/fastafiles/")

##Produce a phyloseq object
physeq_ancombc<-qza_to_phyloseq(
  features="17-samples-table-no-mitochondria-no-chloroplast-no-contaminants-sepp-filtered.qza",
  tree="insertion-tree-silva-128.qza",
  "taxonomy-silva-138-515-806-nb.qza",
  metadata = "metadata_new.tsv"
)
physeq_ancombc
phylum_data = aggregate_taxa(physeq_ancombc, "Phylum")
ancom_bc <- ancombc(phyloseq=phylum_data, formula = "category",
                    p_adj_method = "bonferroni", zero_cut = 0.90, lib_cut = 1000, 
                    group = "category", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
                    max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE)

ancom_res_df <- data.frame(
  Species = row.names(ancom_bc$res$beta),
  beta = unlist(ancom_bc$res$beta),
  se = unlist(ancom_bc$res$se),
  W = unlist(ancom_bc$res$W),
  p_val = unlist(ancom_bc$res$p_val),
  q_val = unlist(ancom_bc$res$q_val),
  diff_abn = unlist(ancom_bc$res$diff_abn))

bonferroni_ancom <- ancom_res_df %>%
  dplyr::filter(q_val < 0.05)
dim(bonferroni_ancom)
bonferroni_ancom
write_csv(bonferroni_ancom, "ancombc_output_phylum_level.csv")
##
res = ancom_bc$res
##
samp_frac = ancom_bc$samp_frac
samp_frac[is.na(samp_frac)] = 0 # Replace NA with 0
log_obs_abn = log(abundances(phylum_data) + 1) # Add pesudo-count (1) to avoid taking the log of 0
log_obs_abn_adj = t(t(log_obs_abn) - samp_frac) # Adjust the log observed abundances
round(log_obs_abn_adj[, 1:17], 2) %>% 
  data.table(caption = "Bias-adjusted log observed abundances")
###
df_fig1 = data.frame(res$beta * res$diff_abn, check.names = FALSE) %>% 
  rownames_to_column("taxon_id")
df_fig2 = data.frame(res$se * res$diff_abn, check.names = FALSE) %>% 
  rownames_to_column("taxon_id")
colnames(df_fig2)[-1] = paste0(colnames(df_fig2)[-1], "SD")
df_fig = df_fig1 %>% left_join(df_fig2, by = "taxon_id") %>%
  transmute(taxon_id, categoryPb, categoryPbSD) %>%
  filter(categoryPb != 0) %>% arrange(desc(categoryPb)) %>%
  mutate(group = ifelse(categoryPb > 0, "g1", "g2"))
df_fig$taxon_id = factor(df_fig$taxon_id, levels = df_fig$taxon_id)

p1 = ggplot(data = df_fig, 
           aes(x = taxon_id, y = categoryPb, fill = group, color = group)) + 
  geom_bar(stat = "identity", width = 0.7, 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = categoryPb - categoryPbSD, ymax = categoryPb + categoryPbSD), width = 0.2,
                position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Fomitopsis betulina vs Fomes fomentarius \n log fold change") + 
  theme_bw() + 
  theme(legend.position = "none",
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(face = "italic", size = 6, colour = "black",angle = 70, hjust = 1),
        axis.title.y = element_text(face = "bold.italic", size = 6)
)
p1
ggsave("phylum_ancom_watefall.tiff", width = 3.5, height = 4.5, compression="lzw", dpi = 600)

###

class_data = aggregate_taxa(physeq_ancombc, "Class")
ancom_bc <- ancombc(phyloseq=class_data, formula = "category",
                    p_adj_method = "bonferroni", zero_cut = 0.90, lib_cut = 1000, 
                    group = "category", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
                    max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE)

ancom_res_df <- data.frame(
  Species = row.names(ancom_bc$res$beta),
  beta = unlist(ancom_bc$res$beta),
  se = unlist(ancom_bc$res$se),
  W = unlist(ancom_bc$res$W),
  p_val = unlist(ancom_bc$res$p_val),
  q_val = unlist(ancom_bc$res$q_val),
  diff_abn = unlist(ancom_bc$res$diff_abn))

bonferroni_ancom <- ancom_res_df %>%
  dplyr::filter(q_val < 0.05)
dim(bonferroni_ancom)
bonferroni_ancom
write_csv(bonferroni_ancom, "ancombc_output_Class_level.csv")
###

res = ancom_bc$res
##
samp_frac = ancom_bc$samp_frac
samp_frac[is.na(samp_frac)] = 0 # Replace NA with 0
log_obs_abn = log(abundances(phylum_data) + 1) # Add pesudo-count (1) to avoid taking the log of 0
log_obs_abn_adj = t(t(log_obs_abn) - samp_frac) # Adjust the log observed abundances
round(log_obs_abn_adj[, 1:17], 2) %>% 
  data.table(caption = "Bias-adjusted log observed abundances")
###
df_fig1 = data.frame(res$beta * res$diff_abn, check.names = FALSE) %>% 
  rownames_to_column("taxon_id")
df_fig2 = data.frame(res$se * res$diff_abn, check.names = FALSE) %>% 
  rownames_to_column("taxon_id")
colnames(df_fig2)[-1] = paste0(colnames(df_fig2)[-1], "SD")
df_fig = df_fig1 %>% left_join(df_fig2, by = "taxon_id") %>%
  transmute(taxon_id, categoryPb, categoryPbSD) %>%
  filter(categoryPb != 0) %>% arrange(desc(categoryPb)) %>%
  mutate(group = ifelse(categoryPb > 0, "g1", "g2"))
df_fig$taxon_id = factor(df_fig$taxon_id, levels = df_fig$taxon_id)

p2 = ggplot(data = df_fig, 
           aes(x = taxon_id, y = categoryPb, fill = group, color = group)) + 
  geom_bar(stat = "identity", width = 0.7, 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = categoryPb - categoryPbSD, ymax = categoryPb + categoryPbSD), width = 0.2,
                position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Fomitopsis betulina vs Fomes fomentarius \n log fold change") + 
  theme_bw() + 
  theme(legend.position = "none",
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(face = "italic", size = 6, colour = "black",angle = 70, hjust = 1),
        axis.title.y = element_text(face = "bold.italic", size = 6))
p2
ggsave("class_ancom_watefall.tiff", width = 3.5, height = 4.5, compression="lzw", dpi = 600)


###
order_data = aggregate_taxa(physeq_ancombc, "Order")
ancom_bc <- ancombc(phyloseq=order_data, formula = "category",
                    p_adj_method = "bonferroni", zero_cut = 0.90, lib_cut = 1000, 
                    group = "category", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
                    max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE)

ancom_res_df <- data.frame(
  Species = row.names(ancom_bc$res$beta),
  beta = unlist(ancom_bc$res$beta),
  se = unlist(ancom_bc$res$se),
  W = unlist(ancom_bc$res$W),
  p_val = unlist(ancom_bc$res$p_val),
  q_val = unlist(ancom_bc$res$q_val),
  diff_abn = unlist(ancom_bc$res$diff_abn))

bonferroni_ancom <- ancom_res_df %>%
  dplyr::filter(q_val < 0.05)
dim(bonferroni_ancom)
bonferroni_ancom
write_csv(bonferroni_ancom, "ancombc_output_order_level.csv")
###
res = ancom_bc$res
##
samp_frac = ancom_bc$samp_frac
samp_frac[is.na(samp_frac)] = 0 # Replace NA with 0
log_obs_abn = log(abundances(phylum_data) + 1) # Add pesudo-count (1) to avoid taking the log of 0
log_obs_abn_adj = t(t(log_obs_abn) - samp_frac) # Adjust the log observed abundances
round(log_obs_abn_adj[, 1:17], 2) %>% 
  data.table(caption = "Bias-adjusted log observed abundances")
###
df_fig1 = data.frame(res$beta * res$diff_abn, check.names = FALSE) %>% 
  rownames_to_column("taxon_id")
df_fig2 = data.frame(res$se * res$diff_abn, check.names = FALSE) %>% 
  rownames_to_column("taxon_id")
colnames(df_fig2)[-1] = paste0(colnames(df_fig2)[-1], "SD")
df_fig = df_fig1 %>% left_join(df_fig2, by = "taxon_id") %>%
  transmute(taxon_id, categoryPb, categoryPbSD) %>%
  filter(categoryPb != 0) %>% arrange(desc(categoryPb)) %>%
  mutate(group = ifelse(categoryPb > 0, "g1", "g2"))
df_fig$taxon_id = factor(df_fig$taxon_id, levels = df_fig$taxon_id)

p3 = ggplot(data = df_fig, 
           aes(x = taxon_id, y = categoryPb, fill = group, color = group)) + 
  geom_bar(stat = "identity", width = 0.7, 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = categoryPb - categoryPbSD, ymax = categoryPb + categoryPbSD), width = 0.2,
                position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Fomitopsis betulina vs Fomes fomentarius \n log fold change") + 
  theme_bw() + 
  theme(legend.position = "none",
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(face = "italic", size = 6, colour = "black", angle = 70, hjust = 1),
        axis.title.y = element_text(face = "bold.italic", size = 6))
p3
ggsave("order_ancom_watefall.tiff", width = 6, height = 7, compression="lzw", dpi = 600)


###
family_data = aggregate_taxa(physeq_ancombc, "Family")
ancom_bc <- ancombc(phyloseq=family_data, formula = "category",
                    p_adj_method = "bonferroni", zero_cut = 0.90, lib_cut = 1000, 
                    group = "category", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
                    max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE)

ancom_res_df <- data.frame(
  Species = row.names(ancom_bc$res$beta),
  beta = unlist(ancom_bc$res$beta),
  se = unlist(ancom_bc$res$se),
  W = unlist(ancom_bc$res$W),
  p_val = unlist(ancom_bc$res$p_val),
  q_val = unlist(ancom_bc$res$q_val),
  diff_abn = unlist(ancom_bc$res$diff_abn))

bonferroni_ancom <- ancom_res_df %>%
  dplyr::filter(q_val < 0.05)
dim(bonferroni_ancom)
bonferroni_ancom
write_csv(bonferroni_ancom, "ancombc_output_family_level.csv")
###
res = ancom_bc$res
##
samp_frac = ancom_bc$samp_frac
samp_frac[is.na(samp_frac)] = 0 # Replace NA with 0
log_obs_abn = log(abundances(family_data) + 1) # Add pesudo-count (1) to avoid taking the log of 0
log_obs_abn_adj = t(t(log_obs_abn) - samp_frac) # Adjust the log observed abundances
round(log_obs_abn_adj[, 1:17], 2) %>% 
  data.table(caption = "Bias-adjusted log observed abundances")
###
df_fig1 = data.frame(res$beta * res$diff_abn, check.names = FALSE) %>% 
  rownames_to_column("taxon_id")
df_fig1$taxon_id <- str_replace(df_fig1$taxon_id, "d__Bacteria_", "")
df_fig2 = data.frame(res$se * res$diff_abn, check.names = FALSE) %>% 
  rownames_to_column("taxon_id")
df_fig2$taxon_id <- str_replace_all(df_fig2$taxon_id, "d__Bacteria_", "")
colnames(df_fig2)[-1] = paste0(colnames(df_fig2)[-1], "SD")
df_fig = df_fig1 %>% left_join(df_fig2, by = "taxon_id") %>%
  transmute(taxon_id, categoryPb, categoryPbSD) %>%
  filter(categoryPb != 0) %>% arrange(desc(categoryPb)) %>%
  mutate(group = ifelse(categoryPb > 0, "g1", "g2"))
df_fig$taxon_id = factor(df_fig$taxon_id, levels = df_fig$taxon_id)
 
p4 = ggplot(data = df_fig, 
           aes(x = taxon_id, y = categoryPb, fill = group, color = group)) + 
  geom_bar(stat = "identity", width = 0.7, 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = categoryPb - categoryPbSD, ymax = categoryPb + categoryPbSD), width = 0.2,
                position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Fomitopsis betulina vs Fomes fomentarius effect size (log fold change)") + 
  theme_q2r() + 
  coord_flip()+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(face = "italic", size = 10, colour = "black"),
        axis.title.x = element_text(face = "bold.italic", size = 12))
p4
ggsave("family_ancom_watefall.tiff", width = 9.5, height = 7, compression="lzw", dpi = 600)


###
genus_data = aggregate_taxa(physeq_ancombc, "Genus")
ancom_bc <- ancombc(phyloseq=genus_data, formula = "category",
                    p_adj_method = "bonferroni", zero_cut = 0.90, lib_cut = 1000, 
                    group = "category", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
                    max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE)

ancom_res_df <- data.frame(
  Species = row.names(ancom_bc$res$beta),
  beta = unlist(ancom_bc$res$beta),
  se = unlist(ancom_bc$res$se),
  W = unlist(ancom_bc$res$W),
  p_val = unlist(ancom_bc$res$p_val),
  q_val = unlist(ancom_bc$res$q_val),
  diff_abn = unlist(ancom_bc$res$diff_abn))

bonferroni_ancom <- ancom_res_df %>%
  dplyr::filter(q_val < 0.05)
dim(bonferroni_ancom)
bonferroni_ancom
write_csv(bonferroni_ancom, "ancombc_output_genus_level.csv")
###
res = ancom_bc$res
##
samp_frac = ancom_bc$samp_frac
samp_frac[is.na(samp_frac)] = 0 # Replace NA with 0
log_obs_abn = log(abundances(phylum_data) + 1) # Add pesudo-count (1) to avoid taking the log of 0
log_obs_abn_adj = t(t(log_obs_abn) - samp_frac) # Adjust the log observed abundances
round(log_obs_abn_adj[, 1:17], 2) %>% 
  data.table(caption = "Bias-adjusted log observed abundances")
###
df_fig1 = data.frame(res$beta * res$diff_abn, check.names = FALSE) %>% 
  rownames_to_column("taxon_id")
df_fig2 = data.frame(res$se * res$diff_abn, check.names = FALSE) %>% 
  rownames_to_column("taxon_id")
colnames(df_fig2)[-1] = paste0(colnames(df_fig2)[-1], "SD")
df_fig = df_fig1 %>% left_join(df_fig2, by = "taxon_id") %>%
  transmute(taxon_id, categoryPb, categoryPbSD) %>%
  filter(categoryPb != 0) %>% arrange(desc(categoryPb)) %>%
  mutate(group = ifelse(categoryPb > 0, "g1", "g2"))
df_fig$taxon_id = factor(df_fig$taxon_id, levels = df_fig$taxon_id)

p5 = ggplot(data = df_fig, 
           aes(x = taxon_id, y = categoryPb, fill = group, color = group)) + 
  geom_bar(stat = "identity", width = 0.7, 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = categoryPb - categoryPbSD, ymax = categoryPb + categoryPbSD), width = 0.2,
                position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Fomitopsis betulina vs Fomes fomentarius \n log fold change") + 
  theme_bw() + 
  theme(legend.position = "none",
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(face = "italic", size = 6, colour = "black", angle = 80, hjust = 1),
        axis.title.y = element_text(face = "bold.italic", size = 6))
p5
ggsave("genus_ancom_watefall.tiff", width = 6, height = 7, compression="lzw", dpi = 600)
###
species_data = aggregate_taxa(physeq_ancombc, "Species")
ancom_bc <- ancombc(phyloseq=species_data, formula = "category",
                    p_adj_method = "bonferroni", zero_cut = 0.90, lib_cut = 1000, 
                    group = "category", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
                    max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE)

ancom_res_df <- data.frame(
  Species = row.names(ancom_bc$res$beta),
  beta = unlist(ancom_bc$res$beta),
  se = unlist(ancom_bc$res$se),
  W = unlist(ancom_bc$res$W),
  p_val = unlist(ancom_bc$res$p_val),
  q_val = unlist(ancom_bc$res$q_val),
  diff_abn = unlist(ancom_bc$res$diff_abn))

bonferroni_ancom <- ancom_res_df %>%
  dplyr::filter(q_val < 0.05)
dim(bonferroni_ancom)
bonferroni_ancom
write_csv(bonferroni_ancom, "ancombc_output_species_level.csv")
###
res = ancom_bc$res
##
###
samp_frac = ancom_bc$samp_frac
samp_frac[is.na(samp_frac)] = 0 # Replace NA with 0
log_obs_abn = log(abundances(phylum_data) + 1) # Add pesudo-count (1) to avoid taking the log of 0
log_obs_abn_adj = t(t(log_obs_abn) - samp_frac) # Adjust the log observed abundances
round(log_obs_abn_adj[, 1:17], 2) %>% 
  data.table(caption = "Bias-adjusted log observed abundances")
###
df_fig1 = data.frame(res$beta * res$diff_abn, check.names = FALSE) %>% 
  rownames_to_column("taxon_id")
df_fig2 = data.frame(res$se * res$diff_abn, check.names = FALSE) %>% 
  rownames_to_column("taxon_id")
colnames(df_fig2)[-1] = paste0(colnames(df_fig2)[-1], "SD")
  df_fig = df_fig1 %>% left_join(df_fig2, by = "taxon_id") %>%
  transmute(taxon_id, categoryPb, categoryPbSD) %>%
  filter(categoryPb != 0) %>% arrange(desc(categoryPb)) %>%
  mutate(group = ifelse(categoryPb > 0, "g1", "g2"))
df_fig$taxon_id = factor(df_fig$taxon_id, levels = df_fig$taxon_id)

p6 = ggplot(data = df_fig, 
           aes(x = taxon_id, y = categoryPb, fill = group, color = group)) + 
  geom_bar(stat = "identity", width = 0.7, 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = categoryPb - categoryPbSD, ymax = categoryPb + categoryPbSD), width = 0.2,
                position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Fomitopsis betulina vs Fomes fomentarius \n log fold change") + 
  theme_bw() + 
  theme(legend.position = "none",
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(face = "italic", size = 6, colour = "black", angle = 80, hjust = 1),
        axis.title.y = element_text(face = "bold.italic", size = 6))
p6
ggsave("species_ancom_watefall.tiff", width = 6, height = 7, compression="lzw", dpi = 600)
###
suppl_fig <- p1+p2+p3
suppl_fig + plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face="bold"))
ggsave("supplementary_fig_2abc.tiff", width = 9.5, height = 5.8, compression="lzw", dpi = 600)

suppl_figg <- p5+p6
suppl_figg + plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face="bold"))
ggsave("supplementary_fig_3ab.tiff", width = 10.5, height = 8.5, compression="lzw", dpi = 600)

### ANCOM-BC analyses
ancom_bc <- ancombc(phyloseq=physeq_ancombc, formula = "category",
                    p_adj_method = "bonferroni", zero_cut = 0.90, lib_cut = 1000, 
                    group = "category", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
                    max_iter = 100, conserve = TRUE, alpha = 0.05, global = FALSE)

ancom_res_df <- data.frame(
  Species = row.names(ancom_bc$res$beta),
  beta = unlist(ancom_bc$res$beta),
  se = unlist(ancom_bc$res$se),
  W = unlist(ancom_bc$res$W),
  p_val = unlist(ancom_bc$res$p_val),
  q_val = unlist(ancom_bc$res$q_val),
  diff_abn = unlist(ancom_bc$res$diff_abn))

bonferroni_ancom <- ancom_res_df %>%
  dplyr::filter(q_val < 0.05)
dim(bonferroni_ancom)
bonferroni_ancom
write_csv(bonferroni_ancom, "ancombc_output_whole_table_no_agregation.csv")

###

### Figure_03 (New_analysis)
setwd("/home/irshad/Schilling_Project_075/processed/shi7_learning/fastqs_stitched/fastafiles/figures/New_raw_data_figures/")
ffbar <- readxl::read_xlsx("ff_phyla_10_new.xlsx")
ffbar
ffbar1 <- reshape2::melt(ffbar, value.name = "relative", variable.name= "phyla")
ffbar1

colours = c( "#A54657",  "#582630", "#F7EE7F", "#4DAA57","#F1A66A","#F26157", "#679289","#F9ECCC", "#33658A", "#F6AE2D")
ffbar1$Sample <- factor(ffbar1$index,levels=unique(ffbar1$index))

pl1 = ggplot(ffbar1, aes(x = index, fill = phyla, y = relative)) + 
  geom_bar(stat = "identity", colour = "black")+
  theme_classic()+
  theme_q2r()+
  labs_pubr()+
  labs(x = "", y = "Relative Abundance (%)", fill = "Phyla") + 
  scale_fill_manual(values=colours)+
  theme(axis.text.x = element_text(angle = 35, size = 10, colour = "black", vjust =1, hjust = 1, face= "italic"),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(face = "italic"),
        legend.position = "none")

pl1

###
pbbar <- readxl::read_xlsx("pb_phyla_10_new.xlsx")
pbbar
pbbar1 <- reshape2::melt(pbbar, value.name = "relative", variable.name= "phyla")
pbbar1

colours = c( "#A54657",  "#582630", "#F7EE7F", "#4DAA57","#F1A66A","#F26157", "#679289","#F9ECCC", "#33658A", "#F6AE2D")
pbbar1$Sample <- factor(pbbar1$index,levels=unique(pbbar1$index))

pl2 = ggplot(pbbar1, aes(x = index, fill = phyla, y = relative)) + 
  geom_bar(stat = "identity", colour = "black")+
  theme_classic()+
  theme_q2r()+
  labs_pubr()+
  labs(x = "", y = "Relative Abundance (%)", fill = "Phyla") + 
  scale_fill_manual(values=colours)+
  theme(axis.text.x = element_text(angle = 35, size = 10, colour = "black", vjust = 1, hjust = 1, face= "italic"),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(face = "italic"),
        axis.title.y = element_blank())

pl2

patchwork <-(pl1 + pl2) 
patchwork + plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face="bold"))

ggsave("Fig_0333.tiff", width = 8, height =6, compression="lzw", dpi = 600)

###
setwd("/home/irshad/Schilling_Project_075/processed/shi7_learning/fastqs_stitched/fastafiles/figures/New_raw_data_figures/")

all4 <- read_excel("top_04_phyla_12-09-2021.xlsx")
all4
##
means_all <- all4 %>%
  group_by(category) %>%
  summarise_at(.vars = c("Proteobacteria", "Firmicutes", "Acidobacteria", "Actinobacteria"), list(name=mean))
view(means_all)
##
firm<- ggplot(all4, aes(x=category, y= Firmicutes, fill=category))+
  geom_half_boxplot(side = "r")+
  geom_half_dotplot(binwidth = 1.8)+
  theme_classic()+
  ggtitle('Firmicutes' )+
  theme(
    plot.title=element_text(face='bold.italic', size=8))+
  labs(x = "", y = "Relative Abundance (%)") +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size =8, face = "bold"),
    strip.text.x = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 12, face = "bold.italic"),
    legend.position = "none"
    
  )
firm
firm <- firm +
  geom_line(data = tibble(x=c(0.99, 2.3), y=c(76, 76)), 
            aes(x=x, y=y), 
            inherit.aes = FALSE)+
  geom_text(data = tibble(x=c(1.645), y=c(77.5)), 
            aes(x=x, y=y, label="*"), size=5,
            inherit.aes = FALSE)
firm
###
prot<- ggplot(all4, aes(x=category, y= Proteobacteria, fill=category))+
  geom_half_boxplot(side = "r")+
  geom_half_dotplot(binwidth = 1.8)+
  theme_classic()+
  ggtitle('Proteobacteria' )+
  theme(
    plot.title=element_text(face='bold.italic', size=8))+
  labs(x = "", y = "Relative Abundance (%)") +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size =8, face = "bold"),
    strip.text.x = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 12, face = "bold.italic"),
    legend.position = "none"
    
  )
prot
prot <- prot +
  geom_line(data = tibble(x=c(0.9, 2.4), y=c(81, 81)), 
            aes(x=x, y=y), 
            inherit.aes = FALSE)+
  geom_text(data = tibble(x=c(1.65), y=c(83)), 
            aes(x=x, y=y, label="*"), size=5,
            inherit.aes = FALSE)
prot


###
acid<- ggplot(all4, aes(x=category, y= Acidobacteria, fill=category))+
  geom_half_boxplot(side = "r")+
  geom_half_dotplot(binwidth = 0.5)+
  theme_classic()+
  ggtitle('Acidobacteria' )+
  theme(
    plot.title=element_text(face='bold.italic', size=8))+
  labs(x = "", y = "Relative Abundance (%)") +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size =8, face = "bold"),
    strip.text.x = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 6, face = "bold.italic"),
    legend.position = "bottom",
    legend.key.size = unit(0.3, "cm")
    
  )
acid
acid <- acid +
  geom_line(data = tibble(x=c(1, 2.38), y=c(19, 19)), 
            aes(x=x, y=y), 
            inherit.aes = FALSE)+
  geom_text(data = tibble(x=c(1.65), y=c(19.4)), 
            aes(x=x, y=y, label="*"), size=5,
            inherit.aes = FALSE)
acid


###
actino<- ggplot(all4, aes(x=category, y= Actinobacteria, fill=category))+
  geom_half_boxplot(side = "r")+
  geom_half_dotplot(binwidth = 1.3)+
  theme_classic()+
  ggtitle('Actinobacteria' )+
  theme(
    plot.title=element_text(face='bold.italic', size=8))+
  labs(x = "", y = "Relative Abundance (%)") +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size =8, face = "bold"),
    strip.text.x = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 6, face = "bold.italic"),
    legend.position = "bottom",
    legend.key.size = unit(0.3, "cm")
    
  )
actino
actino <- actino +
  geom_line(data = tibble(x=c(1, 2.38), y=c(54, 54)), 
            aes(x=x, y=y), 
            inherit.aes = FALSE)+
  geom_text(data = tibble(x=c(1.65), y=c(57)), 
            aes(x=x, y=y, label="n.s."), size=4,
            inherit.aes = FALSE)
actino

## Use patchwork to combine and annotate all four plots

p<-firm+prot+acid+actino
p + plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face="bold", size = 10))
##
ggsave("Relative_abundance_stats_separate_new.tiff", width = 5.3, height = 5.2, compression="lzw", dpi=600)
##
kt <- kruskal.test(Actinobacteria ~ category, data=all4)
kt$p.value
####
data_fir <- read.csv("firmicutes_new.csv")
data_fir
ggboxplot(data_fir, x="category", y="relative_abundance", color = "category")
res.kruskal <- data_fir %>% kruskal_test(relative_abundance ~ category)
res.kruskal
##effect_size
data_fir %>% kruskal_effsize(relative_abundance ~ category)
##
pwc <-data_fir %>%
  dunn_test(relative_abundance ~ category, p.adjust.method =  "fdr")
pwc
descriptives_firm = data_fir %>% group_by(category) %>% summarise(mean = mean(relative_abundance), sd = sd(relative_abundance), n = n())
descriptives_firm
###
pwc2 <- data_fir %>% 
  wilcox_test(relative_abundance ~ category, p.adjust.method = "fdr")
pwc2

# Visualization: box plots with p-values
pwc3 <- pwc2 %>% add_xy_position(x = "category")
ggboxplot(data_fir, x = "category", y = "relative_abundance", color = "category") +
  stat_pvalue_manual(pwc3, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.kruskal, detailed = TRUE),
    caption = get_pwc_label(pwc3)
  )
###
proteo <- read_csv("proteobact_new.csv")
view(proteo)
ggboxplot(proteo, x="category", y="relative_abundance", color = "category")
res.kruskal <- proteo %>% kruskal_test(relative_abundance ~ category)
res.kruskal
##effect_size
proteo %>% kruskal_effsize(relative_abundance ~ category)
##
pwc <-proteo %>%
  dunn_test(relative_abundance ~ category, p.adjust.method = "fdr")
pwc
descriptives_prot = proteo %>% group_by(category) %>% summarise(mean = mean(relative_abundance), sd = sd(relative_abundance), n = n())
descriptives_prot
###
pwc2 <- proteo %>% 
  wilcox_test(relative_abundance ~ category, p.adjust.method = "bonferroni")
pwc2

# Visualization: box plots with p-values
pwc3 <- pwc2 %>% add_xy_position(x = "category")
ggboxplot(proteo, x = "category", y = "relative_abundance", color = "category") +
  stat_pvalue_manual(pwc3, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.kruskal, detailed = TRUE),
    caption = get_pwc_label(pwc3)
  )
###
actino<- read_csv("actino_new.csv")
view(actino)
ggboxplot(actino, x="category", y="relative_abundance", color = "category")
res.kruskal <- actino %>% kruskal_test(relative_abundance ~ category)
res.kruskal

##effect size
actino %>% kruskal_effsize(relative_abundance ~ category)
##
pwwc1 <-actino %>%
  dunn_test(relative_abundance ~ category, p.adjust.method = "fdr")
pwwc1
descriptives_actino = actino %>% group_by(category) %>% summarise(mean = mean(relative_abundance), sd = sd(relative_abundance), n = n())
descriptives_actino
##
pwwc <- actino %>% 
  wilcox_test(relative_abundance ~ category, p.adjust.method = "bonferroni")
pwwc

##
acido<- read_csv("acido_new.csv")
ggboxplot(acido, x="category", y="relative_abundance", color = "category")
res.kruskal <- acido %>% kruskal_test(relative_abundance ~ category)
res.kruskal
###

##effect size
acido %>% kruskal_effsize(relative_abundance ~ category)
##
pwwwc1 <-acido %>%
  dunn_test(relative_abundance ~ category, p.adjust.method = "fdr")
pwwwc1
descriptives_acido = acido %>% group_by(category) %>% summarise(mean = mean(relative_abundance), sd = sd(relative_abundance), n = n())
descriptives_acido
##
pwwwwc1 <- acido %>% 
  wilcox_test(relative_abundance ~ category, p.adjust.method = "bonferroni")
pwwwwc1

###
library(tidyverse)
library(readxl)
library(ggplot2)
library(patchwork)
setwd("/home/irshad/16S_rRNA_Analysis/Molly_data/")

wph <- read_xlsx("whiterot_pH.xlsx")
wph
shapiro.test(wph$pH)

bph <- read_xlsx("brownrot_pH.xlsx")
bph
shapiro.test(bph$pH)


###
indp_ttest <- read_xlsx("indp_test_pH.xlsx")
indp_ttest
t.test (pH ~ category, var.equal=TRUE, data=indp_ttest)
###
descriptives = indp_ttest %>% group_by(category) %>% summarise(mean = mean(pH), sd = sd(pH), n = n())
descriptives
mean_difference = descriptives[2,2] - descriptives[1,2]
mean_difference
###
ggplot(descriptives, aes(category, mean))+
  geom_col()+
  geom_errorbar(data = descriptives, mapping = aes(x = category, y = mean, ymin = mean - sd, ymax = mean + sd), size=0.5, color="black", width=.1) 


###
ggplot(indp_ttest, aes(category, pH))+
  geom_jitter(position = position_jitter(w=0.05, h=0.01))+
  geom_point(data=descriptives, aes(x = category, y = mean), size=2, color= "red")+
  geom_errorbar(data = descriptives, mapping = aes(x = category, y = mean, ymin = mean - sd, ymax = mean + sd), size=0.5, color="grey", width=.4) 
###
library(ggpubr)
ggstripchart(indp_ttest, x = "category", y = "pH",
             color = "category",
             palette = c("#00AFBB", "#FC4E07"))
###
ph_plot <- ggplot(indp_ttest, aes(category, pH, color=category))+
  geom_boxplot(width=0.5)+
  geom_jitter(position=position_jitter(w=0.05, h=0.01), size=0.6, color="black")+
  labs(x="")+
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_text(colour = "black", size = 6, face = "italic"),
        axis.text.y = element_text(colour = "black", size = 6),
        axis.title.y = element_text(color = "black", size =6))
ph_plot
ph_plot <- ph_plot+
  geom_line(data = tibble(x=c(0.99, 2), y=c(4.6, 4.6)), 
            aes(x=x, y=y), 
            size=0.5,
            inherit.aes = FALSE)+
  geom_text(data = tibble(x=c(1.5), y=c(4.62)), 
            aes(x=x, y=y, label="*"), size=3,
            inherit.aes = FALSE)
ph_plot
#
#ggsave("pH.tiff", width = 4, height =4.5, compression="lzw", dpi = 600)
###
wl <- read_xlsx("white_rot_lignin.xlsx")
wl
shapiro.test(wl$AIM)
bl <- read_xlsx("brown_rot_lignin.xlsx")
bl
shapiro.test(bl$AIM)
lig <- read_xlsx("indp_test_lignin.xlsx")
lig
t.test (Lignin ~ category, var.equal=TRUE, data=lig)
descriptives_lig = lig %>% group_by(category) %>% summarise(mean = mean(Lignin), sd = sd(Lignin), n = n())
descriptives_lig
###
ggplot(lig, aes(category, Lignin))+
  geom_jitter(position = position_jitter(w=0.05, h=0.01))+
  geom_point(data=descriptives_lig, aes(x = category, y = mean), size=2, color= "red")+
  geom_errorbar(data = descriptives_lig, mapping = aes(x = category, y = mean, ymin = mean - sd, ymax = mean + sd), size=0.5, color="grey", width=.4) 
###
###
lignin_plot <- ggplot(lig, aes(category, Lignin, color=category))+
  geom_boxplot(width=0.5)+
  geom_jitter(position=position_jitter(w=0.05, h=0.01), size=0.6, color="black")+
  labs(x="", y= "Lignin (%)")+
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_text(colour = "black", size = 6, face = "italic"),
        axis.text.y = element_text(colour = "black", size = 6),
        axis.title.y = element_text(color = "black", size =6))

lignin_plot
lignin_plot <- lignin_plot+
  geom_line(data = tibble(x=c(0.99, 2), y=c(42, 42)), 
            aes(x=x, y=y), 
            size=0.5,
            inherit.aes = FALSE)+
  geom_text(data = tibble(x=c(1.5), y=c(42.2)), 
            aes(x=x, y=y, label="*"), size=3,
            inherit.aes = FALSE)+
  geom_line(data=tibble(x=c(0.5,2.5), y=c(30.2, 30.2)),
            aes(x=x, y=y),
            inherit.aes = FALSE,
            linetype="twodash",
            size=0.3,
            color="red")+
  geom_text(data = tibble(x=c(1.4), y=c(31.4)), 
            aes(x=x, y=y, label="Brown rot dominated"), size=1.5,
            inherit.aes = FALSE)+
  geom_text(data = tibble(x=c(1.4), y=c(30.6)), 
            aes(x=x, y=y, label=" '0.8 threshold' "), size=1.5,
            inherit.aes = FALSE)+
  geom_text(data = tibble(x=c(1.4), y=c(29.7)), 
            aes(x=x, y=y, label="White rot dominated"), size=1.5,
            inherit.aes = FALSE)

lignin_plot

#ggsave("lignin.tiff", width = 3, height =2.5, compression="lzw", dpi = 600)
###
plots <- ph_plot + lignin_plot
plots <- plots + plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face="bold", size = 6))
plots
ggsave("/home/irshad/Schilling_Project_075/processed/shi7_learning/fastqs_stitched/fastafiles/figures/New_raw_data_figures/Figure_1_pH_lignin.tiff", width = 5, height =3, compression="lzw", dpi = 600)
##
