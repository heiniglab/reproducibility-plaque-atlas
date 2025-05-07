suppressWarnings(library(BayesPrism))
#suppressWarnings(library(TSIclient))
suppressWarnings(library(dplyr))
suppressWarnings(library(SummarizedExperiment))
library(SingleCellExperiment)
#suppressPackageStartupMessages(library(rtracklayer)) # Parse GTF
suppressPackageStartupMessages(library(ggplot2))
#library(gplots) # heatmap.2
library(zellkonverter)
library(compositions)
library(ggsignif)
library(ggpubr)
library(reshape2)
library(tidyr)


bulk.count_attributes = readRDS("/Users/korbinian.traeuble/PhD-local/projects/main_Roche/data/bulk/Maegdefessel/gene_summarized_experiment.RDS")
# calculate those on cluster

theta <- readRDS("/Users/korbinian.traeuble/PhD-local/projects/main_Roche/data/bulk/Maegdefessel/cell_type_deconvolution/sign/rev2_normalECs/theta.rds")
theta.cv <- readRDS("/Users/korbinian.traeuble/PhD-local/projects/main_Roche/data/bulk/Maegdefessel/cell_type_deconvolution/sign/rev2_normalECs/theta.cv.rds")
output_dir <- "/Users/korbinian.traeuble/PhD-local/projects/main_Roche/data/bulk/Maegdefessel/cell_type_deconvolution/sign/rev2_normalECs"
############################

# Abundances

# Filter samples according to Matthias M.
# Samples
samples.rm.pca = c("22L000719", "22L000728", "22L000872", "22L000884", "22L000886", "22L000889", "22L000893")
samples.rm.perc_mapped = c('120658','120659','120661','120663','120664','120665','120666','120680',
                           '120681','120682','120683','120684','120685','130004','22L000714','22L000719',
                           '22L000722','22L000728','22L000736','22L000817','22L000826','22L000837',
                           '22L000852','22L000872','22L000876','22L000884','22L000886','22L000889',
                           '22L000893','22L001767','22L001770','22L001785') # < 75%
samples.rm.pct_duplication = '120655' # > 60%
samples.rm.num_mapped = c('22L000714','22L000719','22L000722','22L000728','22L000736',
                          '22L000817','22L000826','22L000837','22L000852','22L000860',
                          '22L000872','22L000876','22L000884','22L000886','22L000889',
                          '22L000893','22L001767','22L001770','22L001785') # < 10,000,000
samples.rm = unique(c(samples.rm.pca, samples.rm.perc_mapped, samples.rm.pct_duplication, samples.rm.num_mapped))



df.abund = reshape2::melt(data.frame(sample_id = rownames(theta), theta, check.names=FALSE), id.vars = c("sample_id"), variable.name = "cell_type") %>%
  group_by(cell_type) %>%
  mutate(mean_prop = mean(value)) %>%
  arrange(mean_prop) %>%
  ungroup() %>%
  mutate(cell_type = factor(cell_type, levels = unique(cell_type)))


# Remove QC samples
df.abund <- df.abund[!df.abund$sample_id %in% samples.rm, ]
theta <- subset(theta, !rownames(theta) %in% samples.rm)

# Apply CLR

df.abund$clr_value <- as.numeric(clr(df.abund$value))

# stratified abundances

is_diseased <- bulk.count_attributes$is_diseased
is_stable <- bulk.count_attributes$is_stable
Symptomatic <- bulk.count_attributes$Symptomatic
Perc_stenosis <- bulk.count_attributes$Perc_stenosis


df <- data.frame(sample_id = colnames(bulk.count_attributes),
                 is_diseased = is_diseased,
                 is_stable = is_stable,
                 Symptomatic = Symptomatic,
                 Perc_stenosis= Perc_stenosis)

df.abund <- merge(df.abund, df, by = "sample_id")

df.abund$is_diseased <- factor(df.abund$is_diseased, labels = c("Non-Diseased", "Diseased"))
df.abund$Symptomatic <- factor(df.abund$Symptomatic,levels = c("a", "s"), labels = c("Asymptomatic", "Symptomatic"))
df.abund$is_stable <- factor(df.abund$is_stable, levels = c(0, 1), labels = c("Unstable", "Stable"))

# NO CLR
ggplot(df.abund, aes(y = cell_type, x = value, fill = is_diseased)) +
  geom_boxplot() +
  theme_bw() + 
  labs(x = "Proportion", y = "", fill = "Group") +
  scale_x_continuous(labels = scales::percent, expand = expansion(c(0, 0.05), 0)) +
  scale_fill_manual(values = c("Non-Diseased" = "blue", "Diseased" = "red"),
                    labels = c("Non-Diseased" = "Early Lesion", "Diseased" = "Late Lesion"))+
  theme(
    axis.title.x = element_text(size = 16), # X-axis title font size
    axis.text.x = element_text(size = 14),  # X-axis text font size
    axis.text.y = element_text(size = 14),  # Y-axis text font size
    legend.title = element_text(size = 16), # Legend title font size
    legend.text = element_text(size = 14)   # Legend text font size
  )


ggsave(file.path(output_dir, "abundance_disease2.pdf"), height = 8, width = 8, scale = 1, device = pdf)


# WITH CLR

#with pvals diseased
p <- ggplot(df.abund, aes(y = cell_type, x = clr_value, fill = is_diseased)) +
  geom_boxplot() +
  theme_bw() +
  labs(x = "CLR Transformed Proportion", y = "", fill = "Group") +
  scale_x_continuous(expand = expansion(c(0, 0.05), 0)) +
  scale_fill_manual(values = c("Non-Diseased" = "blue", "Diseased" = "red"),
                    labels = c("Non-Diseased" = "Early Lesion", "Diseased" = "Late Lesion")) +
  theme(
    axis.title.x = element_text(size = 16), # X-axis title font size
    axis.text.x = element_text(size = 14),  # X-axis text font size
    axis.text.y = element_text(size = 14),  # Y-axis text font size
    legend.title = element_text(size = 16), # Legend title font size
    legend.text = element_text(size = 14)   # Legend text font size
  )

# Determine the rightmost position for the text, adding extra space for the labels
max_value <- max(df.abund$clr_value, na.rm = TRUE)
label_position <- max_value + (max_value * 0.25) # Adjust the multiplier as needed for spacing

# Calculate p-values
p_values <- df.abund %>%
  group_by(cell_type) %>%
  summarise(p_value = t.test(clr_value ~ is_diseased)$p.value) %>%
  mutate(label = ifelse(p_value < 0.05, paste0("p = ", format(p_value, digits = 2)), "ns"))


# Adjust plot margins to make space for the p-value labels
p <- p + theme(plot.margin = margin(r = 2, unit = "pt")) # Adjust the right margin as needed

# Add the p-values to the plot, positioning them at the calculated label_position
p_final <- p + geom_text(data = p_values, aes(y = cell_type, x = label_position, label = label), hjust = 1, vjust = 0, inherit.aes = FALSE)
print(p_final)

ggsave(file.path(output_dir, "abundance_disease_CLR_pvals_level22.pdf"), plot=p_final,height = 8, width = 12, scale = 1, device=pdf)
