setwd('D:/Data/0_lixia/0)_8e_HTseq_all/0_Paper_Figure_Data_source_2024/0_d10mu2/080_cas13_gRNA_overlap_PLOT')

library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(stringr)

# Set color palette for fc1 and fc2
cols = c("Day7" = "#1F78B4", "Day14" = "#F0027F")  # Green for fc1, Pink for fc2
sample = 'NLS'
xlim = c(-0.2, 1.52) 
xb=seq(0, 1.5, by = 0.25)
#sample = 'K562_NLS'
#sample ='NES'
#sample = 'K562_NES'
# Load the data
data1 = read.table(str_c(sample, '.sup5_gRNA_summary.HEK293T.base_count'), header = TRUE, sep = '\t')
data1$fc1=0-data1$fc1
data1$fc2=0-data1$fc2

min(data1$fc2)
# Reshape the data into long format
data_long <- reshape(data1, 
                     varying = c('number', 'A'), 
                     v.names = 'value', 
                     timevar = 'type', 
                     times = c('number', 'A'), 
                     direction = 'long')

# Modify 'type' column to have values "A" and "I"
data_long$type <- ifelse(data_long$type == "A", "A", "I")

# Filter out rows where 'value' is greater than 10
data_long <- data_long %>% filter(value <= 10)

# Reshape data for fc1 and fc2 comparison
data_long_fc <- data_long %>%
  gather(key = "fc_type", value = "value_fc", fc1, fc2) %>%
  mutate(
    fc_type = factor(fc_type, levels = c("fc1", "fc2")),  # Ensure the right order in the plot
    linetype = ifelse(fc_type == "fc1", "solid", "dashed"),  # Line type for fc1 vs fc2
    shape = ifelse(type == "A", 16, 17),  # Shape for A (solid circle) and I (solid triangle)
    fc_label = ifelse(fc_type == "fc1", "Day7", "Day14")  # Assign labels for legends
  )

# Summarize the data: Calculate mean, standard deviation, standard error, and confidence intervals
df_summary_fc <- data_long_fc %>%
  group_by(fc_type, type, value) %>%
  summarise(
    mean_value = mean(value_fc, na.rm = TRUE),
    sd_value = sd(value_fc, na.rm = TRUE),
    n = n(),
    se_value = sd_value / sqrt(n),
    ei_lower = mean_value - 1 * se_value,  # 1 * SE for the lower bound
    ei_upper = mean_value + 1 * se_value   # 1 * SE for the upper bound
  ) %>%
  # Add the fc_label to the summary data
  left_join(select(data_long_fc, fc_type, fc_label) %>% distinct(), by = "fc_type")

# Perform Mann-Whitney U test for each pair and mark significance
significance_results <- data_long_fc %>%
  group_by(fc_type, value) %>%
  summarise(
    mw_test_result = wilcox.test(value_fc ~ type)$p.value,
    significance = ifelse(mw_test_result < 0.0001, "**", NA)
  )

# Join significance results to the summarized data
df_summary_fc <- left_join(df_summary_fc, significance_results, by = c("fc_type", "value"))

# Create the plot with both fc1 and fc2 together
pdf(str_c(sample, '.fc1_fc2.all.compare.pdf'), wi = 8, he = 8)
ggplot(df_summary_fc, aes(x = value, y = mean_value, color = fc_label, 
                          linetype = fc_type, shape = type)) +
  geom_point(size = 4.5) +  # Add points for each type and value
  geom_line(size = 1) +  # Add lines connecting points
  geom_errorbar(aes(ymin = ei_lower, ymax = ei_upper), width = 0.1) +  # Add error bars
  labs(
    x = 'The number of background- or modified- 
A sites in gRNA', 
    y = 'On-target activity',
    title = ''
  ) +
  theme_bw(base_rect_size = 1.5) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 1, colour = 'black'),
    axis.title = element_text(size = 18, face = 'bold', color = 'black'),
    axis.text.x = element_text(size = 14, face = 'bold', color = 'black'),
    axis.text.y = element_text(size = 14, face = 'bold', color = 'black'),
    axis.ticks = element_line(linetype = 1, color = 'black', linewidth = 1.2),
    axis.ticks.length = unit(0.25, 'cm'),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, face = 'bold', color = 'black'),
    legend.key.width = unit(1, 'cm')
  ) +
  scale_color_manual(values = cols) +  # Color lines based on fc_label (log2FC.D7 and log2FC.D14)
  scale_linetype_manual(values = c("fc1" = "solid", "fc2" = "dashed")) +  # Line types for fc1 and fc2
  scale_shape_manual(values = c("A" = 16, "I" = 17)) +  # Shapes for A and I types
  scale_y_continuous(limits = xlim, breaks = xb) +
  scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 1)) +
  # Add stars for significant results with different colors
  geom_text(data = df_summary_fc %>%
              filter(type == "A" & !is.na(significance) & fc_type == "fc1"),  # Filter to only include type "A" and significance
            aes(label = significance,
                y = ei_upper - 0.1),
            color = "#1F78B4", size = 10
            ) +  # Adjust size of the star markers
  geom_text(data = df_summary_fc %>%
              filter(type == "I" & !is.na(significance) & fc_type == "fc2"),  # Filter to only include type "A" and significance
            aes(label = significance,
                y = ei_upper + 0.05),
            color = "#F0027F", size = 10
  ) +  # Adjust size of the star markers
  # Add the 'n' count for 'A' type only
  geom_text(data = df_summary_fc %>% filter(type == "A" & fc_type == "fc1"),  # Filter for type A only
            aes(label = n, 
                y = ei_upper - 0.15),  # Adjust y position for count text of A
            size = 3, color = "#1F78B4",face='bold') +  # Set color for A counts (blue)
  # Add the 'n' count for 'I' type only
  geom_text(data = df_summary_fc %>% filter(type == "I" & fc_type == "fc2"),  # Filter for type I only
            aes(label = n, 
                y = ei_upper + 0.15),  # Adjust y position for count text of I
            size = 3, color = "#F0027F",face='bold') +  # Set color for I counts (red)
  guides(
    shape = guide_legend(title = "Type", 
                         override.aes = list(size = 3)),
    color = guide_legend(title = "FC Group", 
                         override.aes = list(size = 3)),
    linetype = guide_legend(title = "Line Type", 
                            override.aes = list(size = 1.5))
  )
dev.off()
