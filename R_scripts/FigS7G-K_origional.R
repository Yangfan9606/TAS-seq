setwd('D:/Data/0_lixia/0)_8e_HTseq_all/0_Paper_Figure_Data_source_2024/0_d10mu2/077_RBP_region_ext_summary_PLOT/K562_encode/RMR_select_RBP')

library(ggplot2)
library(dplyr)
library(stringr)

# Example data (replace df1 and df2 with your actual data frames)
data1=read.table('K562_NLS.norRMR.RBP.sum.long',header=F,sep='\t')
data2=read.table('K562_NES.norRMR.RBP.sum.long',header=F,sep='\t')
names(data1)=c('RBP','type','value')
names(data2)=c('RBP','type','value')
cols=c("#7FC97F","#F0027F")

rbp_name='ENCSR001KKZ.PHF6'
df1=subset(data1,RBP==rbp_name)
df2=subset(data2,RBP==rbp_name)
# Define custom order for types
sample_order <- c("ext5_10","ext5_9","ext5_8","ext5_7","ext5_6","ext5_5","ext5_4", 
                  "ext5_3", "ext5_2", "ext5_1", "K562", "ext3_1", "ext3_2", "ext3_3", 
                  "ext3_4", "ext3_5", "ext3_6", "ext3_7", "ext3_8", "ext3_9", "ext3_10")
# Reorder 'type' column based on 'sample_order'
df1$type <- factor(df1$type, levels = sample_order)
df2$type <- factor(df2$type, levels = sample_order)
# Combine df1 and df2 into one data frame for plotting
df_combined <- bind_rows(
  df1 %>% mutate(dataset = "Nuc"),
  df2 %>% mutate(dataset = "Cyt")
)
# Summarize the data (mean, SD, and error interval - EI)
df_summary <- df_combined %>%
  group_by(dataset, type) %>%
  summarise(
    mean_value = mean(value),
    sd_value = sd(value),
    n = n(),
    se_value = sd_value / sqrt(n),
    ei_lower = mean_value - 1 * se_value,  # Use Standard Error (1 * SE)
    ei_upper = mean_value + 1 * se_value   # Use Standard Error (1 * SE)
  )
# 计算每对相同 `type` 下 `df1` 和 `df2` 的显著性
significance_results <- df_combined %>%
  # 只选择每个 type 对应的 df1 和 df2 数据
  group_by(type) %>%
  summarise(
    # Mann-Whitney U 检验在每个 type 下对 df1 和 df2 进行
    mw_test_result = wilcox.test(value ~ dataset)$p.value,
    # 如果 p 值小于 0.05，则标记为星号
    significance = case_when(
      mw_test_result < 0.0001 ~ "****",  # p-value < 0.0001
      mw_test_result < 0.001  ~ "***",   # p-value < 0.001
      mw_test_result < 0.01   ~ "**",    # p-value < 0.01
      mw_test_result < 0.05   ~ "*",     # p-value < 0.05
      TRUE ~ NA_character_             # p-value >= 0.05, no significance
    )
  )
# 合并显著性结果到 df_summary 数据框
df_summary <- df_summary %>%
  left_join(significance_results %>%
              select(type, significance), by = "type")

# 仅选择每个 `type` 中 `mean_value` 最大的行来标记星号
df_summary_with_significance <- df_summary %>%
  group_by(type) %>%
  mutate(
    significance = ifelse(mean_value == max(mean_value), significance, NA)  # 只保留最高的点
  )
#
ymax=max(df_summary_with_significance$mean_value)

pdf(str_c(rbp_name,'.point_line.pdf'),wi=5,he=5)
ggplot(df_summary_with_significance, aes(x = type, y = mean_value, color = dataset, group = dataset)) +
  geom_point(size = 4, aes(color = dataset)) +  # 点使用 dataset 的颜色
  geom_vline(aes(xintercept='K562'),linewidth=1.2,color='#FFC102')+
  geom_line(aes(color = dataset), linetype = "solid") +  # 线使用 dataset 的颜色，设置线型为实线
  geom_errorbar(aes(ymin = ei_lower, ymax = ei_upper), width = 0.2) +  # 绘制误差条（EI）
  #  scale_x_discrete(limits = sample_order,labels = c('-10','-9',	'-8','-7','-6','-5','-4','-3','-2','-1','0','1','2','3','4','5','6','7','8','9','10')) +  # 确保按自定义顺序排列 x 轴   
  scale_x_discrete(limits = sample_order,labels = c('-10','',	'-8','','-6','','-4','','-2','','0','','2','','4','','6','','8','','10')) +
  scale_y_continuous(limits = c(0,0.5),breaks = seq(0,0.5,0.05))+
  scale_color_manual(values = cols) +
  labs(title = "",
       x = "", y = "norRMR") +
  theme_bw(base_rect_size = 1.5) +  # 应用黑白主题
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 24, face = 'bold', color = 'black'),
    axis.text = element_text(size = 16, face = 'bold', color = 'black'),
    legend.title = element_blank(),
    legend.text = element_text(size = 18, face = 'bold', color = 'black'),
    legend.spacing.y = unit(1,'cm'),         ### 调整legend上下距离
    legend.key.height = unit(.8, 'cm'),     ### 调整legend上下高度
    legend.position = c(0.85,0.8)
  ) + 
  # Add shaded region for K562 on the x-axis (light red shade) with no border
  geom_rect(data = df_summary %>% filter(type == "K562"), 
            aes(xmin = which(sample_order == "K562") - 2.5, 
                xmax = which(sample_order == "K562") + 3.5, 
                ymin = -Inf, ymax = Inf), 
            fill = "gray50", alpha = 0.1, color = NA) +  # Remove border by setting color = NA
  # 在每个显著性点上添加红色星号（*），并让星号靠近 ei_upper
  geom_text(data = df_summary_with_significance %>% filter(!is.na(significance) & (type == 'ext5_2' | type == 'ext5_1' | type == 'K562' | type == 'ext3_1' | type == 'ext3_2' | type == 'ext3_3')), 
            aes(label = significance, y = ei_upper + ymax/50), color = "darkblue", size = 10)
dev.off()
##############################################################################
rbp_name='ENCSR181NRW.ZC3H8'
df1=subset(data1,RBP==rbp_name)
df2=subset(data2,RBP==rbp_name)
# Define custom order for types
sample_order <- c("ext5_10","ext5_9","ext5_8","ext5_7","ext5_6","ext5_5","ext5_4", 
                  "ext5_3", "ext5_2", "ext5_1", "K562", "ext3_1", "ext3_2", "ext3_3", 
                  "ext3_4", "ext3_5", "ext3_6", "ext3_7", "ext3_8", "ext3_9", "ext3_10")
# Reorder 'type' column based on 'sample_order'
df1$type <- factor(df1$type, levels = sample_order)
df2$type <- factor(df2$type, levels = sample_order)
# Combine df1 and df2 into one data frame for plotting
df_combined <- bind_rows(
  df1 %>% mutate(dataset = "Nuc"),
  df2 %>% mutate(dataset = "Cyt")
)
# Summarize the data (mean, SD, and error interval - EI)
df_summary <- df_combined %>%
  group_by(dataset, type) %>%
  summarise(
    mean_value = mean(value),
    sd_value = sd(value),
    n = n(),
    se_value = sd_value / sqrt(n),
    #    ei_lower = mean_value - sd_value,  # Error Interval: mean - SD
    #    ei_upper = mean_value + sd_value   # Error Interval: mean + SD
    ei_lower = mean_value - 1 * se_value,  # Use Standard Error (1 * SE)
    ei_upper = mean_value + 1 * se_value   # Use Standard Error (1 * SE)
  )
# 计算每对相同 `type` 下 `df1` 和 `df2` 的显著性
significance_results <- df_combined %>%
  # 只选择每个 type 对应的 df1 和 df2 数据
  group_by(type) %>%
  summarise(
    # Mann-Whitney U 检验在每个 type 下对 df1 和 df2 进行
    mw_test_result = wilcox.test(value ~ dataset)$p.value,
    # 如果 p 值小于 0.05，则标记为星号
    significance = case_when(
      mw_test_result < 0.0001 ~ "****",  # p-value < 0.0001
      mw_test_result < 0.001  ~ "***",   # p-value < 0.001
      mw_test_result < 0.01   ~ "**",    # p-value < 0.01
      mw_test_result < 0.05   ~ "*",     # p-value < 0.05
      TRUE ~ NA_character_             # p-value >= 0.05, no significance
    )
  )
# 合并显著性结果到 df_summary 数据框
df_summary <- df_summary %>%
  left_join(significance_results %>%
              select(type, significance), by = "type")

# 仅选择每个 `type` 中 `mean_value` 最大的行来标记星号
df_summary_with_significance <- df_summary %>%
  group_by(type) %>%
  mutate(
    significance = ifelse(mean_value == max(mean_value), significance, NA)  # 只保留最高的点
  )
#
ymax=max(df_summary_with_significance$mean_value)

pdf(str_c(rbp_name,'.point_line.pdf'),wi=5,he=5)
ggplot(df_summary_with_significance, aes(x = type, y = mean_value, color = dataset, group = dataset)) +
  geom_point(size = 4, aes(color = dataset)) +  # 点使用 dataset 的颜色
  geom_vline(aes(xintercept='K562'),linewidth=1.2,color='#FFC102')+
  geom_line(aes(color = dataset), linetype = "solid") +  # 线使用 dataset 的颜色，设置线型为实线
  geom_errorbar(aes(ymin = ei_lower, ymax = ei_upper), width = 0.2) +  # 绘制误差条（EI）
#  scale_x_discrete(limits = sample_order,labels = c('-10','-9',	'-8','-7','-6','-5','-4','-3','-2','-1','0','1','2','3','4','5','6','7','8','9','10')) +  # 确保按自定义顺序排列 x 轴
  scale_x_discrete(limits = sample_order,labels = c('-10','',	'-8','','-6','','-4','','-2','','0','','2','','4','','6','','8','','10')) +
  scale_y_continuous(limits = c(0,0.5),breaks = seq(0,0.5,0.05))+
  scale_color_manual(values = cols) +
  labs(title = "",
       x = "", y = "norRMR") +
  theme_bw(base_rect_size = 1.5) +  # 应用黑白主题
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 24, face = 'bold', color = 'black'),
    axis.text = element_text(size = 16, face = 'bold', color = 'black'),
    legend.title = element_blank(),
    legend.text = element_text(size = 18, face = 'bold', color = 'black'),
    legend.spacing.y = unit(1,'cm'),         ### 调整legend上下距离
    legend.key.height = unit(.8, 'cm'),     ### 调整legend上下高度
    legend.position = c(0.85,0.8)
  ) + 
  # Add shaded region for K562 on the x-axis (light red shade) with no border
  geom_rect(data = df_summary %>% filter(type == "K562"), 
            aes(xmin = which(sample_order == "K562") - 2.5, 
                xmax = which(sample_order == "K562") + 2.5, 
                ymin = -Inf, ymax = Inf), 
            fill = "gray50", alpha = 0.1, color = NA) +  # Remove border by setting color = NA
  # 在每个显著性点上添加红色星号（*），并让星号靠近 ei_upper
  geom_text(data = df_summary_with_significance %>% filter(!is.na(significance) & (type == 'ext5_2' | type == 'ext5_1' | type == 'K562' | type == 'ext3_1' | type == 'ext3_2')), 
            aes(label = significance, y = ei_upper + ymax/50), color = "darkblue", size = 10)
dev.off()
##############################################################################
rbp_name='ENCSR295OKT.RBM22'
df1=subset(data1,RBP==rbp_name)
df2=subset(data2,RBP==rbp_name)
# Define custom order for types
sample_order <- c("ext5_10","ext5_9","ext5_8","ext5_7","ext5_6","ext5_5","ext5_4", 
                  "ext5_3", "ext5_2", "ext5_1", "K562", "ext3_1", "ext3_2", "ext3_3", 
                  "ext3_4", "ext3_5", "ext3_6", "ext3_7", "ext3_8", "ext3_9", "ext3_10")
# Reorder 'type' column based on 'sample_order'
df1$type <- factor(df1$type, levels = sample_order)
df2$type <- factor(df2$type, levels = sample_order)
# Combine df1 and df2 into one data frame for plotting
df_combined <- bind_rows(
  df1 %>% mutate(dataset = "Nuc"),
  df2 %>% mutate(dataset = "Cyt")
)
# Summarize the data (mean, SD, and error interval - EI)
df_summary <- df_combined %>%
  group_by(dataset, type) %>%
  summarise(
    mean_value = mean(value),
    sd_value = sd(value),
    n = n(),
    se_value = sd_value / sqrt(n),
    #    ei_lower = mean_value - sd_value,  # Error Interval: mean - SD
    #    ei_upper = mean_value + sd_value   # Error Interval: mean + SD
    ei_lower = mean_value - 1 * se_value,  # Use Standard Error (1 * SE)
    ei_upper = mean_value + 1 * se_value   # Use Standard Error (1 * SE)
  )
# 计算每对相同 `type` 下 `df1` 和 `df2` 的显著性
significance_results <- df_combined %>%
  # 只选择每个 type 对应的 df1 和 df2 数据
  group_by(type) %>%
  summarise(
    # Mann-Whitney U 检验在每个 type 下对 df1 和 df2 进行
    mw_test_result = wilcox.test(value ~ dataset)$p.value,
    # 如果 p 值小于 0.05，则标记为星号
    significance = case_when(
      mw_test_result < 0.0001 ~ "****",  # p-value < 0.0001
      mw_test_result < 0.001  ~ "***",   # p-value < 0.001
      mw_test_result < 0.01   ~ "**",    # p-value < 0.01
      mw_test_result < 0.05   ~ "*",     # p-value < 0.05
      TRUE ~ NA_character_             # p-value >= 0.05, no significance
    )
  )
# 合并显著性结果到 df_summary 数据框
df_summary <- df_summary %>%
  left_join(significance_results %>%
              select(type, significance), by = "type")

# 仅选择每个 `type` 中 `mean_value` 最大的行来标记星号
df_summary_with_significance <- df_summary %>%
  group_by(type) %>%
  mutate(
    significance = ifelse(mean_value == max(mean_value), significance, NA)  # 只保留最高的点
  )
#
ymax=max(df_summary_with_significance$mean_value)

pdf(str_c(rbp_name,'.point_line.pdf'),wi=5,he=5)
ggplot(df_summary_with_significance, aes(x = type, y = mean_value, color = dataset, group = dataset)) +
  geom_point(size = 4, aes(color = dataset)) +  # 点使用 dataset 的颜色
  geom_vline(aes(xintercept='K562'),linewidth=1.2,color='#FFC102')+
  geom_line(aes(color = dataset), linetype = "solid") +  # 线使用 dataset 的颜色，设置线型为实线
  geom_errorbar(aes(ymin = ei_lower, ymax = ei_upper), width = 0.2) +  # 绘制误差条（EI）
  #  scale_x_discrete(limits = sample_order,labels = c('-10','-9',	'-8','-7','-6','-5','-4','-3','-2','-1','0','1','2','3','4','5','6','7','8','9','10')) +  # 确保按自定义顺序排列 x 轴   
  scale_x_discrete(limits = sample_order,labels = c('-10','',	'-8','','-6','','-4','','-2','','0','','2','','4','','6','','8','','10')) +
  scale_y_continuous(limits = c(0,0.5),breaks = seq(0,0.5,0.05))+
  scale_color_manual(values = cols) +
  labs(title = "",
       x = "", y = "norRMR") +
  theme_bw(base_rect_size = 1.5) +  # 应用黑白主题
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 24, face = 'bold', color = 'black'),
    axis.text = element_text(size = 16, face = 'bold', color = 'black'),
    legend.title = element_blank(),
    legend.text = element_text(size = 18, face = 'bold', color = 'black'),
    legend.spacing.y = unit(1,'cm'),         ### 调整legend上下距离
    legend.key.height = unit(.8, 'cm'),     ### 调整legend上下高度
    legend.position = c(0.85,0.8)
  ) + 
  # Add shaded region for K562 on the x-axis (light red shade) with no border
  geom_rect(data = df_summary %>% filter(type == "K562"), 
            aes(xmin = which(sample_order == "K562") - 0.5, 
                xmax = which(sample_order == "K562") + 2.5, 
                ymin = -Inf, ymax = Inf), 
            fill = "gray50", alpha = 0.1, color = NA) +  # Remove border by setting color = NA
  # 在每个显著性点上添加红色星号（*），并让星号靠近 ei_upper
  geom_text(data = df_summary_with_significance %>% filter(!is.na(significance) & (type == 'K562' | type == 'ext3_1' | type == 'ext3_2')), 
            aes(label = significance, y = ei_upper + ymax/50), color = "darkblue", size = 10)
dev.off()
##############################################################################
rbp_name='ENCSR657TZB.XRN2'
df1=subset(data1,RBP==rbp_name)
df2=subset(data2,RBP==rbp_name)
# Define custom order for types
sample_order <- c("ext5_10","ext5_9","ext5_8","ext5_7","ext5_6","ext5_5","ext5_4", 
                  "ext5_3", "ext5_2", "ext5_1", "K562", "ext3_1", "ext3_2", "ext3_3", 
                  "ext3_4", "ext3_5", "ext3_6", "ext3_7", "ext3_8", "ext3_9", "ext3_10")
# Reorder 'type' column based on 'sample_order'
df1$type <- factor(df1$type, levels = sample_order)
df2$type <- factor(df2$type, levels = sample_order)
# Combine df1 and df2 into one data frame for plotting
df_combined <- bind_rows(
  df1 %>% mutate(dataset = "Nuc"),
  df2 %>% mutate(dataset = "Cyt")
)
# Summarize the data (mean, SD, and error interval - EI)
df_summary <- df_combined %>%
  group_by(dataset, type) %>%
  summarise(
    mean_value = mean(value),
    sd_value = sd(value),
    n = n(),
    se_value = sd_value / sqrt(n),
    #    ei_lower = mean_value - sd_value,  # Error Interval: mean - SD
    #    ei_upper = mean_value + sd_value   # Error Interval: mean + SD
    ei_lower = mean_value - 1 * se_value,  # Use Standard Error (1 * SE)
    ei_upper = mean_value + 1 * se_value   # Use Standard Error (1 * SE)
  )
# 计算每对相同 `type` 下 `df1` 和 `df2` 的显著性
significance_results <- df_combined %>%
  # 只选择每个 type 对应的 df1 和 df2 数据
  group_by(type) %>%
  summarise(
    # Mann-Whitney U 检验在每个 type 下对 df1 和 df2 进行
    mw_test_result = wilcox.test(value ~ dataset)$p.value,
    # 如果 p 值小于 0.05，则标记为星号
    significance = case_when(
      mw_test_result < 0.0001 ~ "****",  # p-value < 0.0001
      mw_test_result < 0.001  ~ "***",   # p-value < 0.001
      mw_test_result < 0.01   ~ "**",    # p-value < 0.01
      mw_test_result < 0.05   ~ "*",     # p-value < 0.05
      TRUE ~ NA_character_             # p-value >= 0.05, no significance
    )
  )
# 合并显著性结果到 df_summary 数据框
df_summary <- df_summary %>%
  left_join(significance_results %>%
              select(type, significance), by = "type")

# 仅选择每个 `type` 中 `mean_value` 最大的行来标记星号
df_summary_with_significance <- df_summary %>%
  group_by(type) %>%
  mutate(
    significance = ifelse(mean_value == max(mean_value), significance, NA)  # 只保留最高的点
  )
#
ymax=max(df_summary_with_significance$mean_value)

pdf(str_c(rbp_name,'.point_line.pdf'),wi=5,he=5)
ggplot(df_summary_with_significance, aes(x = type, y = mean_value, color = dataset, group = dataset)) +
  geom_point(size = 4, aes(color = dataset)) +  # 点使用 dataset 的颜色
  geom_vline(aes(xintercept='K562'),linewidth=1.2,color='#FFC102')+
  geom_line(aes(color = dataset), linetype = "solid") +  # 线使用 dataset 的颜色，设置线型为实线
  geom_errorbar(aes(ymin = ei_lower, ymax = ei_upper), width = 0.2) +  # 绘制误差条（EI）
  #  scale_x_discrete(limits = sample_order,labels = c('-10','-9',	'-8','-7','-6','-5','-4','-3','-2','-1','0','1','2','3','4','5','6','7','8','9','10')) +  # 确保按自定义顺序排列 x 轴   
  scale_x_discrete(limits = sample_order,labels = c('-10','',	'-8','','-6','','-4','','-2','','0','','2','','4','','6','','8','','10')) +
  scale_y_continuous(limits = c(0,0.5),breaks = seq(0,0.5,0.05))+
  scale_color_manual(values = cols) +
  labs(title = "",
       x = "", y = "norRMR") +
  theme_bw(base_rect_size = 1.5) +  # 应用黑白主题
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 24, face = 'bold', color = 'black'),
    axis.text = element_text(size = 16, face = 'bold', color = 'black'),
    legend.title = element_blank(),
    legend.text = element_text(size = 18, face = 'bold', color = 'black'),
    legend.spacing.y = unit(1,'cm'),         ### 调整legend上下距离
    legend.key.height = unit(.8, 'cm'),     ### 调整legend上下高度
    legend.position = c(0.85,0.8)
  ) + 
  # Add shaded region for K562 on the x-axis (light red shade) with no border
  geom_rect(data = df_summary %>% filter(type == "K562"), 
            aes(xmin = which(sample_order == "K562") - 0.5, 
                xmax = which(sample_order == "K562") + 0.5, 
                ymin = -Inf, ymax = Inf), 
            fill = "gray50", alpha = 0.1, color = NA) +  # Remove border by setting color = NA
  # 在每个显著性点上添加红色星号（*），并让星号靠近 ei_upper
  geom_text(data = df_summary_with_significance %>% filter(!is.na(significance) & (type == 'K562')), 
            aes(label = significance, y = ei_upper + ymax/50), color = "darkblue", size = 10)
dev.off()
###############################################################
rbp_name='ENCSR658IQB.SMNDC1'
df1=subset(data1,RBP==rbp_name)
df2=subset(data2,RBP==rbp_name)
# Define custom order for types
sample_order <- c("ext5_10","ext5_9","ext5_8","ext5_7","ext5_6","ext5_5","ext5_4", 
                  "ext5_3", "ext5_2", "ext5_1", "K562", "ext3_1", "ext3_2", "ext3_3", 
                  "ext3_4", "ext3_5", "ext3_6", "ext3_7", "ext3_8", "ext3_9", "ext3_10")
# Reorder 'type' column based on 'sample_order'
df1$type <- factor(df1$type, levels = sample_order)
df2$type <- factor(df2$type, levels = sample_order)
# Combine df1 and df2 into one data frame for plotting
df_combined <- bind_rows(
  df1 %>% mutate(dataset = "Nuc"),
  df2 %>% mutate(dataset = "Cyt")
)
# Summarize the data (mean, SD, and error interval - EI)
df_summary <- df_combined %>%
  group_by(dataset, type) %>%
  summarise(
    mean_value = mean(value),
    sd_value = sd(value),
    n = n(),
    se_value = sd_value / sqrt(n),
    #    ei_lower = mean_value - sd_value,  # Error Interval: mean - SD
    #    ei_upper = mean_value + sd_value   # Error Interval: mean + SD
    ei_lower = mean_value - 1 * se_value,  # Use Standard Error (1 * SE)
    ei_upper = mean_value + 1 * se_value   # Use Standard Error (1 * SE)
  )
# 计算每对相同 `type` 下 `df1` 和 `df2` 的显著性
significance_results <- df_combined %>%
  # 只选择每个 type 对应的 df1 和 df2 数据
  group_by(type) %>%
  summarise(
    # Mann-Whitney U 检验在每个 type 下对 df1 和 df2 进行
    mw_test_result = wilcox.test(value ~ dataset)$p.value,
    # 如果 p 值小于 0.05，则标记为星号
    significance = case_when(
      mw_test_result < 0.0001 ~ "****",  # p-value < 0.0001
      mw_test_result < 0.001  ~ "***",   # p-value < 0.001
      mw_test_result < 0.01   ~ "**",    # p-value < 0.01
      mw_test_result < 0.05   ~ "*",     # p-value < 0.05
      TRUE ~ NA_character_             # p-value >= 0.05, no significance
    )
  )
# 合并显著性结果到 df_summary 数据框
df_summary <- df_summary %>%
  left_join(significance_results %>%
              select(type, significance), by = "type")

# 仅选择每个 `type` 中 `mean_value` 最大的行来标记星号
df_summary_with_significance <- df_summary %>%
  group_by(type) %>%
  mutate(
    significance = ifelse(mean_value == max(mean_value), significance, NA)  # 只保留最高的点
  )
#
ymax=max(df_summary_with_significance$mean_value)

pdf(str_c(rbp_name,'.point_line.pdf'),wi=5,he=5)
ggplot(df_summary_with_significance, aes(x = type, y = mean_value, color = dataset, group = dataset)) +
  geom_point(size = 4, aes(color = dataset)) +  # 点使用 dataset 的颜色
  geom_line(aes(color = dataset), linetype = "solid") +  # 线使用 dataset 的颜色，设置线型为实线
  geom_vline(aes(xintercept='K562'),linewidth=1.2,color='#FFC102')+
  geom_errorbar(aes(ymin = ei_lower, ymax = ei_upper), width = 0.2) +  # 绘制误差条（EI）
  #  scale_x_discrete(limits = sample_order,labels = c('-10','-9',	'-8','-7','-6','-5','-4','-3','-2','-1','0','1','2','3','4','5','6','7','8','9','10')) +  # 确保按自定义顺序排列 x 轴   
  scale_x_discrete(limits = sample_order,labels = c('-10','',	'-8','','-6','','-4','','-2','','0','','2','','4','','6','','8','','10')) +
  scale_y_continuous(limits = c(0,0.5),breaks = seq(0,0.5,0.05))+
  scale_color_manual(values = cols) +
  labs(title = "",
       x = "", y = "norRMR") +
  theme_bw(base_rect_size = 1.5) +  # 应用黑白主题
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 24, face = 'bold', color = 'black'),
    axis.text = element_text(size = 16, face = 'bold', color = 'black'),
    legend.title = element_blank(),
    legend.text = element_text(size = 18, face = 'bold', color = 'black'),
    legend.spacing.y = unit(1,'cm'),         ### 调整legend上下距离
    legend.key.height = unit(.8, 'cm'),     ### 调整legend上下高度
    legend.position = c(0.85,0.8)
  ) + 
  # Add shaded region for K562 on the x-axis (light red shade) with no border
  geom_rect(data = df_summary %>% filter(type == "K562"), 
            aes(xmin = which(sample_order == "K562") - 4.5, 
                xmax = which(sample_order == "K562") + 1.5, 
                ymin = -Inf, ymax = Inf), 
            fill = "gray50", alpha = 0.1, color = NA) +  # Remove border by setting color = NA
  # 在每个显著性点上添加红色星号（*），并让星号靠近 ei_upper
  geom_text(data = df_summary_with_significance %>% filter(!is.na(significance) & (type=='ext5_4' | type=='ext5_3' | type=='ext5_2' | type=='ext5_1' | type == 'K562' | type == 'ext3_1')), 
            aes(label = significance, y = ei_upper + ymax/50), color = "darkblue", size = 10)
dev.off()
##############################################################################
rbp_name='ENCSR819XBT.AATF'
df1=subset(data1,RBP==rbp_name)
df2=subset(data2,RBP==rbp_name)
# Define custom order for types
sample_order <- c("ext5_10","ext5_9","ext5_8","ext5_7","ext5_6","ext5_5","ext5_4", 
                  "ext5_3", "ext5_2", "ext5_1", "K562", "ext3_1", "ext3_2", "ext3_3", 
                  "ext3_4", "ext3_5", "ext3_6", "ext3_7", "ext3_8", "ext3_9", "ext3_10")
# Reorder 'type' column based on 'sample_order'
df1$type <- factor(df1$type, levels = sample_order)
df2$type <- factor(df2$type, levels = sample_order)
# Combine df1 and df2 into one data frame for plotting
df_combined <- bind_rows(
  df1 %>% mutate(dataset = "Nuc"),
  df2 %>% mutate(dataset = "Cyt")
)
# Summarize the data (mean, SD, and error interval - EI)
df_summary <- df_combined %>%
  group_by(dataset, type) %>%
  summarise(
    mean_value = mean(value),
    sd_value = sd(value),
    n = n(),
    se_value = sd_value / sqrt(n),
    #    ei_lower = mean_value - sd_value,  # Error Interval: mean - SD
    #    ei_upper = mean_value + sd_value   # Error Interval: mean + SD
    ei_lower = mean_value - 1 * se_value,  # Use Standard Error (1 * SE)
    ei_upper = mean_value + 1 * se_value   # Use Standard Error (1 * SE)
  )
# 计算每对相同 `type` 下 `df1` 和 `df2` 的显著性
significance_results <- df_combined %>%
  # 只选择每个 type 对应的 df1 和 df2 数据
  group_by(type) %>%
  summarise(
    # Mann-Whitney U 检验在每个 type 下对 df1 和 df2 进行
    mw_test_result = wilcox.test(value ~ dataset)$p.value,
    # 如果 p 值小于 0.05，则标记为星号
    significance = case_when(
      mw_test_result < 0.0001 ~ "****",  # p-value < 0.0001
      mw_test_result < 0.001  ~ "***",   # p-value < 0.001
      mw_test_result < 0.01   ~ "**",    # p-value < 0.01
      mw_test_result < 0.05   ~ "*",     # p-value < 0.05
      TRUE ~ NA_character_             # p-value >= 0.05, no significance
    )
  )
# 合并显著性结果到 df_summary 数据框
df_summary <- df_summary %>%
  left_join(significance_results %>%
              select(type, significance), by = "type")

# 仅选择每个 `type` 中 `mean_value` 最大的行来标记星号
df_summary_with_significance <- df_summary %>%
  group_by(type) %>%
  mutate(
    significance = ifelse(mean_value == max(mean_value), significance, NA)  # 只保留最高的点
  )
#
ymax=max(df_summary_with_significance$mean_value)

pdf(str_c(rbp_name,'.point_line.pdf'),wi=5,he=5)
ggplot(df_summary_with_significance, aes(x = type, y = mean_value, color = dataset, group = dataset)) +
  geom_point(size = 4, aes(color = dataset)) +  # 点使用 dataset 的颜色
  geom_vline(aes(xintercept='K562'),linewidth=1.2,color='#FFC102')+
  geom_line(aes(color = dataset), linetype = "solid") +  # 线使用 dataset 的颜色，设置线型为实线
  geom_errorbar(aes(ymin = ei_lower, ymax = ei_upper), width = 0.2) +  # 绘制误差条（EI）
  #  scale_x_discrete(limits = sample_order,labels = c('-10','-9',	'-8','-7','-6','-5','-4','-3','-2','-1','0','1','2','3','4','5','6','7','8','9','10')) +  # 确保按自定义顺序排列 x 轴   
  scale_x_discrete(limits = sample_order,labels = c('-10','',	'-8','','-6','','-4','','-2','','0','','2','','4','','6','','8','','10')) +
  scale_y_continuous(limits = c(0,0.5),breaks = seq(0,0.5,0.05))+
  scale_color_manual(values = cols) +
  labs(title = "",
       x = "", y = "norRMR") +
  theme_bw(base_rect_size = 1.5) +  # 应用黑白主题
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 24, face = 'bold', color = 'black'),
    axis.text = element_text(size = 16, face = 'bold', color = 'black'),
    legend.title = element_blank(),
    legend.text = element_text(size = 18, face = 'bold', color = 'black'),
    legend.spacing.y = unit(1,'cm'),         ### 调整legend上下距离
    legend.key.height = unit(.8, 'cm'),     ### 调整legend上下高度
    legend.position = c(0.85,0.8)
  ) + 
  # Add shaded region for K562 on the x-axis (light red shade) with no border
  geom_rect(data = df_summary %>% filter(type == "K562"), 
            aes(xmin = which(sample_order == "K562") - 2.5, 
                xmax = which(sample_order == "K562") + 3.5, 
                ymin = -Inf, ymax = Inf), 
            fill = "gray50", alpha = 0.1, color = NA) +  # Remove border by setting color = NA
  # 在每个显著性点上添加红色星号（*），并让星号靠近 ei_upper
  geom_text(data = df_summary_with_significance %>% filter(!is.na(significance) & (type=='ext5_2' | type=='ext5_1' | type == 'K562' | type == 'ext3_1'| type == 'ext3_2'| type == 'ext3_3')), 
            aes(label = significance, y = ei_upper + ymax/50), color = "darkblue", size = 10)
dev.off()
##############################################################################
rbp_name='ENCSR840DRD.CSTF2T'
df1=subset(data1,RBP==rbp_name)
df2=subset(data2,RBP==rbp_name)
# Define custom order for types
sample_order <- c("ext5_10","ext5_9","ext5_8","ext5_7","ext5_6","ext5_5","ext5_4", 
                  "ext5_3", "ext5_2", "ext5_1", "K562", "ext3_1", "ext3_2", "ext3_3", 
                  "ext3_4", "ext3_5", "ext3_6", "ext3_7", "ext3_8", "ext3_9", "ext3_10")
# Reorder 'type' column based on 'sample_order'
df1$type <- factor(df1$type, levels = sample_order)
df2$type <- factor(df2$type, levels = sample_order)
# Combine df1 and df2 into one data frame for plotting
df_combined <- bind_rows(
  df1 %>% mutate(dataset = "Nuc"),
  df2 %>% mutate(dataset = "Cyt")
)
# Summarize the data (mean, SD, and error interval - EI)
df_summary <- df_combined %>%
  group_by(dataset, type) %>%
  summarise(
    mean_value = mean(value),
    sd_value = sd(value),
    n = n(),
    se_value = sd_value / sqrt(n),
    #    ei_lower = mean_value - sd_value,  # Error Interval: mean - SD
    #    ei_upper = mean_value + sd_value   # Error Interval: mean + SD
    ei_lower = mean_value - 1 * se_value,  # Use Standard Error (1 * SE)
    ei_upper = mean_value + 1 * se_value   # Use Standard Error (1 * SE)
  )
# 计算每对相同 `type` 下 `df1` 和 `df2` 的显著性
significance_results <- df_combined %>%
  # 只选择每个 type 对应的 df1 和 df2 数据
  group_by(type) %>%
  summarise(
    # Mann-Whitney U 检验在每个 type 下对 df1 和 df2 进行
    mw_test_result = wilcox.test(value ~ dataset)$p.value,
    # 如果 p 值小于 0.05，则标记为星号
    significance = case_when(
      mw_test_result < 0.0001 ~ "****",  # p-value < 0.0001
      mw_test_result < 0.001  ~ "***",   # p-value < 0.001
      mw_test_result < 0.01   ~ "**",    # p-value < 0.01
      mw_test_result < 0.05   ~ "*",     # p-value < 0.05
      TRUE ~ NA_character_             # p-value >= 0.05, no significance
    )
  )
# 合并显著性结果到 df_summary 数据框
df_summary <- df_summary %>%
  left_join(significance_results %>%
              select(type, significance), by = "type")

# 仅选择每个 `type` 中 `mean_value` 最大的行来标记星号
df_summary_with_significance <- df_summary %>%
  group_by(type) %>%
  mutate(
    significance = ifelse(mean_value == max(mean_value), significance, NA)  # 只保留最高的点
  )
#
ymax=max(df_summary_with_significance$mean_value)

pdf(str_c(rbp_name,'.point_line.pdf'),wi=5,he=5)
ggplot(df_summary_with_significance, aes(x = type, y = mean_value, color = dataset, group = dataset)) +
  geom_point(size = 4, aes(color = dataset)) +  # 点使用 dataset 的颜色
  geom_vline(aes(xintercept='K562'),linewidth=1.2,color='#FFC102')+
  geom_line(aes(color = dataset), linetype = "solid") +  # 线使用 dataset 的颜色，设置线型为实线
  geom_errorbar(aes(ymin = ei_lower, ymax = ei_upper), width = 0.2) +  # 绘制误差条（EI）
  #  scale_x_discrete(limits = sample_order,labels = c('-10','-9',	'-8','-7','-6','-5','-4','-3','-2','-1','0','1','2','3','4','5','6','7','8','9','10')) +  # 确保按自定义顺序排列 x 轴   
  scale_x_discrete(limits = sample_order,labels = c('-10','',	'-8','','-6','','-4','','-2','','0','','2','','4','','6','','8','','10')) +
  scale_y_continuous(limits = c(0,0.5),breaks = seq(0,0.5,0.05))+
  scale_color_manual(values = cols) +
  labs(title = "",
       x = "", y = "norRMR") +
  theme_bw(base_rect_size = 1.5) +  # 应用黑白主题
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 24, face = 'bold', color = 'black'),
    axis.text = element_text(size = 16, face = 'bold', color = 'black'),
    legend.title = element_blank(),
    legend.text = element_text(size = 18, face = 'bold', color = 'black'),
    legend.spacing.y = unit(1,'cm'),         ### 调整legend上下距离
    legend.key.height = unit(.8, 'cm'),     ### 调整legend上下高度
    legend.position = c(0.85,0.8)
  ) + 
  # Add shaded region for K562 on the x-axis (light red shade) with no border
  geom_rect(data = df_summary %>% filter(type == "K562"), 
            aes(xmin = which(sample_order == "K562") - 0.5, 
                xmax = which(sample_order == "K562") + 0.5, 
                ymin = -Inf, ymax = Inf), 
            fill = "gray50", alpha = 0.1, color = NA) +  # Remove border by setting color = NA
  # 在每个显著性点上添加红色星号（*），并让星号靠近 ei_upper
  geom_text(data = df_summary_with_significance %>% filter(!is.na(significance) & (type=='ext5_2' | type=='ext5_1' | type == 'K562' | type == 'ext3_1'| type == 'ext3_2'| type == 'ext3_3')), 
            aes(label = significance, y = ei_upper + ymax/50), color = "darkblue", size = 10)
dev.off()
##############################################################################
rbp_name='ENCSR947JVR.DGCR8'
df1=subset(data1,RBP==rbp_name)
df2=subset(data2,RBP==rbp_name)
# Define custom order for types
sample_order <- c("ext5_10","ext5_9","ext5_8","ext5_7","ext5_6","ext5_5","ext5_4", 
                  "ext5_3", "ext5_2", "ext5_1", "K562", "ext3_1", "ext3_2", "ext3_3", 
                  "ext3_4", "ext3_5", "ext3_6", "ext3_7", "ext3_8", "ext3_9", "ext3_10")
# Reorder 'type' column based on 'sample_order'
df1$type <- factor(df1$type, levels = sample_order)
df2$type <- factor(df2$type, levels = sample_order)
# Combine df1 and df2 into one data frame for plotting
df_combined <- bind_rows(
  df1 %>% mutate(dataset = "Nuc"),
  df2 %>% mutate(dataset = "Cyt")
)
# Summarize the data (mean, SD, and error interval - EI)
df_summary <- df_combined %>%
  group_by(dataset, type) %>%
  summarise(
    mean_value = mean(value),
    sd_value = sd(value),
    n = n(),
    se_value = sd_value / sqrt(n),
    #    ei_lower = mean_value - sd_value,  # Error Interval: mean - SD
    #    ei_upper = mean_value + sd_value   # Error Interval: mean + SD
    ei_lower = mean_value - 1 * se_value,  # Use Standard Error (1 * SE)
    ei_upper = mean_value + 1 * se_value   # Use Standard Error (1 * SE)
  )
# 计算每对相同 `type` 下 `df1` 和 `df2` 的显著性
significance_results <- df_combined %>%
  # 只选择每个 type 对应的 df1 和 df2 数据
  group_by(type) %>%
  summarise(
    # Mann-Whitney U 检验在每个 type 下对 df1 和 df2 进行
    mw_test_result = wilcox.test(value ~ dataset)$p.value,
    # 如果 p 值小于 0.05，则标记为星号
    significance = case_when(
      mw_test_result < 0.0001 ~ "****",  # p-value < 0.0001
      mw_test_result < 0.001  ~ "***",   # p-value < 0.001
      mw_test_result < 0.01   ~ "**",    # p-value < 0.01
      mw_test_result < 0.05   ~ "*",     # p-value < 0.05
      TRUE ~ NA_character_             # p-value >= 0.05, no significance
    )
  )
# 合并显著性结果到 df_summary 数据框
df_summary <- df_summary %>%
  left_join(significance_results %>%
              select(type, significance), by = "type")

# 仅选择每个 `type` 中 `mean_value` 最大的行来标记星号
df_summary_with_significance <- df_summary %>%
  group_by(type) %>%
  mutate(
    significance = ifelse(mean_value == max(mean_value), significance, NA)  # 只保留最高的点
  )
#
ymax=max(df_summary_with_significance$mean_value)

pdf(str_c(rbp_name,'.point_line.pdf'),wi=5,he=5)
ggplot(df_summary_with_significance, aes(x = type, y = mean_value, color = dataset, group = dataset)) +
  geom_point(size = 4, aes(color = dataset)) +  # 点使用 dataset 的颜色
  geom_line(aes(color = dataset), linetype = "solid") +  # 线使用 dataset 的颜色，设置线型为实线
  geom_errorbar(aes(ymin = ei_lower, ymax = ei_upper), width = 0.2) +  # 绘制误差条（EI）
  #  scale_x_discrete(limits = sample_order,labels = c('-10','-9',	'-8','-7','-6','-5','-4','-3','-2','-1','0','1','2','3','4','5','6','7','8','9','10')) +  # 确保按自定义顺序排列 x 轴   
  scale_x_discrete(limits = sample_order,labels = c('-10','',	'-8','','-6','','-4','','-2','','0','','2','','4','','6','','8','','10')) +
  scale_y_continuous(limits = c(0,0.5),breaks = seq(0,0.5,0.05))+
  scale_color_manual(values = cols) +
  labs(title = "",
       x = "", y = "norRMR") +
  theme_bw(base_rect_size = 1.5) +  # 应用黑白主题
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 24, face = 'bold', color = 'black'),
    axis.text = element_text(size = 16, face = 'bold', color = 'black'),
    legend.title = element_blank(),
    legend.text = element_text(size = 18, face = 'bold', color = 'black'),
    legend.spacing.y = unit(1,'cm'),         ### 调整legend上下距离
    legend.key.height = unit(.8, 'cm'),     ### 调整legend上下高度
    legend.position = c(0.85,0.8)
  ) + 
  # Add shaded region for K562 on the x-axis (light red shade) with no border
  geom_rect(data = df_summary %>% filter(type == "K562"), 
            aes(xmin = which(sample_order == "K562") - 0.5, 
                xmax = which(sample_order == "K562") + 0.5, 
                ymin = -Inf, ymax = Inf), 
            fill = "gray50", alpha = 0.1, color = NA) +  # Remove border by setting color = NA
  # 在每个显著性点上添加红色星号（*），并让星号靠近 ei_upper
  geom_text(data = df_summary_with_significance %>% filter(!is.na(significance) & (type == 'K562')), 
            aes(label = significance, y = ei_upper + ymax/50), color = "darkblue", size = 10)
dev.off()
