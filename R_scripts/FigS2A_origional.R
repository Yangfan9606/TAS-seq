setwd('D:/Data/0_lixia/0)_8e_HTseq_all/0_Paper_Figure_Data_source_2024/0_d10mu2/051_sample_Rate_correlation_PLOT')
library(GGally)
library(ggplot2)
library(dplyr)
library(tidyr)
file1 <- read.table("NES.Padj_005.d10mu2.A",header = F,sep = '\t')
file2 <- read.table("NLS.Padj_005.d10mu2.A",header = F,sep = '\t')
file3 <- read.table("ERM.Padj_005.d10mu2.A",header = F,sep = '\t')

file1 <- file1 %>%
  mutate(ID = paste(V2, V3, V4, sep = "_"))
file2 <- file2 %>%
  mutate(ID = paste(V2, V3, V4, sep = "_"))
file3 <- file3 %>%
  mutate(ID = paste(V2, V3, V4, sep = "_"))

d1=file1[,c(6,9,13)]
d2=file2[,c(6,9,13)]
d3=file3[,c(6,9,13)]

names(d1)=c('NESr1','NESr2','ID')
names(d2)=c('NLSr1','NLSr2','ID')
names(d3)=c('ERMr1','ERMr2','ID')

# 合并数据框，保留所有匹配的ID
merged_data <- full_join(d1, d2, by = "ID") %>%
  full_join(d3, by = "ID") 
head(merged_data)
# 提取感兴趣的列（ID 和 rate）
data <- merged_data %>%
  select(ID,
         Cyt_Repl = NESr1, Cyt_Rep2 = NESr2,
         Nuc_Repl = NLSr1, Nuc_Rep2 = NLSr2,
         Erm_Repl = ERMr1, Erm_Rep2 = ERMr2,
  )

sample_df <- data %>% select(-ID)
head(sample_df)

# 自定义的散点图函数，包含对角线
custom_points_with_diag <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +  # 确保正确传递 data 和 mapping
    geom_point(size = 2, color = '#277fb8', alpha = 0.2)+    # 设置点的大小、颜色、透明度
    geom_smooth(color = 'red', linetype = 'solid',size=1,method = 'lm',se=F)+
    scale_x_continuous(limits = c(0,1),breaks = seq(0,1,0.2))+
    scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.2))+
    theme_classic()+
    theme(
      axis.line = element_line(size=1.2, colour = 'black'),
      axis.ticks = element_line(linetype = 1, color = 'black', linewidth = 1.2),
      axis.ticks.length = unit(0.2, 'cm'),
    )
}
#
# 自定义密度图
custom_density_with_diag <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_density(fill = '#B3DE69', color = 'gray', alpha = 0.7) +
    theme_classic() +
    scale_x_continuous(limits = c(0,1),breaks = seq(0,1,0.2))+
    theme(
      axis.line = element_line(size = 1.2, colour = 'black'),
      axis.ticks = element_line(linetype = 1, color = 'black', linewidth = 1.2),
      axis.ticks.length = unit(0.2, 'cm'),
    )
}
############################################
sample_df <- data %>% select(-ID)

# Matrix of plots
p1 <- ggpairs(sample_df,
              lower = list(continuous = wrap(custom_points_with_diag)),
              upper = list(continuous = wrap('cor', size = 5,color='black')),
              diag = list(continuous = wrap(custom_density_with_diag))
)  
# 获取列名并镜像反转
reversed_columns <- rev(colnames(sample_df))
#
p2 <- ggcorr(sample_df[, reversed_columns], label = TRUE, label_round = 2, low = "white", high = "#FF2052",
             limits = c(0.5, 1), label_size = 8, legend.size = 18, label_color = "black", midpoint = 0.6) +
  scale_x_discrete(limits = reversed_columns) +  # 反转 x 轴顺序
  theme(
    legend.key.size = unit(1, 'cm'),
    legend.key.height = unit(2.5, 'cm'),
    legend.position = 'right'
  )

g2 <- ggplotGrob(p2)
colors <- g2$grobs[[6]]$children[[3]]$gp$fill
p=length(sample_df)
colors = rev(colors)  ###  注意颜色顺序

# Change background color to tiles in the upper triangular matrix of plots 
idx <- 1
for (k1 in 1:(p-1)) {
  for (k2 in (k1+1):p) {
    plt <- getPlot(p1,k1,k2) +
      theme(
        panel.background = element_rect(fill = colors[idx], color='black'),
        panel.grid.major = element_line(color=colors[idx])
        )
    p1 <- putPlot(p1,plt,k1,k2)
    idx <- idx+1
  }
}

pdf('No_CytH_all.8Re.test.head.pdf',wi=16,he=16)
p1
dev.off()
