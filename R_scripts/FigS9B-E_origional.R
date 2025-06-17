library(ggplot2)
library(ggpubr)
library(dplyr)
setwd('D:/Data/0_lixia/0)_8e_HTseq_all/0_Paper_Figure_Data_source_2024/0_Review_CellRepotyMethod/gRNA_/PLOT')
file=read.table('A_number_binding activity-2.txt',header=F,sep='\t')
f1=file$V1
f2=file$V2[which(file$V2!='NA')]
f3=file$V3[which(file$V3!='NA')]
f4=file$V4[which(file$V4!='NA')]
f5=file$V5[which(file$V5!='NA')]

data <- data.frame(
  value = c(f1,f2,f3,f4,f5),
  group = factor(
    rep(
      c("All", 'A1-4_num-0',	'A1-4_num-1(2)',	'A1-4_num-2(2,3)','A1-4_num-3(2,3,4)'),
      times = c(length(f1),length(f2),length(f3),length(f4),length(f5))
    )
  )
)
data$group <- factor(data$group, levels = c("All", 'A1-4_num-0',	'A1-4_num-1(2)',	'A1-4_num-2(2,3)','A1-4_num-3(2,3,4)'))
my_comparisons = list(c("All", 'A1-4_num-0'),c("All", 'A1-4_num-3(2,3,4)'))

data <- data %>%
  group_by(group) %>%
  mutate(n = n()) %>%  # 计算每组的样本数
  ungroup()

p <- ggboxplot(data, x = 'group', y = 'value', fill = 'group',
               color = 'black', 
               size=1.5, width = 0.8,
               bxp.errorbar = T,
               bxp.linewidth=1,
               outlier.shape = 16,
) +
  geom_beeswarm(data=subset(data,n<20),size=4,priority = "density",cex=0.8,shape = 16,color='blue')+
  stat_summary(fun = median, geom = "crossbar", width = 0.2, color = "black") +  # 添加中位数线
####  scale_fill_manual(values = c("NLS_Intron" = Intron_col, "NLS_exon" = Exon_col, "total_NLS_Intron" = Intron_col, "total_NLS_exon" = Exon_col,"K562_NLS_Intron" = Intron_col, "K562_NLS_exon" = Exon_col)) +  # 自定义颜色
  theme_minimal() +  # 使用简约主题
  labs(title = "", 
       x = "Modified As", 
       y = expression(bold(italic(Ka)*" (1/"*italic(Kd)*")"))) +  # 添加标签
  theme_bw()+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1.5),
        axis.title = element_text(size = 24, face = 'bold', color = 'black'),
        axis.text.y = element_text(size = 18, face = 'bold', color = 'black'),
        axis.text.x = element_text(angle = 45,hjust = 1,size = 20, face = 'bold', color = 'black'),
        axis.line = element_line(linetype = 1, color = 'black', linewidth = 1),
        axis.ticks = element_line(linetype = 1, color = 'black', linewidth = 1),
        axis.ticks.length = unit(0.25, 'cm'),
        legend.text = element_text(size = 20, face = 'bold', color = 'black'),
        legend.title = element_text(size = 24, color = 'black'),
        plot.margin = margin(0,1,0.2,0.2, unit = "cm"),
  )+  # 隐藏图例
  scale_y_continuous(limits = c(0,0.4),breaks = seq(0,0.35,0.05))+
  scale_x_discrete(labels = c(
    "All" = expression(bold('All RNA targets')),
    "A1-4_num-0"= expression(bold(None)),
    'A1-4_num-1(2)'= expression(bold(A[18])),
    'A1-4_num-2(2,3)'= expression(bold(A[18]+A[16])),
    'A1-4_num-3(2,3,4)'= expression(bold(A[18]+A[16]+A[14]))
  )
  )+
  stat_compare_means(
    comparisons = my_comparisons, 
    label.y = c(0.35,0.38),
#    label.y = c(0.39),
    method = "wilcox.test", tip.length = 0.01, bracket.size = 1, step.increase = 0, size = 12, vjust = 0.5,
    label='p.signif'
    )

pdf('A_number_binding activity-2.violin.pdf',wi=6,he=8)
print(p)
dev.off()
