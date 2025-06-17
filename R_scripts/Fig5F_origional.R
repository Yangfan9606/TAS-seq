setwd('D:/Data/0_lixia/0)_8e_HTseq_all/0_Paper_Figure_Data_source_2024/0_d10mu2/062_variation_site_on_gene_level_PLOT')
f1=read.table('var1.txt',header = T,sep='\t')
f2=read.table('var2.txt',header = T,sep='\t')
f3=read.table('var3.txt',header = T,sep='\t')
f4=read.table('var4.txt',header = T,sep='\t')
f5=read.table('var5.txt',header = T,sep='\t')
d1=f1$Nuc.Cyt.Var..HEK293T.
d2=f2$Cyt.Erm.Var..HEK293T.
d3=f3$Nuc.Cyt.Var..K562.
d4=f4$Nuc_Cell.line.Var
d5=f5$Cyt_Cell.line.Var
data=data.frame(
  value=c(d1,d2,d3,d4,d5),
  group=rep(c('Nuc/Cyt (HEK293T)','Cyt/Erm (HEK293T)','Nuc/Cyt (K562)','Nuc Cell_line','Cyt Cell_line'),times=c(nrow(f1),nrow(f2),nrow(f3),nrow(f4),nrow(f5)))
)
data$group <- factor(data$group, levels = c('Nuc/Cyt (HEK293T)','Cyt/Erm (HEK293T)','Nuc/Cyt (K562)','Nuc Cell_line','Cyt Cell_line'))
data$value=log2(data$value)
#

p <- ggplot(data, aes(x = group, y = value, fill = group)) +
  geom_violin(trim = F) +
  geom_boxplot(width = 0.1, color = "black", outlier.shape = NA, fill='white',linewidth=1.2) +  # 在中间添加箱线图，隐藏离群点
  ###  stat_summary(fun = mean, geom = "point", color = "black", size = 2) +  # 添加平均数点
  stat_summary(fun = median, geom = "crossbar", width = 0.1, color = "black") +  # 添加中位数线
  scale_fill_manual(values =c("#D7D4D5", "#949091", "#8DBC9D", "#E08692", "#9397C6"))+
  theme_minimal() +  # 使用简约主题
  labs(title = "", x = "", y = "Variation sites number \non gene level (log2)") +  # 添加标签
  theme_bw()+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1.5),
        axis.title = element_text(size = 20, face = 'bold', color = 'black'),
        axis.text.y = element_text(size = 16, face = 'bold', color = 'black'),
        axis.text.x = element_text(,size = 16, angle = 45, vjust = .5,face = 'bold', color = 'black'),
        axis.line = element_line(linetype = 1, color = 'black', linewidth = 1),
        axis.ticks = element_line(linetype = 1, color = 'black', linewidth = 1),
        axis.ticks.length = unit(0.25, 'cm'),
        legend.text = element_text(size = 20, face = 'bold', color = 'black'),
        legend.title = element_text(size = 24, color = 'black'),
        plot.margin = margin(0,1,0.2,0.2, unit = "cm"),
  )+  
  scale_x_discrete(labels=c('Nuc/Cyt Var\n(HEK293T)','Cyt/Erm Var\n(HEK293T)','Nuc/Cyt Var\n(K562)','Nuc_Cell line Var','Cyt_Cell line Var'))
#  scale_y_continuous(limits = c(0,50),breaks = seq(0,50,10))

pdf('Variation_sites_gene_level.pdf',wi=6,he=6)
p
dev.off()
