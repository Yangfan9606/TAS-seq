setwd('D:/Data/0_lixia/0)_8e_HTseq_all/0_Paper_Figure_Data_source_2024/000_Manuscript/task/0_Download_DATA/6_overlap_boxplot_with 8e_PLOT')
library(ggplot2)
library(ggpubr)
library(stringr)
library(RColorBrewer)
library(patchwork)
cols <- brewer.pal(9, "Set1")
cols = brewer.pal(10,"Paired")[1:10]
#pal <- colorRampPalette(cols)
#colors = pal(1)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors,rownames(qual_col_pals)))
col=col_vector[9:13]
##################
sample='NES'
#name='2014_K562_DMS'
#marker=c('*','****','****')
name='2016_icSHAPE_HS'
marker=c('****','****','****')
col=c('#66C2A3','#8DA1CC','#FC8C62','#FEF6C8')
name='2019_icSHAPE_cy'
marker=c('****','****','****')
col=c('#EBF4E5','#BADCEC','#7CB0D4','#4476B2')
#name='2019_icSHAPE_np'
#name='2021_icSHAPE_HS'
#name='2021_K562_icSHAPE_HS'
#name='2022_293T_smartSHAPE'
#name='DMSmap'
#marker=c('****','****','****')
#col=c('#BFEDC7','#BFE6E9','#A1EEE1','#1DCAAC')


data1=read.table(str_c(sample,'.',name,".overlap.RMR"),sep='\t',header=T)

data1=data1[,c(1,2)]
names(data1)=c('Group','value')
min = 0
max = quantile(data1$value,.75) + sd(data1$value)*3
top=quantile(data1$value,.95)
### 2014_K562_DMS
#p_pos=c(top-top*0.15,top,top+top*0.15)
#m_pos=c(top-top*0.05,top+top*0.1,top+top*0.25)
#pdf(str_c(sample,'.',name,".overlap.RMR.pdf"),w=5,h=7)
###  2016_icSHAPE_HS 
p_pos=c(top-top*0.15,top,top+top*0.15)
m_pos=c(top-top*0.025,top+top*0.125,top+top*0.275)
pdf(str_c(sample,'.',name,".overlap.RMR.pdf"),w=5,h=7)
ylab_m=1
### 2019_icSHAPE_cy
pdf(str_c(sample,'.',name,".overlap.RMR.pdf"),w=5,h=8)
p_pos=c(top-top*0.24,top-top*0.11,top+top*0.02)
m_pos=c(top-top*0.1,top+top*0.03,top+top*0.16)
ylab_m=0.25
#######  only DMSmap use
#pdf(str_c(sample,'.',name,".overlap.RMR.pdf"),w=5,h=8)
#max = quantile(data1$value,.75) + sd(data1$value)*3
#max=0.5139107 
#top=quantile(data1$value,.95)
#p_pos=c(0.293694, 0.336149, 0.378604)
#m_pos=c(0.314167, 0.356622, 0.399077)
#ylab_m = 0.35
#####################################

order=as.character(as.data.frame(table(data1$Group))$Var1)
my_comparisons = list(c('0 ~ 0.25','0.25 ~ 0.5'),c('0 ~ 0.25','0.5 ~ 0.75'), c('0 ~ 0.25','0.75 ~ 1'))

p <- ggboxplot(data1, x = 'Group', y = 'value',
               color = 'black', 
	       fill = 'Group',
               palette = col,
               size=1.5, width = 0.5,
              outlier.shape = NA,
              bxp.errorbar = T,
              bxp.linewidth=1,
#             bxp.errorbar.width=0.25,
#             label.select = list(top.up = 10, top.down = 4),
              xlab = '',ylab='Delta rate'
)+
  scale_x_discrete(limits = order)+
  theme_classic()+
  theme(
    axis.line = element_line(size=1.2, colour = 'black'),
    axis.title.y = element_text(size=18,face = 'bold',color='black'),
    axis.text.x = element_text(size=12,face = 'bold',color='black'),
    axis.text.y = element_text(size=14,face = 'bold',color='black'),
#   legend.text = element_text(size=16,face = 'bold',color='black'),
    legend.position = 'none',
    axis.ticks = element_line(linetype = 1, color = 'black', linewidth = 1.5),
    axis.ticks.length = unit(0.25, 'cm'),
  )+
  scale_y_continuous(name = 'RMR', limits = c(0, max), breaks=seq(0,ylab_m,by=0.05))+
  annotate('text', x = 1:length(table(data1$Group)), y = 0,
           label = paste('n=',table(data1$Group)[order]), col = 'black', vjust = 1.5,fontface = "bold"
           )+
  stat_compare_means(comparisons = my_comparisons,
                     label.y=p_pos,   #  注意设置标签的y坐标位置
                     method='wilcox.test',tip.length = 0.005, bracket.size = 1.2,step.increase = 0,
                     label = 'p.format',size=5,vjust = 0,fontface = "bold"
                     )+
  geom_signif(
    y_position = m_pos, xmin = c(1, 1 , 1), xmax = c(2, 3, 4),
    annotation = marker,  linetype = 'blank', tip_length = 0, textsize = 5,fontface = "bold"
    )
p
dev.off()
###########################
