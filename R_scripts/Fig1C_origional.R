library(ggpubr)
library(stringr)
library(plyr)
library(scales)
###############################
data1=read.table('NES.Padj_005.d10mu2.A',header=F,sep='\t')
data2=read.table('NES.Padj_005.d10mu2.A',header=F,sep='\t')

d1=data1[,6]
d2=data2[,9]

d1=(d1)
d2=(d2)

df <- data.frame(d1, d2)
colnames(df) <- c('Rep1', 'Rep2')

pdf("NES.R1R2.density2d.pdf", w=7, h=6)

name1='Rep1'
name2='Rep2'

ggscatter(df, x = 'Rep1', y = 'Rep2',
          alpha = 0,
          add = 'none',
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = 'pearson', label.x = 0.05, label.sep = '\n',size = 6, face = 'bold', color='black')
)+
  geom_hex(bins = 100,)+
  scale_fill_gradientn(
    colours=c("lightgray","#451077","#721F81","#9F2F7F","#CD4071","#F1605D","#FD9567","#FEC98D" ,"#FCFDBF"
    )
  )+
  theme_bw(base_rect_size = 1.5)+
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 20, face = 'bold',color='black'),
    axis.text = element_text(size = 16, face = 'bold',color='black'),
    axis.ticks = element_line(linewidth = 1.5),
    axis.ticks.length = unit(0.2,'cm'),
    plot.margin = margin(0.2,0.5,0.2,0.2, unit = 'cm'),
    legend.title = element_text(size = 16, face = 'bold',color='black'),
    legend.text = element_text(size = 14, face = 'bold',color='black'),
    #    legend.position = 'none'
  )+
  scale_x_continuous(name = name1,limits=c(0,1),breaks=seq(0,1,by=0.1),expand = c(0,0.01))+
  scale_y_continuous(name = name2,limits=c(0,1),breaks=seq(0,1,by=0.1),expand = c(0.01,0))

dev.off()
