library(ggplot2)
library(ggpubr)
library(stringr)
library(RColorBrewer)
library(patchwork)
cols <- brewer.pal(9, "Set1")
cols = brewer.pal(10,"Paired")[1:10]
#pal <- colorRampPalette(cols)
#colors = pal(3)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors,rownames(qual_col_pals)))
# '#FC8D62' '#8DA0CB' '#E78AC3' '#A6D854' '#FFD92F' '#E5C494' '#B3B3B3' '#8DD3C7' '#FFFFB3' '#BEBADA' '#FB8072'
#  '#7FC97F' '#1B9E77' '#1F78B4' '#7570B3' '#BF5B17' '#FDC086' '#A6CEE3'
#
col=col_vector[c(5,6,9)]
col=c('#62CE0C','#F0027F','#1F78B4')
###############################

data1=read.table(file="NES.r1r2.chr.site.unique.gene",header=F,sep='\t')
data2=read.table(file="NLS.r1r2.chr.site.unique.gene",header=F,sep='\t')
data3=read.table(file="ERM.r1r2.chr.site.unique.gene",header=F,sep='\t')
d1=data1[,12]
d2=data2[,12]
d3=data3[,12]
groups_order=c("Cyt","Erm","Nuc")
data <- data.frame(
		value=c(d1,d3,d2),
		group= factor(rep(groups_order,times=c(length(d1),length(d3),length(d2))),levels = groups_order)
)

# col=col_vector[1:4]
pdf("NES.NLS.ERM.r1r2.chr.site.unique.gene.pdf",w=6,h=6)
ggplot(data, aes(x = value, color = group)) +
  geom_density(alpha = 1,linewidth=1,show_guide=FALSE) +  # Density plot with transparency
  stat_density(aes(x=value, colour=group),geom='line',position='identity',size=1.1)+  
  labs(title = '',
       x = 'RMR',
       y = 'Density') +  # Title and axis labels
  theme_minimal() +  # Minimal theme for a clean look
  scale_fill_brewer(palette = 'Set2') +  # Color palette
  theme_bw(base_rect_size = 1.5) +  
  theme(
    panel.border = element_rect(color = "black", linewidth = 1.5),
    axis.title = element_text(size = 18, face = 'bold',color='black'),
    axis.text = element_text(size = 14, face = 'bold',color='black'),
    legend.text = element_text(size=20,face = 'bold',color='black'),
    legend.position = c(0.8,0.8),
    legend.title = element_blank(),
    legend.ticks = element_line(size=1.5),
    axis.ticks = element_line(size=1.5),
    axis.ticks.length = unit(0.2,'cm'),
    panel.grid = element_blank(),
  )+
  scale_y_continuous()+
  scale_x_continuous(limits=c(0,1),)+
  scale_color_manual(values=col)

dev.off()
