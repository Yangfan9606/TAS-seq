library(RColorBrewer)
library(gridExtra)
library(stringr)
library(ggplot2)
library(dplyr)
library(reshape2)

setwd('D:/Data/0_lixia/0)_8e_HTseq_all/0_Paper_Figure_Data_source_2024/0_single_cell/Final_PLOT_of_cosine_and_R2-HEAT_map')

name='BUB3'
name='MIF'  # 更换name 记得修改 坐标轴刻度

data <- read.table(str_c(name,'.SC_NES_1_10.adj_R2.txt'),header=F,sep = '\t')
####################################################################
heatmap_data <- data[, 1:10]
heatmap_data_matrix <- as.matrix(heatmap_data[, 1:10])
hm=melt(heatmap_data_matrix,value.name='value')
colnames(hm)=c('base','cell','value')
m <- tidyr::pivot_wider(hm, names_from = 'cell', values_from = 'value')
m <- as.matrix(m[, -1]) # -1 to omit categories from matrix
# Cluster based on euclidean distance

clust <- hclust(dist(t(m),method='euclidean'))
clust$order
clust <- hclust(dist(t(m),method='maximum'))
clust$order
clust <- hclust(dist(t(m),method='manhattan'))
clust$order
clust <- hclust(dist(t(m),method='canberra'))
clust$order
clust <- hclust(dist(t(m),method='binary'))
clust$order
clust <- hclust(dist(t(m),method='minkowski'))
clust$order
clust <- hclust(dist(t(m),method='euclidean'),method = 'complete')
######################
R2_values <- data[, 11]
R2_data <- data.frame(Base = 1:length(R2_values), R2 = R2_values)
##################################################################
heatmap_grob <- ggplot(hm, aes(cell,base, fill=value)) +
  geom_tile() +
  scale_x_discrete(limits = colnames(m)[clust$order],expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  scale_fill_gradientn(
    colors = colorRampPalette(c('#6a99d0', '#F7F7F7', 'red'))(100),
    limits = c(0,1),
    name = "RMR",
    guide = guide_colorbar(
      title.position = "left", 
      title.theme = element_text(
        angle = 90,
        hjust = 0.5,   # 水平居中
        )
      ),
    ) +
  theme_classic()+
  labs(title = str_c('Modification on ',name,' mRNA'),y='Modified A sites',x='Cells')+
  theme(
    axis.text.x = element_blank(),
    legend.position = c(-0.15,0),
    legend.justification = 'bottom',
    legend.direction = 'vertical',
    legend.title = element_text(size=30, vjust = 1.5, hjust = .1,face = 'bold', color = 'black'),
    legend.text = element_text(size=24, face = 'bold', color = 'black'),
    axis.title = element_text(size = 30, face = 'bold', color = 'black'),
    title = element_text(size=22,hjust = 0.5, face = 'bold', color = 'black'),
    plot.margin = unit(c(0.25, 0.1, 0.2, 3.1), 'cm'),
  )
heatmap_grob
#
R2_plot <- ggplot(R2_data, aes(x = Base, y = R2)) +
  geom_line() +
  geom_hline(yintercept = 0.5, linetype = 'dashed', color = 'red',linewidth = 1) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 18, face = 'bold', color = 'black'),
    title = element_text(size=16,face = 'bold', color = 'black'),
    axis.ticks = element_line(linetype = 1, color = 'black', linewidth = 1.2),
    axis.ticks.length = unit(0.2, 'cm'),
    plot.margin = unit(c(0.2, 0, -0.1, 0), 'cm'),
  )+
  scale_y_continuous(expand = c(0,0),limits = c(-0.1,1.1),breaks = seq(0,1,0.5))+   # MIF
  scale_x_continuous(expand = c(0,0), position = 'top',breaks = seq(0,200,20)) +  # MIF
#  scale_y_continuous(expand = c(0,0),limits = c(-0.2,1.2),breaks = seq(0,1,0.5))+   # BUB3
#  scale_x_continuous(expand = c(0,0), position = 'top', breaks = seq(0,100,20))+ # BUB3
  labs(title = expression(bold(paste('Homogeneity (adj ',italic(R)^2, ")"))), x = '', y = '')+
  coord_flip()

#
pdf(str_c(name,'.R2.pdf'),wi=10,he=8)
grid.arrange(heatmap_grob, R2_plot, ncol = 2, widths = c(8, 3))
dev.off()
