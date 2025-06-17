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
col_vector
##################
setwd('D:/Data/0_lixia/0)_8e_HTseq_all/0_Paper_Figure_Data_source_2024/0_single_cell/gene_number')
pdf("SC_gene.number.box.pdf",w=6,h=6)
data1=read.table(file="SC_gene.number.NES.txt",sep='\t',header=T)
col=col_vector[1:11]
data1=data1[,c(2,3)]
names(data1)=c('Group','value')
data1$value=data1$value/1000

max=round(max(data1$value)+1)
min=round(min(data1$value)-1)
#####################################
order=as.character(as.data.frame(table(data1$Group))$Var1)
p <- ggboxplot(data1, x = 'Group', y = 'value',
               color = 'Group', 
               #palette = col,
               palette=colorRampPalette(c("#42808E", "#42808E"))(11), 
#               palette = 'rickandmorty', 
#ggsci scientific journal palettes, e.g.: 'npg', 'aaas', 'lancet', 'jco', 'ucscgb', 'uchicago', 'simpsons' and 'rickandmorty'.
               size=1.5,
#             add = 'median_iqr',
              outlier.shape = NA,
              bxp.errorbar = T,
              bxp.linewidth=1.5,
#             bxp.errorbar.width=0.25,
#             label.select = list(top.up = 10, top.down = 4),
              xlab = 'Read number (Million)',ylab='No. of genes (X1000)'
)+
  scale_x_discrete(limits = order,
                   labels = c("1", "2", "3", "4", "5", "7.5", "10", "12.5", "15", "17.5", "20" ))+
  theme(
    axis.title = element_text(size=18,face = 'bold',color='black'),
    axis.text = element_text(size=16,face = 'bold',color='black'),
    legend.position = 'none',
  )+
  scale_y_continuous(limits = c(min,14),breaks = seq(min,14,1))+
  annotate('text', size = 4,
           x = 1:length(table(data1$Group)),
           y = 5.5,
           label = paste('n=',table(data1$Group)[order]),
           col = 'black',
           vjust = 0)
p

dev.off()
