setwd('D:/Data/0_lixia/0)_8e_HTseq_all/0_Paper_Figure_Data_source_2024/0_Review_CellRepotyMethod/DNA_mutation')
library(ggplot2)
library(stringr)
name='K562_DNA.NLS_CTL.d10m2.delta_rate.no0.Fisher.p'

data=read.table(name,header=T,sep='\t')
head(data)
data$qvalue=1

d1=data[,4]
d2=data[,7]

df=as.data.frame(data$PValue)
names(df)=c('p')
nrow(df)

scientific_10 <- function(x) {   
  parse(text=gsub("e\\+*", " %*% 10^", scales::scientific_format()(x))) 
  }

p=ggplot(df, aes(x = p)) + 
  geom_histogram(bins = 50, fill = "skyblue", color = "black") +
  geom_vline(xintercept = 0.05, color = "red", linetype = "dashed") +
  labs(title = "", x = "p-value", y = "Count")+
  theme_bw(base_rect_size = 1.5) +
  scale_y_continuous(limits = c(0,600000),breaks = seq(0,600000,100000),expand = c(0.01,0.01),labels = scientific_10,) +
  theme(
    axis.title = element_text(size = 20, face = 'bold',color='black'),
    axis.text = element_text(size = 14),
    legend.text = element_text(size=18,face = 'bold',color='black'),
    legend.position = c(0.8,0.8),
    legend.title = element_blank(),
    legend.ticks = element_line(size=1.5),
    axis.ticks = element_line(size=1.5),
    axis.ticks.length = unit(0.2,'cm'),
    panel.grid = element_blank(),
  )+
  scale_x_continuous(limits = c(0,1.05),breaks = seq(0,1,0.1),expand = c(0,0))
pdf(str_c(name,'.hist.pdf'),h=8,w=8)
p
dev.off()

qf=as.data.frame(data$qvalue)
names(qf)=c('q')
nrow(qf)

p1=ggplot(qf, aes(x = q)) + 
  geom_histogram(bins = 50, fill = "skyblue", color = "black") +
  geom_vline(xintercept = 0.05, color = "red", linetype = "dashed") +
  labs(title = "", x = "p-adjust", y = "Count")+
  theme_bw(base_rect_size = 1.5) +
  scale_y_continuous(limits = c(0,900000),breaks = seq(0,900000,100000),expand = c(0.01,0.01),labels = scientific_10,) +
  theme(
    axis.title = element_text(size = 24, face = 'bold',color='black'),
    axis.text = element_text(size = 16),
    legend.text = element_text(size=18,face = 'bold',color='black'),
    legend.position = c(0.8,0.8),
    legend.title = element_blank(),
    legend.ticks = element_line(size=1.5),
    axis.ticks = element_line(size=1.5),
    axis.ticks.length = unit(0.2,'cm'),
    panel.grid = element_blank(),
  )+
  scale_x_continuous(limits = c(0,1.05),breaks = seq(0,1,0.1),expand = c(0,0))

pdf(str_c(name,'.padj.hist.pdf'),h=5,w=8)
p1
dev.off()
