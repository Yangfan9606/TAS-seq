setwd('D:/Data/0_lixia/0)_8e_HTseq_all/0_Paper_Figure_Data_source_2024/000_Manuscript/task/0_Download_DATA/1_longest.transcript.A50.A25_rate.gini.gene.avg')
   
library(GGally)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(ggpubr)
name='all.A.sites_n50.A25.gini.gene'
file <- read.table(str_c(name),header = T,sep = '\t')

sample_df=file[,2:13]

#####   var_names <- paste0("`", colnames(sample_df), "`")
names(sample_df) = c('DMS-seq (K562, Who)','icSHAPE (HEK293T, Who, 2016)','DMS-Mapseq (HEK293T, Who)','icSHAPE (HEK293, Cyt)','icSHAPE (HEK293, Nuc)','icSHAPE (HEK293, Who, 2021)','icSHAPE (K562, Who)','smartSHAPE (HEK293T, Who)','TAS-seq (HEK293T, Cyt)','TAS-seq (HEK293T, Nuc)','TAS-seq (K562, Cyt)','TAS-seq (K562, Nuc)')

TAS_df=sample_df$`TAS-seq (HEK293T, Cyt)`
ic_df=sample_df$`icSHAPE (HEK293, Cyt)`
DMS_df=sample_df$`DMS-Mapseq (HEK293T, Who)`
TAS_K562_Cyt=sample_df$`TAS-seq (K562, Cyt)`
DMS_K562_Wh=sample_df$`DMS-seq (K562, Who)`
ic_Wh=sample_df$`icSHAPE (HEK293T, Who, 2016)`

########################################################################
cor_data1=data.frame(TAS_df,ic_df)
colnames(cor_data1) <- c('T','ic')
nrow(cor_data1[!(is.na(cor_data1$T)) & !(is.na(cor_data1$ic)), ])
pdf("TAS_NES.icSHAPE.cor.pdf",w=6,h=6)
ggscatter(cor_data1, x = 'ic', y = 'T',
          color = '#277fb8', size = 2, # Points color, shape and size
          alpha = 0.2,
          # shape = 21, 
          add = 'reg.line',  # Add regressin line
          add.params = list(color = 'darkblue', fill = 'lightgray', size = 1.5), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = 'pearson', label.x = 0.05, label.sep = '\n',size = 6, face = 'bold', color='black')
)+
  theme_bw(base_rect_size = 1.5)+
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 20, face = 'bold',color='black'),
    axis.text = element_text(size = 16, face = 'bold',color='black'),
    axis.ticks = element_line(linewidth = 1.5),
    axis.ticks.length = unit(0.2,'cm'),
    plot.margin = margin(0.2,0.5,0.2,0.2, unit = 'cm')
  )+
  scale_x_continuous(name = str_c('icSHAPE (HEK293, Cyt)','\n','Gini index (A)'),limits=c(0,1),breaks=seq(0,1,by=0.1),expand = c(0,0.01))+
  scale_y_continuous(name = str_c('TAS-seq (HEK293T, Cyt)','\n','Gini index (A)'),limits=c(0,1),breaks=seq(0,1,by=0.1),expand = c(0.01,0))
dev.off()

#############################################################
cor_data1=data.frame(TAS_df,DMS_df)
colnames(cor_data1) <- c('T','dm')
nrow(cor_data1[!(is.na(cor_data1$T)) & !(is.na(cor_data1$dm)), ])
length(cor_data1$T[is.na(cor_data1$T)==F])
length(cor_data1$dm[is.na(cor_data1$dm)==F])
nrow(cor_data1)
pdf("TAS_NES.DMSmap.cor.pdf",w=6,h=6)
ggscatter(cor_data1, x = 'dm', y = 'T',
          color = '#8BC9B4', size = 2, # Points color, shape and size
          alpha = 0.2,
          # shape = 21, 
          add = 'reg.line',  # Add regressin line
          add.params = list(color = 'darkblue', fill = 'lightgray', size = 1.5), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = 'pearson', label.x = 0.05, label.sep = '\n',size = 6, face = 'bold', color='black')
)+
  theme_bw(base_rect_size = 1.5)+
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 20, face = 'bold',color='black'),
    axis.text = element_text(size = 16, face = 'bold',color='black'),
    axis.ticks = element_line(linewidth = 1.5),
    axis.ticks.length = unit(0.2,'cm'),
    plot.margin = margin(0.2,0.5,0.2,0.2, unit = 'cm')
  )+
  scale_x_continuous(name = str_c('DMS-Mapseq (HEK293T, Who)','\n','Gini index (A)'),limits=c(0,1),breaks=seq(0,1,by=0.1),expand = c(0,0.01))+
  scale_y_continuous(name = str_c('TAS-seq (HEK293T, Cyt)','\n','Gini index (A)'),limits=c(0,1),breaks=seq(0,1,by=0.1),expand = c(0.01,0))
dev.off()

#############################################################
cor_data1=data.frame(TAS_K562_Cyt,DMS_K562_Wh)
colnames(cor_data1) <- c('T','dm')
nrow(cor_data1[!(is.na(cor_data1$T)) & !(is.na(cor_data1$dm)), ])

pdf("TAS_K562_NES.DMS_K562_Who.cor.pdf",w=6,h=6)
ggscatter(cor_data1, x = 'dm', y = 'T',
          color = '#FECC4F', size = 2, # Points color, shape and size
          alpha = 0.2,
          # shape = 21, 
          add = 'reg.line',  # Add regressin line
          add.params = list(color = 'darkblue', fill = 'lightgray', size = 1.5), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = 'pearson', label.x = 0.05, label.sep = '\n',size = 6, face = 'bold', color='black')
)+
  theme_bw(base_rect_size = 1.5)+
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 20, face = 'bold',color='black'),
    axis.text = element_text(size = 16, face = 'bold',color='black'),
    axis.ticks = element_line(linewidth = 1.5),
    axis.ticks.length = unit(0.2,'cm'),
    plot.margin = margin(0.2,0.5,0.2,0.2, unit = 'cm')
  )+
  scale_x_continuous(name = str_c('DMS-seq (K562, Who)','\n','Gini index (A)'),limits=c(0,1),breaks=seq(0,1,by=0.1),expand = c(0,0.01))+
  scale_y_continuous(name = str_c('TAS-seq (K562, Cyt)','\n','Gini index (A)'),limits=c(0,1),breaks=seq(0,1,by=0.1),expand = c(0.01,0))
dev.off()

#############################################################
cor_data1=data.frame(ic_df,DMS_df)
colnames(cor_data1) <- c('T','dm')
nrow(cor_data1[!(is.na(cor_data1$T)) & !(is.na(cor_data1$dm)), ])

pdf("icSHAPE_cy.DMSmap_HEK293T_Who.cor.pdf",w=6,h=6)
ggscatter(cor_data1, x = 'T', y = 'dm',
          color = '#FECC4F', size = 2, # Points color, shape and size
          alpha = 0.2,
          # shape = 21, 
          add = 'reg.line',  # Add regressin line
          add.params = list(color = 'darkblue', fill = 'lightgray', size = 1.5), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = 'pearson', label.x = 0.05, label.sep = '\n',size = 6, face = 'bold', color='black')
)+
  theme_bw(base_rect_size = 1.5)+
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 20, face = 'bold',color='black'),
    axis.text = element_text(size = 16, face = 'bold',color='black'),
    axis.ticks = element_line(linewidth = 1.5),
    axis.ticks.length = unit(0.2,'cm'),
    plot.margin = margin(0.2,0.5,0.2,0.2, unit = 'cm')
  )+
  scale_x_continuous(name = 'icSHAPE (HEK293, Cyt)',limits=c(0,1),breaks=seq(0,1,by=0.1),expand = c(0,0.01))+
  scale_y_continuous(name = 'DMS-Mapseq (HEK293T, Who)',limits=c(0,1),breaks=seq(0,1,by=0.1),expand = c(0.01,0))
dev.off()

########################################################################
cor_data1=data.frame(ic_df,DMS_df)
colnames(cor_data1) <- c('T','ic')
nrow(cor_data1[!(is.na(cor_data1$T)) & !(is.na(cor_data1$ic)), ])
pdf("icSHAPE_DMSmap.cor.pdf",w=6,h=6)
ggscatter(cor_data1, x = 'T', y = 'ic',
          color = '#277fb8', size = 2, # Points color, shape and size
          alpha = 0.2,
          # shape = 21, 
          add = 'reg.line',  # Add regressin line
          add.params = list(color = 'darkblue', fill = 'lightgray', size = 1.5), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = 'pearson', label.x = 0.05, label.sep = '\n',size = 6, face = 'bold', color='black')
)+
  theme_bw(base_rect_size = 1.5)+
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 20, face = 'bold',color='black'),
    axis.text = element_text(size = 16, face = 'bold',color='black'),
    axis.ticks = element_line(linewidth = 1.5),
    axis.ticks.length = unit(0.2,'cm'),
    plot.margin = margin(0.2,0.5,0.2,0.2, unit = 'cm')
  )+
  scale_x_continuous(name = 'icSHAPE (HEK293, Cyt)',limits=c(0,1),breaks=seq(0,1,by=0.1),expand = c(0,0.01))+
  scale_y_continuous(name = 'DMS-Mapseq (HEK293T, Who)',limits=c(0,1),breaks=seq(0,1,by=0.1),expand = c(0.01,0))
dev.off()
##############

cor_data1=data.frame(ic_Wh,DMS_df)
colnames(cor_data1) <- c('T','ic')
nrow(cor_data1[!(is.na(cor_data1$T)) & !(is.na(cor_data1$ic)), ])
pdf("icSHAPE_Who_DMSmap.cor.pdf",w=6,h=6)
ggscatter(cor_data1, x = 'T', y = 'ic',
          color = '#fb67ae', size = 2, # Points color, shape and size
          alpha = 0.2,
          # shape = 21, 
          add = 'reg.line',  # Add regressin line
          add.params = list(color = 'darkblue', fill = 'lightgray', size = 1.5), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = 'pearson', label.x = 0.05, label.sep = '\n',size = 6, face = 'bold', color='black')
)+
  theme_bw(base_rect_size = 1.5)+
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 20, face = 'bold',color='black'),
    axis.text = element_text(size = 16, face = 'bold',color='black'),
    axis.ticks = element_line(linewidth = 1.5),
    axis.ticks.length = unit(0.2,'cm'),
    plot.margin = margin(0.2,0.5,0.2,0.2, unit = 'cm')
  )+
  scale_x_continuous(name = str_c('icSHAPE (HEK293T, Who, 2016)','\n','Gini index (A)'),limits=c(0,1),breaks=seq(0,1,by=0.1),expand = c(0,0.01))+
  scale_y_continuous(name = str_c('DMS-Mapseq (HEK293T, Who)','\n','Gini index (A)'),limits=c(0,1),breaks=seq(0,1,by=0.1),expand = c(0.01,0))
dev.off()
