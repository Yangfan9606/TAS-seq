library(ComplexUpset)
library(ggplot2)

setwd('D:/Data/0_lixia/0)_8e_HTseq_all/0_Paper_Figure_Data_source_2024/0_d10mu2/063_gene_overlap_PLOT')
d1=read.table('NES.NLS.variation.sites',header = F,sep = '\t')
d2=read.table('NES.ERM.variation.sites',header = F,sep = '\t')
d3=read.table('K562_NES.K562_NLS.variation.sites',header = F,sep = '\t')
d4=read.table('NLS.K562_NLS.variation.sites',header = F,sep = '\t')
d5=read.table('NES.K562_NES.variation.sites',header = F,sep = '\t')

d1n='Nuc/Cyt Var (HEK293T)'
d2n='Cyt/Erm Var (HEK293T)'
d3n='Nuc/Cyt Var (K562)'
d4n='Nuc_Cell line Var'
d5n='Cyt_Cell line Var'

df <- list(
  'Nuc/Cyt Var (HEK293T)' = d1$V1,
  'Cyt/Erm Var (HEK293T)' = d2$V1,
  'Nuc/Cyt Var (K562)' = d3$V1,
  'Nuc_Cell line Var' = d4$V1,
  'Cyt_Cell line Var' = d5$V1
)

all_elements <- unique(unlist(df))
# 创建一个空的数据框
upset_df <- data.frame(
  element = all_elements,
  stringsAsFactors = FALSE
)
# 对每个元素检查是否存在于每个类别中
for (name in names(df)) {
  upset_df[[name]] <- upset_df$element %in% df[[name]]
}
# 打印结果
print(upset_df)

pdf('Sites_overlap.colored.pdf',wi=12,he=6)

###########################
upset(
  upset_df,
  names(df),
  name='',
  width_ratio=0.3,
  height_ratio = 0.7,
  #########################     
  base_annotations=list( #overlap number plot
    'Intersection size'=intersection_size(
      counts=F,
      mapping=aes(fill='bars_color')
    ) + theme_classic(base_line_size = 1)+
      theme(
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=12,face = 'bold',color = 'black'),
        axis.title.y = element_text(size=16,face = 'bold',color = 'black'),
      )+
 #     scale_fill_manual(values=c('bars_color'='#999999'), guide='none')+
      scale_fill_manual(values=c('bars_color'='#1F78B4'), guide='none')+
      scale_y_continuous(name='Variation sites\noverlap number',expand = c(0,0),breaks = seq(0,80000,20000),limits = c(0,80000))
  ), 
  ##########################    
  set_sizes=( #bar plot
    upset_set_size(position='right')
    + theme_classic(base_line_size = 1)
    + theme(
      axis.ticks.y = element_blank(), 
      axis.line.y =  element_blank(),
      axis.title.y = element_blank(),
      axis.text.y=element_blank(),
      axis.text.x=element_text(size=12,face = 'bold',color = 'black'),
      axis.title.x = element_text(size=14,face = 'bold',color = 'black'),
      plot.margin = margin(0,1,0,0, unit = "cm"), # t = 2, r = 2, b = 2, l = 2, unit = "pt" t、r、b、l分别表示上、右、下、左四侧边距；unit为间距单位，可以使用pt、cm、in等。
    )
    +scale_y_continuous(name='Variation sites number (x1000)',expand = c(0,0),labels = function(x) x / 1000,breaks = seq(0,150000,50000),limits = c(0,150000))
    #   +scale_y_reverse(name='Gene number',expand = c(0,0))  ### 对应 upset_set_size(position='right')
  ),
  #########################               
  themes=upset_modify_themes( #整个图的 theme
    list(
      'intersections_matrix'=
        theme_bw()+theme(
          axis.title.y = element_text(size=16,face = 'bold',color = 'black'),
          axis.title.x = element_text(size=16,face = 'bold',color = 'black'),
          axis.text.y=element_text(size=12,face = 'bold',color = 'black'),
          axis.text.x=element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks = element_line(size=1.5),
          axis.ticks.length = unit(0.2,'cm'),
          panel.border = element_rect(color = "black", linewidth = 1),
        )
    )
  ),
  ########################### 
  matrix=(#  调节点的 类型 大小，链接线段的长度
    intersection_matrix(
      geom=geom_point(shape='circle filled', size=5),
      segment=geom_segment(linetype='solid',linewidth = 1.5),
      ###      outline_color=list(active='black',inactive='grey70')
    )
    +scale_color_manual(
      values=c('TRUE'='#e15759', 'FALSE'='grey'),
      na.value='transparent',
      guide='none'  
    )
  ),
  # 'a'='#7FC97F','b'='#1F78B4', 'c'='#F0027F','d'='#6A3D9A','e'='#BF5B17'
  #  "#A6CEE3","#E78AC3", "#8DD3C7", "#BEAED4"  "#6A3D9A" "#BF5B17" "#FDC086"
  #c("red", "green", "blue", "purple", "orange", "brown", "pink", "cyan", "yellow", "black")
  ##############################
  queries=list(  ## 填充 bar 和点 颜色
    upset_query(set=d1n, fill='#D79C9C'),
    upset_query(set=d2n, fill='#BF5B17'),
    upset_query(set=d3n, fill='#FDC086'),
    upset_query(set=d4n, fill='#7FC97F'),
    upset_query(set=d5n, fill='#BEAED4')
)
dev.off()
