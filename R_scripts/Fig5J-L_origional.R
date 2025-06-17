library(ggplot2)
library(ggplotify)
library(stringr)

setwd('D:/Data/0_lixia/0)_8e_HTseq_all/0_Paper_Figure_Data_source_2024/0_d10mu2/065_GO_PLOT')
#name='NC_var_K562'
#name='C_var_cell_line'
#name='CE_var_293T'
#name='NC_var_293T'

data=read.table(str_c(name,'.txt'),header = T,sep='\t')
names(data)=c('Term','Count','fold','logP')

log10_name=expression(bold(paste('-',log[10], '(Padj)')))
legend_name = 'Fold Enrichment'

pdf(str_c(name,'.GO.pdf'),wi=12,he=8)
ggplot(data, aes(x = reorder(Term, logP), y = logP, size = Count, color = fold)) +
  geom_point() +
  coord_flip() +
  scale_color_gradient(name = legend_name, low = "#6a99d0", high = "red") +
  scale_size(name = "Count",range = c(5, 15)) +
  labs(title = "", 
       x = "", 
       y = log10_name) +
  theme_bw(base_line_size = 1) +
  theme(
    panel.border = element_rect(color = "black", linewidth = 1.5),
    axis.text = element_text(size = 24,face = 'bold', color = 'black'),
    axis.title = element_text(size = 24,face = 'bold', color = 'black'),
    axis.ticks = element_line(linetype = 1, color = 'black', linewidth = 1.5),
    axis.ticks.length = unit(0.2, 'cm'),
    legend.text = element_text(size = 24,face = 'bold', color = 'black',vjust = .5),
    legend.title = element_text(size = 24,face = 'bold', color = 'black'),
    legend.key.size = unit(1, 'cm'),
    
    legend.position = "top",                 # 图例整体在顶部
    legend.justification = "right",           # 左对齐
    legend.spacing.x = unit(1.5, "cm"),      # 调整图例间距
    
  )+
  guides(
    # 通过 order 参数明确指定优先级
    size = guide_legend(order = 1),     # 先显示 Count
    color = guide_colorbar(order = 2)   # 后显示 Fold Enrichment
  )
dev.off()
