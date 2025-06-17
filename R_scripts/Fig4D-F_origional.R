setwd('D:/Data/0_lixia/0)_8e_HTseq_all/0_Paper_Figure_Data_source_2024/0_Review_CellRepotyMethod/ERM_spatial_specificity')

library(pROC)
library(ggplot2)
library(RColorBrewer)

# 定义不同 True False list
true_symbol='True_list.ERM_640.2'
false_symbol='p_s_t.N.overlap.no640.no_uniprot_Secreted'

true_symbol='True_list.ERM_640'
true_symbol='True_list.ERM_640_a30'

true_ensg='True_list.ERM_640.ENSG'
false_ensg='False_list.no_ERM_640.new.2.ENSG'
datasets=list()
# 定义多个数据集的信息（可扩展添加）
datasets <- list(
#  list(
#    name = "loRNA",            # 数据集2名称（用于图例）
#    true_list = true_symbol,
#    false_list = false_symbol,
#    file = "loRNA.txt",
#    gene_col = 2,
#    rate_col = 4
#  ),
  # 可继续添加更多数据集...
  list(
    name = "ERM enriched RNA\n(Erm/Cyt)\n",            # 数据集2名称（用于图例）
    true_list = true_symbol,
    false_list = false_symbol,
    file = "ERM_vs_NES.deseq2.padj005",
    gene_col = 1,
    rate_col = 2
  )
#  list(
#    name = "APEX_ERM",            # 数据集1名称（用于图例）
#    true_list = true_ensg,
#    false_list = false_ensg,
#    file = "APEX_ERM.txt",
#    gene_col = 1,                   # 基因名称所在列
#    rate_col = 5                    # 突变数值所在列
#  )
)

#########  注意颜色数量要和线条数量对应
cols <- brewer.pal(9, "Set1")[1:length(datasets)]
#cols = brewer.pal(10,"Paired")[1:10]
#pal <- colorRampPalette(cols)
#colors = pal(1)
#qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors,rownames(qual_col_pals)))
#cols=col_vector[9:13]

specificity_points=list()

# ================== 数据读取与处理 ==================
roc_list <- list()  # 存储所有ROC对象

for (dataset in datasets) {
  # 读取当前数据集的基因列表和突变数据
  true_genes <- read.table(dataset$true_list, header = FALSE, stringsAsFactors = FALSE)[, 1]
  false_genes <- read.table(dataset$false_list, header = FALSE, stringsAsFactors = FALSE)[, 1]
  
  mutation_data <- read.table(dataset$file, header = F, sep = "\t", stringsAsFactors = FALSE)
  gene_names <- mutation_data[, dataset$gene_col]
  pred_scores <- mutation_data[, dataset$rate_col]
  
  # 筛选属于正类或负类的基因
  keep_rows <- gene_names %in% c(true_genes, false_genes)
  gene_names_filtered <- gene_names[keep_rows]
  pred_scores_filtered <- pred_scores[keep_rows]
  
  # 生成真实标签
  true_labels <- ifelse(gene_names_filtered %in% true_genes, 1, 0)
  
  # 计算 ROC 和 Youden's Index
  roc_obj <- roc(response = true_labels, predictor = pred_scores_filtered)
  best_threshold <- coords(roc_obj, x = "best", best.method = "youden", ret = c("specificity","sensitivity",'threshold'))
  
  # 保存 ROC 对象和最佳点
  roc_list[[dataset$name]] <- list(
    roc = roc_obj,
    auc = auc(roc_obj),
    specificity = best_threshold$specificity,
    sensitivity = best_threshold$sensitivity,
    name = dataset$name
  )
  
  # 收集最佳点数据
  specificity_points <- rbind(
    specificity_points,
    data.frame(
      FPR = 1 - best_threshold$specificity,
      TPR = best_threshold$sensitivity,
      Specificity = round(best_threshold$specificity, 2),
      Dataset = dataset$name,
      Threshold=best_threshold$threshold,
      Sensitivity=round(best_threshold$sensitivity, 2)
    )
  )
  
}

# ================== 合并数据并绘图 ==================
# 生成统一的绘图数据框
roc_data <- data.frame()
for (dataset_name in names(roc_list)) {
  obj <- roc_list[[dataset_name]]
  temp_df <- data.frame(
    FPR = 1 - obj$roc$specificities,
    TPR = obj$roc$sensitivities,
    Dataset = paste0(obj$name, "AUC = ", round(obj$auc, 2), "")
  )
  roc_data <- rbind(roc_data, temp_df)
}

# ================== 指定阈值竖线颜色 ==================
# ROC 画图顺序 默认按照 字母排序，注意可能出错
best_cutoff=data.frame(name=sort(unique(roc_data$Dataset)),val=cols)[order(specificity_points$Dataset),]
best_cutoff$FPR=specificity_points$FPR

# ================== 指定cutoff显示文字内容和颜色 ==================
specificity_points$label <- sprintf(
  "Log[2]~(Erm/Cyt) == %.2f", 
  round(specificity_points$Threshold, 2)
)

# 绘制多曲线ROC图
p=ggplot(roc_data, aes(x = FPR, y = TPR, color = Dataset)) +
  geom_line(size = 1.2) +
  geom_vline(data=best_cutoff,aes(xintercept = FPR),colour = best_cutoff$val, linetype='dashed',linewidth = 1)+

  geom_text(
    data = specificity_points, aes(x = FPR, y = TPR, label = label), size = 6, hjust = -0.02, vjust = -7.5,
    parse = TRUE,  # 关键：解析字符串为公式
    fontface='bold',
    vjust = -0.5,  # 垂直调整位置
    color = best_cutoff$val # 标签颜色
  )+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40",linewidth=1) +
  labs(
    title = "",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)",
    color = "Dataset"
  ) +
#  scale_color_brewer(palette = cols) +  # 使用鲜艳的颜色区分曲线
  scale_color_manual(values = cols)+
  theme_bw(base_rect_size = 2,base_size = 12) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 20, face = 'bold',color='black'),
    axis.text = element_text(size = 18, face = 'bold', color = 'black'),
    axis.ticks = element_line(linewidth=1.5, color = 'black'),
    axis.ticks.length = unit(0.2,'cm'),
    legend.position = c(0.65, 0.2),       # 调整图例位置
    legend.text = element_text(size=18,face = 'bold',color='black',hjust=0.5),
    legend.title = element_blank(),
    legend.key = element_blank(),
#    legend.background = element_rect(fill = "white", color = "white")
legend.background = element_blank(),
plot.margin = margin(0,0.6,0.2,0.2, unit = "cm"), # t = 2, r = 2, b = 2, l = 2, unit = "pt" t、r、b、l分别表示上、右、下、左四侧边距；unit为间距单位，可以使用pt、cm、in等。
  )+
  scale_x_continuous(expand = c(0.001,0.001))+
  scale_y_continuous(expand = c(0.001,0.001))+
  guides(
    color = guide_legend(
      override.aes = list(
        linewidth = 2,   # 增大图例中的线条宽度
#        color = cols  # 修改图例中的线条颜色
color=NA
      )
    )
  )

specificity_points
pdf('ERM_NES.ROC.pdf',he=6,wi=6)
p
dev.off()
