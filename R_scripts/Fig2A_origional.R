library(pROC)
library(ggplot2)

files <- list('ROC1_rRNA_Probe1.txt','ROC2_rRNA.txt','ROC4_Probe1.txt')     ### 可多个文件一起
file_name = list('Structure-known RNAS + RNA probes','Structure-known RNAs','RNA probes')    ### 指定文件对应name

true_ids <- read.table('true.txt', header = F)
false_ids <- read.table('false.txt', header = F)
files <- list('ROC7_rRNA2.txt','ROC7_rRNA3.txt','ROC7_rRNA4.txt','ROC7_rRNA5.txt','ROC7_rRNA6.txt')     ### 可多个文件一起
file_name = list('ROC7_rRNA2','ROC7_rRNA3','ROC7_rRNA4','ROC7_rRNA5','ROC7_rRNA6')    ### 指定文件对应name
true_ids <- read.table('roc7_T.txt', header = F)
false_ids <- read.table('roc7_F.txt', header = F)

#
files <- list('ROC7_rRNA6.txt')     ### 可多个文件一起
file_name = list('ROC7_rRNA6')    ### 指定文件对应name
true_ids <- read.table('roc7_rna6_T.txt', header = F)
false_ids <- read.table('roc7_rna6_F.txt', header = F)
#
files <- list('ROC7_rRNA9.txt','ROC7_rRNA10.txt')     ### 可多个文件一起
file_name = list('ROC7_rRNA9','ROC7_rRNA10')    ### 指定文件对应name
true_ids <- read.table('roc7_rna9_T.txt', header = F)
false_ids <- read.table('roc7_rna9_F.txt', header = F)
#
outname='roc7_rna9_ROC'   # outname.roc.pdf
pdf_h=8     # 高
pdf_w=8    #  宽
###################################
names(true_ids)=c('True')
names(false_ids)=c('False')
all_roc_data <- data.frame()
all_auc <- data.frame()
#
for (i in seq_along(files)) {
  file_data <- read.table(files[[i]],header=T,sep='\t') ### 读取当前文件
  ids <- file_data[, 1]
  values <- file_data[, 2]
  labels <- ifelse(ids %in% true_ids$True, 1,            ### 生成实际的标签 
                   ifelse(ids %in% false_ids$False, 0, NA)) ### 如果ID在True list中，标签为1；如果在False list中，标签为0
  roc_data <- data.frame(labels, values)
  roc_data <- na.omit(roc_data)            ### 确保没有NA的标签值
  if (length(unique(roc_data$labels)) == 2) { ### 确保标签有两个不同的值
    roc_curve <- roc(roc_data$labels, roc_data$values) ### 计算ROC曲线
    auc = auc(roc_curve) ### AUC 计算
    ###### 提取ROC曲线的数据点并添加到总数据框中
    roc_data_points <- data.frame(
      specificity = rev(roc_curve$specificities),
      sensitivity = rev(roc_curve$sensitivities),
      model = rep(file_name[[i]],times=length(roc_curve$specificities))  # 添加模型标签
    )
    all_auc = c(all_auc,auc)
    all_roc_data <- rbind(all_roc_data, roc_data_points) ###  收集 AUC数值
  } else {
    stop(paste('File', files[[i]], 'lables only have 1 level (only True or False)，Please check the data'))
  }
}
all_auc=as.numeric(all_auc)   
# 颜色
custom_colors <- c('red', 'orange', 'blue','green', 'purple', 'brown', 'pink', 'cyan', 'yellow', 'black')
#########  设置AUC数值标签
model_names <- unique(file_name)
legend_labels <- paste0(' (AUC = ', round(all_auc, 2), ')')  
all_roc_data$model <- factor(all_roc_data$model,levels = unique(all_roc_data$model), labels = paste0(model_names,legend_labels))

########## 绘制ROC曲线
p=ggplot(all_roc_data, aes(x = 1 - specificity, y = sensitivity, color = model)) +
  geom_line(size = 1.2) +
  #  geom_smooth(size = 1.2,se = F,)+
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'gray50',linewidth=1) +
  ggtitle('') +
  scale_x_continuous(name = '1 - Specificity',expand = c(0,0),limits = c(-0.01,1.01))+
  scale_y_continuous(name = 'Sensitivity',expand = c(0,0),limits = c(-0.01,1.01))+
  scale_color_manual(values = custom_colors) +  # 使用自定义颜色
  theme_bw(base_rect_size = 2) +
  guides(color = guide_legend(override.aes = list(linewidth = 2))) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 24, face = 'bold', color = 'black'),
    axis.text = element_text(size = 20, face = 'bold', color = 'black'),
    axis.line = element_line(linetype = 1, color = 'black', linewidth = 1),
    axis.ticks = element_line(linetype = 1, color = 'black', linewidth = 1.2),
    axis.ticks.length = unit(0.25, 'cm'),
    legend.title = element_blank(),
    legend.position = c(0.58, 0.2),       ### 调整legend位置
    legend.text = element_text(size = 14, face = 'bold', color = 'black'),
    legend.spacing.y = unit(1,'cm'),         ### 调整legend上下距离
    legend.key.height = unit(.8, 'cm'),     ### 调整legend上下高度
    plot.margin = margin(0,1,0.2,0.2, unit = "cm"), # t = 2, r = 2, b = 2, l = 2, unit = "pt" t、r、b、l分别表示上、右、下、左四侧边距；unit为间距单位，可以使用pt、cm、in等。
  )

pdf(paste0(outname,'.roc.pdf'),wi=pdf_w,he=pdf_h)
p
dev.off()
