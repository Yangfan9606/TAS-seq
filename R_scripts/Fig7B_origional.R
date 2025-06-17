library(ggplot2)
library(dplyr)
library(ggrepel)
setwd('D:/Data/0_lixia/0)_8e_HTseq_all/0_Paper_Figure_Data_source_2024/0_single_cell/Final_PLOT_of_cosine_and_R2-HEAT_map')

# Load the data
all_genes <- read.table('cosine_distance.txt',header=T,sep='\t')
selected_genes = read.table('selected_gene.txt',header=F,sep='\t')
names(all_genes)=c('gene_name','cosine_distance')
names(selected_genes)=c('gene_name')

# Calculate percentiles for the red dashed lines
percentile_25 <- quantile(all_genes$cosine_distance, 0.25)
percentile_75 <- quantile(all_genes$cosine_distance, 0.75)
percentile_25
percentile_75
# Merge the data to get only selected genes
selected_data <- all_genes %>% 
  filter(gene_name %in% selected_genes$gene_name)

# Create the density plot and calculate the density
density_data <- density(all_genes$cosine_distance)
density_df <- data.frame(x = density_data$x, y = density_data$y)

# Get the y-values for the selected genes from the density
selected_data <- selected_data %>%
  rowwise() %>%
  mutate(y = approx(density_df$x, density_df$y, xout = cosine_distance)$y)
# Create the density plot
p=ggplot() +
  geom_density(data = all_genes, aes(x = cosine_distance), fill = "grey", alpha = 0.5) +
  labs(x = "Deviation of cosine distance", y = "Density") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 24, face = 'bold', color = 'black'),
    axis.text = element_text(size = 20, face = 'bold', color = 'black'),
    axis.line = element_line(linetype = 1, color = 'black', linewidth = 1),
    axis.ticks = element_line(linetype = 1, color = 'black', linewidth = 1),
    axis.ticks.length = unit(0.25, 'cm'),
    #    legend.text = element_text(size = 20, face = 'bold', color = 'black'),
    #    legend.title = element_text(size = 24, color = 'black')
  )+
  annotate("text", x = percentile_25/4, y = max(density_df$y) + 1, label = "Homogeneous\n    transcripts", size=5, hjust = 0) +
  geom_segment(aes(x=0,y=max(density_df$y)+max(density_df$y)/10,xend=percentile_25,yend=max(density_df$y)+max(density_df$y)/10), linetype = 1,size=1.5,color='#66A61E') +
  annotate("text", x = percentile_75+percentile_25/4, y = max(density_df$y) + 1, label = "Heterogeneous\n    transcripts", size=5, hjust = 0)+
  geom_segment(aes(x=percentile_75,y=max(density_df$y)+max(density_df$y)/10,xend=percentile_75+percentile_25,yend=max(density_df$y)+max(density_df$y)/10), linetype = 1,size=1.5,color='#E41A1C') +
  scale_x_continuous(limits = c(0,0.6),breaks = seq(0,0.6,0.1))+
  scale_y_continuous(limits =c(0,max(density_df$y) + 2),expand = c(0,0),breaks = seq(0,max(density_df$y)+1,1))+
  geom_vline(xintercept = percentile_25, color = "#6a99d0", linetype = "dashed",linewidth=1.1) +
  geom_vline(xintercept = percentile_75, color = "#6a99d0", linetype = "dashed",linewidth=1.1) +
  geom_text_repel(data = selected_data, aes(x = cosine_distance, y = y, label = gene_name), nudge_y = 1, segment.size = 0.5, point.padding = 0.5, min.segment.length = 0,size = 5) 
  ##  geom_text_repel(data = subset(selected_data, gene_name == "18S rRNA"), aes(x = cosine_distance, y = y, label = gene_name), color = "red", nudge_y = 1, segment.size = 0.5, point.padding = 0.5, min.segment.length = 0) +
p
pdf('cosine_distance.pdf',he=6,wi=8)
p
dev.off()
