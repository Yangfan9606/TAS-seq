setwd('D:/Data/0_lixia/0)_8e_HTseq_all/0_Paper_Figure_Data_source_2024/0_Review_CellRepotyMethod/splicing/plot')
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)

# Read the data
data <- read.table('splicing_5SS_1-90.txt', header = T, sep = '\t')

# Set color palette for fc1 and fc2
cols <- c("Spliced" = "#1F78B4", "Unspliced" = "#FF7F00")

# Reshape the data into long format
data_long <- reshape(data, 
                     varying = c('Spliced','Unspliced'), 
                     v.names = 'value', 
                     timevar = 'type', 
                     times = c('Spliced','Unspliced'), 
                     direction = 'long')

data_long$Position=data_long$Position+1
data_long$value=as.numeric(data_long$value)
# Create the plot
p=ggplot(data = data_long, aes(x = Position, y = value,group = type,colour = type)) +
  geom_line(linewidth=1.2)+
  geom_point(data = data_long, aes(x = Position, y = value, group = type, color = type), size = 0.5, stroke = 1.8)+
  scale_color_manual(values = cols)+ # Set colors 
  theme_classic(base_rect_size = 1.5) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 1, colour = 'black'),
    axis.title = element_text(size = 24, face = 'bold', color = 'black'),
    axis.text.x = element_text(size = 18, face = 'bold', color = 'black'),
    axis.text.y = element_text(size = 18, face = 'bold', color = 'black'),
    axis.ticks = element_line(linetype = 1, color = 'black', linewidth = 1.2),
    axis.ticks.length = unit(0.25, 'cm'),
    legend.title = element_blank(),
    legend.text = element_text(size = 24, face = 'bold', color = 'black'),
    legend.key.width = unit(2, 'cm'),
    legend.position = c(0.8,0.9)
  )+ 
  labs(x = "Relative Position (nt)", y = "Modification rate") +  # Labels for axes
  scale_x_continuous(limits = c(-90,0), breaks = seq(-90, 0, by = 10),expand = c(0.01,0.5))+
  scale_y_continuous(limits = c(0,0.031), breaks = seq(0, 0.05, 0.005),expand = c(0,0))

pdf('splicing_5SS_1-90.pdf',wi=12,he=8)
p
dev.off()
