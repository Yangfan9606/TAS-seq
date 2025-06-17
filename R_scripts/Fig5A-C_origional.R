setwd('D:/Data/0_lixia/0)_8e_HTseq_all/0_Paper_Figure_Data_source_2024/0_d10mu2/078_half_life_RMR_density_PLOT/20241210_select')
library(ggpubr)
library(stringr)
library(plyr)

data1 <- read.table('K562_4sU_half_life.K562_NES.out', header = FALSE, sep = '\t')

d1 <- data1[, 3]
d2 <- data1[, 2]
xname='Modification density\nin TAS-seq (K562, Cyt)'
yname=expression(bold("RNA half-life in K562"~log[10]~"(hours)"))

df <- data.frame(Rep1 = d1, Rep2 = d2)
n_points <- nrow(df)
cor_res <- cor.test(df$Rep1, df$Rep2, method = "pearson")  # Pearson correlation test
cor_value <- cor_res$estimate  # Pearson correlation coefficient
p_value <- cor_res$p.value  # P-value from the test
cor_value

pdf('K562_4sU_half_life.K562_NES.cor.pdf', width = 6, height = 6)
ggscatter(df, x = 'Rep1', y = 'Rep2',
          color = '#8ECBB7', size = 2, alpha = 0.2,
          add = 'reg.line', 
          add.params = list(color = 'darkblue', fill = 'lightgray', size = 1.5),
#          conf.int = TRUE, 
#          cor.coef = TRUE,
          cor.coeff.args = list(method = 'pearson', label.x = max(d1)/10, label.sep = '\n', size = 6, face = 'bold', color = 'black')
          ) +
  theme_bw(base_rect_size = 2) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 20, face = 'bold', color = 'black'),
    axis.text = element_text(size = 16, face = 'bold', color = 'black'),
    axis.ticks = element_line(linewidth = 1.5),
    axis.ticks.length = unit(0.2, 'cm'),
    plot.margin = margin(0.2, 0.5, 0.2, 0.2, unit = 'cm')
  ) +
  scale_x_continuous(name = xname, limits = c(0, 0.28), breaks = seq(0, 0.25, by = 0.05), expand = c(0, 0.0025)) +
  scale_y_continuous(name = yname, limits = c(-0.65, 1.3), breaks=seq(-0.5,1.5,0.2),expand = c(0.01, 0))
#  annotate("text", x = max(d1)/10, y = max(d2) * 0.7, 
#           label = paste("n =", n_points),
#           size = 6, face = "bold", color = "black") 
dev.off()  # Close the PDF device
############################################
data1 <- read.table('K562_4sU_half_life.K562_NLS.out', header = FALSE, sep = '\t')

d1 <- data1[, 3]
d2 <- data1[, 2]
xname='Modification density\nin TAS-seq (K562, Nuc)'
yname= expression(bold("RNA half-life in K562"~log[10]~"(hours)"))

df <- data.frame(Rep1 = d1, Rep2 = d2)
n_points <- nrow(df)
cor_res <- cor.test(df$Rep1, df$Rep2, method = "pearson")  # Pearson correlation test
cor_value <- cor_res$estimate  # Pearson correlation coefficient
p_value <- cor_res$p.value  # P-value from the test
cor_value

pdf('K562_4sU_half_life.K562_NLS.cor.pdf', width = 6, height = 6)
ggscatter(df, x = 'Rep1', y = 'Rep2',
          color = '#277fb8', size = 2, alpha = 0.2,
          add = 'reg.line', 
          add.params = list(color = 'darkblue', fill = 'lightgray', size = 1.5),
          conf.int = TRUE, 
          #cor.coef = TRUE,
#          cor.coeff.args = list(method = 'pearson', label.x = max(d1)/10, label.sep = '\n', size = 6, face = 'bold', color = 'black')
) +
  theme_bw(base_rect_size = 2) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 20, face = 'bold', color = 'black'),
    axis.text = element_text(size = 16, face = 'bold', color = 'black'),
    axis.ticks = element_line(linewidth = 1.5),
    axis.ticks.length = unit(0.2, 'cm'),
    plot.margin = margin(0.2, 0.5, 0.2, 0.2, unit = 'cm')
  ) +
  scale_x_continuous(name = xname, limits = c(0, 0.4), breaks = seq(0, 0.4, by = 0.05), expand = c(0, 0.0025)) +
  scale_y_continuous(name = yname, limits = c(-0.65, 1.3), breaks = seq(-0.5, 1.5, by = 0.2), expand = c(0.01, 0))
#  annotate("text", x = max(d1)/10, y = max(d2) * 0.7, 
#           label = paste("n =", n_points),
#           size = 6, face = "bold", color = "black") 
dev.off()
