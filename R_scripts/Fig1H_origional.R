library(ggplot2)
library("scales")
library(gridExtra)
library(grid)
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors,rownames(qual_col_pals)))
col=c('#62CE0C','#1F78B4','#F0027F','#999999','#FC8D62',"#B15928","#6A3D9A","#FB8072","#A6CEE3")

text_5U <- textGrob("5'UTR", gp=gpar(fontsize=18, fontface="bold"))
text_CDS <- textGrob("CDS", gp=gpar(fontsize=18, fontface="bold"))
text_3U <- textGrob("3'UTR", gp=gpar(fontsize=18, fontface="bold"))

m6a.dist <- read.delim ("../NES.dist.measures.txt", header = T)
trx_len <- m6a.dist$utr5_size + m6a.dist$cds_size + m6a.dist$utr3_size # Determine transcript length
temp <- data.frame(m6a.dist$gene_name, m6a.dist$refseqID, trx_len)
colnames(temp) <- c("gene_name", "gid", "trx_len") 
temp.df <- temp[order(temp$gene_name, temp$gid, -temp$trx_len),]
temp.df <- temp[!duplicated(temp$gene_name),]
# limit m6a data to one transcript per gene (longest)
m6a.dist <- m6a.dist[m6a.dist$refseqID %in% temp.df$gid,]

utr5.SF <- median(m6a.dist$utr5_size, na.rm = T)/median(m6a.dist$cds_size, na.rm = T)
utr3.SF <- median(m6a.dist$utr3_size, na.rm = T)/median(m6a.dist$cds_size, na.rm = T)
# assign the regions to new dataframes
utr5.m6a.dist <- m6a.dist[m6a.dist$rel_location < 1, ]
cds.m6a.dist <- m6a.dist [m6a.dist$rel_location < 2 & m6a.dist$rel_location >= 1, ]
utr3.m6a.dist <- m6a.dist[m6a.dist$rel_location >= 2, ]
xrange=1-utr5.SF-0.2
yrange=3.3
utr5.m6a.dist$rel_location <- rescale(utr5.m6a.dist$rel_location, to = c(xrange, 1), from = c(0,1))
utr3.m6a.dist$rel_location <- rescale(utr3.m6a.dist$rel_location, to = c(2, yrange), from = c(2,3.5))
m6a.metagene.coord <- c(utr5.m6a.dist$rel_location, cds.m6a.dist$rel_location, utr3.m6a.dist$rel_location)
##########################
m6a.dist1 <- read.delim ("../NES.TA.dist.measures.txt", header = T)
trx_len1 <- m6a.dist1$utr5_size + m6a.dist1$cds_size + m6a.dist1$utr3_size
temp1 <- data.frame(m6a.dist1$gene_name, m6a.dist1$refseqID, trx_len1)
colnames(temp1) <- c("gene_name", "gid", "trx_len")
temp1.df <- temp1[order(temp1$gene_name, temp1$gid, -temp1$trx_len),]
temp1.df <- temp1[!duplicated(temp1$gene_name),]
m6a.dist1 <- m6a.dist1[m6a.dist1$refseqID %in% temp1.df$gid,]
utr5.m6a.dist1 <- m6a.dist1[m6a.dist1$rel_location < 1, ]
cds.m6a.dist1 <- m6a.dist1 [m6a.dist1$rel_location < 2 & m6a.dist1$rel_location >= 1, ]
utr3.m6a.dist1 <- m6a.dist1[m6a.dist1$rel_location >= 2, ]
xrange=1-utr5.SF-0.2
yrange=3.3
utr5.m6a.dist1$rel_location <- rescale(utr5.m6a.dist1$rel_location, to = c(xrange, 1), from = c(0,1))
utr3.m6a.dist1$rel_location <- rescale(utr3.m6a.dist1$rel_location, to = c(2, yrange), from = c(2,3.5))
m6a.metagene.coord1 <- c(utr5.m6a.dist1$rel_location, cds.m6a.dist1$rel_location, utr3.m6a.dist1$rel_location)
############################
m6a.dist <- read.delim ("../NES.AA.dist.measures.txt", header = T)
trx_len <- m6a.dist$utr5_size + m6a.dist$cds_size + m6a.dist$utr3_size # Determine transcript length
temp <- data.frame(m6a.dist$gene_name, m6a.dist$refseqID, trx_len)
colnames(temp) <- c("gene_name", "gid", "trx_len") 
temp.df <- temp[order(temp$gene_name, temp$gid, -temp$trx_len),]
temp.df <- temp[!duplicated(temp$gene_name),]
m6a.dist <- m6a.dist[m6a.dist$refseqID %in% temp.df$gid,]
utr5.SF <- median(m6a.dist$utr5_size, na.rm = T)/median(m6a.dist$cds_size, na.rm = T)
utr3.SF <- median(m6a.dist$utr3_size, na.rm = T)/median(m6a.dist$cds_size, na.rm = T)
utr5.m6a.dist <- m6a.dist[m6a.dist$rel_location < 1, ]
cds.m6a.dist <- m6a.dist [m6a.dist$rel_location < 2 & m6a.dist$rel_location >= 1, ]
utr3.m6a.dist <- m6a.dist[m6a.dist$rel_location >= 2, ]
xrange=1-utr5.SF-0.2
yrange=3.3
utr5.m6a.dist$rel_location <- rescale(utr5.m6a.dist$rel_location, to = c(xrange, 1), from = c(0,1))
utr3.m6a.dist$rel_location <- rescale(utr3.m6a.dist$rel_location, to = c(2, yrange), from = c(2,3.5))
m6a.metagene.coord2 <- c(utr5.m6a.dist$rel_location, cds.m6a.dist$rel_location, utr3.m6a.dist$rel_location)
#####################
m6a.dist <- read.delim ("../NES.CA.dist.measures.txt", header = T)
trx_len <- m6a.dist$utr5_size + m6a.dist$cds_size + m6a.dist$utr3_size # Determine transcript length
temp <- data.frame(m6a.dist$gene_name, m6a.dist$refseqID, trx_len)
colnames(temp) <- c("gene_name", "gid", "trx_len") 
temp.df <- temp[order(temp$gene_name, temp$gid, -temp$trx_len),]
temp.df <- temp[!duplicated(temp$gene_name),]
m6a.dist <- m6a.dist[m6a.dist$refseqID %in% temp.df$gid,]
utr5.SF <- median(m6a.dist$utr5_size, na.rm = T)/median(m6a.dist$cds_size, na.rm = T)
utr3.SF <- median(m6a.dist$utr3_size, na.rm = T)/median(m6a.dist$cds_size, na.rm = T)
utr5.m6a.dist <- m6a.dist[m6a.dist$rel_location < 1, ]
cds.m6a.dist <- m6a.dist [m6a.dist$rel_location < 2 & m6a.dist$rel_location >= 1, ]
utr3.m6a.dist <- m6a.dist[m6a.dist$rel_location >= 2, ]
xrange=1-utr5.SF-0.2
yrange=3.3
utr5.m6a.dist$rel_location <- rescale(utr5.m6a.dist$rel_location, to = c(xrange, 1), from = c(0,1))
utr3.m6a.dist$rel_location <- rescale(utr3.m6a.dist$rel_location, to = c(2, yrange), from = c(2,3.5))
m6a.metagene.coord3 <- c(utr5.m6a.dist$rel_location, cds.m6a.dist$rel_location, utr3.m6a.dist$rel_location)
##########################
#####################
m6a.dist <- read.delim ("../NES.GA.dist.measures.txt", header = T)
trx_len <- m6a.dist$utr5_size + m6a.dist$cds_size + m6a.dist$utr3_size # Determine transcript length
temp <- data.frame(m6a.dist$gene_name, m6a.dist$refseqID, trx_len)
colnames(temp) <- c("gene_name", "gid", "trx_len") 
temp.df <- temp[order(temp$gene_name, temp$gid, -temp$trx_len),]
temp.df <- temp[!duplicated(temp$gene_name),]
m6a.dist <- m6a.dist[m6a.dist$refseqID %in% temp.df$gid,]
utr5.SF <- median(m6a.dist$utr5_size, na.rm = T)/median(m6a.dist$cds_size, na.rm = T)
utr3.SF <- median(m6a.dist$utr3_size, na.rm = T)/median(m6a.dist$cds_size, na.rm = T)
utr5.m6a.dist <- m6a.dist[m6a.dist$rel_location < 1, ]
cds.m6a.dist <- m6a.dist [m6a.dist$rel_location < 2 & m6a.dist$rel_location >= 1, ]
utr3.m6a.dist <- m6a.dist[m6a.dist$rel_location >= 2, ]
xrange=1-utr5.SF-0.2
yrange=3.3
utr5.m6a.dist$rel_location <- rescale(utr5.m6a.dist$rel_location, to = c(xrange, 1), from = c(0,1))
utr3.m6a.dist$rel_location <- rescale(utr3.m6a.dist$rel_location, to = c(2, yrange), from = c(2,3.5))
m6a.metagene.coord4 <- c(utr5.m6a.dist$rel_location, cds.m6a.dist$rel_location, utr3.m6a.dist$rel_location)
##########################
m6a.dist1 <- read.delim ("../NES.annot.transcripts.UA.sites.dist.measures.txt", header = T)
trx_len1 <- m6a.dist1$utr5_size + m6a.dist1$cds_size + m6a.dist1$utr3_size
temp1 <- data.frame(m6a.dist1$gene_name, m6a.dist1$refseqID, trx_len1)
colnames(temp1) <- c("gene_name", "gid", "trx_len")
temp1.df <- temp1[order(temp1$gene_name, temp1$gid, -temp1$trx_len),]
temp1.df <- temp1[!duplicated(temp1$gene_name),]
m6a.dist1 <- m6a.dist1[m6a.dist1$refseqID %in% temp1.df$gid,]
utr5.m6a.dist1 <- m6a.dist1[m6a.dist1$rel_location < 1, ]
cds.m6a.dist1 <- m6a.dist1 [m6a.dist1$rel_location < 2 & m6a.dist1$rel_location >= 1, ]
utr3.m6a.dist1 <- m6a.dist1[m6a.dist1$rel_location >= 2, ]
xrange=1-utr5.SF-0.2
yrange=3.3
utr5.m6a.dist1$rel_location <- rescale(utr5.m6a.dist1$rel_location, to = c(xrange, 1), from = c(0,1))
utr3.m6a.dist1$rel_location <- rescale(utr3.m6a.dist1$rel_location, to = c(2, yrange), from = c(2,3.5))
m6a.metagene.coord5 <- c(utr5.m6a.dist1$rel_location, cds.m6a.dist1$rel_location, utr3.m6a.dist1$rel_location)
##########################
m6a.dist1 <- read.delim ("../NES.annot.transcripts.AA.sites.dist.measures.txt", header = T)
trx_len1 <- m6a.dist1$utr5_size + m6a.dist1$cds_size + m6a.dist1$utr3_size
temp1 <- data.frame(m6a.dist1$gene_name, m6a.dist1$refseqID, trx_len1)
colnames(temp1) <- c("gene_name", "gid", "trx_len")
temp1.df <- temp1[order(temp1$gene_name, temp1$gid, -temp1$trx_len),]
temp1.df <- temp1[!duplicated(temp1$gene_name),]
m6a.dist1 <- m6a.dist1[m6a.dist1$refseqID %in% temp1.df$gid,]
utr5.m6a.dist1 <- m6a.dist1[m6a.dist1$rel_location < 1, ]
cds.m6a.dist1 <- m6a.dist1 [m6a.dist1$rel_location < 2 & m6a.dist1$rel_location >= 1, ]
utr3.m6a.dist1 <- m6a.dist1[m6a.dist1$rel_location >= 2, ]
xrange=1-utr5.SF-0.2
yrange=3.3
utr5.m6a.dist1$rel_location <- rescale(utr5.m6a.dist1$rel_location, to = c(xrange, 1), from = c(0,1))
utr3.m6a.dist1$rel_location <- rescale(utr3.m6a.dist1$rel_location, to = c(2, yrange), from = c(2,3.5))
m6a.metagene.coord6 <- c(utr5.m6a.dist1$rel_location, cds.m6a.dist1$rel_location, utr3.m6a.dist1$rel_location)
##########################
m6a.dist1 <- read.delim ("../NES.annot.transcripts.CA.sites.dist.measures.txt", header = T)
trx_len1 <- m6a.dist1$utr5_size + m6a.dist1$cds_size + m6a.dist1$utr3_size
temp1 <- data.frame(m6a.dist1$gene_name, m6a.dist1$refseqID, trx_len1)
colnames(temp1) <- c("gene_name", "gid", "trx_len")
temp1.df <- temp1[order(temp1$gene_name, temp1$gid, -temp1$trx_len),]
temp1.df <- temp1[!duplicated(temp1$gene_name),]
m6a.dist1 <- m6a.dist1[m6a.dist1$refseqID %in% temp1.df$gid,]
utr5.m6a.dist1 <- m6a.dist1[m6a.dist1$rel_location < 1, ]
cds.m6a.dist1 <- m6a.dist1 [m6a.dist1$rel_location < 2 & m6a.dist1$rel_location >= 1, ]
utr3.m6a.dist1 <- m6a.dist1[m6a.dist1$rel_location >= 2, ]
xrange=1-utr5.SF-0.2
yrange=3.3
utr5.m6a.dist1$rel_location <- rescale(utr5.m6a.dist1$rel_location, to = c(xrange, 1), from = c(0,1))
utr3.m6a.dist1$rel_location <- rescale(utr3.m6a.dist1$rel_location, to = c(2, yrange), from = c(2,3.5))
m6a.metagene.coord7 <- c(utr5.m6a.dist1$rel_location, cds.m6a.dist1$rel_location, utr3.m6a.dist1$rel_location)
##########################
m6a.dist1 <- read.delim ("../NES.annot.transcripts.GA.sites.dist.measures.txt", header = T)
trx_len1 <- m6a.dist1$utr5_size + m6a.dist1$cds_size + m6a.dist1$utr3_size
temp1 <- data.frame(m6a.dist1$gene_name, m6a.dist1$refseqID, trx_len1)
colnames(temp1) <- c("gene_name", "gid", "trx_len")
temp1.df <- temp1[order(temp1$gene_name, temp1$gid, -temp1$trx_len),]
temp1.df <- temp1[!duplicated(temp1$gene_name),]
m6a.dist1 <- m6a.dist1[m6a.dist1$refseqID %in% temp1.df$gid,]
utr5.m6a.dist1 <- m6a.dist1[m6a.dist1$rel_location < 1, ]
cds.m6a.dist1 <- m6a.dist1 [m6a.dist1$rel_location < 2 & m6a.dist1$rel_location >= 1, ]
utr3.m6a.dist1 <- m6a.dist1[m6a.dist1$rel_location >= 2, ]
xrange=1-utr5.SF-0.2
yrange=3.3
utr5.m6a.dist1$rel_location <- rescale(utr5.m6a.dist1$rel_location, to = c(xrange, 1), from = c(0,1))
utr3.m6a.dist1$rel_location <- rescale(utr3.m6a.dist1$rel_location, to = c(2, yrange), from = c(2,3.5))
m6a.metagene.coord8 <- c(utr5.m6a.dist1$rel_location, cds.m6a.dist1$rel_location, utr3.m6a.dist1$rel_location)
#####################
groups_order <- c("Cyt", 'Cyt.UA','Cyt.AA','Cyt.CA','Cyt.GA','background.UA','background.AA','background.CA','background.GA')
groups_order <- c("Cyt", 'Cyt.UA','background.UA')
df=data.frame(
		value = c(m6a.metagene.coord, m6a.metagene.coord1, m6a.metagene.coord5),
		group = factor(rep(groups_order, times = c(length(m6a.metagene.coord),length(m6a.metagene.coord1),length(m6a.metagene.coord5))),
		levels = groups_order)
	       )

pdf('NES.UA_split.NES.pdf',wi=12,he=8)
ymax=1
############################
ggplot(df,aes(x = value, color = group)) +
  geom_density(size=1.5,linewidth=1.2,show.legend=FALSE)+
  stat_density(aes(x=value, colour=group),geom="line",position="identity",size=1.1)+
  geom_vline(xintercept = 1:2, col = "red",linetype="dashed", linewidth=1) +
  annotate('rect',xmin=1,xmax=2,ymin=-0.05,ymax=0,color="black",alpha=0.8,size=1,fill="grey50")+
  annotate('rect',xmin=0.6,xmax=1,ymin=-0.015,ymax=-0.035,alpha=1,fill="black")+
  annotate('rect',xmin=2,xmax=3,ymin=-0.015,ymax=-0.035,alpha=1,fill="black")+
  theme_bw(base_rect_size = 1) + 
  theme(
	panel.grid = element_blank(),
        axis.title = element_text(size=28,face = "bold",color='black'),
        axis.text = element_text(size=24,color='black',face = "bold"),
        axis.text.x = element_blank(),
	axis.ticks.x = element_blank(),
	axis.ticks.y = element_line(linewidth=1.5),
	axis.ticks.length.y = unit(0.2,'cm'),
	legend.text = element_text(size=28,face = "bold",color='black'),
	legend.position = 'top',
	legend.ticks = element_line(linewidth=2),
	legend.title = element_blank(),
        )+
  scale_x_continuous(name="")+
  annotation_custom(text_5U,xmin=1.05,xmax=0.5,ymin=-0.07,ymax=-0.07)+
  annotation_custom(text_CDS,xmin=1.5,xmax=1.5,ymin=-0.07,ymax=-0.07)+
  annotation_custom(text_3U,xmin=2.5,xmax=2.5,ymin=-0.07,ymax=-0.07)+
  coord_cartesian(ylim=c(-0.05,0.75), clip="off")+
  scale_color_manual(values=col)+
  scale_y_continuous(name="Density",limits = c(-0.06,ymax),expand = c(0,0),breaks = seq(0,ymax,0.1))+
  coord_cartesian(ylim=c(-0.09,ymax))
#  guides(color = guide_legend(nrow = 2))
dev.off()
