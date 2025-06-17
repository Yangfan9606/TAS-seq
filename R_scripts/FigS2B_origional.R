setwd('D:/Data/0_lixia/0)_8e_HTseq_all/0_Paper_Figure_Data_source_2024/0_d10mu2/050_RPKM_correlation_PLOT')

data = read.table("Re.11_samples.8e_HTseq.all.FPKM",sep="\t",header=T)
a = data[,c(2,3,4,5,6,7,8,9,12,13,14,15,16,17)]
names(a)=c("Ctl_Rep1(PolyA+)","Ctl_Rep2(PolyA+)","Cyt_Rep1(PolyA+)","Cyt_Rep2(PolyA+)","Nuc_Rep1(PolyA+)","Nuc_Rep2(PolyA+)","Erm_Rep1(PolyA+)","Erm_Rep2(PolyA+)",'Ctl_Rep1(Ribo-)','Ctl_Rep2(Ribo-)','Cyt_Rep1(Ribo-)','Cyt_Rep2(Ribo-)','Nuc_Rep1(Ribo-)','Nuc_Rep2(Ribo-)')

library("corrplot")
a=log10(a+1)
matrix = cor(a)
#pdf("8e_HTseq.11_samples.correlation.pdf",width = 16,height = 16)
pdf("No_CytH_8e_HTseq.11_samples.correlation.pdf",width = 16,height = 16)

corrplot(matrix,
		is.corr = FALSE,
#		diag = FALSE,
		addrect = 2, rect.col = 'black',
		rect.lwd = 5,
		method="color", col = colorRampPalette(c("#318ce7","#ff2052"))(5),
#		col.lim = c(0.7, 1),
		order = "hclust",  cl.pos = 'b', cl.cex=1.5, addgrid.col = 'white',
#		tl.pos = 'n',
		tl.col = "black", tl.cex=1.5, 
		bg = "white", addCoef.col = "black",
		cex.main=1.1,mar=c(0, 0, 0, 0), number.cex = 1.5,
		)
dev.off()
