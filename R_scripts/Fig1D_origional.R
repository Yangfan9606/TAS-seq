##############################
# part 1
##############################
library(VennDiagram)
Af = read.table("NES.r1r2.chr.site.unique.gene.temp", header = F, sep='\t')
colnames(Af)=c('V1')
A=Af$V1
Bf = read.table("NLS.r1r2.chr.site.unique.gene.temp", header = F, sep='\t')
colnames(Bf)=c('V1')
B=Bf$V1
Cf = read.table("ERM.r1r2.chr.site.unique.gene.temp", header = F, sep='\t')
colnames(Cf)=c('V1')
C=Cf$V1
# intersect() 取交集
Length_A=nrow(Af)
Length_B=nrow(Bf)
Length_C=nrow(Cf)
Length_AB=length(intersect(A,B))
Length_AC=length(intersect(A,C))
Length_BC=length(intersect(B,C))
Length_ABC=length(intersect(intersect(A,B),C))
pdf('NES_NLS_ERM.sites.overlap.VennDiagram.pdf', wi=7,he=7)
draw.triple.venn(
area1=Length_A, area2=Length_B, area3=Length_C,
n12=Length_AB, n23=Length_BC, n13=Length_AC,
n123=Length_ABC,
category=c('Cyt','Nuc','Erm'),
lwd=1,		#边框线宽度
alpha = 0.75,		#透明度
fill=c('#99C4B2','#FFC799','#FE999A'),	#填充色
margin = 0.05,		#画图边际距离
print.mode=c('raw','percent'),
sigdigs = 4,		#percent mode下保留有效数字个数
fontfamily = 'ArialMT',
cex = 1.5,		#标签字体大小
cat.cex = 1.5,		#category 即标签的字体大小
cat.pos=c(0,0,180),
)
dev.off()

##############################
# part 2
##############################
library(eulerr)

venn_data <- c(
  Cyt = 380078,      # 只属于C的元素
  Nuc = 127356,      # 只属于N的元素
  Erm = 10798,       # 只属于E的元素
  "Cyt&Nuc" = 171214,  # 属于C和N但不属于E的元素
  "Cyt&Erm" = 84608,    # 属于C和E但不属于N的元素
  "Nuc&Erm" = 4258,    # 属于N和E但不属于C的元素
  "Cyt&Nuc&Erm" = 272035 # 属于C、N和E的元素
)

# 绘制Venn图c("#1957FF", "#FF750C", "#FF220C")
pdf('NES_NLS_ERM.venn.pdf',he=8,wi=8)
par(oma=c(2,6,2,6), mar=c(6,6,6,2))
p=plot(euler(venn_data, shape = "ellipse"),
#     quantities = list(type = c("counts","percent"), cex=2),
     fill = c("#1F78B4","#F0027F", "#BF5B17"),
#     labels = list(cex=5),
     alpha = 0.4)
p
dev.off()
