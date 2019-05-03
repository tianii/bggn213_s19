# class 5 R graphics and plots 
#?read.table
weight <- read.table("bimm143_05_rstats/weight_chart.txt", header = T)
plot(weight, typ="o", pch=15, cex=1.5, lwd=2, ylim=c(2,10), xlab="Age (months)", ylab="Weight(kg)", main="FAT BABIES > EVERYONE ELSE")
read.table("bimm143_05_rstats/feature_counts.txt", header = T, sep = "\t")


read.table("bimm143_05_rstats/feature_counts.txt", header = T, sep = "\t")
par(mar=c(5, 13, 5.5, 3.5))

par()$mar
#?barplot
# use the par function to change parameters of the graph including the margins
mf <- read.delim("bimm143_05_rstats/male_female_counts.txt", header = T)
barplot(mf$Count, names.arg = mf$Sample, las=2, 
        col = rainbow(10), ylim=c(0,20))
#coloring by value
genes <- read.delim("bimm143_05_rstats/up_down_expression.txt", header = T)
table(genes$State)
plot(genes$Condition1, genes$Condition2, col=genes$State)
palette()
levels(genes$State)
palette(c("pink", "grey", "orange"))
plot(genes$Condition1, genes$Condition2, 
     col=genes$State, xlab="Expression condition 1", 
     ylab="Expression condition 2")
#coloring by point density 
meth <- read.delim("bimm143_05_rstats/expression_methylation.txt")
plot(meth$gene.meth, meth$expression)
dcols <- densCols(meth$gene.meth,  meth$expression)
plot(meth$gene.meth, meth$expression, col = dcols, pch=20)
inds <- meth$expression > 0 
plot(meth$gene.meth[inds], meth$expression[inds], col=dcols, pch=20)
dcols.custom <- densCols(meth$gene.meth[inds], meth$expression[inds],
                         colramp = colorRampPalette(c("pink",
                                                      "orange",
                                                      "purple",
                                                      "red")) )
plot(meth$gene.meth[inds], meth$expression[inds], 
     col = dcols.custom, pch = 20)

