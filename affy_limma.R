required_Packages = c("affy", "limma","hgu133plus2.db","sva","ggplot2")

if(!all(required_Packages %in% installed.packages())){
  source("https://bioconductor.org/biocLite.R")
  biocLite(setdiff(required_Packages, installed.packages()))
}

require(affy)
require(limma)
require(sva)
require(hgu133plus2hsrefseqcdf)
require(ggplot2)

rm(list = ls())

#need GSEnum
GSEnum <- 'GSE12056'
designname <- 'A2_GSE12056_CREB1_KO.txt'
#need design_path
design_path <- paste0('microarray/design_matrix/',designname)
#need data_path
file_path <- paste0('microarray/data/',GSEnum,'/')
#need result_path
result_path <- paste0('microarray/result/',designname)

# args = commandArgs(T)
# file_path = args[1] #data file path
# design_path = args[2] #design matrix path
# result_path = args[3] #result matrix path

# read design
design_mat = read.table(design_path, sep="\t", header=F)
gsm = as.vector(design_mat[,1])
grouplist = as.vector(design_mat[,2])

# read files in and process data
celFiles <- list.celfiles(path = file_path, full.names = T)
data.affy <- ReadAffy(filenames = celFiles)
data.affy@cdfName = "hgu133plus2hsrefseqcdf"
data.rma <- rma(data.affy)
data.expr <- exprs(data.rma)
colnames(data.expr) <- sapply(sapply(colnames(data.expr),strsplit,'_',USE.NAMES = F),head,1,USE.NAMES = F)

# batcheffect
data.affy@protocolData@data
# groupName <- colnames(data.expr)
# batchIndex <- c(1,2,2,3,3,1)
# Condition <- c(1,1,1,2,2,2)
# batchInfo <- data.frame(groupName, batchIndex, Condition)
# cbmod <- model.matrix(~Condition, data = batchInfo)
# data.expr <- ComBat(dat = data.expr, batch = batchIndex, mean.only = TRUE, mod = cbmod)

#get exprSet
exprSet <- data.expr[,gsm]

#do limma
design <- model.matrix(~0+factor(grouplist))
colnames(design) = levels(factor(grouplist))
rownames(design) = colnames(exprSet)
fit <- lmFit(exprSet,design)
contrast.matrix <- makeContrasts(paste0(rev(unique(grouplist)),collapse = '-'), levels = design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)
tempoutput <- topTable(fit, number = Inf, lfc = 0, p.value = 1)

#Get gene symbol
rownames(tempoutput) <- sapply(sapply(rownames(tempoutput),strsplit,'[.]',USE.NAMES = F),head,1,USE.NAMES = F)
tempoutput$SYMBOL <- select(org.Hs.eg.db, keys=rownames(tempoutput), columns="SYMBOL", keytype="REFSEQ")[,2]
output <- na.omit(tempoutput[!duplicated(tempoutput$SYMBOL),])
rownames(output) <- output$SYMBOL
output <- output[,-7]
colnames(output) <- c("log2FoldChange","AveExpr","t","P.Value","padj","B")

sigout <- subset(output, abs(log2FoldChange) > 1 & padj < 0.05)
p <- ggplot(output,aes(log2FoldChange,-log10(padj))) + 
geom_point(aes(color=ifelse((log2FoldChange > 1 | log2FoldChange < -1) & padj < 0.05, "sign","non-sign")), cex = 1, alpha = 0.3) + 
scale_colour_manual("significance", values = c("steelblue","red")) +
labs(title = designname, x = "log2 Fold Change", y = "-log10 adj.P.Val")+
geom_text(data = sigout, aes(x = log2FoldChange, y=-log10(padj), label = rownames(sigout)), cex = 2, vjust = "inward", hjust = "inward")
ggsave(paste0(result_path,".png"), p, width = 7, height = 7)

write.table(output,result_path,sep = '\t',quote = F)

