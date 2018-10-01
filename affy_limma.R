required_Packages = c("affy", "limma","hgu133plus2.db","sva")

if(!all(required_Packages %in% installed.packages())){
  source("https://bioconductor.org/biocLite.R")
  biocLite(setdiff(required_Packages, installed.packages()))
}

library(affy)
library(limma)
library(hgu133plus2.db)
library(sva)


rm(list = ls())
#need GSEnum
GSEnum <- 'GSE51447'
designname <- 'A2_GSE51447_CREB1_KD.txt'
#need design_path
design_path <- paste0('microarray/design_matrix/',designname)
#need data_path
file_path <- paste0('microarray/data/',GSEnum,'/')
#need result_path
result_path <- paste0('microarray/result/',designname)

# read design
design_mat = read.table(design_path, sep="\t", header=F)
gsm = as.vector(design_mat[,1])
grouplist = as.vector(design_mat[,2])

# read files in and process data
celFiles <- list.celfiles(path = file_path, full.names = T)
data.affy <- ReadAffy(filenames = celFiles)
data.rma <- rma(data.affy)
data.expr <- exprs(data.rma)
colnames(data.expr) <- sapply(colnames(data.expr),substr,1,10,USE.NAMES = F)

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
contrast.matrix <- makeContrasts(paste0(unique(grouplist),collapse = '-'), levels = design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)
output <- topTable(fit, number = Inf, lfc = log2(1.5), p.value = 0.05)

#Get gene symbol
Annot <- data.frame(SYMBOL=sapply(contents(hgu133plus2SYMBOL), paste, collapse=", "))
cpltedata <- merge(output, Annot, by.x=0, by.y=0, all = TRUE)
output <- na.omit(cpltedata)
output <- output[order(output$adj.P.Val),]
write.table(output,result_path,sep = '\t',quote = F)



