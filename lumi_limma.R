required_Packages = c("lumi", "limma", "GEOquery")

if(!all(required_Packages %in% installed.packages())){
	source("https://bioconductor.org/biocLite.R")
	biocLite(setdiff(required_Packages, installed.packages()))
}

require(lumi)
require(limma)
require(GEOquery)

GSE <- getGEO('GSE30622', destdir=".",getGPL = F)
lumi.N.Q <- lumiExpresso(GSE[[1]])
exprSet <- exprs(lumi.N.Q)[,c("GSM759649","GSM759650","GSM759653","GSM759654")]
saveRDS(exprSet,file = 'expr.rds')
ID <- IlluminaID2nuID(rownames(exprSet),species = "Human")
rownames(exprSet) <- ID[,4]

#rownames(exprSet) <- sapply(sapply(rownames(exprSet),strsplit,"[.]"),head,1)

grouplist <- c('ctrl','ctrl','treat','treat')
design <- model.matrix(~0+factor(grouplist))
colnames(design) = levels(factor(grouplist))
rownames(design) = colnames(exprSet)
fit <- lmFit(exprSet,design)
contrast.matrix <- makeContrasts(ctrl - treat, levels = design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)
output <- topTable(fit, number = Inf, lfc = log2(1.5), p.value = 0.05)
write.table(output,'Result_limma.txt',sep = '\t',quote = F)
