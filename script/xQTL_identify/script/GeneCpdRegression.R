## GeneCpdRegression.R ##
## Zhagn Fei <zhangfei-123@foxmail.com> ##
## 2018-05-14 ##

workingdir <- "."
setwd(workingdir)

## might use this method further, but I can't install data.table package in R 
#library(data.table)

#library(ggplot2)
options(stringsAsFactors = F)
Args <- commandArgs(T)

fl_cpd_normal <- Args[1]
fl_cpd_drought <- Args[2]
fl_gene_normal <- Args[3]
fl_gene_drought <- Args[4]
fl_candidateGene <- Args[5]
prefix <- Args[6]

candidateGene <- read.table(fl_candidateGene, header = T, sep = "\t", check.names = F, comment.char = "", quote = "")
cpd_normal <- read.table(fl_cpd_normal, header = T, row.names = 1, sep = "\t", check.names = F)
cpd_drought <- read.table(fl_cpd_drought, header = T, row.names = 1, sep = "\t", check.names = F)
gene.normal <- read.table(fl_gene_normal, header = T, row.names = 1, sep = "\t", check.names = F)
gene.drought <- read.table(fl_gene_drought, header = T, row.names = 1, sep = "\t", check.names = F)

# candidateGene <- read.table("./test_100kb_CandidateGene.txt", header = T, sep = "\t", check.names = F, comment.char = "", quote = "")
# cpd_normal <- read.table("../data/GC_DN_Normal_GWAS.txt", header = T, row.names = 1, sep = "\t", check.names = F)
# cpd_drought <- read.table("../data/GC_DN_Drought_GWAS.txt", header = T, row.names = 1, sep = "\t", check.names = F)
# gene.normal <- read.table("../data/ephenoMatrixC", header = T, row.names = 1, sep = "\t", check.names = F)
# gene.drought <- read.table("../data/ephenoMatrixD", header = T, row.names = 1, sep = "\t", check.names = F)

candidateGene <- candidateGene[which(candidateGene[ , 16] != "-"), ]
gene_normal <- as.data.frame(t(gene.normal[ ,colnames(gene.normal) %in% candidateGene[ , 16]]))
gene_drought <- as.data.frame(t(gene.drought[ , colnames(gene.drought) %in% candidateGene[ , 16]]))
material_intersect <- intersect(colnames(cpd_normal), colnames(gene_normal))

fl_pdf <- paste(prefix, "_plot.pdf", sep = "")
fl_text <- paste(prefix, "_gene_cpd_regression.txt", sep = "")

pdf(fl_pdf, height = 6, width = 12)
#pdf("./test2.pdf", height = 6, width = 12)
met <- matrix(c(1, 2), 1, 2)
layout(met)
result <- data.frame(cpd = candidateGene[ , 1], gene_v3 = candidateGene[ , 15], gene_v4 = candidateGene[ , 16], 
                     R2.normal = rep(NA, dim(candidateGene)[1]), p.normal = rep(NA, dim(candidateGene)[1]), 
                     R2.drought = rep(NA, dim(candidateGene)[1]), p.drought = rep(NA, dim(candidateGene)[1]), 
                     pathway = candidateGene[ , 13], pathwayAnno = candidateGene[ , 14], anno = candidateGene[ , 30])

for(i in 1: dim(candidateGene)[1]){
    meta.normal <- cpd_normal[rownames(cpd_normal) %in% candidateGene[i, 1], colnames(cpd_normal) %in% material_intersect]
    meta.normal <- as.numeric(meta.normal[ , order(colnames(meta.normal))])
    
    gene.normal <- gene_normal[rownames(gene_normal) %in% candidateGene[i, 16], colnames(gene_normal) %in% material_intersect]
    gene.normal <- as.numeric(gene.normal[ , order(colnames(gene.normal))])
    
    meta.drought <- cpd_drought[rownames(cpd_drought) %in% candidateGene[i, 1], colnames(cpd_drought) %in% material_intersect]
    meta.drought <- as.numeric(meta.drought[ , order(colnames(meta.drought))])
    
    gene.drought <- gene_drought[rownames(gene_drought) %in% candidateGene[i, 16], colnames(gene_drought) %in% material_intersect]
    gene.drought <- as.numeric(gene.drought[ , order(colnames(gene.drought))])
    
    data.normal <- data.frame(meta = meta.normal, gene = gene.normal)
    data.drought <- data.frame(meta = meta.drought, gene = gene.drought)
    fit.normal <- lm(meta ~ gene, data = data.normal)
    fit.drought <- lm(meta ~ gene, data = data.drought)
    summary.normal <- summary(fit.normal)
    summary.drought <- summary(fit.drought)
    if(dim(summary.normal$coefficients)[1] > 1){
        result[i, c(4, 5)] <-  c(summary.normal$r.squared, summary.normal$coefficients[2, 4])
        formula <- paste('y = ', as.numeric(round(fit.normal$coefficients[2], 4)), 'x + ', as.numeric(round(fit.normal$coefficients[1], 4)), sep = "")
        plot(x = c(0.95*min(data.normal$gene) - 0.5, 1.05*max(data.normal$gene) + 0.5), y = c(0.95*min(data.normal$meta, na.rm = T), 1.05*max(data.normal$meta, na.rm = T)), 
             xlim = c(0.95*min(data.normal$gene) - 0.5, 1.05*max(data.normal$gene) + 0.5), ylim = c(0.95*min(data.normal$meta, na.rm = T) - 0.5, 1.05*max(data.normal$meta, na.rm = T) + 0.5), 
             type = "n", xlab = paste("gene: ", result$gene_v3[i], sep = ""), ylab = paste("cpd: ", result$cpd[i], sep = ""), main = "Normal")
        x.length = 1.05* max(data.normal$gene) - 0.95*min(data.normal$gene) + 1
        y.length = 1.05*max(data.normal$meta, na.rm = T) - 0.95*min(data.normal$meta, na.rm = T) + 1
        points(meta ~ gene, data = data.normal)
        abline(fit.normal)
        text(0.95*min(data.normal$gene) - 0.5 + 0.85*x.length, 0.95*min(data.normal$meta, na.rm = T) - 0.5 + 0.95*y.length, formula, cex = 0.8)
        text(0.95*min(data.normal$gene) - 0.5 + 0.85*x.length, 0.95*min(data.normal$meta, na.rm = T) - 0.5 + 0.90*y.length, paste("p = ", round(result[i, 5], 4), sep = ""), cex = 0.8)
        text(0.95*min(data.normal$gene) - 0.5 + 0.85*x.length, 0.95*min(data.normal$meta, na.rm = T) - 0.5 + 0.85*y.length, paste("|R| = ", round(sqrt(result[i, 4]), 4), sep = ""), cex = 0.8)
    }else{
        plot(x = c(0.95*min(data.normal$gene) - 0.5, 1.05*max(data.normal$gene) + 0.5), y = c(0.95*min(data.normal$meta, na.rm = T), 1.05*max(data.normal$meta, na.rm = T)), 
             xlim = c(0.95*min(data.normal$gene) - 0.5, 1.05*max(data.normal$gene) + 0.5), ylim = c(0.95*min(data.normal$meta, na.rm = T) - 0.5, 1.05*max(data.normal$meta, na.rm = T) + 0.5), 
             type = "n", xlab = paste("gene: ", result$gene_v3[i], sep = ""), ylab = paste("cpd: ", result$cpd[i], sep = ""), main = "Normal")
        x.length = 1.05* max(data.normal$gene) - 0.95*min(data.normal$gene) + 1
        y.length = 1.05*max(data.normal$meta, na.rm = T) - 0.95*min(data.normal$meta, na.rm = T) + 1
        points(meta ~ gene, data = data.normal)
    }
    if(dim(summary.drought$coefficients)[1] > 1){
        result[i, c(6, 7)] <-  c(summary.drought$r.squared, summary.drought$coefficients[2, 4])
        formula <- paste('y = ', as.numeric(round(fit.drought$coefficients[2], 4)), 'x + ', as.numeric(round(fit.drought$coefficients[1], 4)), sep = "")
        plot(x = c(0.95*min(data.drought$gene) - 0.5, 1.05*max(data.drought$gene) + 0.5), y = c(0.95*min(data.drought$meta, na.rm = T), 1.05*max(data.drought$meta, na.rm = T)), 
             xlim = c(0.95*min(data.drought$gene) - 0.5, 1.05*max(data.drought$gene) + 0.5), ylim = c(0.95*min(data.drought$meta, na.rm = T) - 0.5, 1.05*max(data.drought$meta, na.rm = T) + 0.5), 
             type = "n", xlab = paste("gene: ", result$gene_v3[i], sep = ""), ylab = paste("cpd: ", result$cpd[i], sep = ""), main = "Drought")
        x.length = 1.05* max(data.drought$gene) - 0.95*min(data.drought$gene) + 1
        y.length = 1.05*max(data.drought$meta, na.rm = T) - 0.95*min(data.drought$meta, na.rm = T) + 1
        points(meta ~ gene, data = data.drought)
        abline(fit.drought)
        text(0.95*min(data.drought$gene) - 0.5 + 0.85*x.length, 0.95*min(data.drought$meta, na.rm = T) - 0.5 + 0.95*y.length, formula, cex = 0.8)
        text(0.95*min(data.drought$gene) - 0.5 + 0.85*x.length, 0.95*min(data.drought$meta, na.rm = T) - 0.5 + 0.90*y.length, paste("p = ", round(result[i, 5], 4), sep = ""), cex = 0.8)
        text(0.95*min(data.drought$gene) - 0.5 + 0.85*x.length, 0.95*min(data.drought$meta, na.rm = T) - 0.5 + 0.85*y.length, paste("|R| = ", round(sqrt(result[i, 4]), 4), sep = ""), cex = 0.8)
    }else{
        plot(x = c(0.95*min(data.drought$gene) - 0.5, 1.05*max(data.drought$gene) + 0.5), y = c(0.95*min(data.drought$meta, na.rm = T), 1.05*max(data.drought$meta, na.rm = T)), 
             xlim = c(0.95*min(data.drought$gene) - 0.5, 1.05*max(data.drought$gene) + 0.5), ylim = c(0.95*min(data.drought$meta, na.rm = T) - 0.5, 1.05*max(data.drought$meta, na.rm = T) + 0.5), 
             type = "n", xlab = paste("gene: ", result$gene_v3[i], sep = ""), ylab = paste("cpd: ", result$cpd[i], sep = ""), main = "drought")
        x.length = 1.05* max(data.drought$gene) - 0.95*min(data.drought$gene) + 1
        y.length = 1.05*max(data.drought$meta, na.rm = T) - 0.95*min(data.drought$meta, na.rm = T) + 1
        points(meta ~ gene, data = data.drought)
    }
    
}
dev.off()
write.table(result, fl_text, quote = F, sep = "\t", row.names = F)
