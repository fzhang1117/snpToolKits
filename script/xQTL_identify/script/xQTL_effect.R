## xQTL_effect.R ##
## zhang fei <zhangfei-123@foxmail.com> ##
## 2018-05-16 ##

Args <- commandArgs(T)
my.gpintersec <- function(genotype, phenotype){
    line.intersect <- intersect(phenotype$line, genotype$line)
    genotype <- genotype[genotype$line %in% line.intersect, -1]
    phenotype <- phenotype[phenotype$line %in% line.intersect, -1]
    result <- cbind(phenotype, genotype)
    result <- as.data.frame(result)
    return(result)
}

fl_hmp <- Args[1]
fl_xQTL <- Args[2]
fl_pheno_drought <- Args[3]
fl_pheno_normal <- Args[4]
prefix <- Args[5]

genotype.hmp <- read.table(fl_hmp, header = T, quote = "", sep = "\t", comment.char = "", row.names = 1, check.names = F, na.strings = "N")[ , c(-4:-10)]
xQTL.summary <- read.table(fl_xQTL, header = T, sep = "\t", check.names = F)
pheno.drought <- read.table(fl_pheno_drought, header = T, quote = "", sep = "\t", check.names = F, row.names = 1)
pheno.normal <- read.table(fl_pheno_normal, header = T, quote = "", sep = "\t", check.names = F, row.names = 1)

Effect <- rep(NA, dim(xQTL.summary)[1])
for(i in 1: dim(xQTL.summary)[1]){
    sigsnp <- t(genotype.hmp[genotype.hmp$chrom == xQTL.summary[i, 4] & genotype.hmp$pos >= xQTL.summary[i, 5] & genotype.hmp$pos <= xQTL.summary[i, 6], ])[c(-1, -2, -3) , ]
    if(is.null(dim(sigsnp)) == F){
    sigsnp <- sigsnp[order(rownames(sigsnp)), ]
    sigsnp <- as.data.frame(cbind(rownames(sigsnp), sigsnp))
    colnames(sigsnp)[1] <- c("line")
    }else{
        sigsnp <- data.frame(line = names(sigsnp), phenotype = sigsnp)
        sigsnp <- sigsnp[order(rownames(sigsnp)), ]
        #sigsnp <- as.data.frame(cbind(rownames(sigsnp), sigsnp))
        colnames(sigsnp)[1] <- c("line")
    }
    if(xQTL.summary[i, 2] == 'drought'){
        pheno <- data.frame(line = colnames(pheno.drought), trait = as.numeric(pheno.drought[rownames(pheno.drought) %in% xQTL.summary[i, 1], ]))
        pheno <- pheno[order(pheno$line), ]
    }else if(xQTL.summary[i, 2] == 'normal'){
        pheno <- data.frame(line = colnames(pheno.normal), trait = as.numeric(pheno.normal[rownames(pheno.normal) %in% xQTL.summary[i, 1], ]))
        pheno <- pheno[order(pheno$line), ]
    }
    data.temp <- my.gpintersec(sigsnp, pheno)
    fit.temp <- lm(phenotype ~., data = data.temp)
    Effect[i] <- summary(fit.temp)$r.squared
}
result <- cbind(xQTL.summary[, c(1, 2, 3, 4, 5, 6, 8, 10)], Effect)
fl_out = paste(prefix, "_Effect.txt", sep = "")
write.table(result, fl_out, sep = '\t', quote = F, row.names = F)
