## meta_analysis.R ##
## Script to perform meta analysis ##
## zhang fei <zhangfei-123@foxmail.com> ##
## 2018-05-29 ##

library(VennDiagram)
Args <- commandArgs(T)

literExtract <- c(4, 5)
my.cpdMetaAnalysis <- function(cpdMapping, Extract = literExtract){
    cpdMapping <- cpdMapping[ , c(2, 4)]
    cpdMapping <- cpdMapping[!(cpdMapping[ , 2] %in% Extract), ]
    cpdMapping <- cpdMapping[!duplicated(cpdMapping), ]
    literFactor <- cpdMapping[ , 2][!duplicated(cpdMapping[ , 2])]
    cpdMapping <- na.omit(cpdMapping)
    venn.list <- list()
    for(i in 1: length(literFactor)){
        venn.list[[i]] <- cpdMapping[cpdMapping[ , 2] %in% literFactor[i], 1]
    }
    names(venn.list) <- literFactor
    return(venn.list)
}

my.geneMetaAnalysis <- function(cpdMapping, Extract = literExtract){
    cpdMapping <- cpdMapping[ , c(3, 4)]
    cpdMapping <- cpdMapping[!(cpdMapping[ , 2] %in% Extract), ]
    cpdMapping <- cpdMapping[!duplicated(cpdMapping), ]
    literFactor <- cpdMapping[ , 2][!duplicated(cpdMapping[ , 2])]
    cpdMapping <- na.omit(cpdMapping)
    venn.list <- list()
    for(i in 1: length(literFactor)){
        entry <- cpdMapping[cpdMapping[ , 2] %in% literFactor[i], 1]
        entry <- entry[!duplicated(entry)]
        venn.list[[i]] <- entry
    }
    names(venn.list) <- literFactor
    return(venn.list)
}

my.cgpairMetaAnalysis <- function(cpdMapping, Extract = literExtract){
    cpdMapping_temp <- cpdMapping[ , c(2, 3, 4)]
    cpdMapping <- t(apply(cpdMapping_temp, 1, function(x){c(paste(x[1], x[2], sep = "-"), x[3])}))
    cpdMapping <- cpdMapping[!(cpdMapping[ , 2] %in% Extract), ]
    cpdMapping <- cpdMapping[!duplicated(cpdMapping), ]
    literFactor <- cpdMapping[ , 2][!duplicated(cpdMapping[ , 2])]
    cpdMapping <- na.omit(cpdMapping)
    venn.list <- list()
    for(i in 1: length(literFactor)){
        entry <- cpdMapping[cpdMapping[ , 2] %in% literFactor[i], 1]
        entry <- entry[!duplicated(entry)]
        venn.list[[i]] <- entry
    }
    names(venn.list) <- literFactor
    return(venn.list)
}

## candidate gene list from literture constucted by cpdMapping.py 
fl_candidate <- Args[1]
print(fl_candidate)
print(getwd())
## our candidate gene list constructed by cpdMapping.py
fl_myGene <- Args[2]
print(fl_myGene)

## prefix of output 
prefix <- Args[3]

cpdMapping_candidate <- read.table(fl_candidate, header = F, sep = "\t", quote = "", comment.char = "", skipNul = T)
#cpdMapping_candidate <- read.table("../../result/meta_5.00_normal_myGene.txt", header = F, sep = "\t", quote = "")
print(dim(cpdMapping_candidate))
cpdMapping_myGene <- read.table(fl_myGene, header = F, sep = "\t", quote = "", comment.char = "", skipNul = T)
print(dim(cpdMapping_myGene))
cpdMapping_merge <- rbind(cpdMapping_candidate, cpdMapping_myGene)


venn.cpd <- my.cpdMetaAnalysis(cpdMapping = cpdMapping_merge)
venn.diagram(venn.cpd, file = paste(prefix, "_venn_cpd.png", sep = ""))
venn.gene <- my.geneMetaAnalysis(cpdMapping = cpdMapping_merge)
venn.diagram(venn.gene, file = paste(prefix, "_venn_gene.png", sep = ""))
venn.cgpair <- my.cgpairMetaAnalysis(cpdMapping = cpdMapping_merge)
venn.diagram(venn.cgpair, file = paste(prefix, "_venn_cgpair.png", sep = ""))
