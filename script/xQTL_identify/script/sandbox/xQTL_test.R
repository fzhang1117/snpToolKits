## xQTL_merge.R ##
## merge the xQTLs in a trait which phenotype is greater than the threshold ##
## zhang fei <zhangfei-123@foxmail.com> ##
## 2018-04-14 ##

options(stringsAsFactors = F)
Args <- commandArgs(T)

fl_xQTL <- Args[1]
fl_ld <- Args[2]
Prefix <- Args[3]
ld_threshold = as.numeric(Args[4])
singlePoint = as.logical(Args[5])

workingdir <- "."
setwd(workingdir)

xQTL_summary <- read.table(fl_xQTL, header = T, sep = "\t")
xQTL_ld <- read.table(fl_ld, header = T, sep = "\t")

trait_group <- xQTL_summary$trait
trait_group <- trait_group[!duplicated(trait_group)]

xQTL_rest <- data.frame()
for(trait in trait_group){
    xQTL_subset <- xQTL_summary[xQTL_summary$trait == trait, ]
    ld_subset <- xQTL_ld[xQTL_ld$trait == trait & xQTL_ld$R2 >= ld_threshold, ]
    remove <- c()
	if(dim(ld_subset)[1] != 0){
		for(i in 1: dim(ld_subset)[1]){
			ld_pair <- ld_subset[i, c(2, 3)]
			xQTL_pair <- xQTL_subset[xQTL_subset$xQTL %in% as.character(ld_pair), ]
			xQTL_remove <- xQTL_pair[xQTL_pair$leadp == max(xQTL_pair$leadp), 2]
			remove <- append(remove, xQTL_remove)
		}
	}
    remove <- remove[!duplicated(remove)]
    xQTL_subset <- xQTL_subset[!xQTL_subset$xQTL %in% remove, ]
    xQTL_rest <- rbind(xQTL_rest, xQTL_subset)
}

## metabolite filter ##
metabolite_filter <- c("GC_047", "GC_075", "GC_080", "GC_082", "GC_088")
if(singlePoint == T){
	xQTL_rest = xQTL_rest[xQTL_rest$snp_number != 1 & !(xQTL_rest$trait %in% metabolite_filter), ]
}else if(singlePoint == F){
	xQTL_rest = xQTL_rest[!(xQTL_rest$trait %in% metabolite_filter), ]
}

fl_output <- paste(Prefix, "_QTLfilter_summary.txt", sep = "")
write.table(xQTL_rest, fl_output, quote = F, sep = "\t", row.names = F)
