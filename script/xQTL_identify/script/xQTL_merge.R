## xQTL_merge.R ##
## merge the xQTLs in a trait which phenotype is greater than the threshold ##
## zhang fei <zhangfei-123@foxmail.com> ##
## 2018-04-14 ##
## update: 2018-05-09: add function in drought and normal condition filter 

options(stringsAsFactors = F)
Args <- commandArgs(T)

fl_xQTL <- Args[1]
fl_ld <- Args[2]
Prefix <- Args[3]
ld_threshold = as.numeric(Args[4])
singlePoint = as.logical(Args[5])
# 
# ld_threshold <- 0.1
# singlePoint <- T

workingdir <- "."
setwd(workingdir)

# xQTL_summary <- read.table("mlm_GC_6.09_20kb_0.1_filter_xQTL_summary.txt", header = T, sep = "\t")
# xQTL_ld <- read.table("mlm_GC_6.09_20kb_0.1_filter_xQTL_ld.txt", header = T, sep = "\t")

xQTL_summary <- read.table(fl_xQTL, header = T, sep = "\t")
xQTL_ld <- read.table(fl_ld, header = T, sep = "\t")

metabolite_filter <- c("GC_047", "GC_075", "GC_080", "GC_082", "GC_088")

if(singlePoint == T){
    xQTL_summary <- xQTL_summary[xQTL_summary$snp_number > 1 & !(xQTL_summary$trait %in% metabolite_filter), ]
    xQTL_list <- xQTL_summary$xQTL
    xQTL_list <- sort(xQTL_list[!duplicated(xQTL_list)])
    xQTL_ld <- xQTL_ld[(xQTL_ld$xQTL1 %in% xQTL_list & xQTL_ld$xQTL2 %in% xQTL_list), ]
}else{
    xQTL_summary <- xQTL_summary[!(xQTL_summary$trait %in% metabolite_filter), ]
    xQTL_list <- xQTL_summary$xQTL
    xQTL_list <- xQTL_list[!duplicated(xQTL_list)]
    xQTL_ld <- xQTL_ld[(xQTL_ld$xQTL1 %in% xQTL_list & xQTL_ld$xQTL2 %in% xQTL_list), ]
}

xQTL_drought <- xQTL_summary[xQTL_summary$condition == 'drought', ]
xQTL_normal <- xQTL_summary[xQTL_summary$condition == 'normal', ]

trait_group <- xQTL_summary$trait
trait_group <- trait_group[!duplicated(trait_group)]

xQTL_rest <- data.frame()
remove.drought <- NULL
remove.normal <- NULL
# xQTL_drought_rest <- NULL
# xQTL_normal_rest <- NULL

for(trait in trait_group){
    xQTL_subset <- xQTL_summary[xQTL_summary$trait == trait, ]
    ld_subset_drought <- xQTL_ld[xQTL_ld$condition == 'drought' & xQTL_ld$trait == trait & xQTL_ld$R2 >= ld_threshold, ]
    ld_subset_normal <- xQTL_ld[xQTL_ld$condition == 'normal' & xQTL_ld$trait == trait & xQTL_ld$R2 >= ld_threshold, ]
    xQTL_subset_drought <- xQTL_subset[xQTL_subset$condition == 'drought', ]
    xQTL_subset_normal <- xQTL_subset[xQTL_subset$condition == 'normal', ]

#     if(dim(ld_subset)[1] != 0){
#         ld_subset_drought <- ld_subset_drought[ld_subset$xQTL1 %in% xQTL_subset_drought$xQTL & ld_subset$xQTL2 %in% xQTL_subset_drought$xQTL, ]
#         ld_subset_normal <- ld_subset_normal[ld_subset$xQTL1 %in% xQTL_subset_normal$xQTL & ld_subset$xQTL2 %in% xQTL_subset_normal$xQTL, ]
#     }
    
#### need to re check, bug fix ####  
    if(dim(ld_subset_drought)[1] > 1){
        for(i in 1: dim(ld_subset_drought)[1]){
            xQTL_pair <- c(as.character(ld_subset_drought$xQTL1[i]), as.character(ld_subset_drought$xQTL2[i]))
            leadp_pair <- xQTL_subset_drought[xQTL_subset_drought$xQTL %in% xQTL_pair, 10]
            if(length(leadp_pair) > 1){
                data_temp <- data.frame(trait = rep(trait, 2), xQTL = xQTL_pair, leadp = leadp_pair)
                xQTL_remove <- data_temp[data_temp$leadp == max(data_temp$leadp), c(1, 2)]
                remove.drought <- rbind(remove.drought, xQTL_remove)    
            }
        }
    }
    if(dim(ld_subset_normal)[1] > 1){
        for(i in 1: dim(ld_subset_normal)[1]){
            xQTL_pair <- c(as.character(ld_subset_normal$xQTL1[i]), as.character(ld_subset_normal$xQTL2[i]))
            leadp_pair <- xQTL_subset_normal[xQTL_subset_normal$xQTL %in% xQTL_pair, 10]
            if(length(leadp_pair) > 1){
                data_temp <- data.frame(trait = rep(trait, 2), xQTL = xQTL_pair, leadp = leadp_pair)
                xQTL_remove <- data_temp[data_temp$leadp == max(data_temp$leadp), c(1, 2)]
                remove.normal <- rbind(remove.normal, xQTL_remove)    
            }
        }
    }
    remove.drought <- remove.drought[!duplicated(remove.drought), ]
    remove.normal <- remove.normal[!duplicated(remove.normal), ]
#     xQTL_drought_rest <- rbind(xQTL_drought_rest, xQTL_drought[!(xQTL_drought$trait %in% trait & (xQTL_drought$xQTL %in% remove.drought)), ])
#     xQTL_normal_rest <- rbind(xQTL_normal_rest, xQTL_normal[!(xQTL_normal$trait %in% trait & xQTL_normal$xQTL %in% remove.normal), ])
}

## metabolite filter ##
xQTL_drought_rest <- xQTL_drought[!(xQTL_drought$trait %in% remove.drought$trait & xQTL_drought$xQTL %in% remove.drought$xQTL), ]
xQTL_normal_rest <- xQTL_normal[!(xQTL_normal$trait %in% remove.normal$trait & xQTL_normal$xQTL %in% remove.normal$xQTL), ]

xQTL_rest <- rbind(xQTL_drought_rest, xQTL_normal_rest)

fl_output <- paste(Prefix, "_final_summary.txt", sep = "")
write.table(xQTL_rest, fl_output, quote = F, sep = "\t", row.names = F)
