## Mahatten_plot_extract.R ##
## Mahatten plot for given region (usually in same chromosome) ##
## The input file should in tassel format ##
## Zhang Fei <zhangfei-123@foxmail.com> ##
## 2018-04-03 ##

listdir = list.files(path = "./plot_file/", pattern = ".txt")

for(fl in listdir){
    my_file = read.table(paste("./plot_file/", fl, sep = ""), header = T, sep = "\t")
    trait <- gsub('.txt', '', fl, perl = T)
    range <- max(my_file[ , 4]) - min(my_file[ , 4])
    png(paste('./plot/',trait, ".png", sep = ""), height = 600, 1200)
    plot(c(min(my_file[ , 4]) - 0.2*range, max(my_file[ , 4]) + 0.2*range), c(0, max(-log10(my_file[ , 7])) + 1), type = "n", xaxs = "i", yaxs = "i", xlab = "Location", ylab = "-log10(P)", axes = T, main = trait)
    points(my_file[ , 4], -log10(my_file[ , 7]), col = 'midnightblue', pch = 16, cex = 2)
	abline(h = 5)
    dev.off()
}
