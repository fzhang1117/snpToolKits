## LD_heatmap.R ##
## Zhang Fei ##
## 2018-04-12 ##

workingdir = "."
setwd(workingdir)

library(ggplot2)
library(reshape2)

ld_data <- read.table("./ld_output.txt", header = T, sep = "\t", check.names = F)
#ld_drawing <- ld_data[ , c(1, 2, 7, 8, 14)]
snp.1 <- paste('chr', ld_data[ , 1], '.S_', ld_data[ , 2], sep = "")
snp.2 <- paste('chr', ld_data[ , 7], '.S_', ld_data[ , 8], sep = "")

ld_drawing <- data.frame(snp.1 = snp.1, snp.2 = snp.2, r.square = ld_data[ , 14])
p <- ggplot(data = ld_drawing, aes(x = snp.1, y = snp.2, fill = r.square)) + geom_tile(color = 'white')
p <- p + scale_fill_gradient2(low = 'white', mid = 'grey', high = 'black', midpoint = 0, space = 'Lab') 
#p <- p + scale_y_discrete(breaks = NULL) + theme(panel.grid = element_blank(), panel.background = element_blank()) + labs(x = "", y = "", title = "")
p <- p + theme(panel.grid = element_blank(), panel.background = element_blank()) + labs(x = "", y = "", title = "") + scale_y_discrete(position = 'right') 
p <- p + theme(axis.text.x = element_text(angle = 90))
p <- p + theme(legend.title = element_blank())

ggsave("ld_plot.png", p, width = 6, height = 6, dpi = 1200)


ggplot(mtcars,aes(mpg,cyl))+ geom_point()+ scale_x_continuous(position="top")
