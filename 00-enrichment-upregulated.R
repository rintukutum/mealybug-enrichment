rm(list=ls())
dir.create('./data',showWarnings = FALSE)
upGO <- read.csv(
  './data/upregulated-ortho-GO.csv',
  stringsAsFactors = FALSE
)
######### remove NA and nothing
idx.rm <- grep('N/A',upGO$GO.Annotation)
upGO.ok <- upGO[-idx.rm,]
idx.rm <- which(nchar(upGO.ok$GO.Annotation) == 0)
upGO.ok <- upGO.ok[-idx.rm,]
########
source('./funcs-room.R')
upGO.forward <- list()
for(i in 1:nrow(upGO.ok)){
  upGO.forward[[i]] <- getGO(x = upGO.ok$GO.Annotation[i])
}
names(upGO.forward) <- upGO.ok$ids
upGO.final <- plyr::ldply(
  upGO.forward
)
####
tb_ <- table(upGO.final$GO.TERM)
names(tb_) <- c(
  'Cellular Component',
  'Molecular Function',
  'Biological Process'
)
org.mar <- par()$mar
mod.mar <- c(5.1,12.1,4.1,2.1)
dir.create('./figures',showWarnings = FALSE)
png('./figures/00-upregulated-GO-terms.png',
    width = 1250,
    height = 800,
    res = 200)
par(mar = mod.mar)
bp <- barplot(
  tb_,
  horiz = TRUE,
  las = 2,
  col = rev(c('#a05a2cff',
          '#d38d5fff',
          '#e9c6afff')),
  main = 'UP-regulated genes statistics on GO terms'
)
text(y = bp[,1],x = 3,labels = tb_,adj = 0,cex = 2,
     col = rev(c('white','white','black'))
    )
dev.off()
par(mar = org.mar)
#################
tb_ <- data.frame(table(upGO.final[,c('description','GO.TERM')]))
idx.rm <- which(tb_$Freq == 0)
tb_ <- tb_[-idx.rm,]

df_ <- tb_[tb_$Freq >= 3,]
library(ggplot2)
p <- ggplot(
  df_,
  aes(x = GO.TERM, y = description)
) +
  geom_tile(aes(fill=Freq)) + 
  scale_fill_gradient(
  low = '#afe9c6ff',
  high = '#217844ff'
) +
  scale_x_discrete(
    labels = c(
      'Cellular Component',
      'Molecular Function',
      'Biological Process'
    )
  ) +
  ggtitle('Up-regulated\nGO processes') +
  theme(axis.text.x = element_text(angle=45,hjust = 1))
png('./figures/00-upregulated-GO-processes.png',
    width = 1200,
    height = 1400,
    res = 200)
p
dev.off()

write.csv(
  upGO.final,
  file = './data/up-GO-final.csv',
  row.names = FALSE
)