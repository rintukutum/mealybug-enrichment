rm(list=ls())
dir.create('./data',showWarnings = FALSE)
downGO <- read.csv(
  './data/downregulated-ortho-GO.csv',
  stringsAsFactors = FALSE
)
######### remove NA and nothing
idx.rm <- grep('N/A',downGO$GO.Annotation)
downGO.ok <- downGO[-idx.rm,]
idx.rm <- which(nchar(downGO.ok$GO.Annotation) == 0)
downGO.ok <- downGO.ok[-idx.rm,]
########
source('./funcs-room.R')
downGO.forward <- list()
for(i in 1:nrow(downGO.ok)){
  downGO.forward[[i]] <- getGO(x = downGO.ok$GO.Annotation[i])
}
names(downGO.forward) <- downGO.ok$ids
downGO.final <- plyr::ldply(
  downGO.forward
)
####
tb_ <- table(downGO.final$GO.TERM)
names(tb_) <- c(
  'Cellular Component',
  'Molecular Function',
  'Biological Process'
)
org.mar <- par()$mar
mod.mar <- c(5.1,12.1,4.1,2.1)
dir.create('./figures',showWarnings = FALSE)
png('./figures/01-downregulated-GO-terms.png',
    width = 1300,
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
  main = 'DOWN-regulated genes statistics on GO terms'
)
text(y = bp[,1],x = 3,labels = tb_,adj = 0,cex = 2,
     col = rev(c('white','white','black'))
)
dev.off()
par(mar = org.mar)
######################
tb_ <- data.frame(table(downGO.final[,c('description','GO.TERM')]))
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
    low = '#ffaaaaff',
    high = '#d40000ff'
  ) +
  scale_x_discrete(
    labels = c(
      'Cellular Component',
      'Molecular Function',
      'Biological Process'
    )
  ) +
  ggtitle('Down-regulated\nGO processes') +
  theme(axis.text.x = element_text(angle=45,hjust = 1))
png('./figures/01-downregulated-GO-processes.png',
    width = 800,
    height = 1400,
    res = 200)
p
dev.off()

write.csv(
  downGO.final,
  file = './data/down-GO-final.csv',
  row.names = FALSE
)