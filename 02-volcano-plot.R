rm(list=ls())
limma.output <- read.csv(
  './data/limma-output.csv',
  stringsAsFactors = FALSE
)
library(ggplot2)
idx <- abs(limma.output$logFC) <= 975
table(idx)
p_ <- limma.output[idx,]
sig_ <- p_$adj.P.Val < 0.1 
cols <- ifelse(
  sig_,yes = '#c8373764',no = '#80808032'
)
p <- ggplot(
  p_,
  aes(
    x = logFC,
    y = -log10(P.Value)
  )
) +
  geom_point(col=cols)
p

idx <- abs(limma.output$logFC) <= 10
table(idx)
p_ <- limma.output[idx,]
sig_ <- p_$adj.P.Val < 0.1 
cols <- ifelse(
  sig_,yes = '#c8373764',no = '#80808032'
)
p__ <- ggplot(
  p_,
  aes(
    x = logFC,
    y = -log10(P.Value)
  )
) +
  geom_point(col=cols)
p__


png('./figures/02-volcano-major.png',
    width = 1000,
    height = 1000,
    res = 200)
p + ggtitle('Differential expression between \nMale vs Female (Male/Female)')
dev.off()
png('./figures/02-volcano-minor.png',
    width = 1000,
    height = 1000,
    res = 200) 
p__ + ggtitle('Differential expression between \nMale vs Female (Male/Female)')
dev.off()