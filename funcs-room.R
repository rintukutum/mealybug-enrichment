getGO <- function(x){
  out <- strsplit(
    x,
    #upGO$GO.Annotation[1],
    split = '; '
  )[[1]]
  idxGO <- grep('^GO',out)
  listGO <- list()
  for(i in 1:length(idxGO)){
    begin <- idxGO[i]
    if(i != length(idxGO)){
      end <- idxGO[i+1]-1
    }else{
      end <- length(out)
    }
    OUT <- out[begin:end]
    dfOUT <- data.frame(
      matrix(data = NA,nrow = 1,ncol = 4),
      stringsAsFactors = FALSE
    )
    colnames(dfOUT) <- c('GO.ID','GO.TERM','description','source')
    dfOUT$GO.ID <- OUT[1]
    dfOUT$GO.TERM <- strsplit(OUT[2],split = '\\:')[[1]][1]
    dfOUT$description <- strsplit(OUT[2],split = '\\:')[[1]][2]
    dfOUT$source <- OUT[3]
    listGO[[i]] <- dfOUT
  }
  outDF <- plyr::ldply(
    listGO
  )
  return(outDF)
}
