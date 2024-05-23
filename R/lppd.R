#install.packages("FKSUM")
shhh <- suppressPackageStartupMessages
shhh(library(FKSUM, include.only = "fk_density"))

meanLppd <- function(y.new, y.test){
  
  dims   <- dim(y.test)
  n.test <- dims[1]
  d      <- dims[2]
  
  out <- lapply(1:d, function(j)
    log(FKSUM::fk_density(y.new[,1], x_eval = y.test[,1])$y))
  
  sum(unlist(out)) / n.test
  
}
