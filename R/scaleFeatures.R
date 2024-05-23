train.testScale <- function(y.train, y.test){
  
  # Scale train data
  y.tr.sc <- scale(y.train)
  
  # Make sure train and test data are scaled in the same way
  y.te.sc  <- scale(y.test, 
                    center = attr(y.tr.sc, "scaled:center"),
                    scale  = attr(y.tr.sc, "scaled:scale"))
  
  attr(y.tr.sc, "scaled:center") <- attr(y.te.sc,  "scaled:center") <- NULL
  attr(y.tr.sc, "scaled:scale")  <- attr(y.te.sc,  "scaled:scale")  <- NULL 

  list(y.train.scaled = y.tr.sc, y.test.scaled = y.te.sc)
  
}
