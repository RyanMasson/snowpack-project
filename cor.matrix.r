cor.matrix<-function(x,data=NA){
  # panel.hist function adds the histogram 
  panel.hist <- function(x, ...)
  {
      usr <- par("usr"); on.exit(par(usr))
      par(usr = c(usr[1:2], 0, 1.5) )
      h <- hist(x, plot = FALSE)
      breaks <- h$breaks; nB <- length(breaks)
      y <- h$counts; y <- y/max(y)
      rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
      box()
  }
  panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
  {
      usr <- par("usr"); on.exit(par(usr))
      par(usr = c(0, 1, 0, 1))
      r <- abs(cor(x, y,method="spearman"))
      txt <- format(c(r, 0.123456789), digits=digits)[1]
      txt <- paste(prefix, txt, sep="")
      if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
      text(0.5, 0.5, txt, cex = cex * r)
  }

  if (class(x)=="formula"){
    x<-model.frame(x,data=data)
  }

  pairs(x,lower.panel=panel.smooth,upper.panel=panel.cor,diag.panel=panel.hist,
	cex.labels = 1, font.labels=2)

}
