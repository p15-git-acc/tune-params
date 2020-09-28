# Rscript plotfit.R table.txt

t <- read.table(commandArgs(trailingOnly=TRUE)[1], header=TRUE)
f <- function(y)
{
	x <- "logT";
	plotname = paste("plot", y, ".png", sep="")
	fit <- lm(t[[y]]~poly(t[[x]], degree=3, raw=TRUE))
	png(filename=plotname, width=1280, height=960)
	plot(t[[x]], t[[y]], xlab=x, ylab=y)
	lines(predict(fit)~t[[x]])
	dev.off()
}
invisible(lapply(c("logJ", "K", "grid", "interp", "h", "H"), f))
