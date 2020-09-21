# Rscript fit.R table.txt

options(digits=5)
t <- read.table(commandArgs(trailingOnly=TRUE)[1], header=TRUE)
f <- function(y)
{
	fit <- lm(t[[y]]~poly(t$logT, degree=3, raw=TRUE))
	u <- as.numeric(coefficients(fit))
	cat("d", y, " = ", sep="")
	cat(u[1])
	cat(" + ")
	cat(u[2])
	cat("*x + ")
	cat(u[3])
	cat("*x2 + ")
	cat(u[4])
	cat("*x3;\n")
}
invisible(lapply(c("logJ", "K", "grid", "interp", "h", "H"), f))
