# Rscript fit.R table.txt

t <- read.table(commandArgs(trailingOnly=TRUE)[1], header=TRUE)
f <- function(y)
{
	fit <- lm(t[[y]]~poly(t$logT, degree=2, raw=TRUE))
	u <- as.numeric(coefficients(fit))
	cat(c(y, " = ", u[1], " + ", u[2], "*x + ", u[3], "*x*x;\n"), sep="")
}
invisible(lapply(c("logJ", "K", "grid", "interp", "h", "H"), f))
