\name{biv.binom}
\title{Marginal Bivariate Binomial Regression}
\usage{
biv.binom(freq, marg1=~1, marg2=~1, interaction=~1, pmarg1=1,
pmarg2=1,pinteraction=1, print.level=0, typsiz=abs(p),
ndigit=10, gradtol=0.00001, stepmax=10*sqrt(p\%*\%p),
steptol=0.00001, iterlim=100, fscale=1)
}
\arguments{
\item{freq}{A four-column matrix containing K 2x2 frequency tables.}
\item{marg1}{The model formula for the first margin.}
\item{marg2}{The model formula for the second margin.}
\item{interaction}{The model formula for the interaction.}
\item{pmarg1}{Initial parameter estimates for the first margin regression.}
\item{pmarg2}{Initial parameter estimates for the second margin regression.}
\item{pinteraction}{Initial parameter estimates for the interaction regression.}
\item{other}{Arguments for nlm.}
}
\description{
\code{biv.binom} fits (logit) linear regression models to a marginal
bivariate binomial distribution. The covariates must be of length K,
that is the number of 2x2 tables.
}
\value{
A list of class \code{bivbinom} is returned.
}
\author{J.K. Lindsey}}
\examples{
# 5 2x2 tables
Freq <- matrix(rpois(20,10),ncol=4)
x <- c(6,8,10,12,14)
print(z <- biv.binom(Freq,marg1=~x,marg2=~x,inter=~x,pmarg1=c(-2,0.08),
	pmarg2=c(-2,0.1),pinter=c(3,0)))
}
