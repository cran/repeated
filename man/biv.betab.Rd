\name{biv.betab}
\title{Bivariate Beta-binomial Regression}
\usage{biv.betab(freq, x=NULL, p, depend=T, print.level=0,
typsiz=abs(p), ndigit=10, gradtol=0.00001, stepmax=10*sqrt(p\%*\%p),
steptol=0.00001, iterlim=100, fscale=1)
}
\arguments{
\item{freq}{A matrix containing four columns corresponding to 00, 01,
10, and 11 responses.}
\item{x}{A matrix of explanatory variables, containing pairs of
columns, one for each response, and the same number of rows as freq.}
\item{p}{Initial parameter estimates: intercept, dependence (if depend
is TRUE, and one for each pair of columns of x.}
\item{depend}{If FALSE, the independence (logistic) model is fitted.}
\item{other}{Arguments for nlm.}
}
\description{
\code{biv.betab} fits dependent (logit) linear regression models to a
bivariate beta-binomial distribution.
}
\value{
A list of class \code{bivbetab} is returned.
}
\author{J.K. Lindsey}
\examples{
y <- matrix(  c( 2, 1, 1,13,
		 4, 1, 3, 5,
		 3, 3, 1, 4,
		15, 8, 1, 6),ncol=4,byrow=T)
first <- c(0,0,1,1)
second <- c(0,1,0,1)
self <- cbind(first,second)
other <- cbind(second,first)
biv.betab(y,cbind(self,other),p=c(-1,2,1,1))
# independence
biv.betab(y,cbind(self,other),p=c(-1,1,1),dep=F)
}
