\name{nbkal}
\title{Negative Binomial Model with Kalman Update}
\usage{
nbkal(response, times, mu, preg, pdepend, kalman=TRUE,
	print.level=0, ndigit=10, gradtol=0.00001, steptol=0.00001,
	fscale=1, iterlim=100, typsiz=abs(p), stepmax=10*sqrt(p\%*\%p))
}
\alias{nbkal}
\alias{coefficients.nbkal}
\alias{deviance.nbkal}
\alias{fitted.nbkal}
\alias{residuals.nbkal}
\alias{print.nbkal}
\arguments{
\item{response}{A list of two column matrices with counts and
corresponding times for each individual, one matrix or dataframe of
counts, or an object of class, response (created by
\code{\link{restovec}}) or repeated (created by \code{\link{rmna}}).}
\item{times}{When response is a matrix, a vector of possibly unequally
spaced times when they are the same for all individuals or a matrix of
times. Not necessary if equally spaced. Ignored if response has class,
response or repeated.}
\item{mu}{The mean function.}
\item{preg}{The initial parameter estimates for the mean function.}
\item{pdepend}{The estimates for the dependence parameters, either one or
three.}
\item{kalman}{If TRUE, fits the kalman update model, otherwise, a
standard negative binomial distribution.}
\item{others}{Arguments controlling \code{\link{nlm}}.}
}
\description{
\code{nbkal} fits a negative binomial regression with Kalman update
over time. The variance is proportional to the mean function, whereas,
for \code{\link{kalcount}} with exponential intensity, it is a quadratic
function of the mean.

Marginal and individual profiles can be plotted using
\code{\link{profile}} and \code{\link{iprofile}} and
residuals with \code{\link{plot.residuals}}.
}
\references{
Lambert, P. (1996) Applied Statistics 45, 31-38.

Lambert, P. (1996) Biometrics 52, 50-55.
}
\value{
A list of classes \code{nbkal} and \code{recursive} is returned.
}
\author{J.K. Lindsey}
\seealso{
\code{\link{gar}}, \code{\link{gnlmm}}, \code{\link{gnlr}},
\code{\link{iprofile}} \code{\link{kalcount}}, \code{\link{profile}}
\code{\link{read.list}}, \code{\link{rmna}}, \code{\link{restovec}},
\code{\link{tcctomat}}, \code{\link{tvctomat}}.
}
\examples{
y <- matrix(rnbinom(20,5,0.5), ncol=5)
times <- matrix(rep(seq(10,50,by=10),4), ncol=5, byrow=T)
y0 <- matrix(rep(rnbinom(5,5,0.5),4), ncol=5, byrow=T)
mu <- function(p) p[1]*log(y0)+(times<30)*p[2]*
	(times-30)+(times>30)*p[3]*(times-30)
nbkal(y, preg=c(1.3,0.008,-0.05), times=times, pdep=1.2, mu=mu)
}
