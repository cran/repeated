\name{hidden}
\title{Hidden Markov Chain Models}
\usage{
hidden(response, totals=NULL, distribution="Bernoulli", pgamma,
	cmu=NULL, tvmu=NULL, pcmu=NULL, ptvmu=NULL, pshape=NULL,
	pfamily=NULL, delta=1, fscale=1, print.level=0, ndigit=10,
	gradtol=0.00001, steptol=0.00001, fscale=1, iterlim=100,
	typsiz=abs(p), stepmax=10*sqrt(p\%*\%p))
}
\arguments{
\item{response}{A list of two or three column matrices with counts or
category indicators, times, and possibly totals (if the distribution
is binomial), for each individual, one matrix or dataframe of counts,
or an object of class, response (created by \code{\link{restovec}}) or
repeated (created by \code{\link{rmna}}). If there is only one series,
a vector of responses may be supplied instead.}
\item{totals}{If response is a matrix, a corresponding matrix of
totals if the distribution is binomial. Ignored if response has class,
response or repeated.}
\item{distribution}{Bernoulli, Poisson, multinomial, binomial,
exponential, beta binomial, negative binomial, normal, inverse Gauss,
logistic, gamma, Weibull, Cauchy, Laplace, Levy, Pareto, gen(eralized)
gamma, gen(eralized) logistic, Hjorth, Burr, gen(eralized) Weibull,
gen(eralized) extreme value, gen(eralized) inverse Gauss, or power
exponential.}
\item{pgamma}{A square mxm matrix of initial estimates of the hidden
Markov transition matrix, where m is the number of hidden states.
Rows must sum to one. If the matrix contains zeroes or ones, these are
fixed and not estimated. (Ones cannot appear on the diagonal.)
If a 1x1 matrix or a scalar value of 1 is given, the independence
model is fitted.}
\item{cmu}{A time-constant mean function returning an array with one
row for each individual, one column for each state of the hidden
Markov chain, and, if multinomial, one layer for each category but the
last.}
\item{tvmu}{A time-varying mean function returning an array with one
row for each time point (maximum number of time points for all
individuals if unequal), one column for each state of the hidden
Markov chain, and, if multinomial, one layer for each category but the
last. This is usually a function of time; it is the same for all
individuals.}
\item{pcmu}{Initial estimates of the unknown parameters in
\code{cmu}.}
\item{ptvmu}{Initial estimates of the unknown parameters in
\code{tvmu}.}
\item{pshape}{Initial estimate(s) of the dispersion parameter, for
those distributions having one. This can be one value or a vector with
a different value for each state.}
\item{pfamily}{Initial estimate of the family parameter, for
those distributions having one.}
\item{delta}{Scalar or vector giving the unit of measurement (always
one for discrete data) for each response value, set to unity by
default. For example, if a response is measured to two decimals,
delta=0.01. If the response is transformed, this must be multiplied by
the Jacobian. For example, with a log transformation,
\code{delta=1/response}. Ignored if response has class, response or
repeated.}
\item{others}{Arguments controlling \code{\link{nlm}}.}
}
\description{
\code{hidden} fits a two or more state hidden Markov chain model to
Bernoulli, binomial, Poisson, or categorical (multinomial) data. All
series on different individuals are assumed to start at the same time
point. Time points are equal, discrete steps.

The two mean functions are additive so that interactions between
time-constant and time-varying variables are not possible. Both
functions are on the (generalized) logit scale for the Bernoulli,
binomial, and multinomial distributions and on the log scale for the
Poisson distribution.

See MacDonald, I.L. and Zucchini, W. (1997) Hidden Markov and Other
Models for Discrete-valued Time Series. Chapman and Hall.
}
\value{
A list of class \code{hidden} is returned.
}
\author{J.K. Lindsey}
\seealso{
\code{\link{chidden}}, \code{\link{gar}}, \code{\link{gnlmm}},
\code{\link{kalcount}}, \code{\link{nbkal}}, \code{\link{read.list}},
\code{\link{rmna}}, \code{\link{restovec}}.
}
\examples{
# generate two random Poisson sequences with change-points
y <- rbind(c(rpois(5,1), rpois(15,5)), c(rpois(15,1), rpois(5,5)))
mu <- function(p) array(rep(p[1:2],rep(2,2)), c(2,2))
print(z <- hidden(y,dist="Poisson", cmu=mu, pcmu=c(1,5),
	pgamma=matrix(c(0.9,0.2,0.1,0.8),ncol=2)))
plot(z, nind=1:2)
plot(z, nind=1:2, smooth=T)
}
