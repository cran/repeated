\name{chidden}
\title{Continuous-Time Hidden Markov Chain Models}
\usage{
chidden(response, totals=NULL, times=NULL, distribution="Bernoulli",
	pgamma, cmu=NULL, tvmu=NULL, pcmu=NULL, ptvmu=NULL,
	pshape=NULL, pfamily=NULL, delta=1, print.level=0,ndigit=10,
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
\item{times}{If \code{response} is a matrix, a vector of corresponding
times, when they are the same for all individuals. Ignored if response
has class, response or repeated.}
\item{distribution}{Bernoulli, Poisson, multinomial, binomial,
exponential, beta binomial, negative binomial, normal, inverse Gauss,
logistic, gamma, Weibull, Cauchy, Laplace, Levy, Pareto, gen(eralized)
gamma, gen(eralized) logistic, Hjorth, Burr, gen(eralized) Weibull,
gen(eralized) extreme value, gen(eralized) inverse Gauss, or power
exponential.}
\item{pgamma}{A square mxm matrix of initial estimates of the
continuous-time hidden Markov transition matrix, where m is the number
of hidden states. Rows can either sum to zero or the diagonal elements
can be zero, in which case they will be replaced by minus the sum of the
other values on the rows. If the matrix contains zeroes off diagonal,
these are fixed and not estimated.}
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
\code{chidden} fits a two or more state hidden Markov chain model to
Bernoulli, binomial, Poisson, or categorical (multinomial) data in
continuous time. All series on different individuals are assumed to
start at the same time point. If the time points are equal, discrete
steps, use \code{\link{hidden}}.

The two mean functions are additive so that interactions between
time-constant and time-varying variables are not possible. Both
functions are on the (generalized) logit scale for the Bernoulli,
binomial, and multinomial distributions and on the log scale for the
Poisson distribution.

The algorithm will run more quickly if the most frequently occurring
time step is scaled to be equal to unity.
}
\value{
A list of class \code{hidden} is returned.
}
\author{J.K. Lindsey}
\seealso{
\code{\link{gar}}, \code{\link{gnlmm}}, \code{\link{hidden}},
\code{\link{kalcount}}, \code{\link{mexp}}, \code{\link{nbkal}},
\code{\link{read.list}}, \code{\link{rmna}}, \code{\link{restovec}}.
}
\examples{
# model for one randomly-generated binary series
y <- c(rbinom(10,1,0.1), rbinom(10,1,0.9))
mu <- function(p) array(p, c(1,2))
print(z <- chidden(y, times=1:20, dist="Bernoulli",
	pgamma=matrix(c(-0.1,0.2,0.1,-0.2),ncol=2),
	cmu=mu, pcmu=c(-2,2)))
mexp(z$gamma)
plot(z)
print(z <- chidden(y, times=(1:20)*2, dist="Bernoulli",
	pgamma=matrix(c(-0.05,0.1,0.05,-0.1),ncol=2),
	cmu=mu, pcmu=c(-2,2)))
mexp(z$gamma) \%*\% mexp(z$gamma)
plot(z)
}

