\name{glmm}
\title{Generalized Linear Mixed Models}
\usage{
glmm(formula, family=family, data=list(), weights=NULL, offset=NULL,
	nest, delta=1, maxiter=20, points=10, print.level=0,
	control=glm.control(epsilon=0.0001,maxit=10,trace=FALSE))
}
\arguments{
\item{formula}{A symbolic description of the model to be fitted. If it
contains transformations of the data, a dataframe must be supplied.}
\item{family}{A description of the error distribution and link
function to be used in the model; see \code{\link{family}} for details.}
\item{data}{An optional dataframe containing the variables in
the model.}
\item{weights}{An optional weight vector. If this is used, data must
be supplied in a data.frame.}
\item{offset}{The known component in the linear predictor. If this is
used, data must be supplied in a data.frame. An offset cannot be
specified in the model formula.}
\item{nest}{The variable classifying observations by the unit
(cluster) upon which they were observed.}
\item{delta}{If the response variable has been transformed, this is
the Jacobian of that transformation, so that AICs are comparable.}
\item{maxiter}{The maximum number of iterations of the outer loop for
numerical integration.}
\item{points}{The number of points for Gauss-Hermite integration of
the random effect.}
\item{print.level}{If set equal to 2, the log probabilities are
printed out when the underflow error is given.}
\item{control}{A list of parameters for controlling the fitting
process.}
}
\description{
\code{glmm} fits a generalized linear mixed model with a random
intercept using a normal mixing distribution computed by Gauss-Hermite
integration. For the normal, gamma, and inverse Gaussian
distributions, the deviances supplied are -2 log likelihood, not the
usual \code{\link{glm}} deviance; the degrees of freedom take into
account estimation of the dispersion parameter.

If weights and/or offset are to be used or the formula transforms the
some variables, all of the data must be supplied in a dataframe.
Because the \code{\link{glm}} function is such a hack, if this is not
done, weird error messages will result.

na.omit is not allowed.
}
\value{
\code{glmm} returns a list of class \code{glmm}
}
\author{J.K. Lindsey}
\seealso{
\code{\link{family}}, \code{\link{glm}}, \code{\link{glm.control}},
\code{\link{gnlmm}}, \code{\link{gnlr}}, \code{\link{gnlr3}},
\code{\link{fmr}}
}
\examples{
nest <- gl(5,4)
y <- rpois(20,5+2*codes(nest))
# overdispersion model
glmm(y~1, family=poisson, nest=gl(20,1), points=3)
# clustered model
glmm(y~1, family=poisson, nest=nest, points=3)
}
