\name{catmiss}
\title{Marginal probabilities for categorical repeated measures with
missing data}
\usage{
catmiss(response, frequency, ccov=NULL)
}
\arguments{
\item{response}{A matrix with one column for each of the repeated
measures and one row for each possible combination of responses,
including the missing values, indicated by NAs.}
\item{frequency}{A vector containing the frequencies. Its length must
be a multiple of the number of rows of \code{response}. Responses are
arranged in blocks corresponding to the various possible combinations
of values of the explanatory variables.}
\item{ccov}{An optional matrix containing the explanatory variables
(time-constant covariates) as columns, with one line per block of
responses in \code{frequency}. Thus, the number of rows of response
times the number of rows of \code{ccov} equals the length of
\code{frequency}.}
}
\description{
\code{catmiss} calculates the marginal probabilities of
repeated responses. If there are missing values, it gives both the
complete data estimates and the estimates using all data. It is
useful, for example, when a log linear model is fitted; the resulting
fitted values can be supplied to \code{catmiss} to obtain the
estimates of the marginal probabilities for the model. (Note however
that the standard errors do not take into account the fitting of the
model.)
}
\value{
A matrix with the probabilities and their standard errors is returned.
}
\author{J.K. Lindsey}
\seealso{
\code{\link{glm}}, \code{\link{nordr}}
}
\examples{
y <- rpois(27,15)
r1 <- gl(3,1,27)
r2 <- gl(3,3,27)
r3 <- gl(3,9)
# r1, r2, and r3 are factor variables with 3 indicating missing
# independence model with three binary repeated measures
# with missing values
print(z <- glm(y~r1+r2+r3, family=poisson))
# obtain marginal estimates (no observations with 3 missing values)
resp <- cbind(codes(r1), codes(r2), codes(r3))[1:26,]
resp <- ifelse(resp==3, NA, resp)
catmiss(resp, y[1:26])
}

