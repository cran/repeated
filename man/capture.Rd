\name{capture}
\title{Capture-recapture Model}
\usage{
n <- periods # number of periods
eval(setup)
z <- glm(y~...,family=poisson,weights=pw)
capture(z,n)
}
\arguments{
\item{z}{A Poisson generalized linear model}
\item{n}{The number of repeated observations.}
}
\description{
\code{capture} fits the Cormack capture-recapture model to n sample
periods. Set n to the appropriate value and type eval(setup). This
produces the following variables -
p[i]: logit capture probabilities,
pbd: constant capture probability,
d[i]: death parameters,
b[i]: birth parameters,
pw: prior weights.
Then set up a Poisson model for log linear models and call the
function, \code{capture}.

If there is constant effort, then all estimates are correct.
Otherwise, n[1], p[1], b[1], are correct only if there is no birth in
period 1.  n[s], p[s], are correct only if there is no death in the
last period.  phi[s-1] is correct only if effort is constant in (s-1,
s).  b[s-1] is correct only if n[s] and phi[s-1] both are.
}
\value{
\code{capture} returns a matrix containing the estimates.
}
\author{J.K. Lindsey}
\examples{
y <- c(0,1,0,0,0,1,0,1,0,0,0,1,0,0,0,14,1,1,0,2,1,2,1,16,0,2,0,11,
	2,13,10,0)
n <- 5
eval(setup)
# closed population
print(z0 <- glm(y~p1+p2+p3+p4+p5, family=poisson, weights=pw))
# deaths and emigration only
print(z1 <- update(z0, .~.+d1+d2+d3))
# immigration only
print(z2 <- update(z1, .~.-d1-d2-d3+b2+b3+b4))
# deaths, emigration, and immigration
print(z3 <- update(z2, .~.+d1+d2+d3))
# add trap dependence
print(z4 <- update(z3, .~.+i2+i3))
# constant capture probability over the three middle periods
print(z5 <- glm(y~p1+pbd+p5+d1+d2+d3+b2+b3+b4, family=poisson, weights=pw))
# print out estimates
capture(z5, n)
}