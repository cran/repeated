\name(marg.hom}
\title{Marginal Homogeneity Model}
\usage{
marg.hom(freq,marg1,marg2)
}
\arguments{
\item{freq}{Vector of frequencies}
\item{marg1}{Factor variable for the first margin}
\item{marg2}{Factor variable for the second margin}
}
\description{
\code{marg.hom} fits a marginal homogeneity model to a contingency
table that has two margins of equal size.
}
\value{
A list containing the call, the model, the deviance, the degrees of
freedom, the aic, the fitted values, and the residuals is returned.
}
\author{J.K. Lindsey}
\examples{
# 4x4x2 table in Freq, with margins indexed by Left and Right
Freq <- rpois(32,10)
Left <- gl(4,1,32)
Right <- gl(4,4,32)
marg.hom(Freq, Left, Right)
}
