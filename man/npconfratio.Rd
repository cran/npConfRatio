\name{npconfratio}
\alias{npconfratio}
\title{ Nonparametric Confidence Intervals for the Ratio of Medians }
\description{
  For two sample problems, this function computes confidence intervals
  for the ratio of the medians \code{med(y)/med(x)}. Exact confidence sets
  can be computed for untied as well as for tied samples up to sample
  sizes of 50 each. Additionally, a Hodges-Lehmann point estimator is
  returned for untied samples only. An nonparametric equivalence test can be
  performed. 
}
\usage{
npconfratio(x, y, conf.level=0.95, alternative=c("two.sided", "less", 
"greater"), equiv.test=TRUE, equiv.range=c(0.8, 1.25), pest=TRUE)
}
\arguments{
  \item{x}{numeric vector of data values for the \code{x} sample. }
  \item{y}{numeric vector of data values for the \code{x} sample. }
  \item{conf.level}{confidence level of the interval.}
  \item{alternative}{the alternative hypothesis must be
    one of \code{"two.sided"} (default), \code{"greater"} or
    \code{"less"}.}
  \item{equiv.test}{a logical indicating weather the corresponding 
     equivalence test should be performed.}
  \item{equiv.range}{equivalence range for the equivalence test.}
  \item{pest}{a logical indicating weather a point estimator should be
computed.}
}
\details{
  In some situations it is more appropriate to make statistical inference
  about a relation between two samples instead of the usual
  difference. This procedure computes nonparametric confidence intervals
  for the ratio of the medians \eqn{\varrho = med(y)/med(x)} using the 
  Wilcoxon-Mann-Whitney statistic. Exact confidence intervals are returned
  if both samples contain less than 50 values each, otherwise a normal
  approximation is used. For tied samples the Streitberg \& R\"ohmel shift
  algorithm is used to track the conditional distribution of the
  WMW-statistic. 

  The corresponding equivalence test of the null-hypothesis \eqn{H_0:
  \varrho \not\in [u,l]} against the alternative \eqn{H_1: \varrho \in
  [u,l]} is performed if \code{equiv.test} is \code{TRUE} (\code{equiv.range}
  specifies the interval \eqn{[u,l]}). The test rejects at level
  \eqn{\alpha} if and only if the \eqn{2\alpha} confidence set in included 
  in the equivalence range. 

}
\value{
A list with class \code{"htest"} containing the following components:
  \item{p.value}{the p-value for the test.}
  \item{alternative}{a character string describing the alternative
    hypothesis.}
  \item{method}{the type of procedure applied.}
  \item{data.name}{a character string giving the names of the data.}
  \item{conf.int}{a confidence interval for the ratio of medians.}
  \item{estimate}{a point estimator}

}
\references{
  Torsten Hothorn & Ullrich Munzel (2000),
  \emph{An Exact Nonparametric Confidence Interval for the Ratio of Medians}
  (to appear)

  Bernd Streitberg & Joachim R\"ohmel (1986),
  \emph{Exact Distributions For Permutations and Rank Tests:  
   An Introduction to Some Recently Published Algorithms}
  Statistical Software Letter, 12 No. 1, 10-17.

}
\author{ Torsten Hothorn <Torsten.Hothorn@rzmail.uni-erlangen.de> }

\seealso{ for confidence sets for a shift parameter see 
  \code{\link{wilcox.test}}}

\examples{

# The ASAT values (a liver enzyme, in a control and treatment group)

contr <- c(1.33, 1.78, 1.53, 1.95, 1.83, 1.47, 1.87, 1.55, 2.58, 2.17, 1.97,
1.62, 2.25, 3.53, 2.92, 1.78, 2.22, 2.55, 2.75)
treat <- c(1.53, 1.75, 2.12, 2.83, 2.58, 2.37, 2.92, 2.08, 3.01, 2.67, 2.08,
2.25, 3.08, 2.01, 1.58)

boxplot(contr, treat, range=1.5, names=c("Control", "Compound"),
ylab=expression(paste("ASAT values in ", mu , mol / L)))

npconfratio(contr, treat, alternative="less", pest=F)

# The NEUT values
g1 <- c(16.2, 16.6, 21.2, 21.4, 21.9)
g2 <- c(16.1, 21.9, 23.1, 23.5, 23.8)
g3 <- c(17.5, 20.2, 21.9, 23.1, 24.1)
g4 <- c(19.4, 23.1, 24.1, 30.1, 31)

npconfratio(g1,g2)
npconfratio(g1,g3)
npconfratio(g1,g4)

#  some random samples

x <- rexp(10, 0.25)
y <- rexp(10, 0.5)
rho <- qexp(0.5, 0.5)/qexp(0.5, 0.25)
cat("true parameter: ", rho, "\n")
npconfratio(x,y, equiv.test=F, pest=F)$conf.int

# perform a equivalence test: H_0: rho < 0.8 | rho > 1.25 vs. H_1: 0.8 < rho < 1.25
x <- rexp(100)
y <- rexp(100)
npconfratio(x,y)

}
\keyword{ median, two-sample, ratio }%-- one or more ...
