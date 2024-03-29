% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RRcor.R
\name{RRcor}
\alias{RRcor}
\title{Bivariate correlations including randomized response variables}
\usage{
RRcor(
  x,
  y = NULL,
  models,
  p.list,
  group = NULL,
  bs.n = 0,
  bs.type = c("se.n", "se.p", "pval"),
  nCPU = 1
)
}
\arguments{
\item{x}{a numeric vector, matrix or data frame.}

\item{y}{\code{NULL} (default) or a vector, matrix or data frame with
compatible dimensions to \code{x}.}

\item{models}{a vector defining which RR design is used for each variable.
Must be in the same order as variables appear in \code{x} and \code{y} (by
columns). Available discrete models: \code{Warner}, \code{Kuk}, \code{FR},
\code{Mangat}, \code{UQTknown}, \code{UQTunknown}, \code{Crosswise},
\code{Triangular}, \code{SLD} and \code{direct} (i.e., no randomized
response design). Available continuous models: \code{mix.norm},
\code{mix.exp}.}

\item{p.list}{a \code{list} containing the randomization probabilities of the
RR models defined in \code{models}. Either, all \code{direct}-variables
(i.e., no randomized response) in \code{models} can be excluded in
\code{p.list}; or, if specified, randomization probabilities \code{p} are
ignored for \code{direct}-variables. See \code{\link{RRuni}} for a detailed
specification of p.}

\item{group}{a matrix defining the group membership of each participant
(values 1 and 2) for all multiple group models(\code{SLD},
\code{UQTunknown}). If only one of these models is included in
\code{models}, a vector can be used. For more than one model, each column
should contain one grouping variable}

\item{bs.n}{number of samples used to get bootstrapped standard errors}

\item{bs.type}{to get boostrapped standard errors, use \code{"se.p"} for the
parametric and/or \code{"se.n"} for the nonparametric bootstrap. Use
\code{"pval"} to get p-values from the parametric bootstrap (assuming a
true correlation of zero). Note that \code{bs.n} has to be larger than 0.
The parametric bootstrap is based on the assumption, that the continuous
variable is normally distributed within groups defined by the true state of
the RR variable. For polytomous forced response (FR) designs, the RR
variable is assumed to have equally spaced distances between categories
(i.e., that it is interval scaled)}

\item{nCPU}{only relevant for the bootstrap: either the number of CPU cores
or a cluster initialized via \code{\link[parallel]{makeCluster}}.}
}
\value{
\code{RRcor} returns a list with the following components::

\code{r} estimated correlation matrix

\code{rSE.p}, \code{rSE.n} standard errors from parametric/nonparametric
bootstrap

\code{prob} two-sided p-values from parametric bootstrap

\code{samples.p}, \code{samples.n} sampled correlations from
parametric/nonparametric bootstrap (for the standard errors)
}
\description{
\code{RRcor} calculates bivariate Pearson correlations of variables measured
with or without RR.
}
\details{
Correlations of RR variables are calculated by the method of Fox &
Tracy (1984) by interpreting the variance induced by the RR procedure as
uncorrelated measurement error. Since the error is independent, the
correlation can be corrected to obtain an unbiased estimator.

Note that the continuous RR model \code{mix.norm} with the randomization
parameter \code{p=c(p.truth, mean, SD)} assumes that participants respond
either to the sensitive question with probability \code{p.truth} or otherwise
to a known masking distribution with known mean and SD. The estimated
correlation only depends on the mean and SD and does not require normality.
However, the assumption of normality is used in the parametric bootstrap to
obtain standard errors.
}
\examples{
# generate first RR variable
n <- 1000
p1 <- c(.3, .7)
gData <- RRgen(n, pi = .3, model = "Kuk", p1)

# generate second RR variable
p2 <- c(.8, .5)
t2 <- rbinom(n = n, size = 1, prob = (gData$true + 1) / 2)
temp <- RRgen(model = "UQTknown", p = p2, trueState = t2)
gData$UQTresp <- temp$response
gData$UQTtrue <- temp$true

# generate continuous covariate
gData$cov <- rnorm(n, 0, 4) + gData$UQTtrue + gData$true

# estimate correlations using directly measured / RR variables
cor(gData[, c("true", "cov", "UQTtrue")])
RRcor(
  x = gData[, c("response", "cov", "UQTresp")],
  models = c("Kuk", "d", "UQTknown"), p.list = list(p1, p2)
)
}
\references{
Fox, J. A., & Tracy, P. E. (1984). Measuring associations with
  randomized response. \emph{Social Science Research, 13}, 188-197.
}
\seealso{
\code{vignette('RRreg')} or
  \url{https://www.dwheck.de/vignettes/RRreg.html} for a
  detailed description of the RR models and the appropriate definition of
  \code{p}
}
