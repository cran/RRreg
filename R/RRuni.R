#' Univariate analysis of randomized response data
#'
#' Analyse a data vector \code{response} with a specified RR model (e.g.,
#' \code{Warner}) with known randomization probability \code{p}
#'
#' @param response either vector of responses containing 0='no' and 1='yes' or
#'   name of response variable in \code{data}. For the Forced Response
#'   (\code{FR}) model, response values are integers from 0 to (m-1), where 'm'
#'   is the number of response categories. In Kuk's card playing method
#'   (\code{Kuk}), the observed response variable gives the number of red cards.
#' @param data optional \code{data.frame} containing the response variable
#' @param model defines RR model. Available models: \code{"Warner"},
#'   \code{"UQTknown"}, \code{"UQTunknown"}, \code{"Mangat"}, \code{"FR"},
#'   \code{"Kuk"},\code{"Crosswise"}, \code{"Triangular"}, \code{"CDM"},
#'   \code{"CDMsym"}, \code{"SLD"}, \code{"mix.norm"},
#'   \code{"mix.exp"},\code{"mix.unknown"}, or \code{"custom"}. See argument
#'   \code{p} or type \code{vignette('RRreg')} for detailed specifications.
#' @param p randomization probability (see details or \code{vignette("RRreg")})
#' @param group a group vector of the same length as \code{response} containing
#'   values 1 or 2, only required for two-group models, which specify different
#'   randomization probabilities for two groups, e.g., \code{CDM} or \code{SLD}.
#'   If a data.frame \code{data} is provided, the variable \code{group} is
#'   searched within it.
#' @param MLest whether to use \code{optim} to get ML instead of moment
#'   estimates (only relevant if pi is outside of [0,1])
#' @param Kukrep number of repetitions of Kuk's card-drawing method
#'
#' @details Each RR design \code{model} differs in the definition of the
#' randomization probability \code{p}, which is defined as a single probability
#' for
#' \itemize{
#'   \item \code{"Warner"}: Probability to get sensitive Question
#'   \item \code{"Mangat"}: Prob. for non-carriers to respond truthfully (i.e., with No=0)
#'   \item \code{"Crosswise"}: Probability to respond 'yes' to irrelevant second 
#'     question (coding of responses: 1=['no-no' or 'yes-yes']; 0=['yes-no' or 'no-yes'])
#'   \item \code{"Triangular"}: Probability to respond 'yes' to irrelevant second 
#'       question (coding of responses: 0='no' to both questions (='circle'); 1='yes' 
#'       to at least one question ('triangle'))
#' }
#' and as a two-valued vector of probabilities for
#' \itemize{
#'   \item \code{"Kuk"}: Probability of red cards in first and second set, 
#'         respectively (red=1, black=0);
#'   \item Unrelated Question (\code{"UQTknown"}): Prob. to respond to sensitive 
#'         question and known prevalence of 'yes' responses to unrelated question
#'   \item Unrelated Question (\code{"UQTunknown"}): Prob. to respond to
#'         sensitive question in group 1 and 2, respectively
#'   \item Cheating Detection (\code{"CDM"}): Prob. to be prompted to say yes 
#'         in group 1 and 2, respectively
#'   \item Symmetric CDM (\code{"CDMsym"}): 4-valued vector: Prob. to be 
#'         prompted to say 'yes'/'no' in group 1 and 'yes'/'no' in group 2
#'   \item Stochastic Lie Detector (\code{"SLD"}): Prob. for noncarriers to 
#'         reply with 0='no' in group 1 and 2, respectively
#'   \item Forced Response model (\code{"FR"}): m-valued vector (m=number of 
#'         response categories) with the probabilities of being prompted to select 
#'         response categories 0,1,..,m-1, respectively (requires \code{sum(p)<1})
#'   \item RR as misclassification (\code{"custom"}): a quadratic misclassification 
#'         matrix is specified, where the entry \code{p[i,j]} defines the 
#'         probability of responding i (i-th row) given a true state of j 
#'         (j-th column)) (see \code{\link{getPW}})
#' }
#' For the continuous RR models:
#' \itemize{
#'   \item \code{"mix.norm"}: 3-valued vector - Prob. to respond to sensitive 
#'         question and mean and SD of the masking normal distribution of the 
#'         unrelated question
#'   \item \code{"mix.exp"}: 2-valued vector - Prob. to respond to sensitive 
#'         question and mean of the masking exponential distribution of the 
#'         unrelated question
#'   \item \code{"mix.unknown"}: 2-valued vector - Prob. of responding to 
#'         sensitive question in group 1 and 2, respectively
#' }
#'  
#' @return an \code{RRuni} object, can by analyzed by using \code{\link{summary}}
#' 
#' @seealso \code{vignette('RRreg')} or 
#' \url{https://www.dwheck.de/vignettes/RRreg.html} for a 
#' detailed description of the RR models and the appropriate definition of \code{p}
#' 
#' @examples
#' # Generate responses of 1000 people according to Warner's model
#' # with an underlying true proportion of .3
#' df <- RRgen(n = 1000, pi = .3, model = "Warner", p = .7)
#' head(df)
#'
#' # Analyse univariate data to estimate prevalence 'pi'
#' estimate <- RRuni(response = df$response, model = "Warner", p = .7)
#' summary(estimate)
#'
#' # Generate data in line with the Stochastic Lie Detector
#' # assuming that 90% of the respondents answer truthfully
#' df2 <- RRgen(
#'   n = 1000, pi = .3, model = "SLD", p = c(.2, .8),
#'   complyRates = c(.8, 1), groupRatio = 0.4
#' )
#' estimate2 <- RRuni(
#'   response = df2$response, model = "SLD",
#'   p = c(.2, .8), group = df2$group
#' )
#' summary(estimate2)
#' 
#' @export
RRuni <- function(
    response, 
    data, 
    model,
    p, 
    group = NULL,
    MLest = TRUE, 
    Kukrep = 1
) {
  
  # extract column 'response' from data.frame 'data'
  if (!missing(data)) {
    try({
      data <- as.data.frame(data)
      response <- eval(substitute(response), data, parent.frame())
    })
  }
  n <- length(response)
  model <- match.arg(model, modelnames())
  if (is2group(model) && !missing(data)) {
    try(
      {
        data <- as.data.frame(data)
        group <- eval(substitute(group), data, parent.frame())
      },
      silent = TRUE
    )
  }
  if (!is2group(model)) {
    group <- rep(1, n)
  }
  RRcheck.xpgroup(model, response, p, group, "response")

  res <- switch(model,
    "Warner" = RRuni.Warner(response, p),
    "UQTknown" = RRuni.UQTknown(response, p),
    "Mangat" = RRuni.Mangat(response, p),
    "FR" = RRuni.FR(response, p),
    "custom" = RRuni.custom(response, p, group),
    "Kuk" = RRuni.Kuk(response, p, Kukrep),
    "UQTunknown" = RRuni.UQTunknown(response, p, group),
    "Crosswise" = RRuni.Crosswise(response, p),
    "Triangular" = RRuni.Triangular(response, p),
    "SLD" = RRuni.SLD(response, p, group),
    "CDM" = RRuni.CDM(response, p, group),
    "CDMsym" = RRuni.CDMsym(response, p, group),
    "mix.norm" = RRuni.mix.norm(response, p),
    "mix.exp" = RRuni.mix.exp(response, p),
    "mix.unknown" = RRuni.mix.unknown(response, p, group)
  )

  estNotML <- any(res$pi <= 0 | res$pi >= 1)
  if (is2group(model) && model != "mix.unknown") {
    estNotML <- estNotML || switch(model,
      "SLD" = res$t <= 0 | res$t >= 1,
      "UQTunknown" = res$piUQ <= 0 | res$piUQ >= 1,
      "CDM" = res$gamma <= 0 | res$gamma >= 1,
      "CDMsym" = res$gamma <= 0 | res$gamma >= 1
    )
  }
  
  ncat <- ifelse(length(res$pi) <= 2, 2, length(res$pi))
  if (MLest && estNotML && !(model %in% c("mix.norm", "mix.exp", "mix.unknown"))) {
    
    # starting values
    start <- min(max(RRcheck.param(res$pi), 1e-4), 1 - 1e-4)
    # adjust for categorical responses:
    if (length(res$pi) > 1 && model != "custom") {
      start <- start[-1]
    }
    if (sum(start) > 1) {
      start <- start / sum(start)
    }
    if (is2group(model)) {
      start <- c(start, runif(1))
    }

    # optimization
    oo <- optim(
      par = start, fn = RRuni.ll, 
      model = model, response = response, pp = p, 
      group = group, ncat = ncat, Kukrep = Kukrep,
      method = "L-BFGS-B", lower = 1e-8, upper = 1 - 1e-8, 
      control = list(fnscale = -1), hessian = TRUE
    )
    SE <- rep(NA, length(oo$par))
    try(SE <- sqrt(diag(solve(oo$hessian))))

    res$pi <- oo$par[1]
    res$piSE <- SE[1]
    if (model %in% c("custom", "FR")) {
      res$pi <- c(1 - sum(oo$par[1:(ncat - 1)]), oo$par[1:(ncat - 1)])
      res$piSE <- c(NA, SE[1:(ncat - 1)])
    }
    if (model == "UQTunknown") {
      res$piUQ <- oo$par[2]
      res$piUQSE <- SE[2]
    } else if (model %in% c("CDM", "CDMsym")) {
      res$gamma <- 1 - sum(oo$par)
      # res$gammaSE <- NA
      res$beta <- oo$par[2]
      res$betaSE <- SE[2]
    } else if (model == "SLD") {
      res$t <- oo$par[2]
      res$tSE <- SE[2]
    }
  }
  res$logLik <- NA
  try(res$logLik <- RRuni.ll(par = res$pi, model = model, response = response, 
                             pp = p, group = group, ncat = ncat, Kukrep = Kukrep))

  class(res) <- "RRuni"
  return(res)
}


#' @export
print.RRuni <- function(x, ...) {
  cat("RR model: \n")
  write(x$call, "")
  cat("\nEstimate of pi:\n")
  write(paste0(round(x$pi, 6), " (", round(x$piSE, 6), ") "), "")
  if (x$model == "SLD") {
    cat("\nEstimate of t:\n")
    write(paste0(round(x$t, 6), " (", round(x$tSE, 6), ")"), "")
  } else if (x$model %in% c("CDM", "CDMsym")) {
    cat("\nEstimate of gamma:\n")
    write(paste0(round(x$gamma, 6), " (", round(x$gammaSE, 6), ")"), "")
  } else if (x$model == "UQTunknown") {
    cat("\nEstimate for prevalence of unrelated question:\n")
    write(paste(round(x$piUQ, 6), " (", round(x$piUQSE, 6), ")", sep = ""), "")
  } else if (x$model == "mix.unknown") {
    cat("\nEstimate for unrelated question:\n")
    write(paste(round(x$piUQ, 6), " (", round(x$piUQSE, 6), ")", sep = ""), "")
  }
}


#' @export
summary.RRuni <- function(object, ...) {
  zval <- object$pi / object$piSE
  TAB <- cbind(
    Estimate = object$pi,
    StdErr = object$piSE,
    z = zval,
    "Pr(>|z|)" = pnorm(zval, lower.tail = FALSE)
  )
  if (!object$model %in% c("custom", "FR")) {
    rownames(TAB) <- "pi"
  } else {
    rownames(TAB) <- paste("pi", 0:(length(object$pi) - 1), sep = "")
  }
  if (object$model == "UQTunknown") {
    zval_piUQ <- object$piUQ / object$piUQSE
    TAB <- rbind(
      TAB,
      cbind(object$piUQ, object$piUQSE, zval_piUQ, 
            pnorm(zval_piUQ, lower.tail = FALSE))
    )
    rownames(TAB) <- c("pi", "piUQ")
  }
  if (object$model == "SLD") {
    ifelse
    zval_t <- (object$t - 1) / object$tSE ## teste  t gegen 1 !
    TAB <- rbind(
      TAB,
      cbind(object$t, object$tSE, zval_t, pnorm(zval_t))
    )
    rownames(TAB) <- c("pi", "t")
  }
  if (object$model %in% c("CDM", "CDMsym")) {
    zval_b <- object$beta / object$betaSE
    zval_g <- object$gamma / object$gammaSE
    TAB <- rbind(
      TAB,
      cbind(object$beta, object$betaSE, zval_b, 
            pnorm(zval_b, lower.tail = FALSE)),
      cbind(object$gamma, object$gammaSE, zval_g, 
            pnorm(zval_g, lower.tail = FALSE))
    )
    rownames(TAB) <- c("pi", "beta", "gamma")
  }
  #   if (object$model %in% c("CDMsym")){
  #     zval_g=object$gamma/object$gammaSE
  #     TAB=rbind(TAB,
  # #               cbind(object$beta,object$betaSE,zval_b,pnorm(-abs(zval_b))),
  #               cbind(object$gamma,object$gammaSE,zval_g,pnorm(zval_g,lower.tail=F)))
  #     rownames(TAB)=c("pi","gamma")
  #   }
  if (object$model %in% c("mix.unknown")) {
    zval_g <- object$piUQ / object$piUQSE
    TAB <- rbind(
      TAB,
      cbind(object$piUQ, object$piUQSE, zval_g, 
            pnorm(zval_g, lower.tail = FALSE))
    )
    rownames(TAB) <- c("pi", "piUQ")
  }
  res <- list(
    call = object$call, n = object$n,
    coefficients = TAB, model = object$model
  )
  class(res) <- "summary.RRuni"
  res
}


#' @export
print.summary.RRuni <- function(x, ...) {
  # cat("Call:\n")
  write(x$call, "")
  cat("Sample size: ")
  write(x$n, "")
  cat("\n")
  printCoefmat(round(x$coefficients, 6))
  if (x$model == "SLD") {
    cat("\n(for the parameter t, i.e., the probability of true responding of carriers,
 the test is H0: t=1; H1: t<1 and the one-sided probability value is given)\n")
  }
}
