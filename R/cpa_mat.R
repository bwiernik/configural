#' Conduct criterion profile analysis using a correlation matrix
#'
#' @param formula Regression formula with a single outcome variable on the left-hand side and one or more predictor variables on the right-hand side (e.g., Y ~ X1 + X2).
#' @param cov_mat Correlation matrix containing the variables to be used in the regression.
#' @param n Sample size. Used to compute adjusted R-squared and, if `se_var_mat` is NULL, standard errors. If NULL and `se_var_mat` is specified, effective sample size is computed based on `se_var_mat` (cf. Revelle et al., 2017).
#' @param se_var_mat Optional. The sampling error covariance matrix among the unique elements of \code{cov_mat}. Used to calculate standard errors. If not supplied, the sampling covariance matrix is calculated using \code{n}.
#' @param se_beta_method Method to use to estimate the standard errors of standardized regression (beta) coefficients. Current options include "normal" (use the Jones-Waller, 2015, normal-theory approach) and "lm" (estimate standard errors using conventional regression formulas).
#' @param adjust Method to adjust R-squared for overfitting. See \code{\link{adjust_Rsq}} for details.
#' @param conf_level Confidence level to use for confidence intervals.
#' @param ... Additional arguments.
#'
#' @return An object of class "cpa" containing the criterion pattern vector and CPA variance decomposition
#' @export
#'
#' @encoding UTF-8
#'
#' @references
#' Jones, J. A., & Waller, N. G. (2015).
#' The normal-theory and asymptotic distribution-free (ADF) covariance matrix of standardized regression coefficients: Theoretical extensions and finite sample behavior.
#' _Psychometrika, 80_(2), 365–378. \url{https://doi.org/10/gckfx5}
#'
#' Revelle, W., Condon, D. M., Wilt, J., French, J. A., Brown, A., & Elleman, L. G. (2017).
#' Web- and phone-based data collection using planned missing designs.
#' In N. G. Fielding, R. M. Lee, & G. Blank, _The SAGE Handbook of Online Research Methods_ (pp. 578–594).
#' SAGE Publications. https://doi.org/10.4135/9781473957992.n33
#'
#' Wiernik, B. M., Wilmot, M. P., Davison, M. L., & Ones, D. S. (2019).
#' _Meta-analytic criterion profile analysis._
#' Manuscript submitted for publication.
#'
#' @examples
#' sevar <- cor_covariance_meta(mindfulness$r, mindfulness$n, mindfulness$sevar_r, mindfulness$source)
#' cpa_mat(mindfulness ~ ES + A + C + Ex + O,
#'           cov_mat = mindfulness$rho,
#'           n = NULL,
#'           se_var_mat = sevar,
#'           adjust = "pop")
cpa_mat <- function(formula, cov_mat, n = NULL,
                    se_var_mat = NULL,
                    se_beta_method = c("normal", "lm"),
                    adjust = c("fisher", "pop", "cv"),
                    conf_level = .95,
                    ...) {

  cov_mat <- as.matrix(cov_mat)
  formula <- stats::as.formula(formula)
  y_col <- as.character(formula[[2]])
  y_col <- unique(y_col[y_col %in% rownames(cov_mat)])
  x_col <- as.character(formula)[[3]]
  x_col <- strsplit(x_col, split = "\\s*[+]\\s*")[[1]]
  x_col <- .remove_charmargins(x = x_col)
  x_col <- unique(x_col[x_col %in% rownames(cov_mat)])

  R <- stats::cov2cor(cov_mat)
  p <- length(x_col)
  Rxx <- R[x_col, x_col]
  rxy <- R[x_col, y_col]
  Rinv <- solve(Rxx)
  Id <- diag(p)
  one <- rep(1, p)
  Q <- Id - one %*% t(one) / p

  beta <- Rinv %*% rxy
  R2_total <- t(rxy) %*% beta
  R_total <- sqrt(R2_total)

  bstar <- Q %*% beta
  lev <- (t(rxy) %*% one) / sqrt(t(one) %*% Rxx %*% one)
  lev2 <- lev^2
  pat <- (t(rxy) %*% bstar) / sqrt(t(bstar) %*% Rxx %*% bstar)
  pat2 <- pat^2
  levpat <- (t(one) %*% Rxx %*% bstar) / sqrt(t(bstar) %*% Rxx %*% bstar %*% t(one) %*% Rxx %*% one)
  levpat2 <- levpat^2

  Rxx_cpa <- matrix(c(1, levpat, levpat, 1), nrow = 2)
  beta_cpa <- solve(Rxx_cpa) %*% c(lev, pat)

  if ((is.null(n) & is.null(se_var_mat))) n <- NA

  se_var <- var_error_cpa(Rxx = Rxx, rxy = rxy, n = n, se_var_mat = se_var_mat, adjust = adjust)

  if (is.null(n)) {
    n <- n_effective_R2(R2_total, se_var$r.total.squared, p)
    n_effective <- n
    n_null <- TRUE
  } else {
    n_effective <- NULL
    n_null <- FALSE
  }

  if (is.na(n) | is.infinite(n)) R2_adj <- R2_total else R2_adj <- adjust_Rsq(R2_total, n, p, adjust)
  R_adj <- if (R2_adj < 0) 0 else sqrt(R2_adj)

  pat_adj <- lev * levpat + sqrt(max((1 - levpat2) * (R2_adj - lev2), 0))
  pat2_adj <- pat_adj^2
  beta_cpa_adj <- solve(Rxx_cpa) %*% c(lev, pat_adj)

  delta_lev <- R2_total - pat2
  delta_pat <- R2_total - lev2
  delta_lev_adj <- R2_adj - pat2_adj
  delta_pat_adj <- R2_adj - lev2

  moe <- stats::qnorm((1 - conf_level) / 2, lower.tail = FALSE)
  if (is.infinite(n) | is.na(n) | n_null) {
    moe_beta <- stats::qnorm((1 - conf_level) / 2, lower.tail = FALSE)
  } else {
    moe_beta <- stats::qt((1 - conf_level) / 2, n - p - 1, lower.tail = FALSE)
  }

  se_beta <- sqrt(diag(se_var$beta))
  se_bstar <- sqrt(diag(se_var$bstar))

  beta_mat <- cbind(beta, se_beta, beta - se_beta * moe_beta, beta + se_beta * moe_beta)
  rownames(beta_mat) <- x_col
  colnames(beta_mat) <- c("beta",
                          "SE",
                          paste0((1 - conf_level) / 2 * 100, "%"),
                          paste0((1 - (1 - conf_level) / 2) * 100, "%")
  )

  bstar_mat <- cbind(bstar, se_bstar, bstar - se_bstar * moe_beta, bstar + se_bstar * moe_beta)
  rownames(bstar_mat) <- x_col
  colnames(bstar_mat) <- c("bstar",
                           "SE",
                           paste0((1 - conf_level) / 2 * 100, "%"),
                           paste0((1 - (1 - conf_level) / 2) * 100, "%")
  )

  se_R <- sqrt(unlist(se_var[c("r.level", "r.pattern", "r.total")]))
  se_R2 <- sqrt(unlist(se_var[c("r.level.squared", "r.pattern.squared", "r.total.squared")]))
  se_beta_cpa <- c(sqrt(diag(se_var$beta.cpa)), NA)
  se_delta <- c(sqrt(unlist(se_var[c("delta.r.squared.level", "delta.r.squared.pattern")])), NA)
  summary_mat <- cbind(c(lev, pat, R_total),
                       se_R,
                       c(lev, pat, R_total) - se_R * moe,
                       c(lev, pat, R_total) + se_R * moe,
                       c(lev2, pat2, R2_total),
                       se_R2,
                       c(lev2, pat2, R2_total) - se_R2 * moe,
                       c(lev2, pat2, R2_total) + se_R2 * moe,
                       c(delta_lev, delta_pat, NA),
                       se_delta,
                       c(delta_lev, delta_pat, NA) - se_delta * moe,
                       c(delta_lev, delta_pat, NA) + se_delta * moe,
                       c(beta_cpa, NA)
  )
  rownames(summary_mat) <- c("Level", "Pattern", "Total")
  colnames(summary_mat) <- c("R", "SE",
                             paste0((1 - conf_level) / 2 * 100, "%"),
                             paste0((1 - (1 - conf_level) / 2) * 100, "%"),
                             "R-squared", "SE",
                             paste0((1 - conf_level) / 2 * 100, "%"),
                             paste0((1 - (1 - conf_level) / 2) * 100, "%"),
                             "Delta R-squared", "SE",
                             paste0((1 - conf_level) / 2 * 100, "%"),
                             paste0((1 - (1 - conf_level) / 2) * 100, "%"),
                             "beta"
                             )

  se_R_adj <- sqrt(unlist(se_var[c("r.level", "adjusted.r.pattern", "adjusted.r.total")]))
  se_R2_adj <- sqrt(unlist(se_var[c("r.level.squared", "adjusted.r.pattern.squared", "adjusted.r.total.squared")]))
  se_beta_cpa_adj <- c(sqrt(diag(se_var$adjusted.beta.cpa)), NA)
  se_delta_adj <- c(sqrt(unlist(se_var[c("adjusted.delta.r.squared.level", "adjusted.delta.r.squared.pattern")])), NA)
  adj_summary_mat <- cbind(c(lev, pat_adj, R_adj),
                           se_R_adj,
                           c(lev, pat_adj, R_adj) - se_R_adj * moe,
                           c(lev, pat_adj, R_adj) + se_R_adj * moe,
                           c(lev2, pat2_adj, R2_adj),
                           se_R2_adj,
                           c(lev2, pat2_adj, R2_adj) - se_R2_adj * moe,
                           c(lev2, pat2_adj, R2_adj) + se_R2_adj * moe,
                           c(delta_lev_adj, delta_pat_adj, NA),
                           se_delta_adj,
                           c(delta_lev_adj, delta_pat_adj, NA) - se_delta_adj * moe,
                           c(delta_lev_adj, delta_pat_adj, NA) + se_delta_adj * moe,
                           c(beta_cpa_adj, NA)
  )
  rownames(adj_summary_mat) <- c("Level", "Pattern", "Total")
  colnames(adj_summary_mat) <- c("Adj. R", "SE",
                                 paste0((1 - conf_level) / 2 * 100, "%"),
                                 paste0((1 - (1 - conf_level) / 2) * 100, "%"),
                                 "Adj. R-squared", "SE",
                                 paste0((1 - conf_level) / 2 * 100, "%"),
                                 paste0((1 - (1 - conf_level) / 2) * 100, "%"),
                                 "Adj. Delta R-squared", "SE",
                                 paste0((1 - conf_level) / 2 * 100, "%"),
                                 paste0((1 - (1 - conf_level) / 2) * 100, "%"),
                                 "Adj. beta"
                                 )

  se_levpat <- sqrt(se_var$r.level.pattern)
  levpat_mat <- data.frame(levpat, se_levpat, levpat - se_levpat * moe, levpat + se_levpat * moe)
  colnames(levpat_mat) <- c("R", "SE", paste0((1 - conf_level) / 2 * 100, "%"), paste0((1 - (1 - conf_level) / 2) * 100, "%"))
  rownames(levpat_mat) <- ""

  out <- list(
    beta = beta_mat,
    bstar = bstar_mat,
    fit = summary_mat,
    adjusted_fit = adj_summary_mat,
    r.level.pattern = levpat_mat,
    df.residual = n - p - 1,
    call = match.call(),
    n = n,
    rank = p + 1,
    model = list(R = R, Rxx = Rxx, rxy = rxy),
    vcov = se_var,
    n_eff = n_null
  )

  class(out) <- "cpa"

  return(out)

}

#' @export
#' @keywords internal
#' @exportClass cpa
#' @method print cpa
print.cpa <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {

  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")

  cat("\nCriterion Pattern:\n")
  print(x$bstar, digits = digits)

  cat("\n\nVariance Decomposition:\n")
  print(x$fit, digits = digits)

  cat("\n\nAdjusted Variance Decomposition:\n")
  print(x$adjusted_fit, digits = digits)

  cat("\n\nCorrelation Between Profile Level and Criterion Pattern Similarity:\n")
  print(x$r.level.pattern, digits = digits)


  if (x$n_eff) {
    cat("\nEffective Sample Size:\n")
  } else cat("\nSample Size:\n")
  cat(round(x$n, digits = digits))

  cat("\n\nDegrees of freedom:\n  Level: 1 and", x$n - 2,
      "\n  Pattern:", x$rank - 2, "and", x$n - x$rank + 1,
      "\n  Total:", x$rank - 1, "and", x$n - x$rank)

  invisible(x)
}

#' @export
#' @keywords internal
#' @exportClass cpa
#' @method summary cpa
summary.cpa <- function(object, ...) {
  object
}

#' @export
#' @keywords internal
#' @exportClass cpa
#' @method vcov cpa
vcov.cpa <- function(object, parameter = NULL, ...) {
  if (is.null(parameter)) {
    object$vcov
  } else if (length(parameter) == 1) {
    object$vcov[[parameter]]
  } else {
    object$vcov[parameter]
  }
}

#' @export
#' @keywords internal
#' @exportClass cpa
#' @method coef cpa
coef.cpa <- function(object, parameter = "bstar", ...) {
  get(parameter, object)
}

#' @export
#' @keywords internal
#' @exportClass cpa
#' @method confint cpa
confint.cpa <- function(object, ...) {
  object$bstar[,3:4]
}

#' @export
#' @keywords internal
#' @exportClass cpa
#' @method predict cpa
predict.cpa <- predict(object, newdata, se.fit = FALSE, scale = NULL, df = Inf,
                       interval = c("none", "confidence", "prediction"),
                       level = 0.95, type = c("response", "terms", "profile"),
                       terms = NULL, na.action = na.pass,
                       pred.var = res.var/weights, weights = 1, ...) {

  type <- match.arg(type)

  tt <- terms(object)
  if (!inherits(object, "cpa"))
    warning("calling predict.cpa(<fake-cpa-object>) ...")
  if (missing(newdata) || is.null(newdata)) {
    mm <- X <- model.matrix(object)
    mmDone <- TRUE
    offset <- object$offset
  }
  else {
    Terms <- delete.response(tt)
    m <- model.frame(Terms, newdata, na.action = na.action,
                     xlev = object$xlevels)
    if (!is.null(cl <- attr(Terms, "dataClasses")))
      .checkMFClasses(cl, m)
    X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
    offset <- rep(0, nrow(X))
    if (!is.null(off.num <- attr(tt, "offset")))
      for (i in off.num) offset <- offset + eval(attr(tt,
                                                      "variables")[[i + 1]], newdata)
    if (!is.null(object$call$offset))
      offset <- offset + eval(object$call$offset, newdata)
    mmDone <- FALSE
  }
  n <- length(object$residuals)
  p <- object$rank
  p1 <- seq_len(p)
  piv <- if (p)
    qr.lm(object)$pivot[p1]
  if (p < ncol(X) && !(missing(newdata) || is.null(newdata)))
    warning("prediction from a rank-deficient fit may be misleading")
  beta <- object$coefficients
  predictor <- drop(X[, piv, drop = FALSE] %*% beta[piv])
  if (!is.null(offset))
    predictor <- predictor + offset
  interval <- match.arg(interval)
  if (interval == "prediction") {
    if (missing(newdata))
      warning("predictions on current data refer to _future_ responses\n")
    if (missing(newdata) && missing(weights)) {
      w <- weights.default(object)
      if (!is.null(w)) {
        weights <- w
        warning("assuming prediction variance inversely proportional to weights used for fitting\n")
      }
    }
    if (!missing(newdata) && missing(weights) && !is.null(object$weights) &&
        missing(pred.var))
      warning("Assuming constant prediction variance even though model fit is weighted\n")
    if (inherits(weights, "formula")) {
      if (length(weights) != 2L)
        stop("'weights' as formula should be one-sided")
      d <- if (missing(newdata) || is.null(newdata))
        model.frame(object)
      else newdata
      weights <- eval(weights[[2L]], d, environment(weights))
    }
  }
  type <- match.arg(type)
  if (se.fit || interval != "none") {
    w <- object$weights
    res.var <- if (is.null(scale)) {
      r <- object$residuals
      rss <- sum(if (is.null(w)) r^2 else r^2 * w)
      df <- object$df.residual
      rss/df
    }
    else scale^2
    if (type != "terms") {
      if (p > 0) {
        XRinv <- if (missing(newdata) && is.null(w))
          qr.Q(qr.lm(object))[, p1, drop = FALSE]
        else X[, piv] %*% qr.solve(qr.R(qr.lm(object))[p1,
                                                       p1])
        ip <- drop(XRinv^2 %*% rep(res.var, p))
      }
      else ip <- rep(0, n)
    }
  }
  if (type == "terms") {
    if (!mmDone) {
      mm <- model.matrix(object)
      mmDone <- TRUE
    }
    aa <- attr(mm, "assign")
    ll <- attr(tt, "term.labels")
    hasintercept <- attr(tt, "intercept") > 0L
    if (hasintercept)
      ll <- c("(Intercept)", ll)
    aaa <- factor(aa, labels = ll)
    asgn <- split(order(aa), aaa)
    if (hasintercept) {
      asgn$"(Intercept)" <- NULL
      avx <- colMeans(mm)
      termsconst <- sum(avx[piv] * beta[piv])
    }
    nterms <- length(asgn)
    if (nterms > 0) {
      predictor <- matrix(ncol = nterms, nrow = NROW(X))
      dimnames(predictor) <- list(rownames(X), names(asgn))
      if (se.fit || interval != "none") {
        ip <- matrix(ncol = nterms, nrow = NROW(X))
        dimnames(ip) <- list(rownames(X), names(asgn))
        Rinv <- qr.solve(qr.R(qr.lm(object))[p1, p1])
      }
      if (hasintercept)
        X <- sweep(X, 2L, avx, check.margin = FALSE)
      unpiv <- rep.int(0L, NCOL(X))
      unpiv[piv] <- p1
      for (i in seq.int(1L, nterms, length.out = nterms)) {
        iipiv <- asgn[[i]]
        ii <- unpiv[iipiv]
        iipiv[ii == 0L] <- 0L
        predictor[, i] <- if (any(iipiv > 0L))
          X[, iipiv, drop = FALSE] %*% beta[iipiv]
        else 0
        if (se.fit || interval != "none")
          ip[, i] <- if (any(iipiv > 0L))
            as.matrix(X[, iipiv, drop = FALSE] %*% Rinv[ii,
                                                        , drop = FALSE])^2 %*% rep.int(res.var,
                                                                                       p)
        else 0
      }
      if (!is.null(terms)) {
        predictor <- predictor[, terms, drop = FALSE]
        if (se.fit)
          ip <- ip[, terms, drop = FALSE]
      }
    }
    else {
      predictor <- ip <- matrix(0, n, 0L)
    }
    attr(predictor, "constant") <- if (hasintercept)
      termsconst
    else 0
  }
  if (interval != "none") {
    tfrac <- qt((1 - level)/2, df)
    hwid <- tfrac * switch(interval, confidence = sqrt(ip),
                           prediction = sqrt(ip + pred.var))
    if (type != "terms") {
      predictor <- cbind(predictor, predictor + hwid %o%
                           c(1, -1))
      colnames(predictor) <- c("fit", "lwr",
                               "upr")
    }
    else {
      if (!is.null(terms))
        hwid <- hwid[, terms, drop = FALSE]
      lwr <- predictor + hwid
      upr <- predictor - hwid
    }
  }
  if (se.fit || interval != "none") {
    se <- sqrt(ip)
    if (type == "terms" && !is.null(terms) && !se.fit)
      se <- se[, terms, drop = FALSE]
  }
  if (missing(newdata) && !is.null(na.act <- object$na.action)) {
    predictor <- napredict(na.act, predictor)
    if (se.fit)
      se <- napredict(na.act, se)
  }
  if (type == "terms" && interval != "none") {
    if (missing(newdata) && !is.null(na.act)) {
      lwr <- napredict(na.act, lwr)
      upr <- napredict(na.act, upr)
    }
    list(fit = predictor, se.fit = se, lwr = lwr, upr = upr,
         df = df, residual.scale = sqrt(res.var))
  }
  else if (se.fit)
    list(fit = predictor, se.fit = se, df = df, residual.scale = sqrt(res.var))
  else predictor

}


#' @export
#' @keywords internal
#' @exportClass cpa
#' @method model.matrix cpa
model.matrix.cpa <- function (object, ...)
{
  if (n_match <- match("x", names(object), 0L))
    object[[n_match]]
  else {
    data <- model.frame(object, xlev = object$xlevels, ...)
    if (exists(".GenericCallEnv", inherits = FALSE))
      NextMethod("model.matrix", data = data, contrasts.arg = object$contrasts)
    else {
      dots <- list(...)
      dots$data <- dots$contrasts.arg <- NULL
      do.call("model.matrix.default", c(list(object = object,
                                             data = data, contrasts.arg = object$contrasts),
                                        dots))
    }
  }
}



