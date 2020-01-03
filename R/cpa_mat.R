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
  x_col <- stringr::str_split(x_col, pattern = "[+]")[[1]]
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

  if (is.null(n)) n <- n_effective_R2(R2_total, se_var$r.total, p)

  if (is.na(n)) R2_adj <- R2_total else R2_adj <- adjust_Rsq(R2_total, n, p, adjust)
  R_adj <- if (R2_adj < 0) 0 else sqrt(R2_adj)

  pat_adj <- lev * levpat + sqrt(max((1 - levpat2) * (R2_adj - lev2), 0))
  pat2_adj <- pat_adj^2
  beta_cpa_adj <- solve(Rxx_cpa) %*% c(lev, pat_adj)

  delta_lev <- R2_total - pat2
  delta_pat <- R2_total - lev2
  delta_lev_adj <- R2_adj - pat2_adj
  delta_pat_adj <- R2_adj - lev2

  if (is.infinite(n)) {
    moe <- stats::qnorm((1 - conf_level) / 2, lower.tail = FALSE)
  } else {
    moe <- stats::qt((1 - conf_level) / 2, n - p - 1, lower.tail = FALSE)
  }

  se_beta <- sqrt(diag(se_var$beta))
  se_bstar <- sqrt(diag(se_var$bstar))

  beta_mat <- cbind(beta, se_beta, beta - se_beta * moe, beta + se_beta * moe)
  rownames(beta_mat) <- x_col
  colnames(beta_mat) <- c("beta",
                          "SE",
                          paste0((1 - conf_level) / 2 * 100, "%"),
                          paste0((1 - (1 - conf_level) / 2) * 100, "%")
  )

  bstar_mat <- cbind(bstar, se_bstar, bstar - se_bstar * moe, bstar + se_bstar * moe)
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
                       (c(lev, pat, R_total) - sign(c(lev, pat, R_total)) * se_R * moe)^2,
                       (c(lev, pat, R_total) + sign(c(lev, pat, R_total)) * se_R * moe)^2,
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
                           (c(lev, pat_adj, R_adj) - se_R_adj * moe)^2,
                           (c(lev, pat_adj, R_adj) + se_R_adj * moe)^2,
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

  levpat_mat <- cbind(levpat)
  rownames(levpat_mat) <- ""
  colnames(levpat_mat) <- c("r.pattern.level")

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
    vcov = se_var
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
