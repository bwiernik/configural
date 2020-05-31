#' Conduct criterion profile analysis using a correlation matrix
#'
#' @param formula Regression formula with a single outcome variable on the left-hand side and one or more predictor variables on the right-hand side (e.g., Y ~ X1 + X2).
#' @param cov_mat Correlation matrix containing the variables to be used in the regression.
#' @param n Sample size. Used to compute adjusted R-squared and, if `se_var_mat` is NULL, standard errors. If NULL and `se_var_mat` is specified, effective sample size is computed based on `se_var_mat` (cf. Revelle et al., 2017).
#' @param se_var_mat Optional. The sampling error covariance matrix among the unique elements of `cov_mat`. Used to calculate standard errors. If not supplied, the sampling covariance matrix is calculated using `n`.
#' @param se_beta_method Method to use to estimate the standard errors of standardized regression (beta) coefficients. Current options include "normal" (use the Jones-Waller, 2015, normal-theory approach) and "lm" (estimate standard errors using conventional regression formulas).
#' @param adjust Method to adjust R-squared for overfitting. See [adjust_Rsq()] for details.
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
#' _Psychometrika, 80_(2), 365–378. <https://doi.org/10.1007/s11336-013-9380-y>
#'
#' Revelle, W., Condon, D. M., Wilt, J., French, J. A., Brown, A., & Elleman, L. G. (2017).
#' Web- and phone-based data collection using planned missing designs.
#' In N. G. Fielding, R. M. Lee, & G. Blank, _The SAGE Handbook of Online Research Methods_ (pp. 578–594).
#' SAGE Publications. <https://doi.org/10.4135/9781473957992.n33>
#'
#' Wiernik, B. M., Wilmot, M. P., Davison, M. L., & Ones, D. S. (2019).
#' Meta-analytic criterion profile analysis.
#' _Psychological Methods_ <https://doi.org/10.1037/met0000305>
#'
#' @examples
#' sevar <- cor_covariance_meta(mindfulness$r, mindfulness$n, mindfulness$sevar_r, mindfulness$source)
#' cpa_mat(mindfulness ~ ES + A + C + Ex + O,
#'           cov_mat = mindfulness$r,
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

#' Compute CPA level and pattern scores for a set of data
#'
#' @param cpa_mod A model returned from [cpa_mat()] (a model of class `"cpa"`)
#' @param newdata A data frame or matrix containing columns with the same names as
#' the predictors in `cpa_mod`.
#' @param cpa_names Character vector of length 2 giving the variable names to assign
#' to the CPA score columns.
#' @param augment Should be CPA score columns be added to `newdata` (`TRUE`, default)
#' or returned alone (`FALSE`)?
#' @param scale Logical. Should the variables in `newdata` be scaled (standardized)?
#' @param scale_center If `scale` is `TRUE`, passed to the `center` argument in
#' [base::scale()]. Can be `TRUE` (center columns of `newdata` around the column
#' means), `FALSE` (don't center), or a numeric vector of length equal to the
#' number of predictors in `cpa_mod` containing the values to center around.
#' @param scale_scale If `scale` is `TRUE`, passed to the `scale` argument in
#' [base::scale()]. Can be `TRUE` (scale/standardize columns of `newdata` using
#' the column standard deviations or root mean squares), `FALSE` (don't scale),
#' or a numeric vector of length equal to the number of predictors in `cpa_mod`
#' containing the values to scale by. See [base::scale()] for details.
#'
#' @return
#' A data frame containing the CPA score variables.
#'
#' @export
#'
#' @importFrom stats coef cov
#'
#' @examples
#' sevar <- cor_covariance_meta(mindfulness$r, mindfulness$n, mindfulness$sevar_r, mindfulness$source)
#' cpa_mod <- cpa_mat(mindfulness ~ ES + A + C + Ex + O,
#'                    cov_mat = mindfulness$r,
#'                    n = NULL,
#'                    se_var_mat = sevar,
#'                    adjust = "pop")
#'
#' nd <- data.frame(ES = c(4.2, 3.2, 3.4, 4.2, 3.8, 4.0, 5.6, 2.8, 3.4, 2.8),
#'                  A  = c(4.0, 4.2, 3.8, 4.6, 4.0, 4.6, 4.6, 2.6, 3.6, 5.4),
#'                  C  = c(2.8, 4.0, 4.0, 3.0, 4.4, 5.6, 4.4, 3.4, 4.0, 5.6),
#'                  Ex = c(3.8, 5.0, 4.2, 3.6, 4.8, 5.6, 4.2, 2.4, 3.4, 4.8),
#'                  O  = c(3.0, 4.0, 4.8, 3.2, 3.6, 5.0, 5.4, 4.2, 5.0, 5.2)
#'                  )
#'
#' nd_cpa <- cpa_scores(cpa_mod, nd, augment = FALSE)
#' nd_augment <- cpa_scores(cpa_mod, nd, augment = FALSE)
cpa_scores <- function(cpa_mod, newdata = NULL, augment = TRUE,
                       cpa_names = c("cpa_lev", "cpa_pat"),
                       scale = FALSE, scale_center = TRUE, scale_scale = TRUE) {
  if (!is.logical(scale)) {
    stop("`scale` must be either TRUE or FALSE.")
  }
  if (!inherits(cpa_mod, "cpa")) {
    stop("`cpa_mod` must be an object of class 'cpa'.")
  }

  bstar <- coef(cpa_mod)[,"bstar"]
  pred_names <- names(bstar)

  data_pred <- as.data.frame(newdata[,pred_names])
  if (scale) {
    if (!is.logical(scale_center)) {
      if (is.numeric(scale_center)) {
        if (!(length(scale_center) == 1 | length(scale_center) == length(bstar))) {
          stop("`scale_center` must be of length 1 or the same length as the number of predictors in `cpa_mod`.")
        }
      } else {
        stop("`scale_center` must be logical or numeric.")
      }
    }
    if (!is.logical(scale_scale)) {
      if (is.numeric(scale_scale)) {
        if (!(length(scale_scale) == 1 | length(scale_scale) == length(bstar))) {
          stop("`scale_scale` must be of length 1 or the same length as the number of predictors in `cpa_mod`.")
        }
      } else {
        stop("`scale_scale` must be logical or numeric.")
      }
    }
    data_pred <- mapply(scale, data_pred, scale_center, scale_scale)
  }

  lev <- rowMeans(data_pred)
  pat <- apply(data_pred, 1, cov, y = bstar)
  cpa_scores <- data.frame(lev, pat)
  colnames(cpa_scores) <- cpa_names

  if (augment) {
    data_return <- cbind(newdata, cpa_scores)
  } else {
    data_return <- cpa_scores
  }
  if (is.matrix(newdata)) {
    data_return <- as.matrix(data_return)
  } else {
    class(data_return) <- class(newdata)
  }

  return(data_return)
}
