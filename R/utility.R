#' Adjust a regression model R-squared for overfitting
#'
#' Estimate shrinkage for regression models
#'
#' @param Rsq Observed model R-squared
#' @param n Sample size
#' @param p Number of predictors
#' @param adjust Which adjustment to apply. Options are "fisher" for the Adjusted R-squared method used in [stats::lm()], "pop" for the positive-part Pratt estimator of the population R-squared, and "cv" for the Browne/positive-part Pratt estimator of the cross-validity R-squared. Based on Shieh (2008), these are the estimators for the population and cross-validity R-squared values that show the least bias with a minimal increase in computational complexity.
#'
#' @return An adjusted R-squared value.
#' @export
#'
#' @encoding UTF-8
#'
#' @references
#' Shieh, G. (2008).
#' Improved shrinkage estimation of squared multiple correlation coefficient and squared cross-validity coefficient.
#' _Organizational Research Methods, 11_(2), 387–407. \doi{10.1177/1094428106292901}
#'
#' @examples
#' adjust_Rsq(.55, 100, 6, adjust = "pop")
adjust_Rsq <- function(Rsq, n, p, adjust = c("fisher", "pop", "cv")) {
  adjust <- match.arg(adjust)
  if (adjust == "fisher") {
    return(1 - (1 - Rsq) * ((n - 1) / (n - p - 1)))
  } else {
    R2_P <- max(1 - ((n - 3) / (n - p - 1)) * (1 - Rsq) * (1 + ((2 * (1 - Rsq)) / (n - p - 2.3))), 0)
    if (adjust == "pop") {
      return(R2_P)
    } else {
      R4_B <- R2_P^2 - (2 * p * (1 - R2_P)^2) / ((n - 1) * (n - p - 1))
      R2_CV <- max(((n - p - 3) * R4_B + R2_P) / ((n - 2 * p - 2) * R2_P + p) , 0)
      return(R2_CV)
    }
  }
}

#' Find the harmonic mean of a vector, matrix, or columns of a data.frame
#'
#' The harmonic mean is merely the reciprocal of the arithmetic mean of the reciprocals.
#'
#' @param x A vector, matrix, or data.frame
#' @param na.rm Logical. If `TRUE`, remove `NA` values before processing
#' @param zero Logical, If `TRUE`, if there are any zeroes, return 0, else, return the harmonic mean of the non-zero elements
#'
#' @return The harmonic mean of x
#'
#' @author Adapted from `psych::harmonic.mean()` by William Revelle
#' @export
#'
#' @encoding UTF-8
#'
#' @examples
#' harmonic_mean(1:10)
harmonic_mean <- function(x, na.rm = TRUE, zero = TRUE) {
  if (!zero) {
    x[x == 0] <- NA
  }
  if (is.null(nrow(x))) {
    1 / mean(1 / x, na.rm = na.rm)
  }
  else {
    1 / (apply(1 / x, 2, mean, na.rm = na.rm))
  }
}

#' Calculate a transition matrix for a symmetric matrix
#'
#' The transition matrix extracts the lower triangular elements from a vectorized symmetric matrix (Nel, 1985).
#'
#' @param p The number of columns in a matrix
#'
#' @author Based on internal functions from the \pkg{fungible} package by Niels Waller
#'
#' @export
#'
#' @encoding UTF-8
#'
#' @keywords internal
#'
#' @references
#' Nel, D. G. (1985).
#' A matrix derivation of the asymptotic covariance matrix of sample correlation coefficients.
#' _Linear Algebra and Its Applications, 67_, 137–145. \doi{10.1016/0024-3795(85)90191-0}
#'
#' @examples
#' transition(5)
transition <- function(p) {
  M <- matrix(nrow = p, ncol = p)
  M[lower.tri(M, diag = TRUE)] <- seq(p * (p + 1) / 2)
  M[upper.tri(M, diag = FALSE)] <- t(M)[upper.tri(M, diag = FALSE)]
  Dp <- outer(c(M), unique(c(M)), FUN = function(x, y) as.numeric(x == y))

  p1 <- p2 <- p
  rows <- rep(1, p)
  for (i in 2:p) {
    rows[i] <- rows[i] + p1
    p1 <- p1 + (p2 - 1)
    p2 <- p2 - 1
  }

  Kpc <- (solve(t(Dp) %*% Dp) %*% t(Dp))[-rows, ]

  return(Kpc)

}


#' Calculate the asymptotic sampling covariance matrix for the unique elements of a correlation matrix
#'
#' @param r A correlation matrix
#' @param n The sample size
#'
#' @return The asymptotic sampling covariance matrix
#' @export
#'
#' @encoding UTF-8
#'
#' @author Based on an internal function from the \pkg{fungible} package by Niels Waller
#'
#' @references
#' Nel, D. G. (1985).
#' A matrix derivation of the asymptotic covariance matrix of sample correlation coefficients.
#' _Linear Algebra and Its Applications, 67_, 137–145. \doi{10.1016/0024-3795(85)90191-0}
#'
#' @examples
#' cor_covariance(matrix(c(1, .2, .3, .2, 1, .3, .3, .3, 1), ncol = 3), 100)
cor_covariance <- function(r, n) {
  p <- ncol(r)
  id <- diag(p)

  Ms <-
    (matrix(
      c(rep(
        c(rep(
          c(1, rep(
            0, times = p * p + (p - 1) )
          ), times = p - 1 ),
          1, rep(0, times = p) ), times = p - 1),
        rep(c(1, rep(0, times = p * p + (p - 1))), times = p - 1),
        1),
      nrow = p^2) + diag(p^2) ) / 2

  Md <-
    diag(replace(rep(0, p^2), seq(1, (p^2), by = (p + 1)), 1))

  Psi <-
    0.5 * (4 * Ms %*% (r %x% r) %*% Ms -
             2 * (r %x% r) %*% Md %*% (id %x% r + r %x% id) -
             2 * (id %x% r + r %x% id) %*% Md %*% (r %x% r) +
             (id %x% r + r %x% id) %*% Md %*% (r %x% r) %*%
             Md %*% (id %x% r + r %x% id))

  Kpc <- transition(p)
  out <- (Kpc %*% Psi %*% t(Kpc)) / (n - 3)
  rownames(out) <- colnames(out) <- cor_labels(colnames(r))
  return(out)
}

#' Generate labels for correlations from a vector of variable names
#'
#' This function returns a vector of labels for the unique correlations between pairs of variables from a supplied vector of variable names
#'
#' @param var_names A character vector of variable names
#'
#' @return A vector of correlation labels
#' @export
#'
#' @encoding UTF-8
#'
#' @examples
#' cor_labels(colnames(mindfulness$r))
cor_labels <- function(var_names) {
  return(vechs(t(outer(var_names, var_names, paste, sep = "-"))))
}

#' Estimate the asymptotic sampling covariance matrix for the unique elements of
#' a meta-analytic correlation matrix
#'
#' @param r A meta-analytic matrix of observed correlations (can be full or lower-triangular).
#' @param n A matrix of total sample sizes for the meta-analytic correlations in `r` (can be full or lower-triangular).
#' @param sevar A matrix of estimated sampling error variances for the meta-analytic correlations in `r` (can be full or lower-triangular).
#' @param source A matrix indicating the sources of the meta-analytic correlations in `r` (can be full or lower-triangular). Used to estimate overlapping sample size for correlations when `n_overlap == NULL`.
#' @param rho A meta-analytic matrix of corrected correlations (can be full or lower-triangular).
#' @param sevar_rho A matrix of estimated sampling error variances for the meta-analytic corrected correlations in `rho` (can be full or lower-triangular).
#' @param n_overlap A matrix indicating the overlapping sample size for the unique (lower triangular) values in `r` (can be full or lower-triangular). Values must be arranged in the order returned by `cor_labels(colnames(R))`.
#'
#' @details If both `source` and `n_overlap` are `NULL`, it is assumed that all meta-analytic correlations come from the the same source.
#'
#' @return The estimated asymptotic sampling covariance matrix
#' @export
#'
#' @encoding UTF-8
#'
#' @references
#' Nel, D. G. (1985).
#' A matrix derivation of the asymptotic covariance matrix of sample correlation coefficients.
#' _Linear Algebra and Its Applications, 67_, 137–145. \doi{10.1016/0024-3795(85)90191-0}
#'
#' Wiernik, B. M. (2018).
#' _Accounting for dependency in meta-analytic structural equations modeling: A flexible alternative to generalized least squares and two-stage structural equations modeling._
#' Unpublished manuscript.
#'
#' @examples
#' cor_covariance_meta(r = mindfulness$r, n = mindfulness$n,
#'                     sevar = mindfulness$sevar_r, source = mindfulness$source)
cor_covariance_meta <- function(r, n, sevar, source = NULL, rho = NULL, sevar_rho = NULL, n_overlap = NULL) {
  if (is.null(colnames(r))) colnames(r) <- 1:ncol(r)
  if (is.null(rownames(r))) rownames(r) <- 1:nrow(r)
  if (!all(colnames(r) == rownames(r))) stop("Row names and column names of `r` must be the same")
  var_names <- colnames(r)
  cor_names <- cor_labels(var_names)

  if (is.null(colnames(n))) if (!all(colnames(n) == var_names)) stop("Column names of `n` and `r` must be the same")
  if (is.null(rownames(n))) if (!all(rownames(n) == var_names)) stop("Row names of `n` and `r` must be the same")

  if (is.null(colnames(sevar))) if (!all(colnames(sevar) == var_names)) stop("Column names of `sevar` and `r` must be the same")
  if (is.null(rownames(sevar))) if (!all(rownames(sevar) == var_names)) stop("Row names of `sevar` and `r` must be the same")

  if (!is.null(source)) {
    if (is.null(colnames(source))) if (!all(colnames(source) == var_names)) stop("Column names of `source` and `r` must be the same")
    if (is.null(rownames(source))) if (!all(rownames(source) == var_names)) stop("Row names of `source` and `r` must be the same")
  }

  if (!is.null(n_overlap)) {
    if (!is.null(colnames(n_overlap))) if (!all(colnames(n_overlap) == cor_names)) stop("Column names of `n_overlap` must be identical to the result of:\n  `cor_labels(colnames(r))`")
    if (!is.null(rownames(n_overlap))) if (!all(rownames(n_overlap) == cor_names)) stop("Row names of `n_overlap` must be identical to the result of:\n  `cor_labels(colnames(r))`")
  }

  r <- vechs2full(vechs(r))
  n <- vechs(n)
  sevar <- vechs(sevar)
  if (!is.null(rho)) {
    A <- vechs(r) / vechs(rho)
  } else if (!is.null(sevar_rho)) {
    A <- sqrt(sevar / vech(sevar_rho))
  }


  if (is.null(n_overlap)) {
    n_overlap <- outer(n, n, FUN = function(X, Y) apply(cbind(X, Y), 1, min))
    if (!is.null(source)) {
      same_source <- outer(vechs(source), vechs(source), `==`)
    } else same_source <- matrix(1, nrow(r), ncol(r))
  } else same_source <- matrix(1, nrow(r), ncol(r))

  n_overlap <- n_overlap * same_source
  n_product <- outer(n, n)

  p <- ncol(r)
  id <- diag(p)

  Ms <-
    (matrix(
      c(rep(
        c(rep(
          c(1, rep(
            0, times = p * p + (p - 1) )
          ), times = p - 1 ),
          1, rep(0, times = p) ), times = p - 1),
        rep(c(1, rep(0, times = p * p + (p - 1))), times = p - 1),
        1),
      nrow = p^2) + diag(p^2)) / 2

  Md <- diag(replace(rep(0, p^2), seq(1, (p^2), by = (p + 1)), 1))

  Psi <-
    0.5 * (4 * Ms %*% (r %x% r) %*% Ms -
             2 * (r %x% r) %*% Md %*% (id %x% r + r %x% id) -
             2 * (id %x% r + r %x% id) %*% Md %*% (r %x% r) +
             (id %x% r + r %x% id) %*% Md %*% (r %x% r) %*%
             Md %*% (id %x% r + r %x% id))

  Kpc <- transition(p)
  out <- (Kpc %*% Psi %*% t(Kpc)) * (n_overlap / n_product)
  diag(out) <- sevar
  if (exists("A")) out <- out / outer(A, A)
  dimnames(out) <- list(cor_names, cor_names)

  if (any(base::eigen(out)$values < 0)) {
    warning(
      paste("Estiamted sampling covariance matrix is not positive definite.",
            "Consider smoothing the matrix with:",
            "Matrix::nearPD(se_var_mat, keepDiag = TRUE, base.matrix = TRUE, ensureSymmetry = TRUE)$mat or",
            "Matrix::nearPD(se_var_mat, base.matrix = TRUE, ensureSymmetry = TRUE)$mat",
            sep = "\n"
      )
    )
  }

  out
}

#' Quadratic form matrix product
#'
#' Calculate the quadratic form \deqn{Q=x^{\prime}Ax}{Q = t(x) \%\*\% A \%\*\% x}
#'
#' @param A A square matrix
#' @param x A vector or matrix
#'
#' @return The quadratic product
#' @export
#'
#' @encoding UTF-8
#'
#' @examples
#' diag(5) %&% 1:5
`%&%` <- function(A, x) {
  t(x) %*% A %*% x
}

.remove_charmargins <- function(x) {
  needs_trim <- nchar(x) > 1 & substr(x, nchar(x), nchar(x)) == " "
  while (any(needs_trim)) {
    x[needs_trim] <- substr(x[needs_trim], 1, nchar(x[needs_trim]) - 1)
    needs_trim <- nchar(x) > 1 & substr(x, nchar(x), nchar(x)) == " "
  }

  needs_trim <- nchar(x) > 1 & substr(x, 1, 1) == " "
  while (any(needs_trim)) {
    x[needs_trim] <- substr(x[needs_trim], 2, nchar(x[needs_trim]))
    needs_trim <- nchar(x) > 1 & substr(x, 1, 1) == " "
  }

  x
}

#' @title Weighted descriptive statistics for a vector of numbers
#'
#' @description
#' Compute the weighted mean and variance of a vector of numeric values. If no weights are supplied, defaults to computing the unweighted mean and the unweighted maximum-likelihood variance.
#'
#' @param x Vector of values to be analyzed.
#' @param wt Weights associated with the values in x.
#' @param unbiased Logical scalar determining whether variance should be unbiased (TRUE) or maximum-likelihood (FALSE).
#' @param df_type Character scalar determining whether the degrees of freedom for unbiased estimates should be based on numbers of cases ("count"; default) or sums of weights ("sum_wts").
#'
#' @return A weighted mean and variance if weights are supplied or an unweighted mean and variance if weights are not supplied.
#' @export
#'
#' @details
#' The weighted mean is computed as
#' \deqn{\bar{x}_{w}=\frac{\Sigma_{i=1}^{k}x_{i}w_{i}}{\Sigma_{i=1}^{k}w_{i}}}{sum(x * wt) / sum(wt)}
#' where \emph{x} is a numeric vector and \emph{w} is a vector of weights.
#'
#' The weighted variance is computed as
#' \deqn{var_{w}(x)=\frac{\Sigma_{i=1}^{k}\left(x_{i}-\bar{x}_{w}\right)^{2}w_{i}}{\Sigma_{i=1}^{k}w_{i}}}{var(x) = sum((x - sum(x * wt) / sum(wt))^2 * wt) / sum(wt)}
#' and the unbiased weighted variance is estimated by multiplying \eqn{var_{w}(x)}{var(x)} by \eqn{\frac{k}{k-1}}{k/(k-1)}.
#'
#' @author Jeffrey A. Dahlke
#'
#' @keywords univar
#'
#' @encoding UTF-8
#'
#' @examples
#' wt_dist(x = c(.1, .3, .5), wt = c(100, 200, 300))
wt_dist <- function(x, wt = rep(1, length(x)), unbiased = TRUE, df_type = c("count", "sum_wts")) {
  df_type <- match.arg(arg = df_type, choices = c("count", "sum_wts"))
  if (length(x) != length(wt)) stop("Lengths of x and wt differ")
  x[is.na(wt)] <- NA
  wt[is.na(x)] <- NA
  mean_x <- sum(as.numeric(x * wt), na.rm = TRUE) / sum(as.numeric(wt), na.rm = TRUE)
  var_x <- sum(as.numeric((x - mean_x)^2 * wt), na.rm = TRUE) / sum(as.numeric(wt), na.rm = TRUE)
  if (unbiased) {
    if (length(x) == 1) {
      var_x <- 0
    } else {
      if (df_type == "count") {
        var_x <- var_x * length(x) / (length(x) - 1)
      } else if (df_type == "sum_wts") {
        var_x <- var_x * sum(wt, na.rm = TRUE) / (sum(wt, na.rm = TRUE) - 1)
      }
    }
  }
  c(mean = mean_x, var = var_x)
}


#' @rdname wt_dist
#' @export
#' @examples
#' wt_mean(x = c(.1, .3, .5), wt = c(100, 200, 300))
wt_mean <- function(x, wt = rep(1, length(x))) {
  if (length(x) != length(wt)) stop("Lengths of x and wt differ")
  x[is.na(wt)] <- NA
  wt[is.na(x)] <- NA
  sum(as.numeric(x * wt), na.rm = TRUE) / sum(as.numeric(wt), na.rm = TRUE)
}

#' @rdname wt_dist
#' @export
#' @examples
#' wt_var(x = c(.1, .3, .5), wt = c(100, 200, 300))
wt_var <- function(x, wt = rep(1, length(x)), unbiased = TRUE, df_type = c("count", "sum_wts")) {
  df_type <- match.arg(arg = df_type, choices = c("count", "sum_wts"))
  if (length(x) != length(wt)) stop("Lengths of x and wt differ")
  x[is.na(wt)] <- NA
  wt[is.na(x)] <- NA
  mean_x <- sum(as.numeric(x * wt), na.rm = TRUE) / sum(as.numeric(wt), na.rm = TRUE)
  var_x <- sum(as.numeric((x - mean_x)^2 * wt), na.rm = TRUE) / sum(as.numeric(wt), na.rm = TRUE)
  if (unbiased) {
    if (length(x) == 1) {
      var_x <- 0
    } else {
      if (df_type == "count") {
        var_x <- var_x * length(x) / (length(x) - 1)
      } else if (df_type == "sum_wts") {
        var_x <- var_x * sum(wt, na.rm = TRUE) / (sum(wt, na.rm = TRUE) - 1)
      }
    }
  }
  var_x
}




#' Compute weighted covariances
#'
#' Compute the weighted covariance among variables in a matrix or between the variables in two separate matrices/vectors.
#'
#' @param x Vector or matrix of x variables.
#' @param y Vector or matrix of y variables
#' @param wt Vector of weights
#' @param as_cor Logical scalar that determines whether the covariances should be standardized (TRUE) or unstandardized (FALSE).
#' @param use Method for handling missing values. "everything" uses all values and does not account for missingness, "listwise" uses only complete cases, and "pairwise" uses pairwise deletion.
#' @param unbiased Logical scalar determining whether variance should be unbiased (TRUE) or maximum-likelihood (FALSE).
#' @param df_type Character scalar determining whether the degrees of freedom for unbiased estimates should be based on numbers of cases (n - 1; "count"; default) or squared sums of weights (1 - sum(w^2); "sum_wts").
#'
#' @return Scalar, vector, or matrix of covariances.
#' @export
#'
#' @author Jeffrey A. Dahlke
#'
#' @encoding UTF-8
#'
#' @examples
#' wt_cov(x = c(1, 0, 2), y = c(1, 2, 3), wt = c(1, 2, 2), as_cor = FALSE, use = "everything")
#' wt_cov(x = c(1, 0, 2), y = c(1, 2, 3), wt = c(1, 2, 2), as_cor = TRUE, use = "everything")
#' wt_cov(x = cbind(c(1, 0, 2), c(1, 2, 3)), wt = c(1, 2, 2), as_cor = FALSE, use = "everything")
#' wt_cov(x = cbind(c(1, 0, 2), c(1, 2, 3)), wt = c(1, 2, 2), as_cor = TRUE, use = "everything")
#' wt_cov(x = cbind(c(1, 0, 2, NA), c(1, 2, 3, 3)),
#'        wt = c(1, 2, 2, 1), as_cor = FALSE, use = "listwise")
#' wt_cov(x = cbind(c(1, 0, 2, NA), c(1, 2, 3, 3)),
#'        wt = c(1, 2, 2, 1), as_cor = TRUE, use = "listwise")
wt_cov <- function(x, y = NULL, wt = NULL, as_cor = FALSE,
                   use = c("everything", "listwise", "pairwise"),
                   unbiased = TRUE, df_type = c("count", "sum_wts")){

  use <- match.arg(arg = use, choices = c("everything", "listwise", "pairwise"))
  df_type <- match.arg(arg = df_type, choices = c("count", "sum_wts"))

  if (is.null(x)) {
    if (is.null(y)) {
      stop("x and y cannot both be NULL")
    } else {
      x <- y
      y <- NULL
    }
  }
  x <- as.matrix(x)
  if (is.null(wt)) wt <- rep(1 / nrow(x), nrow(x))

  if (use != "everything") {
    if (use == "listwise") {
      use_x <- apply(!is.na(x), 1, all) & !is.na(wt)
      x[!use_x,] <- wt[!use_x] <- NA
    } else {
      use_x <- TRUE
    }
  } else {
    use_x <- TRUE
  }

  if (use != "everything") {
    mean_x <- apply(x, 2, function(x) wt_mean(x = x, wt = wt))
  } else {
    mean_x <- (wt %*% x) / sum(wt)
  }

  x <- t(t(x) - drop(mean_x))

  if (is.null(y)) {
    y <- x
  } else {
    y <- as.matrix(y)

    if (use != "everything") {
      if (use == "complete") {
        use_y <- apply(!is.na(y), 1, all) & !is.na(wt)
        y[!use_y,] <- NA
      }
    }
    if (use != "everything") {
      mean_y <- apply(y, 2, function(x) wt_mean(x = x, wt = wt))
    } else {
      mean_y <- (wt %*% y) / sum(wt)
    }

    y <- t(t(y) - drop(mean_y))
  }

  if (use == "everything") {
    out <- t(y) %*% (x * drop(wt)) / sum(wt)
  } else {
    out <- matrix(NA, nrow = ncol(y), ncol = ncol(x))
    for (i in 1:nrow(out)) {
      for (j in 1:ncol(out)) {
        out[i,j] <- wt_mean(x = y[,i] * x[,j], wt = wt)
      }
    }
  }

  if (unbiased) {
    if (df_type == "count") {
      n <- t(!is.na(y)) %*% !is.na(x)
      out <- out * n / (n - 1)
    } else if (df_type == "sum_wts") {
      wt_mat_x <- matrix(wt, nrow = nrow(x), ncol = ncol(x))
      wt_mat_y <- matrix(wt, nrow = nrow(x), ncol = ncol(x))
      wt_mat_x[is.na(wt_mat_x)] <- wt_mat_y[is.na(wt_mat_y)] <- 0
      out <- out / (1 - t(wt_mat_y) %*% wt_mat_x)
    }
  }
  rownames(out) <- colnames(x)
  colnames(out) <- colnames(y)
  if (as_cor) out <- stats::cov2cor(out)
  if (length(out) == 1) {
    drop(out)
  } else out
}

#' @rdname wt_cov
#' @export
wt_cor <- function(x, y = NULL, wt = NULL, use = "everything"){
  wt_cov(x = x, y = y, wt = wt, as_cor = TRUE, use = use)
}

#' Vectorize a matrix
#'
#' `cvec` returns the column-wise vectorization of an input matrix (stacking
#' the columns on one another). `rvec` returns the row-wise vectorization of an
#' input matrix (concatenating the rows after each other). `vech` returns the
#' column-wise half-vectorization of an input matrix (stacking the lower
#' triangular elements of the matrix, including the diagonal). `vechs` returns
#' the strict column-wise half-vectorization of an input matrix (stacking the
#' lower triangular elements of the matrix, excluding the diagonal). All
#' functions return the output as a vector.
#'
#' @param x A matrix
#'
#' @encoding UTF-8
#'
#' @return A vector of values
#'
#' @author Based on functions from the the \pkg{OpenMx} package
#'
#' @export
#'
#' @examples
#' cvec(matrix(1:9, 3, 3))
#' rvec(matrix(1:9, 3, 3))
#' vech(matrix(1:9, 3, 3))
#' vechs(matrix(1:9, 3, 3))
#' vechs(matrix(1:12, 3, 4))
vech <- function(x) {
  return(x[lower.tri(x, diag = TRUE)])
}

#' @rdname vech
#' @export
vechs <- function(x) {
  return(x[lower.tri(x, diag = FALSE)])
}

#' @rdname vech
#' @export
cvec <- function(x) {
  x <- unlist(x)
  return(as.vector(matrix(x, length(x), 1)))
}

#' @rdname vech
#' @export
rvec <- function(x) {
  x <- t(x)
  return(as.vector(matrix(x, length(x), 1)))
}

#' Inverse vectorize a matrix
#'
#' These functions return the symmetric matrix that produces the given
#' half-vectorization result.
#'
#' @details The input consists of a vector of the elements in the lower triangle
#' of the resulting matrix (for `vech2full`, including the elements along the diagonal
#' of the matrix, as a column vector), filled column-wise. For `vechs2full`,
#' the diagonal values are filled as 1 by default, alternative values can be
#' specified using the `diag` argument. The inverse half-vectorization takes a
#' vector and reconstructs a symmetric matrix such that vech2full(vech(x)) is
#' identical to x if x is symmetric.
#'
#' @param x A vector
#' @param diagonal A value or vector of values to enter on the diagonal for `vechs2full` (default = 1)
#'
#' @return A symmetric matrix
#'
#' @encoding UTF-8
#'
#' @author Based on functions from the the \pkg{OpenMx} package
#'
#' @export
#'
#' @examples
#' vech2full(c(1, 2, 3, 5, 6, 9))
#' vechs2full(c(2, 3, 6), diagonal = 0)
vech2full <- function(x) {
  if (is.matrix(x)) x <- as.vector(x)
  if (is.vector(x)) {
    dimension <- length(x)
  }
  else {
    stop("Input to vech2full must be either a matrix or a vector.")
  }
  k <- sqrt(2 * dimension + 0.25) - 0.5
  ret <- matrix(0, nrow = k, ncol = k)
  if (nrow(ret) != k) {
    stop("Incorrect number of elements in vector to construct a matrix from a half-vectorization.")
  }
  ret[lower.tri(ret, diag = TRUE)] <- as.vector(x)
  ret[upper.tri(ret)] <- t(ret)[upper.tri(ret)]
  return(ret)
}

#' @rdname vech2full
#' @export
vechs2full <- function(x, diagonal = 1) {
  if (is.matrix(x)) x <- as.vector(x)
  if (is.vector(x)) {
    dimension <- length(x)
  }
  else {
    stop("Input to the function vechs2full must be either a matrix or a vector.")
  }
  k <- sqrt(2 * dimension + 0.25) + 0.5
  ret <- matrix(0, nrow = k, ncol = k)
  if (nrow(ret) != k) {
    stop("Incorrect number of elements in vector to construct a matrix from a strict half-vectorization.")
  }
  ret[lower.tri(ret, diag = FALSE)] <- as.vector(x)
  ret[upper.tri(ret)] <- t(ret)[upper.tri(ret)]
  diag(ret) <- diagonal
  return(ret)
}
