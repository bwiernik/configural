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
#' @references
#' Shieh, G. (2008).
#' Improved shrinkage estimation of squared multiple correlation coefficient and squared cross-validity coefficient.
#' _Organizational Research Methods, 11_(2), 387–407. <https://doi.org/10/bcwqf3>
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

#' Calculate a transition matric for a symmetric matrix
#'
#' The transition matrix extracts the lower triangular elements from a vectorized symmetric matrix (Nel, 1985).
#'
#' @param p The number of columns in a matrix
#'
#' @author Based on internal functions from the \pkg{fungible} package by Niels Waller
#'
#' @export
#'
#' @keywords internal
#'
#' @importFrom dplyr %>%
#'
#' @references
#' Nel, D. G. (1985).
#' A matrix derivation of the asymptotic covariance matrix of sample correlation coefficients.
#' _Linear Algebra and Its Applications, 67_, 137–145. <https://doi.org/10/c75jmg>
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

  Kpc <- solve(t(Dp) %*% Dp) %*% t(Dp) %>% `[`(-rows, )

  return(Kpc)

}


#' Calculate the asymptotic sampling covariance matrix for the unique elements of a correlation matrix
#'
#' @param R A correlation matrix
#' @param n The sample size
#'
#' @return The asymptotic sampling covariance matrix
#' @export
#'
#' @author Based on an internal function from the \pkg{fungible} package by Niels Waller
#'
#' @references
#' Nel, D. G. (1985).
#' A matrix derivation of the asymptotic covariance matrix of sample correlation coefficients.
#' _Linear Algebra and Its Applications, 67_, 137–145. <https://doi.org/10/c75jmg>
#'
#' @examples
#' cor_covariance(matrix(c(1, .2, .3, .2, 1, .3, .3, .3, 1), ncol = 3), 100)
cor_covariance <- function(R, n) {
  p <- ncol(R)
  id <- diag(p)

  Ms <-
    matrix(
      c(rep(
        c(rep(
          c(1, rep(
            0, times = p * p + (p - 1) )
          ), times = p - 1 ),
          1, rep(0, times = p) ), times = p - 1),
        rep(c(1, rep(0, times = p * p + (p - 1))), times = p - 1),
        1),
      nrow = p^2) %>%
    `+`(diag(p^2)) / 2

  Md <-
    rep(0, p^2) %>%
    replace(seq(1, (p^2), by = (p + 1)), 1) %>%
    diag()

  Psi <-
    0.5 * (4 * Ms %*% (R %x% R) %*% Ms -
             2 * (R %x% R) %*% Md %*% (id %x% R + R %x% id) -
             2 * (id %x% R + R %x% id) %*% Md %*% (R %x% R) +
             (id %x% R + R %x% id) %*% Md %*% (R %x% R) %*%
             Md %*% (id %x% R + R %x% id))

  Kpc <- transition(p)
  (Kpc %*% Psi %*% t(Kpc)) / (n - 3)
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
#' @examples
#' diag(5) %&% 1:5
`%&%` <- function(A, x) {
  t(x) %*% A %*% x
}

.remove_charmargins <- function(x) {
  needs_trim <- nchar(x) > 1 & substr(x, nchar(x), nchar(x)) == " "
  while(any(needs_trim)){
    x[needs_trim] <- substr(x[needs_trim], 1, nchar(x[needs_trim]) - 1)
    needs_trim <- nchar(x) > 1 & substr(x, nchar(x), nchar(x)) == " "
  }

  needs_trim <- nchar(x) > 1 & substr(x, 1, 1) == " "
  while(any(needs_trim)){
    x[needs_trim] <- substr(x[needs_trim], 2, nchar(x[needs_trim]))
    needs_trim <- nchar(x) > 1 & substr(x, 1, 1) == " "
  }

  x
}
