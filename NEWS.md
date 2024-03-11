# configural 0.1.5 (2024-03-11)

- Added warnings when estimated or supplied sampling error covariance matrix is non-positive definite
- Bug fix in `cor_covariance_meta()` to accommodate matrices with missing values both above and below diagonal
- Added `complete_matrix()` function to make a matrix symmetric by averaging with its transpose

# configural 0.1.4 (2021-01-18)

- Bug fix in `print.fungible_extrema()`

# configural 0.1.2 (2020-05-31)

- Added `cpa_scores()` function to compute profile and level scores for new data
- Additional confidence intervals for CPA parameters in `cpa_mat()`
- `cor_covariance()` now assigns dimnames to output matrix
- Bug fixes in `var_error_cpa()`

# configural 0.1.0 (2019-01-20)

- configural has officially launched for public beta!
- Please see the "configural-package" entry in the configural manual for an overview of the package and its applications.
