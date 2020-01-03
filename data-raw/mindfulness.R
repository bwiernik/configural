mindfulness_mat <- matrix(NA, ncol = 6, nrow = 6, dimnames = list(c("ES", "A", "C", "Ex", "O", "mindfulness"), c("ES", "A", "C", "Ex", "O", "mindfulness")))
mindfulness <- list(r = mindfulness_mat,
                    rho = mindfulness_mat,
                    n = mindfulness_mat,
                    sevar_r = mindfulness_mat,
                    sevar_rho = mindfulness_mat,
                    source = mindfulness_mat)
mindfulness$r <-
  psychmeta::reshape_vec2mat(c(0.24, 0.27, 0.22, 0.07, 0.47, 0.32, 0.16, 0.15,
                               0.26, 0.15, 0.09, 0.34, 0.26, 0.17, 0.15),
                             var_names = c("ES", "A", "C", "Ex", "O", "mindfulness")
                             )
# Replace with computations
mindfulness$rho <-
  psychmeta::reshape_vec2mat(c(0.31, 0.33, 0.27, 0.09, 0.566429529, 0.41, 0.2,
                               0.19, 0.320465449, 0.19, 0.12, 0.414361549, 0.33,
                               0.204878766, 0.187009015),
                             var_names = c("ES", "A", "C", "Ex", "O", "mindfulness")
                             )
mindfulness$n[lower.tri(mindfulness_mat)] <-
  c(79610, 84256, 92111, 65095, 6185, 76306, 75274, 61538, 6145, 74154, 62258,
    5976, 71206, 6910, 9441
    )
# Replace with computations
mindfulness$sevar_r[lower.tri(mindfulness_mat)] <-
  c(0.000239521, 0.000174096, 0.000121327, 0.000166234, 0.000441, 0.000228481,
    0.000279114, 0.000114189, 0.000784, 0.000164103, 0.000243919, 4e-04,
    0.000161006, 0.000576, 0.000484
    )
# Replace with computations
mindfulness$sevar_rho[lower.tri(mindfulness_mat)] <-
  c(0.000309, 0.000213, 0.000149, 0.000214, 0.000640523, 0.000293, 0.000349,
    0.000145, 0.001191055, 0.000208, 0.000325, 0.000594102, 0.000204,
    0.000836601, 0.000752295
    )
mindfulness$source[lower.tri(mindfulness_mat)] <-
  c('A', 'A', 'A', 'A', 'B', 'A', 'A', 'A', 'B', 'A', 'A', 'B', 'A', 'B', 'B')

usethis::use_data(mindfulness)

