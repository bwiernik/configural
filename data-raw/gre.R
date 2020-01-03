gre_mat <- matrix(NA, ncol = 4, nrow = 4, dimnames = list(c("V", "Q", "A", "GPA"), c("V", "Q", "A", "GPA")))
gre <- list(r = gre_mat, rho = gre_mat, n = gre_mat, sevar_r = gre_mat, sevar_rho = gre_mat, sevar_r_restricted = gre_mat, source = gre_mat)

gre$rho <- psychmeta::reshape_vec2mat(c(.56, .77, .34, .73, .32, .36), var_names = c("V", "Q", "A", "GPA"))
gre$r <- psychmeta::reshape_vec2mat(c(.56, .77, .23, .73, .21, .24), var_names = c("V", "Q", "A", "GPA"))
gre$r_restricted <- psychmeta::reshape_vec2mat(
  c(psychmeta:::.attenuate_r_bvdrr(.56, ux = .77, uy = .73),
    psychmeta:::.attenuate_r_bvdrr(.77, ux = .77, uy = .74),
    .23,
    psychmeta:::.attenuate_r_bvdrr(.73, ux = .73, uy = .74),
    .21, .24), var_names = c("V", "Q", "A", "GPA")
)
gre$n[lower.tri(diag(4))] <- c(145912, 3895, 14156, 3895, 14425, 1928)
gre$sevar_r[lower.tri(diag(4))] <- c(NA, NA, .14/sqrt(103), NA, .11/sqrt(103), .12/sqrt(20))
gre$sevar_rho[lower.tri(diag(4))] <- gre$sevar_r[lower.tri(diag(4))] / (gre$r[lower.tri(diag(4))] / gre$rho[lower.tri(diag(4))])
gre$sevar_rho[c(2, 3, 7)] <- sqrt(.02^2 / c(7, 2, 2) +
                                    c((1 - gre$rho[c(2, 3, 7)]^2)^2 / (gre$n[c(2, 3, 7)] - c(7, 2, 2)))
)
gre$sevar_r[c(2, 3, 7)] <- sqrt(.02^2 / c(7, 2, 2) +
                                  c((1 - gre$r[c(2, 3, 7)]^2)^2 / (gre$n[c(2, 3, 7)] - c(7, 2, 2)))
)
gre$sevar_r_restricted[c(2, 3, 7)] <- sqrt(.02^2 / c(7, 2, 2) +
                                             c((1 - gre$r_restricted[c(2, 3, 7)]^2)^2 / (gre$n[c(2, 3, 7)] - c(7, 2, 2)))
)
gre$source[lower.tri(diag(4))] <- "A"

usethis::use_data(gre)

