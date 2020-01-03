gre_mat <- matrix(NA, ncol = 4, nrow = 4, dimnames = list(c("V", "Q", "A", "GPA"), c("V", "Q", "A", "GPA")))
gre <- list(r = gre_mat, rho = gre_mat,
            n = gre_mat,
            sevar_r = gre_mat, sevar_rho = gre_mat,
            source = gre_mat)

gre$rho <- psychmeta::reshape_vec2mat(c(.56, .77, .34, .73, .32, .36), var_names = c("V", "Q", "A", "GPA"))
gre$r <- psychmeta::reshape_vec2mat(c(.56, .77, .23, .73, .21, .24), var_names = c("V", "Q", "A", "GPA"))
# N
gre$n[lower.tri(gre_mat)] <- c(145912, 3895, 14156, 3895, 14425, 1928)
# SE var
gre$sevar_r[lower.tri(gre_mat)] <- c(NA, NA, .14^2/103, NA, .11^2/103, .12^2/20)
gre$sevar_rho[lower.tri(gre_mat)] <- gre$sevar_r[lower.tri(gre_mat)] / (gre$r[lower.tri(gre_mat)] / gre$rho[lower.tri(gre_mat)])
gre$sevar_rho[c(2, 3, 7)] <- sqrt(.01^2 / c(7, 2, 2) +
                                    c((1 - gre$rho[c(2, 3, 7)]^2)^2 / (gre$n[c(2, 3, 7)] - c(7, 2, 2)))
)
gre$sevar_r[c(2, 3, 7)] <- sqrt(.01^2 / c(7, 2, 2) +
                                  c((1 - gre$r[c(2, 3, 7)]^2)^2 / (gre$n[c(2, 3, 7)] - c(7, 2, 2)))
)
# Source
gre$source[lower.tri(gre_mat)] <- "A"

usethis::use_data(gre, overwrite = TRUE)

