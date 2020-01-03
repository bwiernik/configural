#' Meta-analytic correlations among Big Five personality traits and trait mindfulness
#'
#' Big Five intercorrelations from Davies et al. (2015). Big Five–Mindfulness
#' correlations from Hanley and Garland (2017). Coefficient alpha for
#' mindfulness measures taken from Giluk (2009).
#'
#' @docType data
#'
#' @encoding UTF-8
#'
#' @usage data(mindfulness)
#'
#' @format list with entries `r` (mean oberseved correlations), `rho` (mean
#' corrected correlations), `n` (sample sizes), `sevar_r` (sampling error
#' variances for mean observed correlations), `sevar_rho` (sampling error
#' variances for mean corrected correlations), and `source` (character labels
#' indicating which meta-analytic correlations came from the same source)
#'
#' @keywords datasets
#'
#' @references
#' Davies, S. E., Connelly, B. L., Ones, D. S., & Birkland, A. S. (2015).
#' The general factor of personality: The “Big One,” a self-evaluative trait, or a methodological gnat that won’t go away?
#' _Personality and Individual Differences, 81_, 13–22. <https://doi.org/10/bc98>
#'
#' Giluk, T. L. (2009).
#' Mindfulness, Big Five personality, and affect: A meta-analysis.
#' _Personality and Individual Differences, 47_(8), 805–811. <https://doi.org/10/d23n2k>
#'
#' Hanley, A. W., & Garland, E. L. (2017).
#' The mindful personality: A meta-analysis from a cybernetic perspective.
#' _Mindfulness, 8_(6), 1456–1470. <https://doi.org/10/gcmvb4>
#'
#' @examples
#' data(mindfulness)
"mindfulness"


#' Meta-analytic correlations among Big Five personality traits and psychological disorders
#'
#' Big Five intercorrelations from Davies et al. (2015). Big Five–psychological
#' disorder correlations from Kotov et al. (2010). Note that there were several
#' duplicate or missing values in the reported data table in the published
#' article. These results are based on corrected data values.
#'
#' @docType data
#'
#' @encoding UTF-8
#'
#' @usage data(disorders)
#'
#' @format list with entries `r` (mean oberseved correlations), `rho` (mean
#' corrected correlations), `n` (sample sizes), `sevar_r` (sampling error
#' variances for mean observed correlations), `sevar_rho` (sampling error
#' variances for mean corrected correlations), and `source` (character labels
#' indicating which meta-analytic correlations came from the same source)
#'
#' @keywords datasets
#'
#' @references
#' Davies, S. E., Connelly, B. L., Ones, D. S., & Birkland, A. S. (2015).
#' The general factor of personality: The “Big One,” a self-evaluative trait, or a methodological gnat that won’t go away?
#' _Personality and Individual Differences, 81_, 13–22. <https://doi.org/10/bc98>
#'
#' Kotov, R., Gamez, W., Schmidt, F., & Watson, D. (2010). Linking “big” personality traits to anxiety, depressive, and substance use disorders: A meta-analysis.
#' _Psychological Bulletin, 136_(5), 768–821. <https://doi.org/10/bxzt2k>
#'
#' @examples
#' data(disorders)
"disorders"


#' Meta-analytic correlations of HRM practices with organizational financial performance
#'
#' Human resource management practice–organizational financial performance
#' correlations from Combs et al. (2006). Intercorrelations among HRM practices
#' from Guest et al. (2004).
#'
#' @docType data
#'
#' @encoding UTF-8
#'
#' @usage data(hrm)
#'
#' @format list with entries `r` (mean oberseved correlations), `rho` (mean
#' corrected correlations), `n` (sample sizes), `sevar_r` (sampling error
#' variances for mean observed correlations), `sevar_rho` (sampling error
#' variances for mean corrected correlations), and `source` (character labels
#' indicating which meta-analytic correlations came from the same source)
#'
#' @keywords datasets
#'
#' @references
#' Combs, J., Liu, Y., Hall, A., & Ketchen, D. (2006).
#' How much do high-performance work practices matter? A meta-analysis of their effects on organizational performance.
#' _Personnel Psychology, 59_(3), 501–528. <https://doi.org/10/cm845n>
#'
#' Guest, D., Conway, N., & Dewe, P. (2004).
#' Using sequential tree analysis to search for ‘bundles’ of HR practices.
#' _Human Resource Management Journal, 14_(1), 79–96. <https://doi.org/10/ffxcfs>

#'
#' @examples
#' data(hrm)
"hrm"


#' Meta-analytic correlations among team processes and team effectiveness
#'
#' Team process intercorrelations and team process–team performance/affect
#' correlations from LePine et al. (2008).
#'
#' Note that LePine et al. (2008) did not report confidence intervals, sampling
#' error variances, or heterogeneity estimates for correlations among team
#' processes; included sampling error variances in this list are based on total
#' sample size only and do not include uncertainty stemmming from any effect
#' size heterogeneity.
#'
#' @docType data
#'
#' @usage data(team)
#'
#' @encoding UTF-8
#'
#' @format list with entries `r` (mean oberseved correlations), `rho` (mean
#' corrected correlations), `n` (sample sizes), `sevar_r` (sampling error
#' variances for mean observed correlations), `sevar_rho` (sampling error
#' variances for mean corrected correlations), and `source` (character labels
#' indicating which meta-analytic correlations came from the same source)
#'
#' @keywords datasets
#'
#' @references
#' LePine, J. A., Piccolo, R. F., Jackson, C. L., Mathieu, J. E., & Saul, J. R. (2008).
#' A meta-analysis of teamwork processes: tests of a multidimensional model and relationships with team effectiveness criteria.
#' _Personnel Psychology, 61_(2), 273–307. <https://doi.org/10/d5nxg9>
#'
#' @examples
#' data(team)
"team"

#' Correlations between study design moderators and effect sizes for prejudice reduction following intergroup contact
#'
#' Correlations among study design moderators and study design
#' moderator–observed prejudice reduction effect sizes from Pettigrew and Tropp
#' (2008). Note that correlations with effect size have been reverse-coded so
#' that a positive correlation indicates that a higher level of the moderator is
#' associated with _larger_ prejudice reduction.
#'
#' @docType data
#'
#' @encoding UTF-8
#'
#' @usage data(prejudice)
#'
#' @format list with entries `r` (observed correlations among moderators) and
#' `k` (number of samples in meta-analysis)
#'
#' @keywords datasets
#'
#' @references
#' Pettigrew, T. F., & Tropp, L. R. (2006).
#' A meta-analytic test of intergroup contact theory.
#' _Journal of Personality and Social Psychology, 90_(5), 751–783. <https://doi.org/10/fmsz26>
#'
#' @examples
#' data(prejudice)
"prejudice"


#' Meta-analytic correlations of job characteristics with performance and satisfaction
#'
#' Self-rated job characteristics intercorrelations and correlations with
#' other-rated job performance and self-rated job satisfaction from Humphrey et
#' al. (2007).
#'
#' @docType data
#'
#' @encoding UTF-8
#'
#' @usage data(jobchar)
#'
#' @format list with entries `r` (mean oberseved correlations), `rho` (mean
#' corrected correlations), `n` (sample sizes), `sevar_r` (sampling error
#' variances for mean observed correlations), `sevar_rho` (sampling error
#' variances for mean corrected correlations), and `source` (character labels
#' indicating which meta-analytic correlations came from the same source)
#'
#' @keywords datasets
#'
#' @references
#' Humphrey, S. E., Nahrgang, J. D., & Morgeson, F. P. (2007).
#' Integrating motivational, social, and contextual work design features: A meta-analytic summary and theoretical extension of the work design literature.
#' _Journal of Applied Psychology, 92_(5), 1332–1356. <https://doi.org/10/dch52t>
#'
#' @examples
#' data(jobchar)
"jobchar"

#' Meta-analytic correlations of Graduate Record Examination subtests with graduate grade point average
#'
#' Correlations between GRE subtests and graduate student GPA from Kuncel et al. (2001).
#'
#' GRE–GPA correlations in `rho` are corrected for direct range restriction on
#' the GRE and unreliability in GPA. Subtest intercorrelations in `rho` are
#' observed correlations computed among applicant norm samples. These values are
#' also used in `r` to permit accurate computation of sampling error in
#' [cor_covariance_meta()]. To conduct analyses with a correlation matrix
#' representing a restricted sample of admitted students, use `r_admitted`. GRE
#' subtest intercorrelations in `r_admitted` are the values in `rho` attenuated
#' using the bivariate direct range restriction formula (Alexander et al., 1987).
#' `sevar_rho` and`sevar_r` for GRE subtest intercorrelations is computed with an assumed
#' \if{latex}{\eqn{SD_\rho}}\ifelse{html}{\out{SD<sub>&rho;</sub>}}{SD_rho} = .02.
#' `sevar_r_admitted` for GRE subtest intercorrelations computed using the estimated
#' values for `r_admitted`, an assumed
#' \if{latex}{\eqn{SD_\rho}}\ifelse{html}{\out{SD<sub>&rho;</sub>}}{SD_rho} = .02,
#' and the N and k values for `rho`.
#'
#' @docType data
#'
#' @encoding UTF-8
#'
#' @usage data(jobchar)
#'
#' @format list with entries `r` (mean oberseved correlations), `rho` (mean
#' corrected correlations), `n` (sample sizes), `sevar_r` (sampling error
#' variances for mean observed correlations), `sevar_rho` (sampling error
#' variances for mean corrected correlations), and `source` (character labels
#' indicating which meta-analytic correlations came from the same source)
#'
#' @keywords datasets
#'
#' @references
#' Alexander, R. A., Carson, K. P., Alliger, G. M., & Carr, L. (1987).
#' Correcting doubly truncated correlations: An improved approximation for correcting the bivariate normal correlation when truncation has occurred on both variables.
#' _Educational and Psychological Measurement, 47_(2), 309–315. <https://doi.org/10.1177/0013164487472002>
#'
#' Kuncel, N. R., Hezlett, S. A., & Ones, D. S. (2001).
#' A comprehensive meta-analysis of the predictive validity of the graduate record examinations: Implications for graduate student selection and performance.
#' _Psychological Bulletin, 127_(1), 162–181. <https://doi.org/10/cffkhc>
#'
#' @examples
#' data(gre)
"gre"
