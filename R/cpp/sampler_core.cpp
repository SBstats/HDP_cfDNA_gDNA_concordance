// sampler_core.cpp
// Rcpp implementations of the three MCMC bottleneck functions:
//   1. compute_allocations_cpp  -- allocation sampling
//   2. update_theta_cpp         -- signature update (3 hierarchical levels)
//   3. compute_loglik_cpp       -- log-likelihood
//
// These are drop-in replacements for the pure-R versions in 01_lib_core.R,
// which are retained as fallbacks.
//
// EQUIVALENCE, PRECISELY STATED. The two implementations are DISTRIBUTIONALLY
// equivalent, NOT stream-identical. Do not expect matching draws from a shared
// seed, and do not treat a mismatch as a bug:
//
//   * update_theta_cpp and compute_loglik_cpp are DETERMINISTIC given their
//     inputs, so these DO agree numerically -- verified to 0.0 (theta, theta_0,
//     pointwise loglik) and ~2e-16 (theta_bar), at inv_temp 1.0 and 0.2.
//
//   * compute_allocations_cpp CANNOT match the R fallback draw-for-draw, for
//     two structural reasons rather than any difference in the math:
//       (i)  R's sample.int(K) consumes K uniforms; the Fisher-Yates loop here
//            draws p-1. The RNG streams desynchronise immediately.
//       (ii) R's sample.int(K, 1, prob=) uses Walker alias / rejection
//            sampling, whereas sample_log_probs() below is an inverse-CDF
//            search. Same seed, same probabilities, different index.
//     Equivalence of the two is therefore established by (a) forcing a
//     degenerate case (K = 1, omega_0 -> 0) in which both are deterministic
//     and agree exactly, and (b) comparing empirical distributions of the
//     allocation counts and sufficient statistics over many paired replicates.
//     Both checks live in R/simulation/run_local_validation.R.
//
// MODEL (see Writeup Draft Sep 14.tex, eq:gdna / eq:cfdna):
//
//   Y^g_{ijt}  ~ sum_k  w^g_{ikt}  N(theta[j,k,i], sigma_k^2)
//   Y^cf_{ijt} ~ (1-w0_{it}) sum_k w^cf_{ikt} N(theta[j,k,i], sigma_k^2)
//                + w0_{it} N(theta_0[j,i], sigma_0^2)
//
// Three properties of this specification drive the layout below.
//   * theta is SHARED across sources and indexed by subject: [p x K x n],
//     linear index j + k*p + i*p*K.
//   * theta_0 is an ESTIMATED background mean [p x n], index j + i*p. It is
//     not fixed at zero; on the VST scale a zero-mean background lies outside
//     the data range and forces w0 -> 0 regardless of the truth.
//   * w0 is per OBSERVATION (length N_obs), not per patient.

#include <Rcpp.h>
using namespace Rcpp;

// --------------------------------------------------------------------------
// Helper: log-sum-exp for two values
// --------------------------------------------------------------------------
inline double log_sum_exp2(double a, double b) {
  double m = std::max(a, b);
  if (!std::isfinite(m)) return -INFINITY;
  return m + std::log(std::exp(a - m) + std::exp(b - m));
}

// --------------------------------------------------------------------------
// Helper: sample one index from unnormalized log-probabilities
// Uses the inverse-CDF method with R's RNG
// --------------------------------------------------------------------------
inline int sample_log_probs(const std::vector<double>& log_p) {
  int n = log_p.size();
  double max_lp = -INFINITY;
  for (int i = 0; i < n; i++) {
    if (log_p[i] > max_lp) max_lp = log_p[i];
  }

  std::vector<double> probs(n);
  double sum_p = 0.0;
  for (int i = 0; i < n; i++) {
    probs[i] = std::exp(log_p[i] - max_lp);
    sum_p += probs[i];
  }

  double u = R::runif(0.0, 1.0) * sum_p;
  double cum = 0.0;
  for (int i = 0; i < n; i++) {
    cum += probs[i];
    if (u <= cum) return i;
  }
  return n - 1; // fallback
}

// --------------------------------------------------------------------------
// 1. compute_allocations_cpp -- COLLAPSED / RAO-BLACKWELLIZED version
//
// theta_vec has dim [p, K, n] (column-major): theta[j,k,i] = theta_vec[j + k*p + i*p*K]
// theta_0_vec has dim [p, n]:                theta_0[j,i]  = theta_0_vec[j + i*p]
// omega0 has length N_obs (per-observation contamination fraction).
//
// Both sources read the SAME theta for a given (j, k, subject), so the
// sufficient statistics for sigma_k^2 and for the signature update pool
// across sources.
//
// Tempering convention (matches writeup Algorithm 1, Step 1):
//   * the Gaussian density factor is raised to 1/T, AND
//   * the Polya-urn running count is divided by T, i.e. the prior factor is
//     (eta_k + n^{(-j)}/T), not (eta_k + n^{(-j)}).
// The anchor eta_k itself is NOT tempered. Tempering the count is what makes
// the collapsed urn the marginal of a tempered Dirichlet-multinomial, so both
// factors see the same likelihood power -- this is the "symmetric" in
// symmetric tempering. At T >> 1 the count term is attenuated and allocations
// are driven mainly by eta, which is the intended flattening of the
// rich-get-richer dynamics during burn-in.
// --------------------------------------------------------------------------
// [[Rcpp::export]]
List compute_allocations_cpp(NumericMatrix Y_g, NumericMatrix Y_cf,
                              NumericVector theta_vec,
                              NumericVector theta_0_vec,
                              NumericVector sigma_k2, double sigma_0_2,
                              NumericVector eta_g,
                              NumericVector eta_cf,
                              NumericVector omega0,
                              IntegerVector patient, IntegerVector time_vec,
                              int K, int p,
                              double temperature = 1.0) {

  int N_obs = Y_g.ncol();
  int n = max(patient);

  // Output matrices
  IntegerMatrix z_g(p, N_obs);
  IntegerMatrix z_cf(p, N_obs);
  IntegerMatrix n_g_counts(K, N_obs);
  IntegerMatrix n_cf_counts(K, N_obs);
  // Background / tumour allocation counts are now PER OBSERVATION, because
  // omega0 is per observation.
  NumericVector N0_vec(N_obs, 0.0);
  NumericVector Nplus_vec(N_obs, 0.0);
  // Residual sum of squares for sigma_0^2, taken about theta_0 (NOT about 0).
  double ss_normal = 0.0;
  NumericVector ss_component(K, 0.0);

  // Precompute
  std::vector<double> log_sigma_k(K);
  std::vector<double> inv_2sigma_k2(K);
  for (int k = 0; k < K; k++) {
    log_sigma_k[k] = 0.5 * std::log(sigma_k2[k]);
    inv_2sigma_k2[k] = 1.0 / (2.0 * sigma_k2[k]);
  }
  double log_sigma_0 = 0.5 * std::log(sigma_0_2);
  double inv_2sigma_0_2 = 1.0 / (2.0 * sigma_0_2);

  double inv_temp = 1.0 / temperature;

  // Working buffer for a random permutation of 1..p (gene scan order).
  std::vector<int> gene_order(p);

  const int pK = p * K;

  for (int obs = 0; obs < N_obs; obs++) {
    int i = patient[obs] - 1;      // 0-indexed subject
    const int theta_off = i * pK;  // subject slice offset into theta_vec
    const int t0_off    = i * p;   // subject slice offset into theta_0_vec

    std::vector<double> n_run_g(K, 0.0);
    std::vector<double> n_run_cf(K, 0.0);

    // Refresh gene scan order (Fisher-Yates with R's RNG).
    for (int j = 0; j < p; j++) gene_order[j] = j;
    for (int j = p - 1; j > 0; j--) {
      int swap_idx = (int)(R::runif(0.0, 1.0) * (j + 1));
      if (swap_idx > j) swap_idx = j;
      std::swap(gene_order[j], gene_order[swap_idx]);
    }

    // --- gDNA allocations: pure mixture, no background component ---
    std::vector<double> lp_g(K);
    for (int idx = 0; idx < p; idx++) {
      int j = gene_order[idx];
      double y_g = Y_g(j, obs);
      for (int k = 0; k < K; k++) {
        double theta_jk = theta_vec[j + k * p + theta_off];
        double resid = y_g - theta_jk;
        double log_prior = std::log(std::max(eta_g[k] + inv_temp * n_run_g[k], 1e-300));
        double log_gauss = (-log_sigma_k[k] - resid * resid * inv_2sigma_k2[k]) * inv_temp;
        lp_g[k] = log_prior + log_gauss;
      }
      int zk = sample_log_probs(lp_g);
      z_g(j, obs) = zk + 1; // 1-indexed for R
      n_run_g[zk] += 1.0;
    }
    for (int k = 0; k < K; k++) n_g_counts(k, obs) = (int) n_run_g[k];

    // --- cfDNA allocations: K tumour components + background at slot 0 ---
    double w0 = omega0[obs];
    double log_w0   = std::log(std::max(w0, 1e-300));
    double log_1mw0 = std::log(std::max(1.0 - w0, 1e-300));

    // Refresh cfDNA gene scan order independently.
    for (int j = 0; j < p; j++) gene_order[j] = j;
    for (int j = p - 1; j > 0; j--) {
      int swap_idx = (int)(R::runif(0.0, 1.0) * (j + 1));
      if (swap_idx > j) swap_idx = j;
      std::swap(gene_order[j], gene_order[swap_idx]);
    }

    std::vector<double> lp_cf(K + 1);
    for (int idx = 0; idx < p; idx++) {
      int j = gene_order[idx];
      double y_cf = Y_cf(j, obs);

      // Background component: residual about the ESTIMATED mean theta_0[j,i].
      double resid0 = y_cf - theta_0_vec[j + t0_off];
      lp_cf[0] = log_w0 + (-log_sigma_0 - resid0 * resid0 * inv_2sigma_0_2) * inv_temp;

      // NORMALISER for the tumour block.
      //
      // The background is a Bernoulli split OUTSIDE the Dirichlet, so its
      // probability must be compared against the tumour block's NORMALISED
      // component probability (eta_k + n_k) / sum_l(eta_l + n_l), not against
      // the raw Polya-urn weight. Omitting this normaliser -- as earlier
      // versions of this sampler did -- hands every tumour component a
      // spurious +log(sum_l(eta_l + n_l)) ~ +log(p) bonus: about +3 log-units
      // at p=40 and +5 at p=310. The background then essentially never wins an
      // allocation and omega_0 collapses to the boundary irrespective of the
      // truth. That artifact, not the biology, is the most likely explanation
      // for the historically "negligible" contamination estimates.
      double urn_tot = 0.0;
      for (int k = 0; k < K; k++) urn_tot += eta_cf[k] + inv_temp * n_run_cf[k];
      double log_urn_tot = std::log(std::max(urn_tot, 1e-300));

      for (int k = 0; k < K; k++) {
        double theta_jk = theta_vec[j + k * p + theta_off];
        double resid = y_cf - theta_jk;
        double log_prior = std::log(std::max(eta_cf[k] + inv_temp * n_run_cf[k], 1e-300))
                           - log_urn_tot;
        double log_gauss = (-log_sigma_k[k] - resid * resid * inv_2sigma_k2[k]) * inv_temp;
        lp_cf[k + 1] = log_1mw0 + log_prior + log_gauss;
      }
      int zk = sample_log_probs(lp_cf); // 0 = background, 1..K = tumour
      z_cf(j, obs) = zk;

      if (zk == 0) {
        N0_vec[obs] += 1.0;
        ss_normal += resid0 * resid0;   // centred residual
      } else {
        n_run_cf[zk - 1] += 1.0;
        Nplus_vec[obs] += 1.0;
      }
    }
    for (int k = 0; k < K; k++) n_cf_counts(k, obs) = (int) n_run_cf[k];
  }

  // Sufficient statistics for sigma_k^2. Both sources' residuals are taken
  // about the same subject-level signature theta[j,k,i].
  for (int k = 0; k < K; k++) {
    double ss_k = 0.0;
    for (int obs = 0; obs < N_obs; obs++) {
      int i = patient[obs] - 1;
      const int theta_off = i * pK;
      for (int j = 0; j < p; j++) {
        if (z_g(j, obs) == k + 1) {
          double resid = Y_g(j, obs) - theta_vec[j + k * p + theta_off];
          ss_k += resid * resid;
        }
        if (z_cf(j, obs) == k + 1) {  // tumour component k (0 = background)
          double resid = Y_cf(j, obs) - theta_vec[j + k * p + theta_off];
          ss_k += resid * resid;
        }
      }
    }
    ss_component[k] = ss_k;
  }

  return List::create(
    Named("n_g") = n_g_counts,
    Named("n_cf") = n_cf_counts,
    Named("N0") = N0_vec,
    Named("Nplus") = Nplus_vec,
    Named("ss_normal") = ss_normal,
    Named("ss_component") = ss_component,
    Named("z_g") = z_g,
    Named("z_cf") = z_cf
  );
}

// --------------------------------------------------------------------------
// 2. update_theta_cpp
//
// Three hierarchical levels, drawn in sequence (writeup Algorithm Step 2a-2c):
//
//   2a. theta[j,k,i] ~ N( mu*, s2* ),
//         s2*^-1 = n_{ijk}/sigma_k^2 + 1/varsigma_k^2
//         mu*    = s2* ( S_{ijk}/sigma_k^2 + theta_bar[j,k]/varsigma_k^2 )
//       Both sources' allocated observations contribute to n_{ijk}, S_{ijk}.
//       When n_{ijk} == 0 the draw reduces to the prior N(theta_bar, varsigma^2)
//       -- i.e. it reverts to the POPULATION value, not to zero.
//
//   2b. theta_bar[j,k] ~ N( mbar*, sbar2* ),
//         sbar2*^-1 = n/varsigma_k^2 + 1/(tau_k^2 lambda_{jk}^2)
//         mbar*     = sbar2* * sum_i theta[j,k,i] / varsigma_k^2
//       Not tempered: conditions on latent parameters, not on data.
//
//   2c. theta_0[j,i] ~ N( mu0*, s02* ),
//         s02*^-1 = N0_{ij}/sigma_0^2 + 1/varsigma_0^2
//         mu0*    = s02* ( S0_{ij}/sigma_0^2 + theta_0_bar[j]/varsigma_0^2 )
//
// inv_temp scales the DATA sufficient statistics only (symmetric tempering).
// --------------------------------------------------------------------------
// [[Rcpp::export]]
List update_theta_cpp(NumericMatrix Y_g, NumericMatrix Y_cf,
                       IntegerMatrix z_g, IntegerMatrix z_cf,
                       NumericVector sigma_k2,
                       NumericVector varsigma_k2,
                       NumericVector theta_bar_vec,
                       NumericVector tau_k2, NumericMatrix lambda2,
                       double sigma_0_2, double varsigma_0_2,
                       NumericVector theta_0_bar,
                       IntegerVector patient,
                       int K, int p, int N_obs, int n,
                       double inv_temp = 1.0) {

  const int pK = p * K;
  NumericVector theta(pK * n, 0.0);
  NumericVector theta_bar_new(pK, 0.0);
  NumericVector theta_0(p * n, 0.0);

  // ---- Accumulate per-(subject, gene, component) sufficient statistics ----
  // Pooled over BOTH sources and all timepoints of that subject.
  std::vector<double> sum_y(pK * (size_t)n, 0.0);
  std::vector<int>    n_ijk(pK * (size_t)n, 0);
  // Background sufficient statistics per (subject, gene).
  std::vector<double> sum_y0(p * (size_t)n, 0.0);
  std::vector<int>    n_ij0(p * (size_t)n, 0);

  for (int obs = 0; obs < N_obs; obs++) {
    int i = patient[obs] - 1;
    const int theta_off = i * pK;
    const int t0_off    = i * p;
    for (int j = 0; j < p; j++) {
      int kg = z_g(j, obs);          // 1..K
      if (kg >= 1) {
        size_t ix = (size_t)(j + (kg - 1) * p + theta_off);
        sum_y[ix] += Y_g(j, obs);
        n_ijk[ix]++;
      }
      int kc = z_cf(j, obs);         // 0 = background, 1..K = tumour
      if (kc == 0) {
        size_t ix0 = (size_t)(j + t0_off);
        sum_y0[ix0] += Y_cf(j, obs);
        n_ij0[ix0]++;
      } else {
        size_t ix = (size_t)(j + (kc - 1) * p + theta_off);
        sum_y[ix] += Y_cf(j, obs);
        n_ijk[ix]++;
      }
    }
  }

  // ---- Step 2a: subject-level signatures ----
  for (int k = 0; k < K; k++) {
    double sig2_k  = sigma_k2[k];
    double vars2_k = varsigma_k2[k];
    double prior_prec_sub = 1.0 / vars2_k;
    for (int i = 0; i < n; i++) {
      const int theta_off = i * pK;
      for (int j = 0; j < p; j++) {
        size_t ix = (size_t)(j + k * p + theta_off);
        double tb = theta_bar_vec[j + k * p];
        double data_prec = ((double)n_ijk[ix] * inv_temp) / sig2_k;
        double post_prec = data_prec + prior_prec_sub;
        double post_var  = 1.0 / post_prec;
        // When n_ijk == 0 this collapses to the prior mean theta_bar.
        double post_mean = post_var * ((sum_y[ix] * inv_temp) / sig2_k
                                       + tb * prior_prec_sub);
        theta[ix] = R::rnorm(post_mean, std::sqrt(post_var));
      }
    }
  }

  // ---- Step 2b: population signatures (prior-side; not tempered) ----
  for (int k = 0; k < K; k++) {
    double vars2_k = varsigma_k2[k];
    double tau2_k  = tau_k2[k];
    for (int j = 0; j < p; j++) {
      double sum_sub = 0.0;
      for (int i = 0; i < n; i++) sum_sub += theta[(size_t)(j + k * p + i * pK)];
      double prior_prec = 1.0 / (tau2_k * lambda2(j, k));
      double post_prec  = (double)n / vars2_k + prior_prec;
      double post_var   = 1.0 / post_prec;
      double post_mean  = post_var * (sum_sub / vars2_k);
      theta_bar_new[j + k * p] = R::rnorm(post_mean, std::sqrt(post_var));
    }
  }

  // ---- Step 2c: subject-level background means ----
  double prior_prec_0 = 1.0 / varsigma_0_2;
  for (int i = 0; i < n; i++) {
    const int t0_off = i * p;
    for (int j = 0; j < p; j++) {
      size_t ix0 = (size_t)(j + t0_off);
      double data_prec = ((double)n_ij0[ix0] * inv_temp) / sigma_0_2;
      double post_prec = data_prec + prior_prec_0;
      double post_var  = 1.0 / post_prec;
      double post_mean = post_var * ((sum_y0[ix0] * inv_temp) / sigma_0_2
                                     + theta_0_bar[j] * prior_prec_0);
      theta_0[ix0] = R::rnorm(post_mean, std::sqrt(post_var));
    }
  }

  theta.attr("dim")         = IntegerVector::create(p, K, n);
  theta_bar_new.attr("dim") = IntegerVector::create(p, K);
  theta_0.attr("dim")       = IntegerVector::create(p, n);

  return List::create(
    Named("theta")     = theta,
    Named("theta_bar") = theta_bar_new,
    Named("theta_0")   = theta_0
  );
}

// --------------------------------------------------------------------------
// 3. compute_loglik_cpp
//
// Observed-data log-likelihood under eq:gdna / eq:cfdna. gDNA is a pure
// mixture; cfDNA adds the background component with mean theta_0[j,i].
//
// Pointwise layout is [2 * N_obs]: entries 0..N_obs-1 are gDNA (summed over
// genes), N_obs..2*N_obs-1 are cfDNA. Unchanged from the previous version so
// that WAIC/LOO consumers need no modification.
// --------------------------------------------------------------------------
// [[Rcpp::export]]
List compute_loglik_cpp(NumericMatrix Y_g, NumericMatrix Y_cf,
                         NumericVector theta_vec,
                         NumericVector theta_0_vec,
                         NumericVector sigma_k2, double sigma_0_2,
                         NumericMatrix omega_g, NumericMatrix omega_cf,
                         NumericVector omega0,
                         IntegerVector patient, IntegerVector time_vec,
                         int K, int p,
                         bool return_pointwise) {

  int N_obs = Y_g.ncol();
  double ll_total = 0.0;

  NumericVector pw;
  if (return_pointwise) {
    pw = NumericVector(2 * N_obs, 0.0);
  }

  std::vector<double> sqrt_sigma_k(K);
  double sqrt_sigma_0 = std::sqrt(sigma_0_2);
  for (int k = 0; k < K; k++) {
    sqrt_sigma_k[k] = std::sqrt(sigma_k2[k]);
  }

  const int pK = p * K;

  for (int obs = 0; obs < N_obs; obs++) {
    int i = patient[obs] - 1;
    const int theta_off = i * pK;
    const int t0_off    = i * p;
    double ll_g = 0.0;
    double ll_cf = 0.0;

    double w0 = omega0[obs];
    double log_w0   = std::log(std::max(w0, 1e-300));
    double log_1mw0 = std::log(std::max(1.0 - w0, 1e-300));

    for (int j = 0; j < p; j++) {
      // gDNA: pure mixture over the subject's own signatures
      double log_mix_g = -INFINITY;
      for (int k = 0; k < K; k++) {
        double theta_jk = theta_vec[j + k * p + theta_off];
        double log_om = std::log(std::max(omega_g(k, obs), 1e-300));
        double ld = R::dnorm(Y_g(j, obs), theta_jk, sqrt_sigma_k[k], 1);
        log_mix_g = log_sum_exp2(log_mix_g, log_om + ld);
      }
      ll_g += log_mix_g;

      // cfDNA: background (estimated mean) + tumour mixture
      double log_normal = log_w0 + R::dnorm(Y_cf(j, obs),
                                            theta_0_vec[j + t0_off],
                                            sqrt_sigma_0, 1);
      double log_tumor = -INFINITY;
      for (int k = 0; k < K; k++) {
        double theta_jk = theta_vec[j + k * p + theta_off];
        double log_om = std::log(std::max(omega_cf(k, obs), 1e-300));
        double ld = R::dnorm(Y_cf(j, obs), theta_jk, sqrt_sigma_k[k], 1);
        log_tumor = log_sum_exp2(log_tumor, log_1mw0 + log_om + ld);
      }
      ll_cf += log_sum_exp2(log_normal, log_tumor);
    }

    ll_total += ll_g + ll_cf;

    if (return_pointwise) {
      pw[obs] = ll_g;
      pw[N_obs + obs] = ll_cf;
    }
  }

  List result;
  result["total"] = ll_total;
  if (return_pointwise) {
    result["pointwise"] = pw;
  }
  return result;
}
