set.seed(123)
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
})

N_ITER  <- 2000
N_POST  <-  500

posterior_prob_A_superior <- function(successA, countA, successB, countB,
                                      alpha0 = 1, beta0 = 1,
                                      n_samples = N_POST) {
  sA <- rbeta(n_samples, alpha0 + successA, beta0 + (countA - successA))
  sB <- rbeta(n_samples, alpha0 + successB, beta0 + (countB - successB))
  mean(sA > sB)
}

mc_se <- function(p_hat, B) sqrt(p_hat * (1 - p_hat) / B)

run_prop_test <- function(assigned, observed, alternative = "two.sided") {
  tab <- table(factor(assigned, levels = c(0, 1)),
               factor(observed, levels = c(0, 1)))
  if (any(rowSums(tab) == 0)) return(1)
  prop.test(tab[, "1"], rowSums(tab),
            alternative = alternative,
            correct = FALSE)$p.value
}

.summarise_oc <- function(reject_null, alloc_a_prop, total_successes,
                          iterations, method = NULL) {
  n_rej      <- sum(reject_null)
  p_hat      <- mean(reject_null)
  alloc_mean <- mean(alloc_a_prop)
  alloc_sd   <- sd(alloc_a_prop)
  alloc_ci   <- alloc_mean + c(-1.96, 1.96) * alloc_sd / sqrt(iterations)
  
  if (!is.null(method) && identical(method, "Fixed 1:1")) {
    alloc_sd <- NA_real_
    alloc_ci <- c(NA_real_, NA_real_)
  }
  
  list(
    typeI_or_power    = p_hat,
    mc_se             = mc_se(p_hat, iterations),
    typeI_or_power_ci = prop.test(n_rej, iterations)$conf.int,
    alloc_mean        = alloc_mean,
    alloc_sd          = alloc_sd,
    alloc_ci          = alloc_ci,
    successes_mean    = mean(total_successes),
    alloc_vals        = alloc_a_prop
  )
}

summarise_run <- function(res, scenario, method) {
  data.frame(
    Scenario       = scenario,
    Method         = method,
    Power_or_TypeI = round(res$typeI_or_power, 3),
    MC_SE          = round(res$mc_se, 4),
    CI_L           = round(res$typeI_or_power_ci[1], 3),
    CI_U           = round(res$typeI_or_power_ci[2], 3),
    Alloc_to_A     = round(res$alloc_mean, 3),
    Alloc_SD       = ifelse(is.na(res$alloc_sd), NA, round(res$alloc_sd, 4)),
    Alloc_CI_L     = ifelse(is.na(res$alloc_ci[1]), NA, round(res$alloc_ci[1], 3)),
    Alloc_CI_U     = ifelse(is.na(res$alloc_ci[2]), NA, round(res$alloc_ci[2], 3)),
    Mean_Successes = round(res$successes_mean, 2),
    stringsAsFactors = FALSE
  )
}

# 1) Fixed 1:1
run_fixed_sim <- function(n_total, pa, pb,
                          iterations = N_ITER,
                          alternative = "two.sided") {
  if (n_total %% 2 != 0) stop("n_total must be even for fixed 1:1.")
  reject_null     <- integer(iterations)
  alloc_a_prop    <- numeric(iterations)
  total_successes <- numeric(iterations)
  
  half <- n_total / 2L
  for (i in seq_len(iterations)) {
    assigned <- sample(c(rep(1L, half), rep(0L, half)))
    p_vec    <- ifelse(assigned == 1L, pa, pb)
    observed <- rbinom(n_total, 1L, p_vec)
    
    reject_null[i]     <- as.integer(
      run_prop_test(assigned, observed, alternative) < 0.05)
    alloc_a_prop[i]    <- mean(assigned)
    total_successes[i] <- sum(observed)
  }
  .summarise_oc(reject_null, alloc_a_prop, total_successes, iterations, "Fixed 1:1")
}

# 2) Simple RAR
run_rar_sim <- function(n_total, pa, pb,
                        iterations  = N_ITER,
                        alternative = "two.sided",
                        p_up = 0.6, p_down = 0.4) {
  reject_null     <- integer(iterations)
  alloc_a_prop    <- numeric(iterations)
  total_successes <- numeric(iterations)
  
  for (i in seq_len(iterations)) {
    assigned <- integer(n_total)
    observed <- integer(n_total)
    sA <- 0L; nA <- 0L; sB <- 0L; nB <- 0L
    
    for (j in seq_len(n_total)) {
      p_A <- if (nA == 0L || nB == 0L) 0.5 else {
        rA <- (sA + 0.5) / (nA + 1)
        rB <- (sB + 0.5) / (nB + 1)
        if (rA > rB) p_up else if (rB > rA) p_down else 0.5
      }
      assigned[j] <- rbinom(1L, 1L, p_A)
      observed[j] <- rbinom(1L, 1L, if (assigned[j] == 1L) pa else pb)
      if (assigned[j] == 1L) { nA <- nA + 1L; sA <- sA + observed[j]
      } else                  { nB <- nB + 1L; sB <- sB + observed[j] }
    }
    
    reject_null[i]     <- as.integer(
      run_prop_test(assigned, observed, alternative) < 0.05)
    alloc_a_prop[i]    <- mean(assigned)
    total_successes[i] <- sum(observed)
  }
  .summarise_oc(reject_null, alloc_a_prop, total_successes, iterations, "Simple RAR")
}

# 3) Bayesian RAR
run_brar_sim <- function(n_total, pa, pb,
                         iterations          = N_ITER,
                         alternative         = "two.sided",
                         alpha0 = 1, beta0 = 1,
                         burn_in             = 10,
                         n_posterior_samples = N_POST,
                         clip_lo = 0.1, clip_hi = 0.9,
                         store_trace = FALSE) {
  reject_null     <- integer(iterations)
  alloc_a_prop    <- numeric(iterations)
  total_successes <- numeric(iterations)
  posterior_trace <- if (store_trace) matrix(NA_real_, iterations, n_total) else NULL
  
  for (i in seq_len(iterations)) {
    assigned <- integer(n_total)
    observed <- integer(n_total)
    sA <- 0L; nA <- 0L; sB <- 0L; nB <- 0L
    
    for (j in seq_len(n_total)) {
      p_A <- if (j <= burn_in || nA == 0L || nB == 0L) 0.5 else {
        p_sup <- posterior_prob_A_superior(sA, nA, sB, nB,
                                           alpha0 = alpha0, beta0 = beta0,
                                           n_samples = n_posterior_samples)
        pmax(clip_lo, pmin(clip_hi, p_sup))
      }
      if (store_trace) posterior_trace[i, j] <- p_A
      
      assigned[j] <- rbinom(1L, 1L, p_A)
      observed[j] <- rbinom(1L, 1L, if (assigned[j] == 1L) pa else pb)
      if (assigned[j] == 1L) { nA <- nA + 1L; sA <- sA + observed[j]
      } else                  { nB <- nB + 1L; sB <- sB + observed[j] }
    }
    
    reject_null[i]     <- as.integer(
      run_prop_test(assigned, observed, alternative) < 0.05)
    alloc_a_prop[i]    <- mean(assigned)
    total_successes[i] <- sum(observed)
  }
  
  oc <- .summarise_oc(reject_null, alloc_a_prop, total_successes,
                      iterations, "Bayesian RAR")
  oc$posterior_trace <- posterior_trace
  oc
}

# 4) Cluster helpers
.beta_params_icc <- function(p, rho) {
  if (rho <= 0) return(list(a = Inf, b = Inf))
  phi <- (1 - rho) / rho
  list(a = p * phi, b = (1 - p) * phi)
}

run_cluster_sim <- function(m_clusters, k_per_clust, pa, pb,
                            icc = 0.05, iterations = N_ITER) {
  bpA <- .beta_params_icc(pa, icc)
  bpB <- .beta_params_icc(pb, icc)
  reject_null  <- integer(iterations)
  alloc_a_prop <- numeric(iterations)
  total_succ   <- numeric(iterations)
  
  for (i in seq_len(iterations)) {
    muA   <- if (icc <= 0) rep(pa, m_clusters) else rbeta(m_clusters, bpA$a, bpA$b)
    muB   <- if (icc <= 0) rep(pb, m_clusters) else rbeta(m_clusters, bpB$a, bpB$b)
    propA <- vapply(muA, function(mu) mean(rbinom(k_per_clust, 1L, mu)), numeric(1))
    propB <- vapply(muB, function(mu) mean(rbinom(k_per_clust, 1L, mu)), numeric(1))
    
    total_succ[i]   <- round(sum(c(propA, propB)) * k_per_clust)
    alloc_a_prop[i] <- 0.5
    p_val <- tryCatch(t.test(propA, propB, var.equal = FALSE)$p.value,
                      error = function(e) NA_real_)
    reject_null[i] <- as.integer(!is.na(p_val) && p_val < 0.05)
  }
  
  deff <- 1 + (k_per_clust - 1) * icc
  oc   <- .summarise_oc(reject_null, alloc_a_prop, total_succ,
                        iterations, "Cluster Fixed")
  oc$design_effect <- round(deff, 3)
  oc$effective_n   <- round((m_clusters * k_per_clust * 2) / deff, 1)
  oc
}

run_cluster_brar_sim <- function(m_clusters, k_per_clust, pa, pb,
                                 icc = 0.05, iterations = N_ITER,
                                 burn_in_clusters     = 4,
                                 n_posterior_samples  = N_POST,
                                 clip_lo = 0.1, clip_hi = 0.9) {
  bpA <- .beta_params_icc(pa, icc)
  bpB <- .beta_params_icc(pb, icc)
  reject_null  <- integer(iterations)
  alloc_a_prop <- numeric(iterations)
  total_succ   <- numeric(iterations)
  
  for (i in seq_len(iterations)) {
    arm_labels <- integer(m_clusters)
    clust_prop <- numeric(m_clusters)
    sA <- 0L; nA <- 0L; sB <- 0L; nB <- 0L
    
    for (cl in seq_len(m_clusters)) {
      p_A <- if (cl <= burn_in_clusters || nA == 0L || nB == 0L) 0.5 else {
        p_sup <- posterior_prob_A_superior(sA, nA, sB, nB,
                                           n_samples = n_posterior_samples)
        pmax(clip_lo, pmin(clip_hi, p_sup))
      }
      arm_labels[cl] <- rbinom(1L, 1L, p_A)
      true_p <- if (arm_labels[cl] == 1L) {
        if (icc <= 0) pa else rbeta(1L, bpA$a, bpA$b)
      } else {
        if (icc <= 0) pb else rbeta(1L, bpB$a, bpB$b)
      }
      y_cl           <- rbinom(1L, size = k_per_clust, prob = true_p)
      clust_prop[cl] <- y_cl / k_per_clust
      
      if (arm_labels[cl] == 1L) { sA <- sA + y_cl; nA <- nA + k_per_clust
      } else                     { sB <- sB + y_cl; nB <- nB + k_per_clust }
    }
    
    propA <- clust_prop[arm_labels == 1L]
    propB <- clust_prop[arm_labels == 0L]
    alloc_a_prop[i] <- mean(arm_labels)
    total_succ[i]   <- sum(round(clust_prop * k_per_clust))
    p_val <- tryCatch(
      if (length(propA) < 2 || length(propB) < 2) NA_real_
      else t.test(propA, propB, var.equal = FALSE)$p.value,
      error = function(e) NA_real_)
    reject_null[i] <- as.integer(!is.na(p_val) && p_val < 0.05)
  }
  .summarise_oc(reject_null, alloc_a_prop, total_succ, iterations, "Cluster-Adapt")
}

# 5) Covariate-adjusted RAR
run_ca_rar_sim <- function(n_total, pa_base, pb_base,
                           beta_x      = 0.5,
                           iterations  = N_ITER,
                           alternative = "two.sided",
                           burn_in     = 20,
                           n_posterior_samples = N_POST,
                           clip_lo = 0.1, clip_hi = 0.9,
                           refit_every = 5) {
  logit     <- function(p) log(p / (1 - p))
  inv_logit <- function(x) 1 / (1 + exp(-x))
  b0_A <- logit(pa_base)
  b0_B <- logit(pb_base)
  
  reject_null     <- integer(iterations)
  alloc_a_prop    <- numeric(iterations)
  total_successes <- numeric(iterations)
  
  for (i in seq_len(iterations)) {
    assigned <- integer(n_total)
    observed <- integer(n_total)
    covar    <- rbinom(n_total, 1L, 0.5)
    sA <- 0L; nA <- 0L; sB <- 0L; nB <- 0L
    cached_p_A <- 0.5
    
    for (j in seq_len(n_total)) {
      xj <- covar[j]
      
      if (j <= burn_in || nA < 5L || nB < 5L) {
        p_A <- 0.5
      } else if ((j - burn_in) %% refit_every == 1L) {
        hist_df <- data.frame(Y   = observed[seq_len(j - 1L)],
                              arm = assigned[seq_len(j - 1L)],
                              X   = covar[seq_len(j - 1L)])
        fit <- tryCatch(glm(Y ~ arm + X, data = hist_df, family = binomial()),
                        error = function(e) NULL)
        if (is.null(fit)) {
          p_A <- cached_p_A
        } else {
          cf <- coef(fit)
          b_int <- unname(cf["(Intercept)"])
          b_arm <- unname(cf["arm"])
          b_x   <- unname(cf["X"])
          if (anyNA(c(b_int, b_arm, b_x))) {
            p_A <- cached_p_A
          } else {
            eta_A   <- b_int + b_arm * 1 + b_x * xj
            eta_B   <- b_int + b_arm * 0 + b_x * xj
            p_adj_A <- inv_logit(eta_A)
            p_adj_B <- inv_logit(eta_B)
            eff_nA  <- max(nA, 1L); eff_nB <- max(nB, 1L)
            eff_sA  <- round(p_adj_A * eff_nA)
            eff_sB  <- round(p_adj_B * eff_nB)
            p_sup   <- posterior_prob_A_superior(eff_sA, eff_nA, eff_sB, eff_nB,
                                                 n_samples = n_posterior_samples)
            p_A <- pmax(clip_lo, pmin(clip_hi, p_sup))
            cached_p_A <- p_A
          }
        }
      } else {
        p_A <- cached_p_A
      }
      
      assigned[j] <- rbinom(1L, 1L, p_A)
      true_p      <- inv_logit(if (assigned[j] == 1L) b0_A + beta_x * xj
                               else                   b0_B + beta_x * xj)
      observed[j] <- rbinom(1L, 1L, true_p)
      if (assigned[j] == 1L) { nA <- nA + 1L; sA <- sA + observed[j]
      } else                  { nB <- nB + 1L; sB <- sB + observed[j] }
    }
    
    alloc_a_prop[i]    <- mean(assigned)
    total_successes[i] <- sum(observed)
    reject_null[i]     <- as.integer(
      run_prop_test(assigned, observed, alternative) < 0.05)
  }
  
  .summarise_oc(reject_null, alloc_a_prop, total_successes, iterations, "CA-RAR")
}

# 6) Main suite
run_suite <- function(n_total    = 100,
                      pb         = 0.5,
                      pa_values  = c(0.5, 0.55, 0.65, 0.8),
                      iterations = N_ITER) {
  out <- vector("list", length(pa_values) * 3L)
  k   <- 0L
  for (pa in pa_values) {
    sc <- sprintf("pa=%.2f, pb=%.2f, n=%d", pa, pb, n_total)
    cat("Running:", sc, "\n")
    fixed <- run_fixed_sim(n_total, pa, pb, iterations)
    rar   <- run_rar_sim  (n_total, pa, pb, iterations)
    brar  <- run_brar_sim (n_total, pa, pb, iterations, store_trace = FALSE)
    k <- k + 1L; out[[k]] <- summarise_run(fixed, sc, "Fixed 1:1")
    k <- k + 1L; out[[k]] <- summarise_run(rar,   sc, "Simple RAR")
    k <- k + 1L; out[[k]] <- summarise_run(brar,  sc, "Bayesian RAR")
  }
  do.call(rbind, out)
}

# Run simulations
cat("\n=== Section 6: Main suite n=100 ===\n")
results_100 <- run_suite(n_total = 100, pb = 0.5)

cat("\n=== Section 6.5: Sensitivity n=500 ===\n")
results_500 <- run_suite(n_total = 500, pb = 0.5)

cat("\n=== Section 7.2: Cluster fixed - varying ICC ===\n")
icc_vals <- c(0.00, 0.02, 0.05, 0.10, 0.20)
cluster_fixed_df <- do.call(rbind, lapply(icc_vals, function(rho) {
  cat(" ICC =", rho, "\n")
  res <- run_cluster_sim(20, 10, 0.65, 0.5, icc = rho)
  data.frame(ICC = rho, Method = "Cluster Fixed",
             Power         = round(res$typeI_or_power, 3),
             MC_SE         = round(res$mc_se, 4),
             CI_L          = round(res$typeI_or_power_ci[1], 3),
             CI_U          = round(res$typeI_or_power_ci[2], 3),
             Alloc_to_A    = round(res$alloc_mean, 3),
             Design_Effect = round(res$design_effect, 3),
             Effective_N   = res$effective_n)
}))

cat("\n=== Section 7.3: Cluster adaptive vs fixed - varying ICC ===\n")
cluster_brar_df <- do.call(rbind, lapply(icc_vals, function(rho) {
  cat(" ICC =", rho, "\n")
  res <- run_cluster_brar_sim(20, 10, 0.65, 0.5, icc = rho)
  data.frame(ICC = rho, Method = "Cluster-Adapt",
             Power         = round(res$typeI_or_power, 3),
             MC_SE         = round(res$mc_se, 4),
             CI_L          = round(res$typeI_or_power_ci[1], 3),
             CI_U          = round(res$typeI_or_power_ci[2], 3),
             Alloc_to_A    = round(res$alloc_mean, 3),
             Design_Effect = NA_real_,
             Effective_N   = NA_real_)
}))
cluster_all_df <- rbind(cluster_fixed_df, cluster_brar_df)

cat("\n=== Section 8.2: CA-RAR vs Fixed 1:1 ===\n")
ca_results <- rbind(
  summarise_run(run_fixed_sim (100, 0.5,  0.5),              "H0 (pa=0.50)", "Fixed 1:1"),
  summarise_run(run_ca_rar_sim(100, 0.5,  0.5, beta_x = 0.5),"H0 (pa=0.50)", "CA-RAR"),
  summarise_run(run_fixed_sim (100, 0.65, 0.5),              "H1 (pa=0.65)", "Fixed 1:1"),
  summarise_run(run_ca_rar_sim(100, 0.65, 0.5, beta_x = 0.5),"H1 (pa=0.65)", "CA-RAR")
)

cat("\n=== Section 8.3: CA-RAR - varying beta_x ===\n")
beta_x_vals     <- c(0.0, 0.3, 0.6, 1.0, 1.5)
fixed_betax_res <- run_fixed_sim(100, 0.65, 0.5)
ca_betax_df <- do.call(rbind, lapply(beta_x_vals, function(bx) {
  cat(" beta_x =", bx, "\n")
  ca_res <- run_ca_rar_sim(100, 0.65, 0.5, beta_x = bx)
  rbind(
    data.frame(beta_x = bx, Method = "Fixed 1:1",
               Power = round(fixed_betax_res$typeI_or_power, 3),
               MC_SE = round(fixed_betax_res$mc_se, 4),
               Alloc = round(fixed_betax_res$alloc_mean, 3)),
    data.frame(beta_x = bx, Method = "CA-RAR",
               Power = round(ca_res$typeI_or_power, 3),
               MC_SE = round(ca_res$mc_se, 4),
               Alloc = round(ca_res$alloc_mean, 3))
  )
}))

cat("\n=== Section 9: Building discussion table ===\n")
.get_suite_res <- function(df, pa_str, method) {
  row <- df[grepl(pa_str, df$Scenario) & df$Method == method, ]
  list(typeI_or_power = row$Power_or_TypeI,
       alloc_mean     = row$Alloc_to_A)
}

cache_carar_H0   <- run_ca_rar_sim(100, 0.5,  0.5, beta_x = 0.5)
cache_carar_H165 <- run_ca_rar_sim(100, 0.65, 0.5, beta_x = 0.5)
cache_carar_H180 <- run_ca_rar_sim(100, 0.8,  0.5, beta_x = 0.5)

make_disc_row <- function(method, scenario, type1, power, alloc) {
  data.frame(Method = method, Scenario = scenario,
             Type_I = type1, Power = power,
             Ethical_Alloc = alloc, stringsAsFactors = FALSE)
}
.r <- function(x) round(x, 3)

discussion_table <- rbind(
  make_disc_row("Fixed 1:1",    "H0",
                .r(.get_suite_res(results_100,"pa=0.50","Fixed 1:1")$typeI_or_power),
                NA, .r(.get_suite_res(results_100,"pa=0.50","Fixed 1:1")$alloc_mean)),
  make_disc_row("Simple RAR",   "H0",
                .r(.get_suite_res(results_100,"pa=0.50","Simple RAR")$typeI_or_power),
                NA, .r(.get_suite_res(results_100,"pa=0.50","Simple RAR")$alloc_mean)),
  make_disc_row("Bayesian RAR", "H0",
                .r(.get_suite_res(results_100,"pa=0.50","Bayesian RAR")$typeI_or_power),
                NA, .r(.get_suite_res(results_100,"pa=0.50","Bayesian RAR")$alloc_mean)),
  make_disc_row("CA-RAR",       "H0",
                .r(cache_carar_H0$typeI_or_power), NA, .r(cache_carar_H0$alloc_mean)),
  make_disc_row("Fixed 1:1",    "H1 moderate (pa=0.65)", NA,
                .r(.get_suite_res(results_100,"pa=0.65","Fixed 1:1")$typeI_or_power),
                .r(.get_suite_res(results_100,"pa=0.65","Fixed 1:1")$alloc_mean)),
  make_disc_row("Simple RAR",   "H1 moderate (pa=0.65)", NA,
                .r(.get_suite_res(results_100,"pa=0.65","Simple RAR")$typeI_or_power),
                .r(.get_suite_res(results_100,"pa=0.65","Simple RAR")$alloc_mean)),
  make_disc_row("Bayesian RAR", "H1 moderate (pa=0.65)", NA,
                .r(.get_suite_res(results_100,"pa=0.65","Bayesian RAR")$typeI_or_power),
                .r(.get_suite_res(results_100,"pa=0.65","Bayesian RAR")$alloc_mean)),
  make_disc_row("CA-RAR",       "H1 moderate (pa=0.65)", NA,
                .r(cache_carar_H165$typeI_or_power), .r(cache_carar_H165$alloc_mean)),
  make_disc_row("Fixed 1:1",    "H1 large (pa=0.80)", NA,
                .r(.get_suite_res(results_100,"pa=0.80","Fixed 1:1")$typeI_or_power),
                .r(.get_suite_res(results_100,"pa=0.80","Fixed 1:1")$alloc_mean)),
  make_disc_row("Simple RAR",   "H1 large (pa=0.80)", NA,
                .r(.get_suite_res(results_100,"pa=0.80","Simple RAR")$typeI_or_power),
                .r(.get_suite_res(results_100,"pa=0.80","Simple RAR")$alloc_mean)),
  make_disc_row("Bayesian RAR", "H1 large (pa=0.80)", NA,
                .r(.get_suite_res(results_100,"pa=0.80","Bayesian RAR")$typeI_or_power),
                .r(.get_suite_res(results_100,"pa=0.80","Bayesian RAR")$alloc_mean)),
  make_disc_row("CA-RAR",       "H1 large (pa=0.80)", NA,
                .r(cache_carar_H180$typeI_or_power), .r(cache_carar_H180$alloc_mean))
)

cat("\n=== Posterior trace for Figure 3 ===\n")
h1_brar_trace <- run_brar_sim(100, 0.8, 0.5, store_trace = TRUE)
df_trace <- data.frame(
  Participant = seq_len(100),
  Median = apply(h1_brar_trace$posterior_trace, 2, median,   na.rm = TRUE),
  Q25    = apply(h1_brar_trace$posterior_trace, 2, quantile, 0.25, na.rm = TRUE),
  Q75    = apply(h1_brar_trace$posterior_trace, 2, quantile, 0.75, na.rm = TRUE)
)

cat("\n=== Allocation data for Figure 2 ===\n")
h1_fixed <- run_fixed_sim(100, 0.8, 0.5)
h1_rar   <- run_rar_sim  (100, 0.8, 0.5)
h1_brar  <- run_brar_sim (100, 0.8, 0.5, store_trace = FALSE)
df_alloc <- data.frame(
  Method = factor(rep(c("Fixed 1:1", "Simple RAR", "Bayesian RAR"), each = N_ITER),
                  levels = c("Fixed 1:1", "Simple RAR", "Bayesian RAR")),
  Proportion = c(h1_fixed$alloc_vals, h1_rar$alloc_vals, h1_brar$alloc_vals)
)

cat("\n--- Table 1: n=100 ---\n");             print(results_100,      row.names = FALSE)
cat("\n--- Table 2: n=500 Sensitivity ---\n"); print(results_500,      row.names = FALSE)
cat("\n--- Table 3: Cluster Fixed ---\n");     print(cluster_fixed_df, row.names = FALSE)
cat("\n--- Table 4: Cluster Adaptive ---\n");  print(cluster_all_df,   row.names = FALSE)
cat("\n--- Table 5: CA-RAR vs Fixed ---\n");   print(ca_results,       row.names = FALSE)
cat("\n--- Table 6: CA-RAR beta_x ---\n");     print(ca_betax_df,      row.names = FALSE)
cat("\n--- Table 7: Discussion ---\n");        print(discussion_table, row.names = FALSE)

# ── Colours ───────────────────────────────────────────────────────────────────
COLS <- c(
  "Fixed 1:1"     = "firebrick",
  "Simple RAR"    = "steelblue",
  "Bayesian RAR"  = "seagreen",
  "CA-RAR"        = "mediumpurple",
  "Cluster Fixed" = "darkorange",
  "Cluster-Adapt" = "goldenrod"
)

BASE_THEME <- theme_minimal(base_size = 16) +
  theme(
    plot.title       = element_text(face = "bold", size = 16),
    plot.subtitle    = element_text(colour = "grey45", size = 13),
    legend.position  = "bottom",
    legend.text      = element_text(size = 13),
    legend.title     = element_text(size = 14),
    axis.text        = element_text(size = 13),
    axis.title       = element_text(size = 14),
    panel.grid.minor = element_blank()
  )

# Figure 1
power_df <- results_100 %>%
  filter(!grepl("pa=0.50", Scenario)) %>%
  mutate(pa     = as.numeric(sub("pa=(.*?),.*", "\\1", Scenario)),
         Method = factor(Method, levels = c("Fixed 1:1", "Simple RAR", "Bayesian RAR")))

fig1 <- ggplot(power_df, aes(pa, Power_or_TypeI, colour = Method, group = Method)) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 1.1) + geom_point(size = 3) +
  geom_errorbar(aes(ymin = CI_L, ymax = CI_U), width = 0.012, alpha = 0.6) +
  scale_colour_manual(values = COLS) +
  scale_y_continuous(labels = percent_format(1), limits = c(0, 1)) +
  BASE_THEME +
  labs(title    = "Figure 1 (Section 6.2): Power Curves by Method (n=100, pb=0.50)",
       subtitle = "Error bars = 95% Monte Carlo CI | dashed = 80% conventional target",
       x = "True success probability of Arm A (pa)",
       y = "Estimated Power", colour = "Method")
ggsave("fig1_power_curves.png", fig1, width = 8, height = 5, dpi = 300)

# Figure 2
fig2 <- ggplot(df_alloc, aes(x = Proportion, fill = Method)) +
  geom_histogram(alpha = 0.55, position = "identity", bins = 35) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey30") +
  scale_fill_manual(values = COLS) +
  scale_x_continuous(labels = percent_format(1)) +
  BASE_THEME +
  labs(title    = "Figure 2 (Section 6.3): Allocation to Arm A Under H\u2081",
       subtitle = "pa=0.80, pb=0.50, n=100 | dashed = 50% equal allocation",
       x = "Proportion allocated to Arm A", y = "Frequency", fill = "Method")
ggsave("fig2_allocation.png", fig2, width = 8, height = 5, dpi = 300)

# Figure 3
fig3 <- ggplot(df_trace, aes(x = Participant)) +
  geom_ribbon(aes(ymin = Q25, ymax = Q75), fill = "seagreen", alpha = 0.25) +
  geom_line(aes(y = Median), colour = "darkgreen", linewidth = 1) +
  geom_hline(yintercept = c(0.5, 0.9),
             linetype = c("dashed", "dotted"), colour = "grey40", linewidth = 0.6) +
  geom_vline(xintercept = 10, linetype = "dotted", colour = "grey55") +
  annotate("text", x = 14, y = 0.32, label = "Burn-in ends",
           colour = "grey45", size = 5) +
  annotate("text", x = 65, y = 0.92, label = "Clip ceiling = 0.90",
           colour = "grey45", size = 5) +
  scale_y_continuous(limits = c(0.28, 1), labels = percent_format(1)) +
  BASE_THEME +
  labs(title    = "Fig 3 (Sec 6.4): Bayesian RAR \u2014 P(Arm A Superior)",
       subtitle = "Median + IQR ribbon | H\u2081: pa=0.80, n=100",
       x = "Participant number", y = "P(\u03b8\u1d2c > \u03b8\u1d2e | data)")
ggsave("fig3_posterior.png", fig3, width = 8, height = 5, dpi = 300)

# Figure 4
.prep_sens <- function(df, n_label) {
  df %>%
    filter(!grepl("pa=0.50", Scenario)) %>%
    mutate(pa      = as.numeric(sub("pa=(.*?),.*", "\\1", Scenario)),
           n_label = n_label,
           Method  = factor(Method, c("Fixed 1:1", "Simple RAR", "Bayesian RAR")))
}
sens_df <- bind_rows(.prep_sens(results_100, "n = 100"),
                     .prep_sens(results_500, "n = 500"))

fig4 <- ggplot(sens_df, aes(pa, Power_or_TypeI,
                            colour   = Method,
                            linetype = n_label,
                            group    = interaction(Method, n_label))) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 1) + geom_point(size = 2.5) +
  scale_colour_manual(values = COLS) +
  scale_linetype_manual(values = c("n = 100" = "solid", "n = 500" = "longdash")) +
  scale_y_continuous(labels = percent_format(1), limits = c(0, 1)) +
  BASE_THEME +
  labs(title    = "Figure 4 (Section 6.5): Sensitivity \u2014 Sample Size Comparison",
       subtitle = "Solid = n=100 | Dashed = n=500 | pb=0.50",
       x = "True success probability of Arm A (pa)", y = "Estimated Power",
       colour = "Method", linetype = "Sample size")
ggsave("fig4_sensitivity.png", fig4, width = 9, height = 5, dpi = 300)

fig5 <- ggplot(cluster_fixed_df, aes(ICC, Power)) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey50") +
  geom_ribbon(aes(ymin = CI_L, ymax = CI_U), fill = "darkorange", alpha = 0.2) +
  geom_line(colour = "darkorange", linewidth = 1.1) +
  geom_point(colour = "darkorange", size = 3) +
  geom_text(aes(label = sprintf("DEFF=%.2f", Design_Effect)),
            vjust = 2, hjust = 0.5, size = 4.5, colour = "grey30") +
  scale_x_continuous(labels = percent_format(1), limits = c(-0.02, 0.23)) +
  scale_y_continuous(labels = percent_format(1), limits = c(0, 1)) +
  BASE_THEME +
  labs(title    = "Fig 5 (Sec 7.2): Cluster Randomisation \u2014 Power vs ICC",
       subtitle = "20 clusters per arm, 10 participants per cluster | pa=0.65, pb=0.50",
       x = "Intra-cluster correlation (ICC, \u03c1)", y = "Estimated Power")
ggsave("fig5_cluster_power.png", fig5, width = 8, height = 5, dpi = 300)
# Figure 6
fig6 <- cluster_all_df %>%
  mutate(Method = factor(Method, levels = c("Cluster Fixed", "Cluster-Adapt"))) %>%
  ggplot(aes(ICC, Power, colour = Method, group = Method)) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 1.1) + geom_point(size = 3) +
  geom_errorbar(aes(ymin = CI_L, ymax = CI_U), width = 0.01, alpha = 0.7) +
  scale_colour_manual(values = c("Cluster Fixed" = "darkorange",
                                 "Cluster-Adapt" = "goldenrod")) +
  scale_x_continuous(labels = percent_format(1)) +
  scale_y_continuous(labels = percent_format(1), limits = c(0, 1)) +
  BASE_THEME +
  labs(title    = "Figure 6 (Section 7.3): Cluster Adaptive vs Fixed \u2014 Power vs ICC",
       subtitle = "20 clusters per arm, 10 per cluster | pa=0.65, pb=0.50",
       x = "ICC (\u03c1)", y = "Estimated Power", colour = "Method")
ggsave("fig6_cluster_adapt.png", fig6, width = 8, height = 5, dpi = 300)

# Figure 7
ca_long <- ca_results %>%
  mutate(Method     = factor(Method, c("Fixed 1:1", "CA-RAR")),
         Hypothesis = ifelse(grepl("H0", Scenario),
                             "H\u2080 (no effect)", "H\u2081 (pa=0.65)"))

f7a <- ggplot(ca_long, aes(Hypothesis, Power_or_TypeI, fill = Method)) +
  geom_col(position = position_dodge(0.7), width = 0.6, alpha = 0.85) +
  geom_errorbar(aes(ymin = CI_L, ymax = CI_U),
                position = position_dodge(0.7), width = 0.25, linewidth = 0.7) +
  geom_hline(yintercept = 0.05, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = c("Fixed 1:1" = "firebrick", "CA-RAR" = "mediumpurple")) +
  scale_y_continuous(labels = percent_format(1), limits = c(0, 1)) +
  BASE_THEME + theme(legend.position = "right") +
  labs(title = "7a: Type I Error & Power", x = NULL, y = "Rate", fill = "Method")

f7b <- ca_long %>% filter(grepl("H1", Scenario)) %>%
  ggplot(aes(Method, Alloc_to_A, fill = Method)) +
  geom_col(width = 0.5, alpha = 0.85) +
  geom_errorbar(aes(ymin = Alloc_CI_L, ymax = Alloc_CI_U), width = 0.2, linewidth = 0.7) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = c("Fixed 1:1" = "firebrick", "CA-RAR" = "mediumpurple")) +
  scale_y_continuous(labels = percent_format(1), limits = c(0, 1)) +
  BASE_THEME + theme(legend.position = "none") +
  labs(title = "7b: Allocation to Superior Arm (H\u2081)",
       x = NULL, y = "Proportion to Arm A")

fig7 <- f7a + f7b + plot_annotation(
  title = "Figure 7 (Section 8.2): Covariate-Adjusted RAR vs Fixed 1:1",
  theme = theme(plot.title = element_text(face = "bold", size = 16))
)
ggsave("fig7_ca_rar.png", fig7, width = 11, height = 5, dpi = 300)

# Figure 8
ca_betax_plot <- ca_betax_df %>%
  mutate(Method = factor(Method, c("Fixed 1:1", "CA-RAR")))

f8a <- ggplot(ca_betax_plot, aes(beta_x, Power, colour = Method, group = Method)) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 1.1) + geom_point(size = 3) +
  scale_colour_manual(values = c("Fixed 1:1" = "firebrick", "CA-RAR" = "mediumpurple")) +
  scale_y_continuous(labels = percent_format(1), limits = c(0, 1)) +
  BASE_THEME + theme(legend.position = "right") +
  labs(title = "8a: Power vs Covariate Effect",
       x = "beta_x (log-odds)", y = "Power", colour = "Method")

f8b <- ggplot(ca_betax_plot, aes(beta_x, Alloc, colour = Method, group = Method)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey40") +
  geom_line(linewidth = 1.1) + geom_point(size = 3) +
  scale_colour_manual(values = c("Fixed 1:1" = "firebrick", "CA-RAR" = "mediumpurple")) +
  scale_y_continuous(labels = percent_format(1), limits = c(0.4, 0.8)) +
  BASE_THEME + theme(legend.position = "none") +
  labs(title = "8b: Allocation vs Covariate Effect",
       x = "beta_x (log-odds)", y = "Proportion to Arm A")

fig8 <- f8a + f8b + plot_annotation(
  title    = paste0("Figure 8 (Section 8.3): CA-RAR \u2014 Varying",
                    " Covariate Effect (n=100, pa_base=0.65)"),
  subtitle = "Power and ethical allocation as covariate confounding grows",
  theme    = theme(plot.title    = element_text(face = "bold", size = 16),
                   plot.subtitle = element_text(colour = "grey45", size = 13))
)
ggsave("fig8_ca_betax.png", fig8, width = 11, height = 5, dpi = 300)

# Figure 9
disc_heat <- discussion_table %>%
  mutate(
    Metric = case_when(
      Scenario == "H0"                    ~ "Type I Error",
      Scenario == "H1 moderate (pa=0.65)" ~ "Power (moderate)",
      Scenario == "H1 large (pa=0.80)"    ~ "Power (large)"
    ),
    Value  = ifelse(!is.na(Power), Power, Type_I),
    Method = factor(Method, levels = c("Fixed 1:1", "Simple RAR",
                                       "Bayesian RAR", "CA-RAR"))
  ) %>%
  filter(!is.na(Value))

disc_long <- bind_rows(
  disc_heat %>% transmute(Method, Metric, Value),
  disc_heat %>%
    filter(Metric != "Type I Error") %>%
    transmute(Method,
              Metric = paste0("Ethical Alloc (", Metric, ")"),
              Value  = Ethical_Alloc)
)

fig9 <- ggplot(disc_long, aes(x = Metric, y = Method, fill = Value)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.3f", Value)),
            size = 5, colour = "white", fontface = "bold") +
  scale_fill_gradient2(low = "firebrick", mid = "yellow", high = "seagreen",
                       midpoint = 0.5, limits = c(0, 1),
                       labels = percent_format(1)) +
  scale_x_discrete(guide = guide_axis(angle = 20)) +
  BASE_THEME +
  theme(legend.position = "right",
        panel.grid      = element_blank(),
        axis.text.x     = element_text(size = 13),
        axis.text.y     = element_text(size = 13)) +
  labs(title    = "Figure 9 (Section 9): Unified Method Comparison",
       subtitle = paste0("Type I error, power, and ethical allocation",
                         " across all methods (n=100)"),
       x = NULL, y = NULL, fill = "Value")
ggsave("fig9_discussion_heatmap.png", fig9, width = 11, height = 5, dpi = 300)

cat("\n=== All 9 figures complete! ===\n")

print(fig1)
print(fig2)
print(fig3)
print(fig4)
print(fig5)
print(fig6)
print(fig7)
print(fig8)
print(fig9)