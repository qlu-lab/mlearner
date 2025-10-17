#!/usr/bin/env Rscript

# =========================
# Load packages
# =========================
require(optparse)
require(mlr3)
require(mlr3learners)
require(sandwich)
require(lmtest)
require(ggplot2)
require(scales)
require(data.table)


# =========================
#  Propensity score
# =========================
# learner ∈ {"constant","lasso","random_forest","tree","svm","logit"}
# returns list(estimates, mlr3_objects)
propensity_score <- function(Z, D, learner = c("constant","lasso","random_forest","tree","svm","logit"),
                             n_threads = 1L) {
  learner <- match.arg(learner)
  stopifnot(is.numeric(D), all(D %in% c(0,1)))
  if (!is.matrix(Z) && !is.data.frame(Z))
    stop("'Z' must be a matrix or data.frame.")
  n <- NROW(Z)

  if (learner == "constant") {
    est <- rep(mean(D), n)
    return(structure(list(
      estimates    = est,
      mlr3_objects = NULL
    ), class = "propensity_score"))
  }

  # ----------------- mlr3 learners -----------------
  lrn <- switch(
    learner,
    "lasso"         = mlr3::lrn("classif.cv_glmnet", s = "lambda.min", alpha = 1),
    "random_forest" = mlr3::lrn("classif.ranger", num.trees = 500),
    "tree"          = mlr3::lrn("classif.rpart"),
    "svm"           = mlr3::lrn("classif.svm"),
    "logit"         = mlr3::lrn("classif.log_reg")
  )
  lrn$predict_type <- "prob"
  invisible(mlr3::set_threads(lrn, n_threads))

  # ----------------- task -----------------
  df <- data.frame(D = factor(D, levels = c(0,1)), Z, check.names = TRUE)
  task <- mlr3::TaskClassif$new(id = "propensity_score", backend = df, target = "D")

  # ----------------- fit & predict -----------------
  lrn$train(task)
  pr   <- lrn$predict(task)$prob
  col1 <- which(colnames(pr) %in% c("1","TRUE","Yes"))
  if (length(col1) != 1)
    stop("Could not locate positive-class column in probabilities.")
  p1 <- pr[, col1]

  structure(list(
    estimates    = p1,
    mlr3_objects = list(task = task, learner = lrn)
  ), class = "propensity_score")
}


# =========================
# Unified learner builder
# =========================
make_learner <- function(name, n_threads = 1L){
  lrn_id  <- paste0("regr.", name)
  learner <- mlr3::lrn(lrn_id, predict_type = "response")
  invisible(mlr3::set_threads(learner, n_threads))
  # sensible defaults
  if (name == "cv_glmnet") {
    learner$param_set$values <- list(alpha = 0.5)
  } else if (name == "nnet") {
    learner$param_set$values <- list(size = 5, decay = 0.1, maxit = 200, MaxNWts = 5e4)
  }
  learner
}

# =========================
# BCA: fit E[Y|D=0,Z] on A_set controls and predict for all
# =========================
proxy_BCA <- function(Z, D, Y, A_set,
                      base_learners, min_variation = 1e-5,
                      n_threads = 1L, verbose = FALSE){
  stopifnot(is.numeric(D), all(D %in% c(0,1)))
  stopifnot(length(Y) == nrow(Z), length(D) == nrow(Z))
  if (!is.matrix(Z) && !is.data.frame(Z)) stop("'Z' must be a matrix or data.frame.")
  if (length(base_learners) < 1) stop("'base_learners' must contain at least one learner.")
  if (length(A_set) < 2) stop("'A_set' is too small to train models.")

  idx_ctrl <- which((seq_len(length(D)) %in% A_set) & D == 0)
  if (length(idx_ctrl) < 2) stop("Not enough control samples in A_set to train E[Y|D=0,Z].")

  df   <- data.frame(Y = Y, Z, check.names = TRUE)
  task <- mlr3::TaskRegr$new(id = "bca.controls", backend = df, target = "Y")

  lrn  <- make_learner(base_learners[[1L]], n_threads)
  lrn$train(task, row_ids = idx_ctrl)
  y0hat <- lrn$predict(task)$response

  if (is.finite(stats::var(y0hat)) && stats::var(y0hat) < min_variation) {
    y0hat <- y0hat + stats::rnorm(length(Y), 0, sqrt(max(stats::var(Y), 1e-8))/20)
    if (verbose) message("[proxy_BCA] Variance too small; added tiny noise.")
  }

  structure(list(
    estimates    = y0hat,
    mlr3_objects = list(Y0_learner = list(task = task, learner = lrn))
  ), class = "proxy_BCA")
}

# =========================
# CATE: estimate E[Y|D=1,Z] - E[Y|D=0,Z]
# =========================
proxy_CATE <- function(Z, D, Y, A_set,
                       base_learners, proxy_BCA = NULL,
                       min_variation = 1e-5, n_threads = 1L, verbose = FALSE){
  stopifnot(is.numeric(D), all(D %in% c(0,1)))
  stopifnot(length(Y) == nrow(Z), length(D) == nrow(Z))
  if (!is.matrix(Z) && !is.data.frame(Z)) stop("'Z' must be a matrix or data.frame.")
  if (length(base_learners) < 1) stop("'base_learners' must contain at least one learner.")
  if (length(A_set) < 2) stop("'A_set' is too small to train models.")

  idx_treat <- which((seq_len(length(D)) %in% A_set) & D == 1)
  idx_ctrl  <- which((seq_len(length(D)) %in% A_set) & D == 0)
  if (length(idx_treat) < 2) stop("Not enough treated in A_set.")
  if (length(idx_ctrl)  < 2) stop("Not enough control in A_set.")

  df   <- data.frame(Y = Y, Z, check.names = TRUE)
  task <- mlr3::TaskRegr$new(id = "cate.task", backend = df, target = "Y")

  lrn1 <- make_learner(base_learners[[1L]], n_threads)
  lrn1$train(task, row_ids = idx_treat)
  y1hat <- lrn1$predict(task)$response

  if (is.null(proxy_BCA)) {
    lrn0 <- make_learner(base_learners[[1L]], n_threads)
    lrn0$train(task, row_ids = idx_ctrl)
    y0hat <- lrn0$predict(task)$response
  } else {
    y0hat <- proxy_BCA
    lrn0  <- "provided"
  }

  cate <- y1hat - y0hat
  if (is.finite(stats::var(cate)) && stats::var(cate) < min_variation) {
    cate <- cate + stats::rnorm(length(Y), 0, sqrt(max(stats::var(Y), 1e-8))/20)
    if (verbose) message("[proxy_CATE] Variance too small; added tiny noise.")
  }

  structure(list(
    estimates    = list(CATE = cate, Y1 = y1hat, Y0 = y0hat),
    mlr3_objects = list(Y1_learner = list(task = task, learner = lrn1),
                        Y0_learner = list(task = task, learner = lrn0))
  ), class = "proxy_CATE")
}

# =========================
# Build X1 design from functions (S/B/p) + covariates/fixed effects
# =========================
get_X1_df <- function(functions_mat, X1_control){
  keep <- intersect(X1_control$funs_Z, colnames(functions_mat))
  X1 <- as.data.frame(functions_mat[, keep, drop = FALSE])
  if (!is.null(X1_control$covariates)) {
    X1 <- cbind(X1, as.data.frame(X1_control$covariates))
  }
  if (!is.null(X1_control$fixed_effects)) {
    X1$fixed.effects <- X1_control$fixed_effects
  }
  X1
}

# =========================
# Sandwich covariance helper
# =========================
get_vcov <- function(x, type = "const"){
  sandwich::vcovHC(x, type = type)
}

# =========================
# BLP 
# =========================
BLP <- function(D, Y, propensity_scores,
                proxy_BCA, proxy_CATE,
                X1_control = list(funs_Z = c("B"), covariates = NULL, fixed_effects = NULL),
                significance_level = 0.05){

  H <- (D - propensity_scores) / (propensity_scores * (1 - propensity_scores))
  X1. <- get_X1_df(cbind(S = proxy_CATE, B = proxy_BCA, p = propensity_scores), X1_control)

  if (is.null(X1_control$fixed_effects)) {
    X1H <- X1. * H
  } else {
    fe  <- X1.$fixed.effects
    Xnr <- X1.[, setdiff(colnames(X1.), "fixed.effects"), drop = FALSE]
    X1H <- data.frame(Xnr * H, model.matrix(~ fe + 0)[, -1, drop = FALSE])
  }
  colnames(X1H) <- paste0(colnames(X1H), ".H")

  X <- data.frame(X1H, beta.2 = proxy_CATE - mean(proxy_CATE))
  fit <- stats::lm(Y * H ~ ., data = X)
  vc  <- get_vcov(fit)
  ct  <- lmtest::coeftest(fit, vcov. = vc)

  pick <- c("(Intercept)", "beta.2")
  keep <- ct[pick, 1:3, drop = FALSE]
  colnames(keep) <- c("Estimate", "Std. Error", "z value")
  z  <- stats::qnorm(1 - significance_level/2)
  ci <- cbind(keep[, "Estimate"] - z * keep[, "Std. Error"],
              keep[, "Estimate"] + z * keep[, "Std. Error"])
  colnames(ci) <- c("CB lower", "CB upper")
  p  <- 2 * stats::pnorm(-abs(keep[, "z value"]))

  generic <- cbind(keep[, "Estimate", drop = FALSE], ci, keep[, "Std. Error", drop = FALSE],
                   keep[, "z value", drop = FALSE],
                   "Pr(<z)" = p/2, "Pr(>z)" = p/2)
  rownames(generic) <- c("beta.1", "beta.2")

  structure(list(
    generic_targets = generic,
    coefficients    = ct,
    lm              = fit
  ), class = "BLP")
}


# =========================
# SATES 
# =========================
SATES <- function(D, Y, propensity_scores,
                  proxy_BCA, proxy_CATE, membership,
                  X1_control = list(funs_Z = c("B"), covariates = NULL, fixed_effects = NULL),
                  significance_level = 0.05){

  groups <- 1 * membership
  K <- ncol(groups)
  H <- (D - propensity_scores) / (propensity_scores * (1 - propensity_scores))

  X1 <- get_X1_df(cbind(S = proxy_CATE, B = proxy_BCA, p = propensity_scores), X1_control)
  if (is.null(X1_control$fixed_effects)) {
    X1H <- X1 * H
  } else {
    fe  <- X1$fixed.effects
    Xnr <- X1[, setdiff(colnames(X1), "fixed.effects"), drop = FALSE]
    X1H <- data.frame(Xnr * H, model.matrix(~ fe + 0)[, -1, drop = FALSE])
  }
  colnames(X1H) <- paste0(colnames(X1H), ".H")

  X <- data.frame(X1H, groups)
  colnames(X) <- c(colnames(X1H), paste0("gamma.", seq_len(K)))

  fit <- stats::lm(stats::as.formula(paste0("I(Y*H) ~ ", paste0(colnames(X), collapse = " + "), " + 0")),
                   data = cbind(Y = Y, X))
  vc  <- get_vcov(fit)
  ct  <- lmtest::coeftest(fit, vcov. = vc)

  gamma_rows <- paste0("gamma.", seq_len(K))
  keep <- ct[gamma_rows, 1:3, drop = FALSE]
  colnames(keep) <- c("Estimate", "Std. Error", "z value")
  z  <- stats::qnorm(1 - significance_level/2)
  ci <- cbind(keep[, "Estimate"] - z * keep[, "Std. Error"],
              keep[, "Estimate"] + z * keep[, "Std. Error"])
  colnames(ci) <- c("CB lower", "CB upper")
  p  <- 2 * stats::pnorm(-abs(keep[, "z value"]))

  generic <- cbind(keep[, "Estimate", drop = FALSE], ci, keep[, "Std. Error", drop = FALSE],
                   keep[, "z value", drop = FALSE],
                   "Pr(<z)" = p/2, "Pr(>z)" = p/2)

  structure(list(
    generic_targets = generic,
    coefficients    = ct,
    lm              = fit
  ), class = "SATES")
}


# =========================
# FIS (per-feature OLS t-stat)
# =========================
FIS <- function(Z, CATE, M_set){
  # Returns a named numeric vector of t-stats (CATE ~ z) for each column of Z.
  Z_M    <- Z[M_set, , drop = FALSE]
  CATE_M <- CATE[M_set]
  p <- ncol(Z_M)
  out <- rep(NA_real_, p)
  names(out) <- colnames(Z_M)
  for (j in seq_len(p)) {
    z <- Z_M[, j]
    if (length(unique(z)) <= 1) next
    cf <- summary(lm(CATE_M ~ z))$coefficients
    out[j] <- if (nrow(cf) >= 2) cf[2, 3] else NA_real_
  }
  out
}


# =========================
# Quantile-based grouping 
# =========================
quantile_group <- function(x, cutoffs){
  K <- length(cutoffs) + 1L
  n <- length(x)
  for (tries in 1:4){
    q  <- stats::quantile(x, cutoffs)
    g  <- matrix(FALSE, n, K)
    for (k in seq_len(K)){
      if (k == 1)      g[,k] <- x <  q[k]
      else if (k == K) g[,k] <- x >= q[k-1]
      else             g[,k] <- (q[k-1] <= x) & (x < q[k])
      if (sum(g[,k]) < 2){ g[,] <- FALSE; break }
    }
    if (any(g)) return(structure(g*1L, type = "quantile_group"))
    x <- x + rnorm(n, 0, 1e-3)
  }
  stop("Unable to form legal quantile groups; adjust cutoffs.")
}


# =========================
# Single fold computation
# =========================
single_fold <- function(Y, D, Z, A_set, M_set,
                        propensity_score, base_learners, quantile_cutoffs){

  bca  <- proxy_BCA(Z, D, Y, A_set, base_learners)
  cate <- proxy_CATE(Z, D, Y, A_set, base_learners, proxy_BCA = bca$estimates)

  pscore_M <- propensity_score[M_set]
  bca_M    <- bca$estimates[M_set]
  cate_M   <- cate$estimates$CATE[M_set]

  X1_cfg_M <- list(funs_Z = c("B"), covariates = NULL, fixed_effects = NULL)

  blp <- BLP(D[M_set], Y[M_set], pscore_M, bca_M, cate_M, X1_cfg_M)
  mem <- quantile_group(cate_M, cutoffs = quantile_cutoffs)
  sat <- SATES(D[M_set], Y[M_set], pscore_M, bca_M, cate_M, mem, X1_cfg_M)
  fim <- FIS(Z, CATE = cate$estimates$CATE, M_set = M_set)

  list(BLP = blp, SATES = sat, FIS = fim)
}


# =========================
# Initialize result containers
# =========================
initializer.for.splits <- function(Z, num_folds, quantile_cutoffs){
  p <- ncol(Z); K <- length(quantile_cutoffs) + 1L
  list(
    BLP    = array(NA_real_, c(2, 7, num_folds),
                   dimnames = list(c("beta.1","beta.2"),
                                   c("Estimate","CB lower","CB upper","Std. Error","z value","Pr(<z)","Pr(>z)"),
                                   NULL)),
    SATES  = array(NA_real_, c(K, 7, num_folds),
                   dimnames = list(paste0("gamma.", 1:K),
                                   c("Estimate","CB lower","CB upper","Std. Error","z value","Pr(<z)","Pr(>z)"),
                                   NULL)),
    FIS = matrix(NA_real_, nrow = num_folds, ncol = p,
                     dimnames = list(NULL, colnames(Z)))
  )
}


# =========================
# Fold-aggregate reporting
# =========================
for.output <- function(result, significance_level = 0.05){
  folds <- dim(result$BLP)[3]
  blp <- matrix(NA_real_, 2, 4, dimnames=list(c("beta.1","beta.2"), c("Estimate","CB lower","CB upper","P value")))
  sat <- matrix(NA_real_, dim(result$SATES)[1], 4, dimnames=list(rownames(result$SATES), c("Estimate","CB lower","CB upper","P value")))

  # BLP
  blp[, "Estimate"]  <- apply(result$BLP[,"Estimate",], 1, mean, na.rm = TRUE)
  zsum               <- apply(result$BLP[,"z value",], 1, function(x) sum(x, na.rm = TRUE)) / sqrt(folds)
  blp[, "P value"]   <- 2*pnorm(-abs(zsum))
  z                  <- stats::qnorm(1 - significance_level/2)
  se_hat             <- abs(blp[, "Estimate"]) / pmax(abs(zsum), 1e-12)
  blp[, "CB lower"]  <- blp[, "Estimate"] - z * se_hat
  blp[, "CB upper"]  <- blp[, "Estimate"] + z * se_hat

  # SATES
  sat[, "Estimate"]  <- apply(result$SATES[,"Estimate",], 1, mean, na.rm = TRUE)
  zsum_g             <- apply(result$SATES[,"z value",], 1, function(x) sum(x, na.rm = TRUE)) / sqrt(folds)
  sat[, "P value"]   <- 2*pnorm(-abs(zsum_g))
  se_hat_g           <- abs(sat[, "Estimate"]) / pmax(abs(zsum_g), 1e-12)
  sat[, "CB lower"]  <- sat[, "Estimate"] - z * se_hat_g
  sat[, "CB upper"]  <- sat[, "Estimate"] + z * se_hat_g

  # FEATIMP/fis: combine per-feature t-stats across folds into two-sided p-values
  featimp_p <- 2 * pnorm(-abs(colSums(result$FIS, na.rm = TRUE) / sqrt(folds)))

  list(BLP = blp, SATES = sat, FIS = featimp_p)
}


# =========================
# Main pipeline
# =========================
M_learner <- function(Y, D, Z, num_folds, seed,
                      learner_propensity_score = "constant",
                      base_learners, quantile_cutoffs, significance_level = 0.05){

  set.seed(seed)
  n <- nrow(Z); p <- ncol(Z)
  fold.id <- sample(rep(seq_len(num_folds), length.out = n))
  parts   <- split(seq_len(n), fold.id)

  score  <- propensity_score(Z, D, learner = learner_propensity_score)
  store  <- initializer.for.splits(Z, num_folds - 1L, quantile_cutoffs)

  for (k in 2:num_folds){
    A_set <- unlist(parts[1:(k-1)], use.names = FALSE)
    M_set <- parts[[k]]

    one <- single_fold(Y, D, Z, A_set, M_set,
                       propensity_score = score$estimates,
                       base_learners    = base_learners,
                       quantile_cutoffs = quantile_cutoffs)

    store$BLP[,,k-1]    <- one$BLP$generic_targets
    store$SATES[,,k-1]  <- one$SATES$generic_targets
    store$FIS[k-1, ] <- one$FIS   # named numeric vector (length = p)
  }

  for.output(store, significance_level)
}


# =========================
# Plotter 
# =========================
plot_sates_png <- function(m_obj,
                           file    = NULL,
                           width   = 8,
                           height  = 6,
                           dpi     = 300,
                           title   = NULL,
                           group_labels = NULL,
                           y_limits = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install ggplot2.")
  if (!requireNamespace("scales", quietly = TRUE))  stop("Please install scales.")
  ggplot2::theme_set(ggplot2::theme_classic())

  has_tables <- !is.null(m_obj$tables)
  BLP   <- if (has_tables) m_obj$tables$BLP   else m_obj$BLP
  SATES <- if (has_tables) m_obj$tables$SATES else m_obj$SATES
  if (is.null(BLP) || is.null(SATES)) stop("BLP and/or SATES tables not found in m_obj.")

  if (!"beta.1" %in% rownames(BLP)) stop("Row 'beta.1' not found in BLP table.")
  req_cols_blp <- c("Estimate","CB lower","CB upper")
  if (!all(req_cols_blp %in% colnames(BLP))) {
    stop("BLP must have columns: ", paste(req_cols_blp, collapse = ", "))
  }
  ate      <- as.numeric(BLP["beta.1","Estimate"])
  ci_lower <- as.numeric(BLP["beta.1","CB lower"])
  ci_upper <- as.numeric(BLP["beta.1","CB upper"])

  req_cols_g <- c("Estimate","CB lower","CB upper")
  if (!all(req_cols_g %in% colnames(SATES))) {
    stop("SATES must have columns: ", paste(req_cols_g, collapse = ", "))
  }
  K <- nrow(SATES)
  default_labels <- if (is.null(group_labels)) {
    if (K == 3) c("Low HTE Score","Medium HTE Score","High HTE Score")
    else paste0("Group ", seq_len(K))
  } else {
    if (length(group_labels) != K) stop("group_labels length must equal number of groups.")
    group_labels
  }

  SATE_df <- data.frame(
    Group    = factor(default_labels, levels = default_labels),
    Estimate = as.numeric(SATES[,"Estimate"]),
    Lower    = as.numeric(SATES[,"CB lower"]),
    Upper    = as.numeric(SATES[,"CB upper"])
  )

  band_df <- data.frame(xmin = 0.5, xmax = K + 0.5, ymin = ci_lower, ymax = ci_upper, label = "ATE 95% CI")
  line_df <- data.frame(y = ate, label = "ATE")

  p <- ggplot2::ggplot(SATE_df, ggplot2::aes(x = Group, y = Estimate)) +
    ggplot2::geom_rect(data = band_df,
                       ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = label),
                       inherit.aes = FALSE, alpha = 0.12, color = NA) +
    ggplot2::geom_hline(data = line_df,
                        ggplot2::aes(yintercept = y, color = label),
                        linetype = "dashed", linewidth = 1) +
    ggplot2::geom_hline(yintercept = 0, color = "grey50", linetype = "dotted", linewidth = 0.6) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = Lower, ymax = Upper), width = 0.15, linewidth = 0.9) +
    ggplot2::geom_point(shape = 21, size = 3.8, stroke = 0.8, color = "white", fill = "#2C7FB8") +
    ggplot2::labs(title = title, x = "Subgroups by Predicted HTE Score", y = "Estimated Treatment Effect") +
    ggplot2::scale_fill_manual(name = "", values = c("ATE 95% CI" = "#2C7FB8")) +
    ggplot2::scale_color_manual(name = "", values = c("ATE" = "#2C7FB8")) +
    { if (is.null(y_limits)) {
        ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 5),
                                    expand = ggplot2::expansion(mult = c(0.03, 0.06)))
      } else {
        ggplot2::scale_y_continuous(limits = y_limits,
                                    breaks = scales::pretty_breaks(n = 5),
                                    expand = ggplot2::expansion(mult = c(0.03, 0.06)))
      }
    } +
    ggplot2::theme(
      plot.title   = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title.x = ggplot2::element_text(size = 14, face = "bold", margin = ggplot2::margin(t = 8)),
      axis.title.y = ggplot2::element_text(size = 14, face = "bold", margin = ggplot2::margin(r = 8)),
      axis.text    = ggplot2::element_text(size = 12),
      axis.line    = ggplot2::element_line(linewidth = 0.6),
      axis.ticks   = ggplot2::element_line(linewidth = 0.6),
      legend.position = "right",
      legend.box.margin = ggplot2::margin(t = -6),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background  = ggplot2::element_rect(fill = "white", colour = NA)
    ) +
    ggplot2::guides(
      fill  = ggplot2::guide_legend(order = 1, override.aes = list(alpha = 0.25, linewidth = 0)),
      color = ggplot2::guide_legend(order = 1, override.aes = list(linetype = "dashed", linewidth = 1))
    )

  if (is.null(file)) file <- file.path(tempdir(), "sates_plot.png")
  ggplot2::ggsave(filename = file, plot = p, width = width, height = height, dpi = dpi)
  invisible(file)
}


# =========================
# CLI main
# =========================
suppressMessages(library(optparse))

options(stringsAsFactors = FALSE)

option_list <- list(
  make_option("--outcome", action = "store", default = NA, type = "character",
              help = "Outcome file path (e.g., Y.txt)"),
  make_option("--treatment", action = "store", default = NA, type = "character",
              help = "Treatment indicator file path (e.g., D.txt)"),
  make_option("--prs", action = "store", default = NA, type = "character",
              help = "PRS/covariates file path (e.g., Z.txt)"),
  make_option("--num_folds", action = "store", default = NA, type = "integer",
              help = "Number of folds for cross-fitting"),
  make_option("--seed", action = "store", default = NA, type = "integer",
              help = "Random seed for reproducibility"),
  make_option("--prop_learner", action = "store", default = NA, type = "character",
              help = "Propensity score learner: constant, lasso, random_forest, tree, svm, or logit"),
  make_option("--base_learners", action = "store", default = NA, type = "character",
              help = "CATE learner(s), comma-separated (e.g., svm or svm,ranger)"),
  make_option("--quantile_cutoffs", action = "store", default = NA, type = "character",
              help = "Subgroup cutoffs, comma-separated decimals (e.g., 0.3333,0.6667)"),
  make_option("--significance_level", action = "store", default = NA, type = "double",
              help = "Significance level for confidence intervals (e.g., 0.05)"),
  make_option("--plot_file", action = "store", default = NULL, type = "character",
              help = "Output path for subgroup treatment effect plot (PNG format)"),
  make_option("--fis_file", action = "store", default = NULL, type = "character",
              help = "Output path for Feature Importance Scores (text file)")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Check required arguments
required <- c("outcome", "treatment", "prs", "num_folds", "seed",
              "prop_learner", "base_learners", "quantile_cutoffs", "significance_level")
missing <- required[sapply(opt[required], is.na)]
if (length(missing) > 0) {
  stop("Missing required arguments: ", paste0("--", missing, collapse = " "))
}

Y_file        <- opt$outcome
D_file        <- opt$treatment
Z_file        <- opt$prs
num_folds     <- opt$num_folds
seed          <- opt$seed
prop_learner  <- opt$prop_learner
base_learners <- strsplit(opt$base_learners, ",")[[1]]
quantiles     <- as.numeric(strsplit(opt$quantile_cutoffs, ",")[[1]])
alpha         <- opt$significance_level
plot_file     <- opt$plot_file
fis_file      <- opt$fis_file

# ---- Load data ----
Y <- as.matrix(fread(Y_file))
D <- as.matrix(fread(D_file))
Z <- as.matrix(fread(Z_file))

# ---- Run M-learner ----
res <- M_learner(
  Y, D, Z,
  num_folds          = num_folds,
  seed               = seed,
  learner_propensity_score = prop_learner,
  base_learners      = base_learners,
  quantile_cutoffs   = quantiles,
  significance_level = alpha
)

# ---- Print results ----
cat("\n===== BLP results =====\n")
print(res$BLP)

cat("\n===== SATES results =====\n")
print(res$SATES)

# ---- Optional outputs ----
if (!is.null(plot_file) && nzchar(plot_file)) {
  plot_sates_png(res, file = plot_file, title = "Sorted Average Treatment Effects (SATEs)")
  cat("\nPlot saved to:", plot_file, "\n")
}
if (!is.null(fis_file) && nzchar(fis_file)) {
  write.table(res$FIS, fis_file, quote = FALSE, row.names = TRUE, col.names = NA)
  cat("\nFIS scores saved to:", fis_file, "\n")
}