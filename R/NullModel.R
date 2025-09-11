NullModel <- function(genofile = NULL, phenofile, outcomeCol, sampleCol,
                              covCol = NULL, PRSCol = NULL, LOCO = TRUE, chr = NULL,
                              use_SPA = FALSE, range=c(-100,100), length.out = 1e4, verbose = FALSE) {

  # --- Part 1: Load and Prepare Phenotype Data ---
  if (!is.null(phenofile)) {
    if (is.character(phenofile)) {
      if (!file.exists(phenofile)) stop("Phenotype file does not exist!")
      use_data <- data.table::fread(phenofile, data.table = FALSE)
    } else {
      use_data <- as.data.frame(phenofile)
    }
    use_data <- na.omit(use_data[, c(sampleCol, outcomeCol, covCol, PRSCol)])
  } else {
    stop("You must provide the phenotype data via the 'phenofile' argument.")
  }

  use_data[[sampleCol]] <- as.character(use_data[[sampleCol]])

  # --- Part 2: Match Samples with Genotype Data (if provided) ---
  if (!is.null(genofile)) {
    if (inherits(genofile, "SeqVarGDSClass")) {
      sample.geno <- SeqArray::seqGetData(genofile, "sample.id")
    } else if (is.character(genofile)) {
      fam_file <- paste0(genofile, ".fam")
      if (!file.exists(fam_file)) stop("PLINK .fam file not found!")
      sample.geno <- data.table::fread(fam_file, data.table = FALSE)$V2
    } else {
      stop("'genofile' must be a GDS object or a path to a PLINK file.")
    }

    sample.pheno <- use_data[, sampleCol]
    common_samples <- intersect(sample.pheno, sample.geno)

    if (length(common_samples) == 0) stop("No common samples found between phenotype and genotype files.")

    use_data <- use_data[use_data[, sampleCol] %in% common_samples, ]
    message(paste(nrow(use_data), "samples remain after matching."))
  }

  # --- Part 3: Dynamic Formula Construction ---
  if (!is.ordered(use_data[[outcomeCol]])) {
    use_data[[outcomeCol]] <- as.ordered(use_data[[outcomeCol]])
  }

  all_covars <- covCol
  if (is.null(PRSCol)) {
    formula_null <- as.formula(paste(outcomeCol, "~", paste(all_covars, collapse = " + ")))
  } else {
    formula_null <- as.formula(paste(outcomeCol, "~", paste(all_covars, collapse = " + "), "+", paste(PRSCol, collapse = " + ")))
  }

  # --- Part 4: Fit the Ordinal Null Model ---
  clm_obj <- ordinal::clm(formula = formula_null, data = use_data, link = "probit", model = TRUE)
  if(verbose) print(summary(clm_obj))

  # --- Part 5: Calculate Residuals and Variance Components ---
  model_data <- clm_obj$model
  kept_row_indices <- as.numeric(rownames(model_data))
  sample_ids <- use_data[[sampleCol]][kept_row_indices]

  alpha_coefs <- clm_obj$beta
  X_mat <- model.matrix(object = formula(clm_obj), data = model_data)
  eta <- as.vector(X_mat[, names(alpha_coefs), drop = FALSE] %*% alpha_coefs)
  thresholds <- c(-Inf, clm_obj$alpha, Inf)
  y_idx <- as.numeric(clm_obj$y)

  lower_b <- thresholds[y_idx] - eta
  upper_b <- thresholds[y_idx + 1] - eta
  prob_interval <- pnorm(upper_b) - pnorm(lower_b)
  prob_interval[prob_interval < 1e-12] <- 1e-12
  residuals <- (dnorm(lower_b) - dnorm(upper_b)) / prob_interval

  term1 <- (lower_b * dnorm(lower_b) - upper_b * dnorm(upper_b)) / prob_interval
  var_y <- 1 + term1 - residuals^2
  var_y[!is.finite(var_y) | var_y < 1e-8] <- 1

  # Var(U) = G'WG - (G'WX) * (X'WX)⁻¹ * (X'WG)
  W_mat <- Diagonal(x = 1 / var_y)
  X_t_W <- crossprod(X_mat, W_mat)
  XWX_mat <- X_t_W %*% X_mat
  XWX_inv <- solve(XWX_mat)
  WX_mat <- t(X_t_W)


  # --- Part 6: Assemble the Final List with LOCO information ---
  base_list <- list(
    residuals = residuals,
    sample_ids = sample_ids,
    W_mat = W_mat,
    X_mat = X_mat,
    WX_mat = WX_mat,
    XWX_inv = XWX_inv,
    formula_null = formula(clm_obj),
    coefficients_null = coef(clm_obj),
    use_data = model_data,
    use_SPA = use_SPA
  )

  if (is.null(LOCO) || isFALSE(LOCO)) {
    fit_null <- c(base_list, list(LOCO = FALSE))
  } else {
    fit_null <- c(base_list, list(LOCO = TRUE, chr = chr))
  }

  if(use_SPA){
    fit_null <- CGF4LatentRes(fit_null = fit_null, range = range, length.out = length.out, verbose = verbose)
    warning("use_SPA=TRUE is not yet implemented. Skipping CGF calculation.")
  }

  message("Ordinal null model fitting complete.")
  return(fit_null)
}
