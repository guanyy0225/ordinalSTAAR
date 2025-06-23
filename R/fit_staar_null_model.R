#' Fit STAAR Null Model
#'
#' This function fits a simplified null model for STAAR, either an ordinal clm
#' model or a generalized linear model depending on the family specified.
#'
#' @param fixed Formula specifying the model equation.
#' @param data Data frame containing the model variables.
#' @param family Family of the model. For ordinal models, use "ordinal".
#' @param ... Additional arguments passed to the fitting functions.
#'
#' @return An object suitable for use in STAAR analysis.
#'
fit_staar_null_model <- function(fixed, data, family = binomial(link = "logit"), ...) {

  is_ordinal <- is.character(family) && family == "ordinal"

  if (is_ordinal) {
    obj_nullmodel <- fit_clm_for_staar(fixed = fixed, data = data, ...)
  } else {
    obj_nullmodel <- stats::glm(formula = fixed, data = data, family = family, ...)
    obj_nullmodel$relatedness <- FALSE
  }

  return(obj_nullmodel)
}

#' Fit STAAR CLM for STAAR
#'
#' This function fits an ordinal clm model specifically for STAAR analysis,
#' calculates predicted probabilities, and assembles the necessary components.
#'
#' @param fixed Formula specifying the model equation.
#' @param data Data frame containing the model variables.
#' @param ... Additional arguments passed to the clm function.
#'
#' @return An object containing the fitted model components for STAAR.
#'
fit_clm_for_staar <- function(fixed, data, ...) {

  # --- Step 1: Fit ordinal clm model ---
  message("Step 1/3: Fitting ordinal clm null model...")
  clm_fit <- tryCatch({
    ordinal::clm(formula = fixed, data = data, model = TRUE, Hess = TRUE, ...)
  }, error = function(e) {
    stop("Error during clm() model fitting: ", e$message)
  })
  message("Model fitting completed.")

  # --- Step 2: Calculate predicted probabilities manually ---
  message("Step 2/3: Calculating predicted probabilities manually...")

  # Extract core components
  beta_coeffs <- clm_fit$beta
  alpha_coeffs <- clm_fit$alpha
  link_obj <- make.link(clm_fit$link)
  link_func_inv <- link_obj$linkinv

  # Manual calculation
  X_matrix <- model.matrix(clm_fit$terms, data = clm_fit$model)
  X_covariates <- X_matrix[, -which(colnames(X_matrix) == "(Intercept)"), drop = FALSE]
  latent_predictor <- as.vector(X_covariates %*% beta_coeffs)

  cum_probs_list <- lapply(alpha_coeffs, function(a) {
    link_func_inv(a - latent_predictor)
  })
  cum_probs_matrix <- do.call(cbind, cum_probs_list)

  cum_probs_padded <- cbind(0, cum_probs_matrix, 1)
  K <- length(clm_fit$y.levels)
  pred_probs <- matrix(NA, nrow = nrow(cum_probs_padded), ncol = K)
  for (k in 1:K) {
    pred_probs[, k] <- cum_probs_padded[, k + 1] - cum_probs_padded[, k]
  }
  pred_probs[pred_probs < 0] <- 0
  pred_probs <- pred_probs / rowSums(pred_probs)
  message("Probability calculation completed.")

  # --- Step 3: Assemble simplified STAAR "pseudo GLM" object ---
  message("Step 3/3: Assembling simplified STAAR null model...")

  # Calculate mean and variance
  categories <- 1:K
  E_y <- as.vector(pred_probs %*% categories)
  E_y_squared <- as.vector(pred_probs %*% (categories^2))
  Var_y <- E_y_squared - (E_y^2)

  # Calculate residuals and weights
  y_factor <- clm_fit$y
  y_numeric <- as.numeric(y_factor)
  residuals <- y_numeric - E_y
  weights <- 1 / (Var_y + 1e-8)

  # Assemble components
  staar_null_obj <- list()

  # STAAR requires these components explicitly
  staar_null_obj$y <- y_numeric
  staar_null_obj$fitted.values <- E_y
  staar_null_obj$weights <- weights
  staar_null_obj$family <- gaussian(link = "identity")
  staar_null_obj$summary <- list(dispersion = 1.0)
  staar_null_obj$relatedness <- FALSE
  staar_null_obj$terms <- clm_fit$terms
  staar_null_obj$model <- clm_fit$model

  # Mask as glm object to handle with model.matrix etc.
  class(staar_null_obj) <- c("glm", "lm")

  message("Simplified STAAR null model created successfully.")
  return(staar_null_obj)
}
