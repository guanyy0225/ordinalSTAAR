library(ordinal)
library(stats)

#' @title Create a STAAR-compatible null model from a clm object
#' @description This function converts a fitted `ordinal::clm` object into a 
#'   null model object compatible with STAAR. It achieves this by setting up the
#'   components for a genetic score test. The ordinal model is transformed into a 
#'   pseudo-linear model using a working variable and weights, where the variance
#'   of the working variable is stabilized. This allows the efficient score test 
#'   implemented in STAAR to be applied to ordinal data.
#'
#' @param clm_obj A fitted object from `ordinal::clm`. Crucially, it must be 
#'   fitted with the `model = TRUE` argument to ensure the model frame is stored.
#'
#' @return A list object with classes "POM", "glm", and "lm". This object mimics 
#'   the structure of a `glm` null model and contains the necessary components 
#'   for STAAR, such as the working dependent variable, weights, and residuals.
#'
create_null_from_clm <- function(clm_obj) {
  
  message("--- Converting 'clm' object to a STAAR-compatible score test framework ---")
  
  # --- Step 1: Extract Null Model Parameters and Calculate Probabilities ---
  message("Step 1: Calculating predictors and probabilities under the null model...")
  
  model_data <- clm_obj$model
  if (is.null(model_data)) {
    stop("The 'clm' object must be fitted with 'model = TRUE'.")
  }
  
  X_matrix <- model.matrix(object = formula(clm_obj), data = model_data)
  beta_names <- names(clm_obj$beta)
  X_no_intercept <- X_matrix[, beta_names, drop = FALSE]
  latent_predictor <- as.vector(X_no_intercept %*% clm_obj$beta)
  
  alpha_coeffs <- clm_obj$alpha
  link_obj <- make.link(clm_obj$link)
  link_func_inv <- link_obj$linkinv
  cum_probs_list <- lapply(alpha_coeffs, function(a) {
    link_func_inv(a - latent_predictor)
  })
  cum_probs_matrix <- do.call(cbind, cum_probs_list)
  
  cum_probs_padded <- cbind(0, cum_probs_matrix, 1)
  K <- length(clm_obj$y.levels)
  pred_probs <- cum_probs_padded[, 2:(K + 1)] - cum_probs_padded[, 1:K]
  
  pred_probs[pred_probs < 1e-12] <- 0
  pred_probs <- pred_probs / rowSums(pred_probs)
  
  # --- Step 2: Calculate Working Variables for the Score Test ---
  message("Step 2: Calculating working residuals and dependent variable...")
  
  y_numeric <- as.numeric(clm_obj$y) 
  categories <- 1:K 
  expected_y <- as.vector(pred_probs %*% categories)
  working_residual <- y_numeric - expected_y
  z_working_variable <- latent_predictor + working_residual
  
  # --- Step 3: Calculate Weights for the Pseudo-Linear Model ---
  message("Step 3: Calculating weights...")
  
  expected_y_squared <- as.vector(pred_probs %*% (categories^2))
  var_y <- expected_y_squared - (expected_y^2)
  model_weights <- 1 / (var_y + 1e-8)
  dispersion_est <- 1.0
  
  # --- Step 4: Assemble the Final STAAR-Compatible Object ---
  message("Step 4: Assembling the final list object...")
  
  stopifnot(
    "Dimension mismatch between design matrix and working variable" = 
      nrow(X_matrix) == length(z_working_variable),
    "Dimension mismatch between design matrix and weights" = 
      nrow(X_matrix) == length(model_weights)
  )
  
  null_obj <- list()
  
  null_obj$y <- as.vector(z_working_variable)
  null_obj$weights <- as.vector(model_weights)
  null_obj$residuals <- as.vector(working_residual)
  null_obj$fitted.values <- as.vector(latent_predictor)
  
  null_obj$family <- gaussian(link = "identity")
  null_obj$dispersion <- dispersion_est
  null_obj$formula <- formula(clm_obj)
  null_obj$terms <- clm_obj$terms
  null_obj$model <- clm_obj$model
  null_obj$relatedness <- FALSE
  
  class(null_obj) <- c("POM", "glm", "lm")
  
  message("--- Conversion complete. Object of class 'POM' created successfully. ---")
  return(null_obj)
}


# #########################################################################
# The main user-facing function now uses the professional-grade converter
# #########################################################################

#' @title Fit a Null Model for STAAR
#' @description A wrapper to fit a null model using either a standard GLM
#'   (for binary/continuous traits) or an ordinal regression model which is then
#'   converted for use with STAAR.
#'
#' @param fixed A formula specifying the null model (e.g., `pheno ~ covar1 + covar2`).
#' @param data A data frame containing the variables in the formula.
#' @param family For GLM, a family object (e.g., `binomial()`). For ordinal
#'   regression, the string `"ordinal"`.
#' @param ... Additional arguments passed to `stats::glm` or `ordinal::clm`.
#'
#' @return A null model object compatible with STAAR.
fit_null_model <- function(fixed, data, family = binomial(link = "logit"), ...){
  
  is_ordinal <- is.character(family) && family == "ordinal"
  
  if (is_ordinal) {
    message("Fitting an ordinal regression model and preparing for STAAR...")
    outcome_var_name <- as.character(fixed[[2]])
    
    if (!is.ordered(data[[outcome_var_name]])) {
      stop(paste("For 'ordinal' family, the outcome variable '", outcome_var_name, "' must be an ordered factor."))
    }
    
    clm_fit <- ordinal::clm(formula = fixed, data = data, model = TRUE, ...)
    obj_nullmodel <- create_null_from_clm(clm_obj = clm_fit)
    
  } else {
    message("Fitting a generalized linear model (glm)...")
    obj_nullmodel <- stats::glm(formula = fixed, data = data, family = family, ...)
    obj_nullmodel$relatedness <- FALSE
  }
  
  return(obj_nullmodel)
}