#' @title Calculate Score Statistics for an Ordinal Model (Final Robust Version)
#' @description This function calculates score statistics for single variants and is the core
#'   computational engine for set-based tests. It includes robust handling of potentially
#'   singular covariance matrices.
#' @keywords internal
Ordinal_exactScore <- function(objNull, G_mat, use_SPA = NULL, SPA_filter = TRUE, SPA_filter_cutoff = 0.05, verbose = FALSE) {

  # Calculate the score vector U = G' * W * residuals
  Score <- as.vector(crossprod(G_mat, objNull$W_mat %*% objNull$residuals))

  # Calculate the main components of the variance-covariance matrix V
  G_prime_W_G <- crossprod(G_mat, objNull$W_mat %*% G_mat)
  G_prime_W_X <- crossprod(G_mat, objNull$WX_mat)
  Correction_Term <- G_prime_W_X %*% objNull$XWX_inv %*% t(G_prime_W_X)
  Var_mat <- G_prime_W_G - Correction_Term

  # --- [CRITICAL FIX] ---
  # Regularize the covariance matrix by adding a small positive value ("nugget")
  # to its diagonal. This ensures the matrix is strictly positive-definite and
  # prevents the variance of aggregate tests (like Burden) from collapsing to zero
  # due to high correlation/collinearity between variants.
  diag(Var_mat) <- diag(Var_mat) + 1e-6
  # --- [END OF FIX] ---

  # Extract the diagonal elements for single-variant tests
  Variance <- diag(Var_mat)

  # A secondary safeguard to ensure no individual variance is zero or negative
  # This is somewhat redundant after the nugget fix but kept for extra safety.
  Variance[Variance <= 1e-8] <- 1e-8

  # Calculate single-variant score test statistics and p-values
  Stest <- Score^2 / Variance
  p_value <- pchisq(Stest, df = 1, lower.tail = FALSE)

  # Assemble the results data frame
  result_df <- data.frame(Score = Score, Variance = Variance, Stest = Stest, Pvalue = p_value)

  # --- Saddlepoint Approximation (SPA) Section ---
  if (is.null(use_SPA)) use_SPA <- objNull$use_SPA

  if (use_SPA) {
    result_df$Pvalue_SPA <- result_df$Pvalue
    SPA_index <- if(SPA_filter) which(p_value < SPA_filter_cutoff) else seq_along(p_value)

    if (length(SPA_index) > 0) {
      if (verbose) print(paste0("Single variants analysis: apply SPA to ", length(SPA_index), " markers"))

      G_tilde <- G_tilde_forSPA_Ordinal(G_mat, objNull)

      p_spa_values <- sapply(SPA_index, function(i) {
        # Pass the robust variance to the SPA function
        single_SPA_Ordinal(G_tilde[, i], Score[i], Variance[i], objNull)
      })
      result_df$Pvalue_SPA[SPA_index] <- p_spa_values
    }
  }

  # --- Final Formatting Section ---
  Est <- result_df$Score / result_df$Variance
  Est_se <- 1 / sqrt(result_df$Variance)

  if (use_SPA) {
    result_df$pvalue_log10 <- -log10(result_df$Pvalue_SPA)
    result_df <- result_df[, c("Pvalue", "Pvalue_SPA", "pvalue_log10", "Score", "Variance", "Stest"), drop = FALSE]
  } else {
    result_df$pvalue_log10 <- -log10(result_df$Pvalue)
    result_df <- result_df[, c("Pvalue", "pvalue_log10", "Stest", "Score", "Variance"), drop = FALSE]
  }
  result_df$Est <- Est
  result_df$Est_se <- Est_se

  # Return all necessary components for downstream aggregate tests
  results <- list(result = result_df, Score = Score, Covariance = Var_mat)
  return(results)
}





genoMatrixPlink = function(Geno, bim_data = bim_data, markerIndex,
                           geno_missing_cutoff = 1, geno_missing_imputation = c("mean", "minor"),
                           min_mac_cutoff = NULL, min_maf_cutoff = 0.01) {

  Geno[Geno == -9 | Geno == 9] = NA
  colnames(Geno) = bim_data$V2[markerIndex]
  Geno_col = ncol(Geno)
  Geno_row = nrow(Geno)

  G_summary = list()
  G_na = list()

  ## calculate MAF MAC missing and flip
  for (p in 1:Geno_col) {
    Geno_p = Geno[, p]

    MAF = mean(Geno_p, na.rm = T)/2
    if (MAF > 0.5) {
      MAF = 1-MAF
      Geno_p = 2-Geno_p
      Geno[, p] = Geno_p
    }
    MAC = sum(Geno_p, na.rm = T)
    missing_index = which(is.na(Geno_p))
    missing_rate = length(missing_index)/Geno_row
    G_summary[[p]] = data.frame(MAF, MAC, missing_rate)
    G_na[[p]] = missing_index
  }

  G_summary = data.table::rbindlist(G_summary)


  ## variants filter
  if (is.null(min_mac_cutoff)) {
    include_index = which(G_summary$missing_rate <= geno_missing_cutoff & G_summary$MAF >= min_maf_cutoff)
  } else {
    maf_mac_cutoff = Geno_row * min_maf_cutoff * 2
    if (maf_mac_cutoff < min_mac_cutoff) {
      include_index = which(G_summary$missing_rate <= geno_missing_cutoff & G_summary$MAF >= min_maf_cutoff)
    } else {
      include_index = which(G_summary$missing_rate <= geno_missing_cutoff & G_summary$MAC >= min_mac_cutoff)
    }
  }

  Geno = Geno[, include_index, drop = FALSE]
  G_summary = G_summary[include_index, ]
  G_na = G_na[include_index]


  ## imputation
  switch (geno_missing_imputation,
          "mean" = {
            for (p in 1:length(include_index)) {
              na_index = G_na[[p]]
              imput_p = G_summary$MAF[p]*2
              Geno[na_index, p] = imput_p
            }
          },
          "minor" = {
            for (p in 1:length(include_index)) {
              na_index = G_na[[p]]
              Geno[na_index, p] = 0
            }
          }
  )

  G_summary = cbind(bim_data[markerIndex[include_index], c(1,4:6,2)], G_summary)
  colnames(G_summary)[1:5] = c("CHR", "POS", "REF1", "REF2", "rsID")

  ## update imputed MAC
  G_summary$MAC = colSums(Geno)

  result = list(Geno = Geno, G_summary = G_summary)

  return(result)
}


genoMatrixGDS = function(Geno, G_info, geno_missing_cutoff = 1, geno_missing_imputation = c("mean", "minor"),
                         min_mac_cutoff = NULL, min_maf_cutoff = 0.01) {

  Geno_col = ncol(Geno)
  Geno_row = nrow(Geno)

  G_summary = list()
  G_na = list()

  ## calculate MAF MAC missing and flip
  for (p in 1:Geno_col) {
    Geno_p = Geno[, p]

    AF = mean(Geno_p, na.rm = T)/2
    ALT_AF = 1-AF
    MAF = AF
    if (MAF > 0.5) {
      MAF = 1-MAF
      Geno_p = 2-Geno_p
      Geno[, p] = Geno_p
    }
    MAC = sum(Geno_p, na.rm = T)
    missing_index = which(is.na(Geno_p))
    missing_rate = length(missing_index)/Geno_row
    G_summary[[p]] = data.frame(MAF, MAC, missing_rate)
    G_na[[p]] = missing_index
  }

  G_summary = data.table::rbindlist(G_summary)
  G_summary = cbind(G_info, G_summary)

  ## variants filter
  if (is.null(min_mac_cutoff)) {
    include_index = which(G_summary$missing_rate <= geno_missing_cutoff & G_summary$MAF >= min_maf_cutoff)
  } else {
    maf_mac_cutoff = Geno_row * min_maf_cutoff * 2
    if (maf_mac_cutoff < min_mac_cutoff) {
      include_index = which(G_summary$missing_rate <= geno_missing_cutoff & G_summary$MAF >= min_maf_cutoff)
    } else {
      include_index = which(G_summary$missing_rate <= geno_missing_cutoff & G_summary$MAC >= min_mac_cutoff)
    }
  }

  Geno = Geno[, include_index, drop = FALSE]
  G_summary = G_summary[include_index, ]
  G_na = G_na[include_index]


  ## imputation
  switch (geno_missing_imputation,
          "mean" = {
            for (p in 1:length(include_index)) {
              na_index = G_na[[p]]
              imput_p = G_summary$MAF[p]*2
              Geno[na_index, p] = imput_p
            }
          },
          "minor" = {
            for (p in 1:length(include_index)) {
              na_index = G_na[[p]]
              Geno[na_index, p] = 0
            }
          }
  )

  ## update imputed MAC
  G_summary$MAC = colSums(Geno)

  result = list(Geno = Geno, G_summary = G_summary, include_index = include_index)

  return(result)
}


genoFlipRV = function(Geno, geno_missing_imputation = c("mean", "minor"),
                      geno_missing_cutoff = 0.1, min_maf_cutoff = 1e-4, rare_maf_cutoff = 0.01, rare_num_cutoff = 2) {

  Geno_col = ncol(Geno)
  Geno_row = nrow(Geno)

  G_summary = list()
  G_na = list()

  ## calculate MAF, MAC, missing, and flip
  for (p in 1:Geno_col) {
    Geno_p = Geno[, p]

    MAF = mean(Geno_p, na.rm = T)/2
    if (MAF > 0.5) {
      MAF = 1-MAF
      Geno_p = 2-Geno_p
      Geno[, p] = Geno_p
    }
    MAC = sum(Geno_p, na.rm = T)
    missing_index = which(is.na(Geno_p))
    missing_rate = length(missing_index)/Geno_row
    G_summary[[p]] = data.frame(MAF, MAC, missing_rate)
    G_na[[p]] = missing_index
  }
  G_summary = data.table::rbindlist(G_summary)


  ## variants filter
  include_index = which(G_summary$missing_rate <= geno_missing_cutoff & G_summary$MAF >= min_maf_cutoff & G_summary$MAF <= rare_maf_cutoff)
  Geno = Geno[, include_index, drop = FALSE]
  G_summary = G_summary[include_index, ]
  G_na = G_na[include_index]


  ## check for rare variants number
  if(length(include_index) < rare_num_cutoff) {
    stop('Number of variants in this gene is less than ', rare_num_cutoff, ", will skip this category...", call. = FALSE)
  }


  ## imputation
  switch (geno_missing_imputation,
          "mean" = {
            for (p in 1:length(include_index)) {
              na_index = G_na[[p]]
              imput_p = G_summary$MAF[p]*2
              Geno[na_index, p] = imput_p
            }
          },
          "minor" = {
            for (p in 1:length(include_index)) {
              na_index = G_na[[p]]
              Geno[na_index, p] = 0
            }
          }
  )

  ## update imputed MAC
  G_summary$MAC = colSums(Geno)

  result = list(Geno = Geno, G_summary = G_summary, include_index = include_index)

  return(result)

}



genoFlip = function(Geno) {

  Geno_col = ncol(Geno)
  Geno_row = nrow(Geno)

  G_summary = list()
  G_na = list()

  ## calculate MAF MAC missing and flip
  for (p in 1:Geno_col) {
    Geno_p = Geno[, p]

    MAF = mean(Geno_p, na.rm = T)/2
    if (MAF > 0.5) {
      MAF = 1-MAF
      Geno_p = 2-Geno_p
      Geno[, p] = Geno_p
    }
    MAC = sum(Geno_p, na.rm = T)
    missing_index = which(is.na(Geno_p))
    missing_rate = length(missing_index)/Geno_row
    G_summary[[p]] = data.frame(MAF, MAC)
    G_na[[p]] = missing_index
  }
  G_summary = data.table::rbindlist(G_summary)

  ## imputation  minor
  for (p in 1:ncol(Geno)) {
    na_index = G_na[[p]]
    Geno[na_index, p] = 0
  }

  ## update imputed MAC
  G_summary$MAC = colSums(Geno)

  result = list(Geno = Geno, G_summary = G_summary)

  return(result)

}


#' @export
rap_load_as = function(file, name) {
  # file <- basename(file)
  data_list <- load(file)
  assign(name, get(data_list[1]), envir = .GlobalEnv)
}


#' @export
argsReshape = function(default_args, args, num_args, log_args) {

  args_list = default_args

  for (arg in args) {

    key_value <- strsplit(arg, "=")[[1]]

    if (length(key_value) == 2) {
      arg_name = sub("--", "", key_value[1])

      if (arg_name %in% names(default_args)) {

        if (arg_name %in% num_args) {

          value <- suppressWarnings(as.numeric(key_value[2]))
          if (is.numeric(value)) {
            args_list[[arg_name]] <- value
          } else {
            warning(paste("Invalid value: ", arg_name, ". A numeric value is required."))
          }

        } else if (arg_name %in% log_args) {

          value <- suppressWarnings(as.logical(key_value[2]))
          if (is.logical(value)) {
            args_list[[arg_name]] <- value
          } else {
            warning(paste("Invalid value: ", arg_name, ". A logical value is required."))
          }

        } else {
          value <- key_value[2]
          args_list[[arg_name]] <- value
        }

      } else {
        warning(paste("Unknown argument: ", arg_name))
      }
    }
  }
  return(args_list)
}











### Following functions are sourced from the GitHub repository SPACox (https://github.com/WenjianBI/SPACox)

CGF4LatentRes <- function(fit_null, range = c(-100, 100), length.out = 1e4, verbose = TRUE) {

  fit_null$use_SPA <- TRUE
  idx0 <- qcauchy(seq_len(length.out) / (length.out + 1))
  idx1 <- idx0 * max(range) / max(idx0)
  res_null <- fit_null$residuals

  if (verbose) print("Start calculating empirical CGF for latent variable residuals...")

  cumul_list <- lapply(idx1, function(t) {
    exp_res <- exp(res_null * t)
    if (any(is.infinite(exp_res))) return(c(t = t, K0 = NA, K1 = NA, K2 = NA))

    M0 <- mean(exp_res)
    M1 <- mean(res_null * exp_res)
    M2 <- mean(res_null^2 * exp_res)

    K0 <- log(M0)
    K1 <- M1 / M0
    K2 <- (M0 * M2 - M1^2) / M0^2

    return(c(t = t, K0 = K0, K1 = K1, K2 = K2))
  })

  cumul <- do.call(rbind, cumul_list)
  cumul <- cumul[complete.cases(cumul), ]

  if (verbose) print("CGF calculation complete. Creating interpolation functions...")

  fit_null$CGF_INFO <- list(
    K_org_emp = approxfun(cumul[, "t"], cumul[, "K0"], rule = 2),
    K_1_emp = approxfun(cumul[, "t"], cumul[, "K1"], rule = 2),
    K_2_emp = approxfun(cumul[, "t"], cumul[, "K2"], rule = 2)
  )

  return(fit_null)
}





#' @title 为 OrdinalSTAAR 计算经过协变量校正的基因型 (G_tilde)
#' @description 本函数计算在有序 Probit 零模型下，经过协变量校正后的基因型
#'   矩阵 G_tilde。它使用了加权最小二乘法的投影，其组件已在 Ordinal_NullModel
#'   中预先计算好。
#'
#' @param G 一个 n x p 的基因型矩阵 (n=样本数, p=变异数)。
#' @param objNull 从 `Ordinal_NullModel` 返回的零模型对象，必须包含
#'   `X_mat`, `W_mat`, `XWX_inv` 等组件。
#'
#' @return 一个 n x p 的、经过协变量校正的基因型矩阵 G_tilde。
#'
G_tilde_forSPA_Ordinal <- function(G, objNull) {

  # --- 1. 输入验证 ---
  required_components <- c("X_mat", "W_mat", "XWX_inv")
  if (!all(required_components %in% names(objNull))) {
    stop("The null model object 'objNull' is missing required components for G_tilde calculation (X_mat, W_mat, XWX_inv).")
  }

  # 确保 G 是一个标准的 matrix，以避免稀疏矩阵的兼容性问题
  if (inherits(G, "sparseMatrix")) {
    G <- as.matrix(G)
  }

  # --- 2. 从 objNull 中解包预计算的矩阵 ---
  X <- objNull$X_mat         # n x k 设计矩阵
  W <- objNull$W_mat         # n x n 对角权重矩阵
  XWX_inv <- objNull$XWX_inv  # k x k 的 (X'WX)^-1 矩阵

  # --- 3. 执行加权投影的矩阵代数运算 ---

  # 计算 X'WG
  # X is (n x k), W is (n x n), G is (n x p)
  # t(X) is (k x n)
  # W %*% G is (n x p)
  # t(X) %*% (W %*% G) gives a (k x p) matrix
  X_t_W_G <- crossprod(X, W %*% G)

  # 计算 X * (X'WX)^-1 * (X'WG)
  # X is (n x k)
  # XWX_inv is (k x k)
  # X_t_W_G is (k x p)
  # The result is an (n x p) matrix, representing the projected part of G
  G_projected <- X %*% (XWX_inv %*% X_t_W_G)

  # --- 4. 计算 G_tilde = G - G_projected ---
  G_tilde <- G - G_projected

  return(G_tilde)
}




#' @title 使用鞍点近似计算序数模型的p值
#' @description 这是为序数模型定制的SPA p值计算核心，使用Lugananni-Rice公式。
#' @keywords internal
GetProb_SPA_Ordinal <- function(Score, objNull, G_w) {

  if (is.null(objNull$CGF_INFO)) stop("CGF information is missing in objNull for SPA.")

  G_w_vec <- as.vector(G_w)

  # 构建CGF的一阶导数方程
  K_q_deriv1 <- function(t, q_obs) {
    sum(G_w_vec * objNull$CGF_INFO$K_1_emp(t * G_w_vec)) - q_obs
  }

  # 求解鞍点 zeta
  saddle_root <- try(uniroot(K_q_deriv1, interval = c(-100, 100), q_obs = Score, extendInt = "yes", tol = 1e-8), silent = TRUE)

  # 如果求解失败，退回使用卡方检验
  if (inherits(saddle_root, "try-error")) {
    variance_approx <- sum(G_w_vec^2 * objNull$CGF_INFO$K_2_emp(0)) # K''(0) is variance
    return(pchisq(Score^2 / variance_approx, df = 1, lower.tail = FALSE))
  }

  zeta <- saddle_root$root
  if (abs(zeta) < 1e-6) { # 如果鞍点接近0，说明正态近似很好
    variance <- sum(G_w_vec^2 * objNull$CGF_INFO$K_2_emp(0))
    return(pchisq(Score^2 / variance, df = 1, lower.tail = FALSE))
  }

  # 计算鞍点处的CGF及其二阶导
  K_q_zeta <- sum(objNull$CGF_INFO$K_org_emp(zeta * G_w_vec))
  K_q_deriv2_zeta <- sum(G_w_vec^2 * objNull$CGF_INFO$K_2_emp(zeta * G_w_vec))

  # 应用Lugananni-Rice公式
  w <- sign(zeta) * sqrt(2 * (zeta * Score - K_q_zeta))
  v <- zeta * sqrt(K_q_deriv2_zeta)

  if (is.na(w) || is.na(v) || abs(w) < 1e-6 || abs(v) < 1e-6) return(1.0)

  # 返回双侧p值
  pval <- pnorm(abs(w), lower.tail = FALSE) * 2
  # 更精确的近似
  # pval <- (pnorm(abs(w), lower.tail=FALSE) + dnorm(abs(w))*(1/abs(v) - 1/abs(w)))*2

  return(min(1, pval))
}






### Following function is sourced from the GitHub repository STAAR (https://github.com/xihaoli/STAAR/blob/master/R/CCT.R)

CCT <- function(pvals, weights=NULL){
  #### check if there is NA
  if(sum(is.na(pvals)) > 0){
    stop("Cannot have NAs in the p-values!")
  }

  #### check if all p-values are between 0 and 1
  if((sum(pvals<0) + sum(pvals>1)) > 0){
    stop("All p-values must be between 0 and 1!")
  }

  #### check if there are p-values that are either exactly 0 or 1.
  is.zero <- (sum(pvals==0)>=1)
  is.one <- (sum(pvals==1)>=1)
  if(is.zero && is.one){
    stop("Cannot have both 0 and 1 p-values!")
  }
  if(is.zero){
    return(0)
  }
  if(is.one){
    warning("There are p-values that are exactly 1!")
    return(1)
  }

  #### check the validity of weights (default: equal weights) and standardize them.
  if(is.null(weights)){
    weights <- rep(1/length(pvals),length(pvals))
  }else if(length(weights)!=length(pvals)){
    stop("The length of weights should be the same as that of the p-values!")
  }else if(sum(weights < 0) > 0){
    stop("All the weights must be positive!")
  }else{
    weights <- weights/sum(weights)
  }

  #### check if there are very small non-zero p-values
  is.small <- (pvals < 1e-16)
  if (sum(is.small) == 0){
    cct.stat <- sum(weights*tan((0.5-pvals)*pi))
  }else{
    cct.stat <- sum((weights[is.small]/pvals[is.small])/pi)
    cct.stat <- cct.stat + sum(weights[!is.small]*tan((0.5-pvals[!is.small])*pi))
  }

  #### check if the test statistic is very large.
  if(cct.stat > 1e+15){
    pval <- (1/cct.stat)/pi
  }else{
    pval <- 1-pcauchy(cct.stat)
  }
  return(pval)
}


### Following function is sourced from the GitHub repository OrdinalSTAAR (https://github.com/Cui-yd/OrdinalSTAAR/blob/main/R/basicFunction.R)

#' @title SPA封装函数
#' @keywords internal
single_SPA_Ordinal <- function(G_w, Score, Variance, objNull) {
  # 直接调用新的核心计算函数
  pval <- GetProb_SPA_Ordinal(Score = Score, objNull = objNull, G_w = G_w)
  return(pval)
}




#' @title 有序模型的 Burden 检验 (支持SPA, 修正版)
#' @keywords internal
OrdinalBurden <- function(Geno, Score, Covariance, weight, objNull, use_SPA = NULL, SPA_filter = TRUE, SPA_filter_cutoff = 0.05, verbose = FALSE) {

  if (is.null(use_SPA)) use_SPA <- objNull$use_SPA

  pval_B <- NULL
  G_tilde <- if(use_SPA) G_tilde_forSPA_Ordinal(Geno, objNull) else NULL

  for (k in 1:ncol(weight)) {
    weight_k <- weight[, k]
    Score_k <- sum(Score * weight_k)
    Variance_k <- as.vector(t(weight_k) %*% Covariance %*% weight_k)

    # --- [关键修正] ---
    # 添加保护，防止方差为0或负数
    if (Variance_k <= 1e-8) {
      Variance_k <- 1e-8
    }
    # --- [修正结束] ---

    pval_k <- pchisq(Score_k^2 / Variance_k, df = 1, lower.tail = FALSE)

    if (use_SPA && (!SPA_filter || (SPA_filter && pval_k < SPA_filter_cutoff))) {
      G_tilde_w_k <- G_tilde %*% weight_k
      # 注意：传递给SPA的方差也应该是修正后的
      pval_k <- single_SPA_Ordinal(G_tilde_w_k, Score_k, Variance_k, objNull)
    }

    pval_B <- c(pval_B, pval_k)
  }

  return(matrix(pval_B, nrow = 2, byrow = TRUE))
}


#' @title 有序模型的 SKAT (与SurvSTAAR完全兼容版, 逻辑重构)
#' @keywords internal
OrdinalSKAT <- function(Geno, Score, Covariance, Pvalue, MAC = NULL, weight_S, weight_B,
                        objNull, use_SPA = NULL, SPA_filter = TRUE, SPA_filter_cutoff = 0.05,
                        combine_ultra_rare = TRUE, ultra_rare_mac_cutoff = 60, verbose = FALSE) {

  # --- 1. 输入验证与初始化 ---
  if (!inherits(Geno, c("matrix", "Matrix"))) stop("基因型必须是一个矩阵！")
  if (is.null(use_SPA)) use_SPA <- objNull$use_SPA
  if (is.null(MAC)) MAC <- colSums(Geno)

  ultra_rare_index <- which(MAC <= ultra_rare_mac_cutoff)
  pval_S <- NULL

  # --- 2. 根据超稀有变异的情况选择分析策略 ---

  # 策略1：不合并，或超稀有变异太少 -> 标准SKAT
  if (!combine_ultra_rare || length(ultra_rare_index) <= 1) {
    if (verbose) print("             在SKAT检验中不合并超稀有变异")
    for (k in 1:ncol(weight_S)) {
      weight_k <- weight_S[, k]
      Qtest <- sum(Score^2 * weight_k^2)
      w_k_mat <- diag(weight_k)
      Cov_weight <- w_k_mat %*% Covariance %*% w_k_mat
      if (sum(abs(Cov_weight)) == 0) { pval_k <- NA } else {
        lambda <- eigen(Cov_weight, only.values = TRUE, symmetric = TRUE)$values
        lambda <- lambda[lambda > 1e-6]
        pval_k <- tryCatch({ p_davies <- CompQuadForm::davies(Qtest, lambda)$Qq; if (is.na(p_davies) || p_davies <= 0 || p_davies > 1) CompQuadForm::liu(Qtest, lambda) else p_davies }, error = function(e) CompQuadForm::liu(Qtest, lambda))
      }
      pval_S <- c(pval_S, pval_k)
    }

    # 策略2：所有变异都是超稀有 -> 直接退化为Burden检验
  } else if (length(ultra_rare_index) == ncol(Geno)) {
    if (verbose) print("             所有变异均为超稀有，SKAT退化为Burden检验")
    pval_S <- OrdinalBurden(Geno, Score, Covariance, weight_B, objNull, use_SPA, SPA_filter, SPA_filter_cutoff, verbose = FALSE)

    # 策略3：混合情况 -> 逻辑重构
  } else {
    if (verbose) print(paste0("             在SKAT检验中合并 ", length(ultra_rare_index), " 个超稀有变异"))

    # 准备普通变异和超稀有变异的子集
    Geno_sub <- Geno[, -ultra_rare_index, drop = FALSE]
    Geno_rare <- Geno[, ultra_rare_index, drop = FALSE]

    weight_S_sub <- weight_S[-ultra_rare_index, , drop = FALSE]
    weight_B_rare <- weight_B[ultra_rare_index, , drop = FALSE]

    for (k in 1:ncol(weight_S)) {
      # 1. 创建“超稀有超级变异”：将超稀有变异按Burden权重加权合并
      weight_B_rare_k <- weight_B_rare[, k]
      G_burden_rare <- as.matrix(Geno_rare %*% weight_B_rare_k)

      # 2. 构建新的基因型矩阵
      G_new <- cbind(G_burden_rare, Geno_sub)

      # 3. 构建新的SKAT权重向量
      # “超级变异”的权重是原来那些变异权重的均值
      weight_S_new_k <- c(mean(weight_S[ultra_rare_index, k]), weight_S_sub[, k])

      # 4. 在新的“组合基因型”上，重新应用标准的SKAT公式
      # Q = (G_new' * w_new * res)' * (G_new' * w_new * res)
      #   = res' * w_new * G_new * G_new' * w_new * res
      # 我们需要的是 Q = U_new' * diag(w_new)^2 * U_new
      # U_new = G_new' * res

      # 重新计算这个新系统的 Score 和 Covariance
      # 这是最稳健的方法
      new_test <- Ordinal_exactScore(objNull, G_new)
      Score_new <- new_test$Score
      Covariance_new <- new_test$Covariance

      # 5. 在这个新的 m+1 系统上进行SKAT检验
      Qtest_k <- sum(Score_new^2 * weight_S_new_k^2)
      w_k_mat <- diag(weight_S_new_k)
      Cov_weight_k <- w_k_mat %*% Covariance_new %*% w_k_mat

      if (sum(abs(Cov_weight_k)) == 0) {
        pval_k <- NA
      } else {
        lambda_k <- eigen(Cov_weight_k, only.values = TRUE, symmetric = TRUE)$values
        lambda_k <- lambda_k[lambda_k > 1e-6]

        pval_k <- tryCatch({ p_davies <- CompQuadForm::davies(Qtest_k, lambda_k)$Qq; if (is.na(p_davies) || p_davies <= 0 || p_davies > 1) CompQuadForm::liu(Qtest_k, lambda_k) else p_davies }, error = function(e) CompQuadForm::liu(Qtest_k, lambda_k))
      }
      pval_S <- c(pval_S, pval_k)
    }
  }

  # --- 3. 格式化输出 ---
  pval_S <- matrix(pval_S, nrow = 2, byrow = TRUE)
  return(pval_S)
}

OrdinalACAT <- function(Geno, Score, Covariance, Pvalue, MAC = NULL,use_SPA = NULL, weight_A, weight_B,
                        objNull, combine_ultra_rare = TRUE, ultra_rare_mac_cutoff = 20, verbose = FALSE) {

  if (!inherits(Geno, c("matrix", "Matrix"))) stop("基因型必须是一个矩阵！")
  if (is.null(MAC)) MAC <- colSums(Geno)

  ultra_rare_index <- which(MAC <= ultra_rare_mac_cutoff)
  pval_A <- NULL

  # 策略1：不合并，或超稀有变异太少
  if (!combine_ultra_rare || length(ultra_rare_index) <= 1) {
    if (verbose) print("             在ACAT检验中应用标准柯西组合")
    for (k in 1:ncol(weight_A)) {
      pval_k <- CCT(Pvalue, weight_A[, k])
      pval_A <- c(pval_A, pval_k)
    }
    # 策略2：所有变异都是超稀有 -> 直接退化为Burden检验
  } else if (length(ultra_rare_index) == ncol(Geno)) {
    if (verbose) print("             所有变异均为超稀有，ACAT退化为Burden检验")
    pval_A <- OrdinalBurden(Geno, Score, Covariance, weight_B, objNull,
                            use_SPA, SPA_filter, SPA_filter_cutoff, verbose = FALSE)
    # 策略3：混合情况 -> 合并超稀有变异
  } else {
    if (verbose) print(paste0("             在ACAT检验中合并 ", length(ultra_rare_index), " 个超稀有变异"))

    # 1. 对超稀有变异子集做Burden检验
    G_ultra_rare <- Geno[, ultra_rare_index, drop = FALSE]
    Score_rare <- Score[ultra_rare_index]
    Covariance_rare <- Covariance[ultra_rare_index, ultra_rare_index, drop = FALSE]
    weight_B_rare <- weight_B[ultra_rare_index, , drop = FALSE]

    Pval_rare_B <- OrdinalBurden(G_ultra_rare, Score = Score_rare, Covariance = Covariance_rare,
                                 weight = weight_B_rare, objNull = objNull, verbose = FALSE)

    # 2. 准备组合检验的输入
    Pvalue_sub <- Pvalue[-ultra_rare_index]
    # 对权重进行调整：超稀有变异的权重合并为1个，其余不变
    weight_A_new <- rbind(colMeans(weight_A[ultra_rare_index, , drop = FALSE]), weight_A[-ultra_rare_index, , drop = FALSE])

    for (k in 1:ncol(weight_A_new)) {
      # 将Burden p值和剩余p值打包
      Pvalue_k <- c(Pval_rare_B[k], Pvalue_sub)
      weight_k <- weight_A_new[, k]

      pval_k <- CCT(Pvalue_k, weight_k)
      pval_A <- c(pval_A, pval_k)
    }
  }

  pval_A <- matrix(pval_A, nrow = 2, byrow = TRUE)
  return(pval_A)
}




OrdinalSTAAR_O = function(Geno, objNull, annotation_rank = NULL, MAC = NULL,
                       use_SPA = NULL, SPA_filter = TRUE, SPA_filter_cutoff = 0.05,
                       weight_A, weight_B, weight_S,
                       combine_ultra_rare = TRUE, ultra_rare_mac_cutoff = 20, verbose = FALSE) {

  if (!inherits(Geno, "matrix") && !inherits(Geno, "Matrix")) stop("Genotype is not a matrix")

  if (ncol(Geno) != nrow(weight_A) | ncol(Geno) != nrow(weight_B) | ncol(Geno) != nrow(weight_S)) stop("Dimensions don't match for genotype and weights")

  if (is.null(use_SPA)) use_SPA = objNull$use_SPA

  if (is.null(MAC)) MAC = colSums(Geno)


  ### individual test
  single_test = Ordinal_exactScore(objNull = objNull, G_mat = Geno, use_SPA = use_SPA, SPA_filter = SPA_filter, SPA_filter_cutoff = SPA_filter_cutoff, verbose = verbose)
  Score = single_test$Score
  Covariance = single_test$Covariance
  if(use_SPA) {
    Pvalue = single_test$result$Pvalue_SPA
  } else {
    Pvalue = single_test$result$Pvalue
  }


  ### set based test
  if (verbose) print(paste0("SKAT test:   begin at ", Sys.time()))
  Pvalue_S = OrdinalSKAT(Geno, Score, Covariance, Pvalue, MAC, weight_S, weight_B, objNull, combine_ultra_rare, ultra_rare_mac_cutoff, verbose)

  if (verbose) print(paste0("ACAT test:   begin at ", Sys.time()))
  Pvalue_A = OrdinalACAT(Geno, Score, Covariance, Pvalue, MAC, weight_A, weight_B, objNull, combine_ultra_rare, ultra_rare_mac_cutoff, verbose)

  if (verbose) print(paste0("Burden test: begin at ", Sys.time()))
  Pvalue_B = OrdinalBurden(Geno, Score, Covariance, weight_B, objNull, verbose)


  ### combine results
  results_pvalue = rbind(Pvalue_S, Pvalue_A, Pvalue_B)
  rownames(results_pvalue) = c("SKAT-(1,25)", "SKAT-(1,1)", "ACAT-(1,25)", "ACAT-(1,1)", "Burden-(1,25)", "Burden-(1,1)")

  if (!is.null(annotation_rank)) {
    colnames(results_pvalue) = c("Beta", colnames(annotation_rank))
  } else {
    colnames(results_pvalue) = "results"
  }


  results_STAAR_O = CCT(c(results_pvalue))
  results_ACAT_O = CCT(c(results_pvalue[, 1]))
  results_STAAR_S = CCT(results_pvalue[1:2, ])
  results_STAAR_A = CCT(results_pvalue[3:4, ])
  results_STAAR_B = CCT(results_pvalue[5:6, ])

  STAAR_S_1_25 = CCT(Pvalue_S[1, ])
  STAAR_S_1_1  = CCT(Pvalue_S[2, ])
  STAAR_A_1_25 = CCT(Pvalue_A[1, ])
  STAAR_A_1_1  = CCT(Pvalue_A[2, ])
  STAAR_B_1_25 = CCT(Pvalue_B[1, ])
  STAAR_B_1_1  = CCT(Pvalue_B[2, ])

  Combined_P = c(STAAR_S_1_25, STAAR_S_1_1, STAAR_A_1_25, STAAR_A_1_1, STAAR_B_1_25, STAAR_B_1_1)
  results_pvalue = cbind(results_pvalue, Combined_P)


  results = list("OrdinalSTAAR_O" = results_STAAR_O,
                 "ACAT_O" = results_ACAT_O,
                 "OrdinalSTAAR_SKAT" = results_STAAR_S,
                 "OrdinalSTAAR_ACAT" = results_STAAR_A,
                 "OrdinalSTAAR_Burden" = results_STAAR_B,
                 "OrdinalSTAAR_pvalue" = results_pvalue)


  return(results)
}






