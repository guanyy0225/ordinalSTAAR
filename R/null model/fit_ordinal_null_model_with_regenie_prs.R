
##########################################################
# fit ordinal null model with regenie binary multitrait
# Zilin Li, Xihao Li, Yuanyuan Guan
# 19/06/2025
##########################################################

# sinfo
# srun -p maths-gpu -w c90 --mem=300G --cpus-per-task=8 --time=24:00:00 --pty bash
# conda activate ukbb
# R

library(data.table)
library(Matrix)
library(stringr)
library(STAAR)
library(GENESIS)
library(GMMAT)
library(STAARpipeline)
library(SeqVarTools)
library(tidyr)
library(dplyr)
library(ordinal)
# library(ordinalSTAAR)

# step1: generate full data ------------------------------------------------

load("/datapool/maths/UKBB/Data/Phenotype/main_survey/bd_676623.Rdata")

# PC1-PC20
pcs <- fread("/datapool/maths/UKBB/Data/PC_GRM/sGRM_and_PCs/output.pca.score")
colnames(pcs)[1:2] <- c("userId", "userID2")
colnames(pcs)[3:22] <- paste0("PC", 1:20)

# Phenotype
pheno <- bd[, c(c(0, 769) + 1)]
colnames(pheno) <- c("userId", "alcohol_intake_frequency")

# Covariates
covariate <- bd[, c(c(0, 22, 10318) + 1)]
colnames(covariate)[1:3] <- c("userId", "sex", "age")

# merge PCs, phenotype and covariates
t1 <- merge(pheno, covariate, by = "userId")
fullDat <- merge(t1, pcs, by = "userId")

# data filter
gId <- get(load("/datapool/maths/UKBB/Data/Phenotype/sample_id/UKB_500K_WGS_sampleid.RData"))
fullDat <- fullDat[
  (fullDat$userId %in% gId) &
    !is.na(fullDat$alcohol_intake_frequency) &
    !is.na(fullDat$PC1) & !is.na(fullDat$age) & !is.na(fullDat$sex),
]

# data transformation
fullDat$age2 <- (fullDat$age)^2
fullDat <- fullDat[fullDat$alcohol_intake_frequency != "Prefer not to answer", ]
fullDat$alcohol_intake_frequency <- ordered(
  factor(
    fullDat$alcohol_intake_frequency,
    levels = c("Never", "Special occasions only",
               "One to three times a month", "Once or twice a week",
               "Three or four times a week", "Daily or almost daily")
  )
)
fullDat$alcohol_ordinal <- as.integer(fullDat$alcohol_intake_frequency)

ncat <- max(fullDat$alcohol_ordinal, na.rm=TRUE) # 6
for(i in 1:(ncat - 1)) {
  alcohol_bin <- paste0("alcohol_bin_", i)
  fullDat[[alcohol_bin]] <- ifelse(fullDat$alcohol_ordinal <= i, 1, 0)
}

# data standardization
cols_to_scale <- c("age", "age2", paste0("PC", 1:10))
fullDat[, (cols_to_scale) := lapply(.SD, scale), .SDcols = cols_to_scale]

save(fullDat, file="/datapool/home/2024102311/UKB/nullmodels/multiclass/alcohol_intake_frequency/binary_regenie_combination/WGS_binary_alcohol_intake_frequency_fullDat.20250601.Rdata")
# load("/datapool/home/2024102311/UKB/nullmodels/multiclass/alcohol_intake_frequency/binary_regenie_combination/WGS_binary_alcohol_intake_frequency_fullDat.20250601.Rdata")
# load("D:/desktop/multiclass/Nullmodel/alcohol_intake_frequency/LOCO/WGS_alcohol_intake_frequency_fullDat.20250619.Rdata")

# step2: generate regenie input files -------------------------------------

# generate phenotype and covariates
pheno_cols <- c("userId", "userId", paste0("alcohol_bin_", 1:5))
pheno_file <- fullDat[, ..pheno_cols]
setnames(pheno_file, old = 1:2, new = c("FID", "IID"))
fwrite(pheno_file,
       file = "/datapool/home/2024102311/UKB/nullmodels/multiclass/alcohol_intake_frequency/binary_regenie_combination/phenotype_alcohol_intake_frequency.20250601.txt",
       sep = "\t",
       row.names = FALSE)

covar_cols <- c("userId", "userId", "sex", "age", "age2", paste0("PC", 1:10))
covar_file <- fullDat[, ..covar_cols]
setnames(covar_file, old = 1:2, new = c("FID", "IID"))
fwrite(covar_file,
       file = "/datapool/home/2024102311/UKB/nullmodels/multiclass/alcohol_intake_frequency/binary_regenie_combination/covariates.20250517.txt",
       sep = "\t",
       row.names = FALSE)

# merge phenotype and covariates
pheno_covar <- merge(pheno_file, covar_file, by = c("FID", "IID"))
pheno_covar$FID_IID <- paste0(pheno_covar$FID, "_", pheno_covar$IID)

# step3: regenie ----------------------------------------------------------

# conda activate regenie_env
# regenie \
# --step 1 \
# --bed /datapool/maths/UKBB/UKB_genotypes_QCed_ZHL/chrall \
# --phenoFile /datapool/home/2024102311/UKB/nullmodels/multiclass/alcohol_intake_frequency/binary_regenie_combination/phenotype_alcohol_intake_frequency.20250601.txt \
# --covarFile /datapool/home/2024102311/UKB/nullmodels/multiclass/alcohol_intake_frequency/binary_regenie_combination/covariates.20250517.txt \
# --covarColList sex,age,age2,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
# --catCovarList sex \
# --bt \
# --bsize 1000 \
# --out /datapool/home/2024102311/UKB/nullmodels/multiclass/alcohol_intake_frequency/binary_regenie_combination \
# --threads 8 \
# --print-prs

# Step 4: Read and Reshape WIDE format PRS files, then Merge

output_path <- "/datapool/home/2024102311/UKB/nullmodels/multiclass/alcohol_intake_frequency/binary_regenie_combination/"

prs_reshaped_list <- list()
ncat <- max(fullDat$alcohol_ordinal, na.rm=TRUE) # 6
ncat <- 6

# Loop through each binary phenotype to read and reshape its PRS file
for (i in 1:(ncat - 1)) {
  # --- 1. Construct file name ---
  prs_file_name <- paste0(output_path, "binary_regenie_combination_", i, ".prs")

  if (file.exists(prs_file_name)) {

    # --- 2. Read the WIDE format PRS file ---
    prs_wide_data <- fread(prs_file_name, header = TRUE, data.table = FALSE)
    colnames(prs_wide_data)[1] <- "chr"

    prs_long_data <- prs_wide_data %>%
      pivot_longer(
        cols = -chr,
        names_to = "FID_IID",
        values_to = paste0("prs_bin_", i)
      ) %>%
      select(FID_IID, !!paste0("prs_bin_", i))

    prs_reshaped_list[[i]] <- prs_long_data

    message(paste("Successfully read and reshaped PRS for alcohol_bin_", i))

  } else {
    warning(paste("PRS file not found:", prs_file_name))
  }
}

# --- Step 5: Merge all reshaped PRS results ---
fullDat <- fullDat %>%
  select(
    userId, alcohol_intake_frequency, sex, age, age2,
    paste0("PC", 1:10)
  )
fullDat$FID_IID <- paste(fullDat$userId, fullDat$userId, sep = "_")

merged_prs_data <- Reduce(
  function(x, y) full_join(x, y, by = "FID_IID"),
  prs_reshaped_list
)

fullDat <- fullDat %>%
  left_join(merged_prs_data, by = "FID_IID") %>%
  select(-FID_IID)

all_pred_cols <- paste0("prs_bin_", 1:(ncat - 1))
all_pred <- fullDat %>% select(all_of(all_pred_cols))

# Handle NAs
na_count <- sum(is.na(all_pred))
if(na_count > 0) {
  warning(paste(na_count, "NA values found in the merged PRS columns. Removing affected rows."))
  fullDat <- fullDat[complete.cases(all_pred), ]
  all_pred <- all_pred[complete.cases(all_pred), ]
}

all_pred_cols <- paste0("prs_bin_", 1:(ncat - 1))
all_probit_pred_cols <- paste0("prs_probit_bin_", 1:(ncat - 1))
fullDat[, all_probit_pred_cols] <- lapply(fullDat[, all_pred_cols, with = FALSE], function(logit_prs) {
  return(qnorm(plogis(logit_prs)))
})




# Step 6: PCA -------------------------------------------------------------
all_probit_pred <- fullDat %>% select(all_of(all_probit_pred_cols))
if (any(!is.finite(as.matrix(all_probit_pred)))) {
  warning("Infinite values found after probit transformation. Removing affected rows.")
  # is.finite() 会对每一列进行判断, rowSums > 0 表示该行至少有一个非有限值
  finite_rows <- rowSums(is.finite(as.matrix(all_probit_pred))) == ncol(all_probit_pred)
  fullDat <- fullDat[finite_rows, ]
  all_probit_pred <- all_probit_pred[finite_rows, ]
}

# 确保在处理完NA/Inf后，数据仍然完整
na_count <- sum(is.na(all_probit_pred))
if(na_count > 0) {
  warning(paste(na_count, "NA values found in the probit PRS columns. Removing affected rows."))
  fullDat <- fullDat[complete.cases(all_probit_pred), ]
  all_probit_pred <- all_probit_pred[complete.cases(all_probit_pred), ]
}


pca_result <- prcomp(all_probit_pred, center = TRUE, scale. = TRUE)
summary(pca_result)


# Choose the number of principal components to keep
keep_num <- 5
pca_components <- as.data.frame(pca_result$x[, 1:keep_num])
colnames(pca_components) <- paste0("prs_pc", 1:keep_num)
# head(pca_components)

data_for_null_model <- bind_cols(fullDat, pca_components)
# 正确的代码
save(data_for_null_model,
     file = "/datapool/home/2024102311/UKB/nullmodels/multiclass/alcohol_intake_frequency/latent_pheno/WGS_alcohol_intake_frequency_fullDat_with_PRS_PCs.20250910.Rdata")


# --- Step 7: Fit Final Ordinal Null Model ---
load("D:\\desktop\\Proj\\Multiclass\\Nullmodel\\alcohol_intake_frequency\\PRS\\WGS_alcohol_intake_frequency_fullDat_with_PRS_PCs.20250619.Rdata")
load("/mnt/project/UKB_500K_WGS_staarpipeline/Multiclass/alcohol_intake_frequency/WGS_alcohol_intake_frequency_fullDat_with_PRS_PCs.20250910.Rdata")

message("Step 7.1: Preparing inputs for the final null model...")

base_covar_cols <- c("sex", "age", "age2", paste0("PC", 1:10))
prs_pc_cols <- paste0("prs_pc", 1:5)

outcome_column <- "alcohol_intake_frequency"
id_column <- "userId"

message(paste("Outcome:", outcome_column))
message("Base Covariates:", paste(base_covar_cols, collapse=", "))
message("PRS Covariates:", paste(prs_pc_cols, collapse=", "))


message("\nStep 7.2: Fitting the final Ordinal Null Model...")

obj.STAAR.UKB.alcohol_intake_frequency <- NullModel(
  phenofile = data_for_null_model,
  outcomeCol = outcome_column,
  sampleCol = id_column,
  covCol = base_covar_cols,
  PRSCol = prs_pc_cols,
  LOCO = FALSE,
  verbose = TRUE
)

message("\nStep 7.3: Saving the final null model object...")
save(
  obj.STAAR.UKB.alcohol_intake_frequency,
  file = "D:/desktop/multiclass/Nullmodel/alcohol_intake_frequency/PRS/obj.STAAR.UKB.alcohol_intake_frequency.20250909.Rdata"
)

message("--- Analysis Complete ---")







# --- 准备工作 ---
# 假设以下变量已经在您的 R 环境中定义好了：
# - data_for_null_model: 包含所有数据的完整数据框
# - outcome_column: 字符变量, "alcohol_intake_frequency"
# - base_covar_cols: 一个包含协变量名称的向量
# - prs_pc_cols: 一个包含PRS主成分名称的向量

# --- 步骤 1: 拟合模型以识别样本 ---
# 注意：这里使用标准模型（firth=FALSE），因为它足以识别出问题样本
# 正如我们之前讨论过的，Firth模型会得到相同的结果。
formula_null <- as.formula(paste(outcome_column, "~", paste(c(base_covar_cols, prs_pc_cols), collapse = " + ")))

clm_obj <- ordinal::clm(
  formula = formula_null,
  data = data_for_null_model,
  link = "probit",
  model = TRUE # model=TRUE 是必需的
)

# --- 步骤 2: 重新计算方差以找到问题样本的索引 ---
# 这段逻辑直接从您的 NullModel 函数中复制而来
model_data <- clm_obj$model
kept_row_indices <- as.numeric(rownames(model_data))
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
var_y_raw <- 1 + term1 - residuals^2

# --- 步骤 3: 识别并分割数据 ---
is_problematic <- !is.finite(var_y_raw) | var_y_raw < 1e-8
original_indices_of_problematic_samples <- kept_row_indices[is_problematic]

# 创建两个数据子集
problematic_data <- data_for_null_model[original_indices_of_problematic_samples, ]
non_problematic_data <- data_for_null_model[-original_indices_of_problematic_samples, ]

# --- 步骤 4: 生成与图片完全一致的输出 ---

# 输出第一部分 (Problematic Samples)
# summary() 函数作用于一个因子(factor)变量时，会产生这种格式的输出
summary(problematic_data[[outcome_column]])

# 输出第二部分 (Non-Problematic Samples)
summary(non_problematic_data[[outcome_column]])
