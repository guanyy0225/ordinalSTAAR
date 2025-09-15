system("dx-restore-folder /STAARpipeline/R_packages/rstudio_workbench_ukbrap_trial.zilinli_iu.2024-11-24T06-56-17.tar.gz")
library(STAAR)
library(STAARpipeline)
library(SeqArray)
library(dplyr)
library(Matrix)

# output_path <- "/UKB_500K_WGS_staarpipeline/Multiclass/alcohol_intake_frequency/"
# output_rdata_name <- paste0("obj.STAAR.UKB.alcohol_intake_frequency",".RData")
# save(obj.STAAR.UKB.alcohol_intake_frequency, file = output_rdata_name)
# dx_upload_command <- paste0("dx upload ", output_rdata_name, " --path ", output_path)
# system(dx_upload_command)

# --- 1. Load Pre-computed Null Model and Gene Info ---
message("--- Loading Pre-computed Null Model and Gene Info ---")

# --- 2. Define Analysis Parameters ---
message("--- Defining Analysis Parameters ---")
Annotation_dir <- "annotation/info/FunctionalAnnotation"
system("dx download -f Annotation_name_catalog.csv")
Annotation_name_catalog <- read.csv("Annotation_name_catalog.csv")
Annotation_name <- c('CADD','LINSIGHT','FATHMM.XF','aPC.EpigeneticActive',
                     'aPC.EpigeneticRepressed','aPC.EpigeneticTranscription',
                     'aPC.Conservation','aPC.LocalDiversity',
                     'aPC.Mappability','aPC.TF','aPC.Protein')
variant_type <- "variant"
geno_missing_imputation <- "mean"
QC_label <- "annotation/filter"
output_path <- "/UKB_500K_WGS_staarpipeline/Multiclass/alcohol_intake_frequency/"

# --- 3. Define Analysis Tasks ---
# analysis_tasks <- data.frame(
#   chr = c(4, 4, 4, 4, 4, 4, 4),
#   gene_name = c("ADH1C", "ADH1C", "ADH1C","ADH1C","ADH1C","ADH1C","ADH1C"),
#   categories = c("plof_ds", "plof", "disruptive_missense", "ptv","synonymous","ptv_ds","missense")
# )

analysis_tasks <- data.frame(
  chr = c(4),
  gene_name = c("ADH1C"),
  categories = c("plof")
)
unique_chromosomes <- unique(analysis_tasks$chr)

# --- 4. Loop and Execute Association Analysis ---
message("--- Starting Main Analysis Loop ---")

for (current_chr in unique_chromosomes) {

  gds_path <- paste0("ukb.500k.wgs.chr", current_chr, ".pass.annotated.gds")
  message(paste0("\n--- Preparing to analyze chromosome: ", current_chr, " ---"))

  if (!file.exists(gds_path)) {
    message(paste0("GDS file not found locally. Downloading from DNAnexus..."))
    gds_dx_path <- paste0("/UKB_500K_WGS_AGDS_uncompressed_newQClabel_HWE/", gds_path)
    system(paste0("dx download -f ", gds_dx_path))
  }

  if (!file.exists(gds_path)) {
    stop(paste0("FATAL ERROR: Failed to download or find GDS file: ", gds_path))
  }

  genofile <- seqOpen(gds_path)
  message(paste0("Successfully opened aGDS file: ", gds_path))

  tasks_for_this_chr <- analysis_tasks[analysis_tasks$chr == current_chr, ]

  for (i in 1:nrow(tasks_for_this_chr)) {

    gene_name <- tasks_for_this_chr$gene_name[i]
    categories <- tasks_for_this_chr$categories[i]

    message(paste0("\n=============================================="))
    message(paste0("Processing Gene: ", gene_name, " (Chr: ", current_chr, "), categories: ", categories))
    message(paste0("==============================================\n"))

    # --- 关键修正：调用正确的顶层函数和参数 ---
    analysis_results <- tryCatch({
      results <- Ordinal_GeneCentricCoding(
        # --- 核心参数 (你已经提供了) ---
        gene_name = gene_name,
        genofile = genofile,
        objNull = obj.STAAR.UKB.alcohol_intake_frequency,
        genes_info = genes_info,
        categories = categories,
        variant_type = variant_type,

        # --- 筛选和编码相关参数 (补全) ---
        rare_maf_cutoff = 0.01,
        rare_num_cutoff = 2,
        geno_missing_cutoff = 0.1,               # <--- 补上这个参数
        geno_missing_imputation = geno_missing_imputation,
        min_maf_cutoff = 0,                      # <--- 补上这个参数
        combine_ultra_rare = TRUE,               # <--- 补上这个参数
        ultra_rare_mac_cutoff = 20,              # <--- 补上这个参数
        rm_long = TRUE,                          # <--- 补上这个参数
        rm_long_cutoff = 3000,                   # <--- 补上这个参数

        # --- 注释和质控相关参数 (你已经提供了) ---
        QC_label = QC_label,
        Annotation_dir = Annotation_dir,
        Annotation_name_catalog = Annotation_name_catalog,
        Use_annotation_weights = TRUE,
        Annotation_name = Annotation_name,

        # --- SPA 相关参数 (补全) ---
        use_SPA = NULL, # 设为 NULL，函数会从 objNull 中自动获取
        SPA_filter = TRUE,                       # <--- 补上这个参数
        SPA_filter_cutoff = 0.05,                # <--- 补上这个参数

        # --- 其他参数 ---
        verbose = FALSE                          # <--- 补上这个参数，可以设为 TRUE 来查看详细输出
      )
    }, error = function(e) {
      message(paste("An error occurred while analyzing categories '", categories, "': ", e$message))
      return(NULL)
    })

    # --- 结果处理 ---

    if (!is.null(analysis_results) && !is.na(analysis_results$OrdinalSTAAR_O)) {
      output_rdata_name <- paste0(gene_name, "_", categories, ".RData")
      message(paste0("Saving results to '", output_rdata_name, "'..."))
      save(analysis_results, file = output_rdata_name)
      dx_upload_command <- paste0("dx upload ", output_rdata_name, " --path ", output_path)
      message(paste0("Uploading file with command: ", dx_upload_command))
      system(dx_upload_command)
    } else {
      message(paste0("No results for categories '", categories, "'. Skipping file save/upload."))
    }
  }

  # seqClose(genofile)
  # system(paste0("rm ", gds_path))
}

message("\nAll analyses completed!")
