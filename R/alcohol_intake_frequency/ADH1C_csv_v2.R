# --- [CORRECTED AND IMPROVED LAPPLY VERSION] ---

# 1. Define a robust function for a SINGLE analysis task.
#    This function explicitly defines all parameters it needs.
run_single_analysis <- function(task_row, genofile, objNull, genes_info, 
                                variant_type, rare_maf_cutoff, rare_num_cutoff,
                                geno_missing_cutoff, geno_missing_imputation,
                                Annotation_dir, Annotation_name_catalog, 
                                Use_annotation_weights, Annotation_name, QC_label,
                                verbose) {
  
  gene_name <- task_row$gene_name
  categories <- task_row$categories
  
  message(paste0("\nProcessing Gene: ", gene_name, ", categories: ", categories))
  
  # Call the main coding function, passing all arguments received.
  analysis_results <- tryCatch({
    Ordinal_GeneCentricCoding(
      gene_name = gene_name,
      categories = categories,
      genofile = genofile,
      objNull = objNull,
      genes_info = genes_info,
      variant_type = variant_type,
      rare_maf_cutoff = rare_maf_cutoff,
      rare_num_cutoff = rare_num_cutoff,
      geno_missing_cutoff = geno_missing_cutoff,
      geno_missing_imputation = geno_missing_imputation,
      Annotation_dir = Annotation_dir,
      Annotation_name_catalog = Annotation_name_catalog,
      Use_annotation_weights = Use_annotation_weights,
      Annotation_name = Annotation_name,
      QC_label = QC_label,
      verbose = verbose
      # We removed hardcoded parameters and now pass them from the function's arguments.
    )
  }, error = function(e) {
    message(paste("An error occurred: ", e$message))
    return(NULL)
  })
  
  # Format and return the result for this single task
  return(format_results_simple(analysis_results))
}


# 2. Convert the tasks data.frame to a list of rows (this part is correct)
tasks_list <- split(tasks_for_this_chr, seq(nrow(tasks_for_this_chr)))

# 3. Use lapply to iterate over the tasks.
#    Now, we pass ALL the fixed parameters that run_single_analysis needs.
all_formatted_results_list <- lapply(tasks_list, run_single_analysis,
                                     # Core objects that change less often
                                     genofile = genofile,
                                     objNull = obj.STAAR.UKB.alcohol_intake_frequency,
                                     genes_info = genes_info,
                                     
                                     # Analysis parameters defined at the top of your script
                                     variant_type = variant_type,
                                     rare_maf_cutoff = 0.01,
                                     rare_num_cutoff = 2,
                                     geno_missing_cutoff = 0.1,
                                     geno_missing_imputation = "mean",
                                     Annotation_dir = Annotation_dir,
                                     Annotation_name_catalog = Annotation_name_catalog,
                                     Use_annotation_weights = FALSE, # Set to FALSE as intended
                                     Annotation_name = Annotation_name,
                                     QC_label = QC_label,
                                     verbose = FALSE
)

# 4. Combine results and save (same as before)
if (length(all_formatted_results_list) > 0) {
  final_summary_table <- dplyr::bind_rows(all_formatted_results_list)
  output_csv_name <- "OrdinalSTAAR_summary_results_lapply.csv"
  write.csv(final_summary_table, file = output_csv_name, row.names = FALSE, quote = FALSE, na = "NA")
  message(paste0("Successfully saved results to '", output_csv_name, "'"))
} else {
  message("No results were generated to save.")
}