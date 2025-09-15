# --- [CORRECTED FORMATTING FUNCTION FOR ANNOTATED RESULTS] ---
# This version correctly extracts p-values by their column names,
# making it robust to whether annotations are used or not.

format_results <- function(results_list) {

  # --- 1. Handle cases with no valid results ---
  if (is.null(results_list) || is.na(results_list$OrdinalSTAAR_O)) {
    return(NULL) # Return NULL for failed runs, dplyr::bind_rows will ignore it.
  }

  # --- 2. Extract components from the results list ---
  gene_info <- results_list$gene_info
  p_matrix <- results_list$OrdinalSTAAR_pvalue

  # --- 3. Create the single-row data.frame ---
  output_df <- data.frame(
    # Basic Info
    `Gene name` = gene_info$hgnc_symbol,
    Chr = gene_info$chr,
    Category = results_list$category,
    `#SNV` = results_list$num_variant,
    cMAC = results_list$cMAC,

    # --- [THE FIX IS HERE] ---
    # Extract p-values by their proper COLUMN NAMES, not by index (1, 2).
    # "Beta" column corresponds to the p-value without annotation weights.

    `SKAT(1,25)` = p_matrix["SKAT-(1,25)", "Beta"],
    `SKAT(1,1)` = p_matrix["SKAT-(1,1)", "Beta"],

    `Burden(1,25)` = p_matrix["Burden-(1,25)", "Beta"],
    `Burden(1,1)` = p_matrix["Burden-(1,1)", "Beta"],

    `ACAT-V(1,25)` = p_matrix["ACAT-(1,25)", "Beta"],
    `ACAT-V(1,1)` = p_matrix["ACAT-(1,1)", "Beta"],

    # Extract combined p-values (across annotations) for each test type.
    # These come from the "Combined_P" column.
    `STAAR-S(1,25)` = p_matrix["SKAT-(1,25)", "Combined_P"],
    `STAAR-S(1,1)`  = p_matrix["SKAT-(1,1)", "Combined_P"],
    `STAAR-B(1,25)` = p_matrix["Burden-(1,25)", "Combined_P"],
    `STAAR-B(1,1)`  = p_matrix["Burden-(1,1)", "Combined_P"],
    `STAAR-A(1,25)` = p_matrix["ACAT-(1,25)", "Combined_P"],
    `STAAR-A(1,1)`  = p_matrix["ACAT-(1,1)", "Combined_P"],

    # Overall omnibus p-values from the main list
    `ACAT-O` = results_list$ACAT_O,
    `STAAR-O` = results_list$OrdinalSTAAR_O,

    check.names = FALSE
  )

  return(output_df)
}
