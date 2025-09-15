Ordinal_GeneCentricCoding <- function(gene_name, genofile, objNull, genes_info, variant_type = NULL,
                                      categories = c("all", "all_ptv", "plof", "plof_ds", "ptv", "ptv_ds", "synonymous", "missense", "disruptive_missense"),
                                      rare_maf_cutoff = 0.01, rare_num_cutoff = 2,
                                      geno_missing_cutoff = 0.1, geno_missing_imputation = c("mean","minor"),
                                      min_maf_cutoff = 0, combine_ultra_rare = TRUE, ultra_rare_mac_cutoff = 20,
                                      QC_label = "annotation/filter", Annotation_dir = "annotation/info/FunctionalAnnotation",
                                      Annotation_name_catalog, Use_annotation_weights = TRUE, Annotation_name = NULL,
                                      use_SPA = FALSE, SPA_filter = TRUE, SPA_filter_cutoff = 0.05,
                                      rm_long = TRUE, rm_long_cutoff = 3000,
                                      instability_variance_cutoff = 10000,
                                      verbose = FALSE) {
  results_list <- list()
  if (categories == "all") {

    results <- Ordinal_coding(gene_name = gene_name, genofile = genofile, objNull = objNull, genes_info = genes_info, variant_type = variant_type,
                              rare_maf_cutoff = rare_maf_cutoff, rare_num_cutoff = rare_num_cutoff, geno_missing_cutoff = geno_missing_cutoff, geno_missing_imputation = geno_missing_imputation,
                              min_maf_cutoff = min_maf_cutoff, combine_ultra_rare = combine_ultra_rare, ultra_rare_mac_cutoff = ultra_rare_mac_cutoff,
                              QC_label = QC_label, Annotation_dir = Annotation_dir,
                              Annotation_name_catalog = Annotation_name_catalog,
                              Use_annotation_weights = Use_annotation_weights,
                              Annotation_name = Annotation_name,
                              use_SPA = use_SPA, SPA_filter = SPA_filter, SPA_filter_cutoff = SPA_filter_cutoff,
                              rm_long = rm_long, rm_long_cutoff = rm_long_cutoff,
                              instability_variance_cutoff = instability_variance_cutoff,
                              verbose = verbose)

  } else if (categories == "all_ptv") {

    results <- Ordinal_coding_incl_ptv(gene_name = gene_name, genofile = genofile, objNull = objNull, genes_info = genes_info, variant_type = variant_type,
                                       rare_maf_cutoff = rare_maf_cutoff, rare_num_cutoff = rare_num_cutoff, geno_missing_cutoff = geno_missing_cutoff, geno_missing_imputation = geno_missing_imputation,
                                       min_maf_cutoff = min_maf_cutoff, combine_ultra_rare = combine_ultra_rare, ultra_rare_mac_cutoff = ultra_rare_mac_cutoff,
                                       QC_label = QC_label, Annotation_dir = Annotation_dir,
                                       Annotation_name_catalog = Annotation_name_catalog,
                                       Use_annotation_weights = Use_annotation_weights,
                                       Annotation_name = Annotation_name,
                                       use_SPA = use_SPA, SPA_filter = SPA_filter, SPA_filter_cutoff = SPA_filter_cutoff,
                                       rm_long = rm_long, rm_long_cutoff = rm_long_cutoff,
                                       instability_variance_cutoff = instability_variance_cutoff,
                                       verbose = verbose)

  } else if (categories == "plof") {

    results <- Ordinal_plof(gene_name = gene_name, genofile = genofile, objNull = objNull, genes_info = genes_info, variant_type = variant_type,
                            rare_maf_cutoff = rare_maf_cutoff, rare_num_cutoff = rare_num_cutoff, geno_missing_cutoff = geno_missing_cutoff, geno_missing_imputation = geno_missing_imputation,
                            min_maf_cutoff = min_maf_cutoff, combine_ultra_rare = combine_ultra_rare, ultra_rare_mac_cutoff = ultra_rare_mac_cutoff,
                            QC_label = QC_label, Annotation_dir = Annotation_dir,
                            Annotation_name_catalog = Annotation_name_catalog,
                            Use_annotation_weights = Use_annotation_weights,
                            Annotation_name = Annotation_name,
                            use_SPA = use_SPA, SPA_filter = SPA_filter, SPA_filter_cutoff = SPA_filter_cutoff,
                            rm_long = rm_long, rm_long_cutoff = rm_long_cutoff,
                            instability_variance_cutoff = instability_variance_cutoff,
                            verbose = verbose)

  } else if (categories == "plof_ds") {

    results <- Ordinal_plof_ds(gene_name = gene_name, genofile = genofile, objNull = objNull, genes_info = genes_info, variant_type = variant_type,
                               rare_maf_cutoff = rare_maf_cutoff, rare_num_cutoff = rare_num_cutoff, geno_missing_cutoff = geno_missing_cutoff, geno_missing_imputation = geno_missing_imputation,
                               min_maf_cutoff = min_maf_cutoff, combine_ultra_rare = combine_ultra_rare, ultra_rare_mac_cutoff = ultra_rare_mac_cutoff,
                               QC_label = QC_label, Annotation_dir = Annotation_dir,
                               Annotation_name_catalog = Annotation_name_catalog,
                               Use_annotation_weights = Use_annotation_weights,
                               Annotation_name = Annotation_name,
                               use_SPA = use_SPA, SPA_filter = SPA_filter, SPA_filter_cutoff = SPA_filter_cutoff,
                               rm_long = rm_long, rm_long_cutoff = rm_long_cutoff,
                               instability_variance_cutoff = instability_variance_cutoff,
                               verbose = verbose)

  } else if (categories == "ptv") {

    results <- Ordinal_ptv(gene_name = gene_name, genofile = genofile, objNull = objNull, genes_info = genes_info, variant_type = variant_type,
                           rare_maf_cutoff = rare_maf_cutoff, rare_num_cutoff = rare_num_cutoff, geno_missing_cutoff = geno_missing_cutoff, geno_missing_imputation = geno_missing_imputation,
                           min_maf_cutoff = min_maf_cutoff, combine_ultra_rare = combine_ultra_rare, ultra_rare_mac_cutoff = ultra_rare_mac_cutoff,
                           QC_label = QC_label, Annotation_dir = Annotation_dir,
                           Annotation_name_catalog = Annotation_name_catalog,
                           Use_annotation_weights = Use_annotation_weights,
                           Annotation_name = Annotation_name,
                           use_SPA = use_SPA, SPA_filter = SPA_filter, SPA_filter_cutoff = SPA_filter_cutoff,
                           rm_long = rm_long, rm_long_cutoff = rm_long_cutoff,
                           instability_variance_cutoff = instability_variance_cutoff,
                           verbose = verbose)

  } else if (categories == "ptv_ds") {

    results <- Ordinal_ptv_ds(gene_name = gene_name, genofile = genofile, objNull = objNull, genes_info = genes_info, variant_type = variant_type,
                              rare_maf_cutoff = rare_maf_cutoff, rare_num_cutoff = rare_num_cutoff, geno_missing_cutoff = geno_missing_cutoff, geno_missing_imputation = geno_missing_imputation,
                              min_maf_cutoff = min_maf_cutoff, combine_ultra_rare = combine_ultra_rare, ultra_rare_mac_cutoff = ultra_rare_mac_cutoff,
                              QC_label = QC_label, Annotation_dir = Annotation_dir,
                              Annotation_name_catalog = Annotation_name_catalog,
                              Use_annotation_weights = Use_annotation_weights,
                              Annotation_name = Annotation_name,
                              use_SPA = use_SPA, SPA_filter = SPA_filter, SPA_filter_cutoff = SPA_filter_cutoff,
                              rm_long = rm_long, rm_long_cutoff = rm_long_cutoff,
                              instability_variance_cutoff = instability_variance_cutoff,
                              verbose = verbose)

  } else if (categories == "synonymous") {

    results <- Ordinal_synonymous(gene_name = gene_name, genofile = genofile, objNull = objNull, genes_info = genes_info, variant_type = variant_type,
                                  rare_maf_cutoff = rare_maf_cutoff, rare_num_cutoff = rare_num_cutoff, geno_missing_cutoff = geno_missing_cutoff, geno_missing_imputation = geno_missing_imputation,
                                  min_maf_cutoff = min_maf_cutoff, combine_ultra_rare = combine_ultra_rare, ultra_rare_mac_cutoff = ultra_rare_mac_cutoff,
                                  QC_label = QC_label, Annotation_dir = Annotation_dir,
                                  Annotation_name_catalog = Annotation_name_catalog,
                                  Use_annotation_weights = Use_annotation_weights,
                                  Annotation_name = Annotation_name,
                                  use_SPA = use_SPA, SPA_filter = SPA_filter, SPA_filter_cutoff = SPA_filter_cutoff,
                                  rm_long = rm_long, rm_long_cutoff = rm_long_cutoff,
                                  instability_variance_cutoff = instability_variance_cutoff,
                                  verbose = verbose)

  } else if (categories == "missense") {

    results <- Ordinal_missense(gene_name = gene_name, genofile = genofile, objNull = objNull, genes_info = genes_info, variant_type = variant_type,
                                rare_maf_cutoff = rare_maf_cutoff, rare_num_cutoff = rare_num_cutoff, geno_missing_cutoff = geno_missing_cutoff, geno_missing_imputation = geno_missing_imputation,
                                min_maf_cutoff = min_maf_cutoff, combine_ultra_rare = combine_ultra_rare, ultra_rare_mac_cutoff = ultra_rare_mac_cutoff,
                                QC_label = QC_label, Annotation_dir = Annotation_dir,
                                Annotation_name_catalog = Annotation_name_catalog,
                                Use_annotation_weights = Use_annotation_weights,
                                Annotation_name = Annotation_name,
                                use_SPA = use_SPA, SPA_filter = SPA_filter, SPA_filter_cutoff = SPA_filter_cutoff,
                                rm_long = rm_long, rm_long_cutoff = rm_long_cutoff,
                                instability_variance_cutoff = instability_variance_cutoff,
                                verbose = verbose)

  } else if (categories == "disruptive_missense") {

    results <- Ordinal_disruptive_missense(gene_name = gene_name, genofile = genofile, objNull = objNull, genes_info = genes_info, variant_type = variant_type,
                                           rare_maf_cutoff = rare_maf_cutoff, rare_num_cutoff = rare_num_cutoff, geno_missing_cutoff = geno_missing_cutoff, geno_missing_imputation = geno_missing_imputation,
                                           min_maf_cutoff = min_maf_cutoff, combine_ultra_rare = combine_ultra_rare, ultra_rare_mac_cutoff = ultra_rare_mac_cutoff,
                                           QC_label = QC_label, Annotation_dir = Annotation_dir,
                                           Annotation_name_catalog = Annotation_name_catalog,
                                           Use_annotation_weights = Use_annotation_weights,
                                           Annotation_name = Annotation_name,
                                           use_SPA = use_SPA, SPA_filter = SPA_filter, SPA_filter_cutoff = SPA_filter_cutoff,
                                           rm_long = rm_long, rm_long_cutoff = rm_long_cutoff,
                                           instability_variance_cutoff = instability_variance_cutoff,
                                           verbose = verbose)

  } else {
    stop(paste("Unknown or unsupported category:", categories), call. = FALSE)
  }

  SeqArray::seqResetFilter(genofile, verbose = FALSE)
  gc()

  return(results)
}
