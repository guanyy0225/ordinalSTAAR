Ordinal_missense <- function(gene_name, genofile, objNull, genes_info, variant_type = NULL,
                                    rare_maf_cutoff = 0.01, rare_num_cutoff = 2,
                                    geno_missing_cutoff = 0.1,
                                    geno_missing_imputation = c("mean","minor"),
                                    min_maf_cutoff = 0, combine_ultra_rare = TRUE, ultra_rare_mac_cutoff = 20,
                                    QC_label = "annotation/filter", Annotation_dir = "annotation/info/FunctionalAnnotation",
                                    Annotation_name_catalog, Use_annotation_weights = TRUE, Annotation_name = NULL,
                                    use_SPA = NULL, SPA_filter = TRUE, SPA_filter_cutoff = 0.05,
                                    rm_long = TRUE, rm_long_cutoff = 3000,
                                    instability_variance_cutoff = 10000,
                                    verbose = FALSE) {

  # --- [INTERNAL HELPER FUNCTION - FULLY ROBUST AND READABLE] ---
  run_sub_analysis <- function(variant_ids_to_analyze, sub_category_name) {

    # --- Initial variant count checks ---
    if (length(variant_ids_to_analyze) < rare_num_cutoff) {
      message(paste0("Variants number of *", sub_category_name, "* is less than ", rare_num_cutoff, ", skipping..."))
      return(list("OrdinalSTAAR_O" = NA))
    }
    if (rm_long && length(variant_ids_to_analyze) > rm_long_cutoff) {
      message(paste0("Variants number of *", sub_category_name, "* is more than ", rm_long_cutoff, ", skipping..."))
      return(list("OrdinalSTAAR_O" = NA))
    }

    seqSetFilter(genofile, variant.id=variant_ids_to_analyze, sample.id=phenotype.id)

    # --- Genotype and Annotation Extraction ---
    Anno.Int.PHRED.sub <- NULL
    Anno.Int.PHRED.sub.name <- NULL
    if(variant_type != "Indel" && Use_annotation_weights) {
      for(k in 1:length(Annotation_name)) {
        if(Annotation_name[k]%in%Annotation_name_catalog$name) {
          Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
          Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))
          if(Annotation_name[k]=="CADD") { Annotation.PHRED[is.na(Annotation.PHRED)] <- 0 }
          if(Annotation_name[k]=="aPC.LocalDiversity") {
            Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
            Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
          }
          Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
        }
      }
      Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
      colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
    }

    id.genotype <- seqGetData(genofile,"sample.id")
    id.genotype.merge <- data.frame(id.genotype, index=seq_along(id.genotype))
    phenotype.id.merge <- data.frame(phenotype.id)
    phenotype.id.merge <- dplyr::left_join(phenotype.id.merge, id.genotype.merge, by=c("phenotype.id"="id.genotype"))
    id.genotype.match <- phenotype.id.merge$index

    Geno <- seqGetData(genofile, "$dosage")[id.genotype.match, , drop=FALSE]

    # --- Filtering, Imputation, AND NA REMOVAL ---
    getGeno = genoFlipRV(Geno=Geno, geno_missing_imputation=geno_missing_imputation, geno_missing_cutoff=geno_missing_cutoff,
                         min_maf_cutoff=min_maf_cutoff, rare_maf_cutoff=rare_maf_cutoff, rare_num_cutoff=rare_num_cutoff)

    Geno = getGeno$Geno
    MAF = getGeno$G_summary$MAF
    MAC = getGeno$G_summary$MAC

    if(!is.null(Anno.Int.PHRED.sub)) {
      Anno.Int.PHRED.sub = Anno.Int.PHRED.sub[getGeno$include_index, , drop = FALSE]
    }

    if (!is.null(Anno.Int.PHRED.sub)) {
      complete_anno_idx <- complete.cases(Anno.Int.PHRED.sub)
      if (sum(!complete_anno_idx) > 0) {
        message(paste0("INFO: In *", sub_category_name, "*, found and removed ", sum(!complete_anno_idx), " variant(s) with missing annotation scores."))
        Geno <- Geno[, complete_anno_idx, drop = FALSE]
        MAF <- MAF[complete_anno_idx]
        MAC <- MAC[complete_anno_idx]
        Anno.Int.PHRED.sub <- Anno.Int.PHRED.sub[complete_anno_idx, , drop = FALSE]
      }
    }

    if (is.null(dim(Geno)) || ncol(Geno) < rare_num_cutoff) {
      message(paste0("After all filtering, variants number of *", sub_category_name, "* is less than ", rare_num_cutoff, ", skipping..."))
      return(list("OrdinalSTAAR_O" = NA))
    }

    # --- Robustness Check (now on NA-free data) ---
    message(paste0("Performing pre-check for numerically unstable variants in *", sub_category_name, "*..."))
    pre_check_stats <- Ordinal_exactScore(objNull = objNull, G_mat = Geno, use_SPA = FALSE)

    unstable_idx <- which(pre_check_stats$result$Variance > instability_variance_cutoff)

    if (length(unstable_idx) > 0) {
      message(paste0("WARNING: Found and removed ", length(unstable_idx), " unstable variant(s) in *", sub_category_name, "*."))

      stable_idx <- setdiff(1:ncol(Geno), unstable_idx)

      Geno <- Geno[, stable_idx, drop = FALSE]
      MAF <- MAF[stable_idx]
      MAC <- MAC[stable_idx]

      if (!is.null(Anno.Int.PHRED.sub)) {
        Anno.Int.PHRED.sub <- Anno.Int.PHRED.sub[stable_idx, , drop = FALSE]
      }

      if (ncol(Geno) < rare_num_cutoff) {
        message(paste0("After removing unstable variants from *", sub_category_name, "*, remaining number is less than ", rare_num_cutoff, ", skipping..."))
        return(list("OrdinalSTAAR_O" = NA))
      }
    } else {
      message(paste0("No unstable variants found in *", sub_category_name, "*."))
    }

    # --- Final Analysis ---
    final_result = try(OrdinalSTAAR(Geno, MAF, MAC, objNull, annotation_phred = Anno.Int.PHRED.sub,
                                    rare_maf_cutoff, rare_num_cutoff, combine_ultra_rare, ultra_rare_mac_cutoff,
                                    use_SPA, SPA_filter, SPA_filter_cutoff, verbose), silent = FALSE)

    if (inherits(final_result, "try-error")) {
      final_result = list("OrdinalSTAAR_O" = NA)
    }

    return(final_result)
  }
  # --- [END OF INTERNAL HELPER FUNCTION] ---


  # --- Part 1: Initial Variant Selection ---
  phenotype.id = objNull$sample_ids
  if(is.null(use_SPA)) use_SPA = objNull$use_SPA

  seqResetFilter(genofile, verbose=FALSE)

  filter <- seqGetData(genofile, QC_label)

  if(variant_type=="variant")
  {
    SNVlist <- filter == "PASS"
  }

  if(variant_type=="SNV")
  {
    SNVlist <- (filter == "PASS") & isSNV(genofile)
  }

  if(variant_type=="Indel")
  {
    SNVlist <- (filter == "PASS") & (!isSNV(genofile))
  }

  position <- as.numeric(seqGetData(genofile, "position"))
  variant.id <- seqGetData(genofile, "variant.id")
  rm(filter); gc()

  kk <- which(genes_info$hgnc_symbol==gene_name)
  gene_info_kk = genes_info[kk, 1:2]

  sub_start_loc <- genes_info[kk,3]
  sub_end_loc <- genes_info[kk,4]

  is.in <- (SNVlist) & (position>=sub_start_loc) & (position<=sub_end_loc)
  variant.id.gene <- variant.id[is.in]

  seqSetFilter(genofile, variant.id=variant.id.gene, sample.id=phenotype.id)

  # --- Part 2: Define Variant Sets for Both Analyses ---
  GENCODE.EXONIC.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.EXONIC.Category")]))
  MetaSVM_pred <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="MetaSVM")]))
  variant.id.gene.current <- seqGetData(genofile, "variant.id")

  is_missense <- (GENCODE.EXONIC.Category=="nonsynonymous SNV")
  is_dmissense <- is_missense & (MetaSVM_pred=="D")

  variant.id.missense <- variant.id.gene.current[is_missense]
  variant.id.dmissense <- variant.id.gene.current[is_dmissense]

  # --- Part 3: Run Both Analyses Using the Helper Function ---
  message("\n--- Analyzing all missense variants ---")
  result.missense <- run_sub_analysis(variant.id.missense, "missense")

  message("\n--- Analyzing disruptive missense variants ---")
  result.dmissense <- run_sub_analysis(variant.id.dmissense, "disruptive missense")

  # --- Part 4: Combine Results ---
  if (!is.na(result.missense[1]) & !is.na(result.dmissense[1])) {
    if(!is.null(result.missense$results_STAAR_O_P) & !is.null(result.dmissense$results_STAAR_O_P))
    {
      p_names <- colnames(result.missense$results_STAAR_O_P)
      p_names <- p_names[p_names != "ACAT-O"]

      result.missense$results_STAAR_O_P <- result.missense$results_STAAR_O_P[, p_names, drop = FALSE]
      result.missense$results_STAAR_O_P <- cbind(result.missense$results_STAAR_O_P, "Disruptive" = result.dmissense$results_STAAR_O_P[, 1])

      combined_p <- apply(result.missense$results_STAAR_O_P, 1, function(p_row) CCT(p_row))

      result.missense$results_STAAR_O['p.value',1] <- CCT(result.missense$results_STAAR_O_P)
      result.missense$results_STAAR_O_P <- rbind(result.missense$results_STAAR_O_P, "STAAR-O" = combined_p)
    }
  }

  seqResetFilter(genofile, verbose=FALSE)

  result = c(list("gene_info" = gene_info_kk, "category" = "missense"), result.missense)

  return(result)
}
