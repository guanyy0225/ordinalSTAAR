OrdinalSTAAR = function(Geno, MAF = NULL, MAC = NULL, objNull, annotation_phred = NULL,
                     rare_maf_cutoff = 0.01, rare_num_cutoff = 2,
                     combine_ultra_rare = TRUE, ultra_rare_mac_cutoff = 20,
                     use_SPA = NULL, SPA_filter = TRUE, SPA_filter_cutoff = 0.05, verbose = FALSE) {

  if (!inherits(Geno, "matrix") && !inherits(Geno, "Matrix")) {
    stop("Genotype is not a matrix!")
  }

  if(inherits(Geno, "sparseMatrix")){
    Geno = as.matrix(Geno)
  }

  if (ncol(Geno) < rare_num_cutoff) {
    stop(paste0("Number of variants in the set is less than ", rare_num_cutoff, ", will skip this category..."), call. = FALSE)
    return(list("results_STAAR_O" = NA))
  }

  annotation_phred <- as.data.frame(annotation_phred)
  if (nrow(annotation_phred) != 0 & ncol(Geno) != nrow(annotation_phred)) {
    stop("Dimensions don't match for genotype and annotation!")
  }

  if (inherits(Geno, "sparseMatrix")) {
    Geno <- as.matrix(Geno)
  }

  if (is.null(MAF) | is.null(MAC)) {

    genotype = genoFlip(Geno = Geno)
    MAF = genotype$G_summary$MAF
    MAC = MAF * nrow(Geno) * 2

    RV_label = as.vector((MAF < rare_maf_cutoff)&(MAF > 0))
    Geno = genotype$Geno[ ,RV_label]

    rm(genotype)

  } else {

    RV_label = as.vector((MAF < rare_maf_cutoff)&(MAF > 0))
    Geno = Geno[ ,RV_label]

  }

  MAF = MAF[RV_label]
  MAC = MAC[RV_label]
  Geno = as(Geno, "CsparseMatrix")
  annotation_phred <- annotation_phred[RV_label,,drop=FALSE]


  if(sum(RV_label) < rare_num_cutoff) {

    stop(paste0("Number of rare variant in the set is less than ", rare_num_cutoff, " will skip this category..."), call. = FALSE)
    return(list("results_STAAR_O" = NA))

  } else {

    if (combine_ultra_rare) {
      ultra_rare_length = length(which(MAC < ultra_rare_mac_cutoff))
      if (verbose) print(paste0("Apply STAAR-O test to ", sum(RV_label), " rare variants, with ", ultra_rare_length, " ultra rare variants"))
    } else {
      if (verbose) print(paste0("Apply STAAR-O test to ", sum(RV_label), " rare variants"))
    }

    annotation_rank <- 1 - 10^(-annotation_phred/10)


    w_1 <- dbeta(MAF, 1, 25)
    w_2 <- dbeta(MAF, 1, 1)

    w_a_1 <- w_1^2/dbeta(MAF,0.5,0.5)^2
    w_a_2 <- w_2^2/dbeta(MAF,0.5,0.5)^2

    if(dim(annotation_phred)[2] == 0){

      ## Burden, SKAT, ACAT-V
      w_B <- w_S <- as.matrix(cbind(w_1, w_2))
      w_A <- as.matrix(cbind(w_a_1, w_a_2))

    }else{

      ## Burden
      w_B = as.matrix(cbind(w_1, annotation_rank*w_1, w_2, annotation_rank*w_2))
      ## SKAT
      w_S = as.matrix(cbind(w_1, sqrt(annotation_rank)*w_1, w_2, sqrt(annotation_rank)*w_2))
      ## ACAT-V
      w_A = as.matrix(cbind(w_a_1, annotation_rank*w_a_1, w_a_2, annotation_rank*w_a_2))

    }


    pvalues <- OrdinalSTAAR_O(Geno = Geno, objNull = objNull, annotation_rank, MAC = MAC,
                           use_SPA = use_SPA, SPA_filter = SPA_filter, SPA_filter_cutoff = SPA_filter_cutoff,
                           weight_A = w_A, weight_B = w_B, weight_S = w_S,
                           combine_ultra_rare, ultra_rare_mac_cutoff, verbose = verbose)
    cMAC <- sum(Geno)

    return(c(pvalues,
             list(num_variant = sum(RV_label),
                  cMAC = cMAC,
                  MAF = MAF)))

  }
}
