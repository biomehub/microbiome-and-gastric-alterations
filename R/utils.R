# Define auxiliary functions ---------------------------------------------------

compute_lowest <- function(.phyloseq, .empty  = c("", "0", 0, NA_character_)) {
  #' Add lca column to tax_table of a phyloseq object
  #'
  #' @description Uses `tax_table`` info to construct addtional lca column
  #' corresponding to the assigned taxonomy at its lowest (non-empty) level.
  #'  
  #' @param .phyloseq A phyloseq object.
  #' 
  #' @details Returns a phyloseq object with additional lowest_taxonomy and lowest_rank columns in the `tax_table` slot.
  #' The lowest_taxonomy column is assumed to be the lowest non-empty taxonomy assigned to each oligotype. 
  
  rank_list <- c("species", "genus", "family", "order", 
                 "class", "phylum", "kingdom")
  tax <- cbind(
    .phyloseq@tax_table@.Data,
    lowest_taxonomy = NA_character_,
    lowest_rank = NA_character_
  )
  for (i in 1:nrow(tax)) {
    for (k in rank_list) {
      taxonomy <- tax[i, k]
      if (!(taxonomy %in% .empty)) {
        tax[i, "lowest_taxonomy"] <- taxonomy
        tax[i, "lowest_rank"] <- k
        break()
      }
    }
  }
  
  phyloseq::tax_table(.phyloseq) <- phyloseq::tax_table(tax)
  return(.phyloseq)
}

tax_glom2 <- function(.phyloseq, .rank) {
  
  if (.rank == "oligotype") {
    return(.phyloseq)
  }
  
  if (.rank  %in% c("lowest", "lowest_taxonomy", "lowest_rank")) {
    .rank <- "lowest_taxonomy"
  }
  
  if (.rank == "lowest_taxonomy" & !("lowest_taxonomy" %in% phyloseq::rank_names(.phyloseq))) {
    .phyloseq <- compute_lowest(.phyloseq = .phyloseq)
  }
  
  .phyloseq <- phyloseq::tax_glom(
    physeq = .phyloseq,
    taxrank = .rank, 
    bad_empty = NA
  )
  
  tax <- .phyloseq@tax_table@.Data
  duplicated_taxa <- table(tax[, .rank]) > 1
  duplicated_taxa_names <- names(duplicated_taxa)[duplicated_taxa]
  ix_duplicated_taxa <- which(tax[,.rank] %in% duplicated_taxa_names)
  
  if (.rank == "lowest_taxonomy") {
    new_names <- paste0(tax[ix_duplicated_taxa, .rank],
                        " (",
                        tax[ix_duplicated_taxa, 'lowest_rank'],
                        ")"
    ) %>% 
      make.unique()
  } else {
    new_names <- make.unique(tax[ix_duplicated_taxa, .rank])
  }
  
  tax[ix_duplicated_taxa, .rank] <- new_names
  
  tax_table(.phyloseq) <- tax
  phyloseq::taxa_names(.phyloseq) <- phyloseq::tax_table(.phyloseq)[, .rank]    
  
  zeros <- phyloseq::taxa_names(.phyloseq)[
    stringr::str_starts(phyloseq::taxa_names(.phyloseq), '0')
  ]
  
  .phyloseq <- phyloseq::merge_taxa(.phyloseq, zeros)
  
  unclassified_merged <- stringr::str_starts(phyloseq::taxa_names(.phyloseq), '0')
  phyloseq::taxa_names(.phyloseq)[unclassified_merged] <- stringr::str_glue(
    "Unclassified {.rank}"
  )
  
  return(.phyloseq)
  
}

taxa_prevalence <- function(.phyloseq) {
  .otutable <- phyloseq::otu_table(.phyloseq, taxa_are_rows = phyloseq::taxa_are_rows(.phyloseq))
  apply(.otutable, 2, function(x) mean(x > 0))
}

f <- \(.string, is_path = TRUE) {
  .string <- stringr::str_glue(.string)
  if (isTRUE(is_path)) {
    dir.create(dirname(.string), recursive = TRUE, showWarnings = FALSE)
  }
  .string
}

get_taxon_couns <- function(.phyloseq, taxon, .rank = NULL) {
  if (isFALSE(is.null(.rank))) {
    .phyloseq <- tax_glom2(
      .phyloseq = .phyloseq,
      .rank = .rank
    )
  }
  .metadata <- data.frame(.phyloseq@sam_data) %>% 
    tibble::as_tibble(rownames = "sample_id")
  .countdata <- tibble::as_tibble(
    phyloseq::get_sample(.phyloseq, i = taxon),
    rownames = "sample_id"
  ) %>% 
    dplyr::mutate(
      taxname = taxon
    )
  
  dplyr::left_join(
    .metadata,
    .countdata,
    by = "sample_id"
  )
}

summary_logistf <- function(object, .print = TRUE, ...) {
  
  if (isTRUE(.print)) {
    print(object$call)
    cat("\nModel fitted by", object$method)
    cat("\nCoefficients:\n")
  }
  
  call <- object$call
  if (!is.null(object$modcontrol$terms.fit)) {
    var.red <- object$var[object$modcontrol$terms.fit, object$modcontrol$terms.fit]
    coefs <- coef(object)[object$modcontrol$terms.fit]
    chi2 <- vector(length = length(object$terms))
    chi2[object$modcontrol$terms.fit] <- qchisq(1 - object$prob[object$modcontrol$terms.fit], 
                                                1)
    chi2[-object$modcontrol$terms.fit] <- 0
  }
  else {
    var.red <- object$var
    coefs <- coef(object)
    chi2 <- qchisq(1 - object$prob, 1)
  }
  out <- cbind(object$coefficients, diag(object$var)^0.5, object$ci.lower, 
               object$ci.upper, chi2, object$prob, ifelse(object$method.ci == 
                                                            "Wald", 1, ifelse(object$method.ci == "-", 3, 2)))
  dimnames(out) <- list(names(object$coefficients), c("coef", 
                                                      "se(coef)", paste(c("lower", "upper"), 1 - object$alpha), 
                                                      "Chisq", "p", "method"))
  LL <- -2 * (object$loglik["null"] - object$loglik["full"])
  wald.z <- tryCatch({
    t(coefs) %*% solve(var.red) %*% coefs
  }, error = function(cond) {
    message("\n Variance-Covariance matrix is singular \n")
    return(NA)
  })
  
  if (isTRUE(.print)) {
    print(out)
    cat("\nMethod: 1-Wald, 2-Profile penalized log-likelihood, 3-None\n")
    
    cat("\nLikelihood ratio test=", LL, " on ", object$df, " df, p=", 
        1 - pchisq(LL, object$df), ", n=", object$n, sep = "")
    cat("\nWald test =", wald.z, "on", object$df, "df, p =", 
        1 - pchisq(wald.z, object$df))
  }
  invisible(out)
}
