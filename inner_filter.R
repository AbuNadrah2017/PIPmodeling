#' Wilcoxon test filter
#' 
#' Simple univariate filter using Wilcoxon (Mann-Whitney) test using the 
#' matrixTests package.
#' 
#' @param y Response vector
#' @param x Matrix of predictors
#' @param force_vars Vector of column names within `x` which are always retained
#'   in the model (i.e. not filtered). Default `NULL` means all predictors will
#'   be passed to `filterFUN`.
#' @param nfilter Number of predictors to return. If `NULL` all predictors with 
#' p values < `p_cutoff` are returned.
#' @param p_cutoff p value cut-off
#' @param rsq_cutoff r^2 cutoff for removing predictors due to collinearity.
#'   Default `NULL` means no collinearity filtering. Predictors are ranked based
#'   on Wilcoxon test. If 2 or more predictors are collinear, the first ranked
#'   predictor by Wilcoxon test is retained, while the other collinear predictors are
#'   removed. See [collinear()].
#' @param type Type of vector returned. Default "index" returns indices,
#' "names" returns predictor names, "full" returns a matrix of p-values.
#' @param exact Logical whether exact or approximate p-value is calculated. 
#' Default is `FALSE` for speed.
#' @param ... Further arguments passed to [matrixTests::row_wilcoxon_twosample]
#' @return Integer vector of indices of filtered parameters (type = "index") or 
#' character vector of names (type = "names") of filtered parameters. If 
#' `type` is `"full"` full output from [matrixTests::row_wilcoxon_twosample] 
#' is returned.
#' @importFrom matrixTests row_wilcoxon_twosample
#' @export
#' 
#' 
#' 

ttest_filter <- function(y,
                         x,
                         force_vars = NULL,
                         nfilter = NULL,
                         p_cutoff = 0.05,
                         rsq_cutoff = NULL,
                         type = c("index", "names", "full")) {
  type <- match.arg(type)
  y <- factor(y)
  indx1 <- as.numeric(y) == 1
  indx2 <- as.numeric(y) == 2
  x <- as.matrix(x)
  if (is.null(colnames(x))) colnames(x) <- seq_len(ncol(x))
  res <- Rfast::ttests(x[indx1, ], x[indx2, ])
  rownames(res) <- colnames(x)
  if (type == "full") return(res)
  filter_end(res[, "pvalue"],
             x, force_vars, nfilter, p_cutoff, rsq_cutoff, type)
}



# wilcoxon_filter <- function(y,
#                             x,
#                             force_vars = NULL,
#                             nfilter = NULL,
#                             p_cutoff = 0.05,
#                             rsq_cutoff = NULL,
#                             type = c("index", "names", "full"),
#                             exact = FALSE,
#                             ...) {
#   
#   type <- match.arg(type)
#   y <- factor(y)
#   indx1 <- as.numeric(y) == 1
#   indx2 <- as.numeric(y) == 2
#   
#   res <- suppressWarnings(
#     
#     matrixTests::row_wilcoxon_twosample(t(x[indx1, ]), t(x[indx2, ]), exact = exact, ...)
#     
#   )
#   
#   
#   if (type == "full") return(res)
#   
#   filter_end(res[, "pvalue"],
#              x, 
#              force_vars, 
#              nfilter, 
#              p_cutoff, 
#              rsq_cutoff, 
#              type)
# }




filter_end <- function(pval, x, 
                       force_vars, 
                       nfilter, 
                       p_cutoff, 
                       rsq_cutoff,
                       type) {
  
  check_vars <- which(!colnames(x) %in% force_vars)
  outp <- pval[check_vars]
  outorder <- order(outp)
  out <- check_vars[outorder]
  outp <- outp[outorder]
  if (!is.null(p_cutoff)) out <- out[outp < p_cutoff]
  if (!is.null(rsq_cutoff)) {
    co <- collinear(x[, out], rsq_cutoff = rsq_cutoff)
    if (length(co) > 0) out <- out[-co]
  }
  if (!is.null(nfilter)) out <- out[1:min(nfilter, length(out))]
  if (length(out) == 0) stop("No predictors left after filtering")
  out <- c(which(colnames(x) %in% force_vars), out)
  switch(type,
         index = out, names = colnames(x)[out])
}


collinear <- function(x, rsq_cutoff = 0.9, verbose = FALSE) {
  rsq <- cor(x, method = 'spearman')^2
  rsq[lower.tri(rsq, diag = TRUE)] <- NA
  combsAboveCutoff <- which(rsq > rsq_cutoff)
  colsToCheck <- ceiling(combsAboveCutoff / nrow(rsq))
  rowsToCheck <- combsAboveCutoff %% nrow(rsq)
  colsToDiscard <- colsToCheck > rowsToCheck
  rowsToDiscard <- !colsToDiscard
  if (any(rowsToDiscard)) warning("Unexpected rows to discard")
  if (verbose) {
    df <- data.frame(keep = rowsToCheck, remove = colsToCheck)
    df <- df[order(df$keep), ]
    remd <- NULL
    for (i in unique(df$keep)) {
      rem <- df$remove[df$keep %in% i]
      rem <- rem[!rem %in% remd]
      if (length(rem) > 0) {
        if (!i %in% df$remove) cat("Keep ") else cat ("Removed ")
        cat(paste0(colnames(x)[i], ", remove "))
        cat(paste(colnames(x)[rem], collapse = ", "))
        remd <- c(remd, rem)
        cat("\n")
      }
    }
  }
  deletecol <- colsToCheck
  deletecol <- unique(deletecol)
  deletecol
}




#' Oversampling and undersampling
#' 
#' Random oversampling of the minority group(s) or undersampling of the majority
#' group to compensate for class imbalance in datasets.
#' 
#' @param y Vector of response outcome as a factor
#' @param x Matrix of predictors
#' @param minor Amount of oversampling of the minority class. If set to `NULL`
#'   then all classes will be oversampled up to the number of samples in the
#'   majority class. To turn off oversampling set `minor = 1`.
#' @param major Amount of undersampling of the majority class
#' @param yminor Optional character value specifying the level in `y` which is
#'   to be oversampled. If `NULL`, this is set automatically to the class with
#'   the smallest sample size.
#' @details
#' `minor` < 1 and `major` > 1 are ignored.
#' @return List containing extended matrix `x` of synthesised data and extended
#'   response vector `y`


randomsample <- function(y, x, minor = NULL, major = 1, yminor = NULL) {
  ytab <- table(y)
  ymajor <- names(ytab)[which.max(ytab)]
  nymajor <- round(max(ytab) * major)
  if (is.null(minor)) {
    # equalise
    yset <- names(ytab)[!names(ytab) %in% ymajor]
    add_samples <- unlist(lapply(yset, function(i) {
      size <- nymajor - ytab[i]
      ind <- which(y == i)
      sample(ind, size, replace = size > length(ind))
    }))
  } else if (minor > 1) {
    # manual
    if (is.null(yminor)) yminor <- names(ytab)[which.min(ytab)]
    size <- round(ytab[yminor] * (minor - 1))
    ind <- which(y == yminor)
    add_samples <- sample(ind, size, replace = (minor >= 2))
  } else add_samples <- NULL
  
  if (major < 1) {
    ind <- which(y == ymajor)
    size <- round(major * length(ind))
    major_samples <- sample(ind, size, replace = size > length(ind))
    rm_samples <- ind[!ind %in% major_samples]
    out <- c(seq_along(y)[-rm_samples], add_samples)
  } else {
    out <- c(seq_along(y), add_samples)
  }
  
  y <- y[out]
  x <- x[out, ]
  rownames(x) <- make.names(rownames(x), unique = TRUE)
  list(y = y, x = x)
}


sample_message <- function(minor = 1, major = 1) {
  if (is.null(minor)) {
    cat("Random oversampling to equalise classes\n")
  } else if (minor > 1) cat("Random oversampling minority x", minor, "\n")
  if (!is.null(major)) {
    if (major < 1) cat("Random undersampling majority x", major, "\n")
  }
}


# 
# 
# 
# df_x = scaled_all_dataDF[, -which(names(scaled_all_dataDF) %in% c('response', clinical_variables_names, TMB_variables_names))]
# 
# 
# 
# args <- list(y = scaled_all_dataDF[,class], 
#              x = df_x,
#              p_cutoff = 0.01, 
#              rsq_cutoff = 0.85, 
#              B = 500,
#              filterFUN = wilcoxon_filter,
#              type = "names")
# 
# fset <- do.call(boot_filter, args)
# 
# 
# filt_xtrain <- scaled_all_dataDF[, c(class, clinical_variables_names, fset, TMB_variables_names)]
# 
# 
# 
# out <- randomsample(y = scaled_all_dataDF[,class], 
#                     x = df_x)
# y2 <- out$y
# x2 <- out$x
# table(y2)
# 
# 
