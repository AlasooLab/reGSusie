#!/usr/bin/env Rscript

options(stringsAsFactors=F)

library(data.table)
library(dplyr)
library(stringr)
library(susieR)
library(argparse)
library(R.utils)
library(MungeSumstats)
library(dotgen)


compute_yty <- function(beta, se, p, R, n, k) {
  # beta and se should be standarized
  beta_s <- beta * sqrt(2 * p * (1 - p))
  se_s <- se * sqrt(2 * p * (1 - p))

  # Y'Y =  Bj^2 (Xj'Xj) + Var(Bj)(Xj'Xj)(N - k)
  XjtXj <- (n - 1) * diag(R)
  yty <- beta_s**2 * XjtXj + se_s**2 * XjtXj * (n - k)

  return(median(yty))
}

summarize.susie.cs <- function(object, orig_vars, R, ..., low_purity_threshold = 0.5) {
  if (is.null(object$sets)) {
    stop("Cannot summarize SuSiE object because credible set information is not available")
  }
  variables <- data.frame(cbind(1:length(object$pip), object$pip, -1, NA, NA, NA))
  colnames(variables) <- c("variable", "variable_prob", "cs", "cs_specific_prob", "low_purity", "lead_r2")
  rownames(variables) <- NULL
  added_vars <- c()
  if (object$null_index > 0) variables <- variables[-object$null_index, ]
  if (!is.null(object$sets$cs)) {
    cs <- data.frame(matrix(NA, length(object$sets$cs), 5))
    colnames(cs) <- c("cs", "cs_log10bf", "cs_avg_r2", "cs_min_r2", "variable")
    for (i in 1:length(object$sets$cs)) {
      if (any(object$sets$cs[[i]] %in% added_vars)) {
        print(
          sprintf("Skipping cs %d as there is an overlap between variants in this cs and previous credible sets", i)
        )
        print("Removed cs variants:")
        print(orig_vars[object$sets$cs[[i]], ], max = length(object$sets$cs[[i]]))
        next
      } else {
        added_vars <- append(added_vars, object$sets$cs[[i]])
      }
      in_cs_idx <- which(variables$variable %in% object$sets$cs[[i]])
      variables$cs[in_cs_idx] <- object$sets$cs_index[[i]]
      variables[in_cs_idx, "cs_specific_prob"] <- object$alpha[object$sets$cs_index[[i]], object$sets$cs[[i]]]
      variables$low_purity[in_cs_idx] <- object$sets$purity$min.abs.corr[i] < low_purity_threshold
      lead_pip_idx <- in_cs_idx[which.max(variables$variable_prob[in_cs_idx])]
      variables$lead_r2 <- R[lead_pip_idx, ]^2

      cs$cs[i] <- object$sets$cs_index[[i]]
      cs$cs_log10bf[i] <- log10(exp(object$lbf[cs$cs[i]]))
      cs$cs_avg_r2[i] <- object$sets$purity$mean.abs.corr[i]^2
      cs$cs_min_r2[i] <- object$sets$purity$min.abs.corr[i]^2
      cs$low_purity[i] <- object$sets$purity$min.abs.corr[i] < low_purity_threshold
      cs$variable[i] <- paste(object$sets$cs[[i]], collapse = ",")
    }
    variables <- variables[order(variables$variable_prob, decreasing = T), ]
  } else {
    cs <- NULL
  }
  return(list(vars = variables, cs = na.omit(cs)))
}


susie_ss_wrapper <-function(df, R, n, L, estimate_residual_variance, var_y = 1, prior_weights = NULL, min_abs_corr = 0.0, low_purity_threshold = 0.5) {
  beta <- df$beta
  se <- df$se
  fitted_bhat <- susie_rss(
    bhat = beta,
    shat = se,
    R = R,
    n = n,
    var_y = var_y,
    L = L,
    prior_weights = prior_weights,
    scaled_prior_variance = 0.1,
    estimate_residual_variance = estimate_residual_variance,
    estimate_prior_variance = TRUE,
    standardize = TRUE,
    check_input = FALSE,
    min_abs_corr = min_abs_corr
  )
  cs_summary <- summarize.susie.cs(fitted_bhat, df, R, low_purity_threshold = low_purity_threshold)
  variables <-
    cs_summary$vars %>%
    rename(prob = variable_prob) %>%
    arrange(variable) %>%
    mutate(
      mean = susie_get_posterior_mean(fitted_bhat),
      sd = susie_get_posterior_sd(fitted_bhat)
    )
  cs <- cs_summary$cs

  cs_summary <- summarize.susie.cs(fitted_bhat, df, R, low_purity_threshold = low_purity_threshold)
  sets_95 <- fitted_bhat$sets
  fitted_bhat$sets <- susieR::susie_get_cs(fitted_bhat, coverage = 0.99, Xcorr = R, min_abs_corr = min_abs_corr)
  cs_summary_99 <- summarize.susie.cs(fitted_bhat, df, R, low_purity_threshold = low_purity_threshold)
  fitted_bhat$sets_99 <- fitted_bhat$sets
  fitted_bhat$sets <- sets_95

  variables_99 <-
    cs_summary_99$vars %>%
    rename(prob = variable_prob) %>%
    arrange(variable) %>%
    mutate(
      mean = susie_get_posterior_mean(fitted_bhat),
      sd = susie_get_posterior_sd(fitted_bhat)
    )

  colnames(variables_99) <- paste0(colnames(variables_99),"_99")

  return(list(
    susie_obj = fitted_bhat,
    variables = variables,
    variables_99 = variables_99,
    cs = cs,
    cs_99 = cs_summary_99$cs
  ))
}

main <- function(args) {
  df <- fread(args$z)
  if (startsWith(as.character(df$chromosome[1]), '0')) {
    df$chromosome <- as.numeric(gsub('0', '', as.character(df$chromosome)))
  }
  n <- args$n_samples
  L <- args$L
  prior_weights <- NULL


  message('Step 1.')
  if (!is.null(args$ld)) {
    R <- read.table(gzfile(args$ld), header=F)
    R <- as.data.frame(R)
    rownames(R) <- df$rsid
    colnames(R) <- df$rsid
    R <- as.matrix(R)
  } else {
    warning('No LD file is provided.')
    require(Matrix)
    R <- Matrix::.sparseDiagonal(nrow(df))
  }
  message(paste(length(colnames(R)), 'in analysis.'))

  
  message('Step 2.')
  yty <- compute_yty(df$beta, df$se, df$maf, R, n, args$n_covariates)
  var_y <- yty / (n - 1)
  
  
  message('Step 3.')
  #t <- try(
  res <- susie_ss_wrapper(
         df, R, n, L, estimate_residual_variance = TRUE,  var_y, prior_weights,
         min_abs_corr = args$min_cs_corr, low_purity_threshold = args$low_purity_threshold)
  #     ))
  #if (inherits(t, 'try-error')) {
  #  message('Error in SusieRSS, parameter estimated_residual_variance should be FALSE.')
  #  res <- susie_ss_wrapper(
  #         df, R, n, L, estimate_residual_variance = FALSE,  var_y, prior_weights,
  #         min_abs_corr = args$min_cs_corr, low_purity_threshold = args$low_purity_threshold
  #     )}
  
  print(str(res))

  
  message('Step 4.')
  susie_obj <- res$susie_obj
  
  lbf_variables <- susie_obj$lbf_variable
  transposed_lbf <- t(lbf_variables)
  transposed_lbf <- data.frame(rsid=row.names(transposed_lbf), transposed_lbf)
  transposed_lbf$old_position <- df$position
  transposed_lbf$A1 <- df$allele1
  transposed_lbf$A2 <- df$allele2
  
  old_rsid <- df$rsid
  old_position <- df$position
  df <- df %>% rename('SNP'='rsid', 'CHR'='chromosome', 'BP'='position', 'A1'='allele1', 'A2'='allele2')
  df$rsid <- old_rsid
  df$old_position <- old_position
  df$CHR <- ifelse(df$CHR == 23, as.character('X'), df$CHR)

  if (args$GRCh == 38) {
    df_hg38 <- df
  } else {
    df_hg38 <- MungeSumstats::liftover(df, convert_ref_genome='hg38', ref_genome='hg19')
  }

  message(paste('Number of variants before liftover:', nrow(df)))
  message(paste('Number of variants after liftover:', nrow(df_hg38)))
  
  df_hg38$molecular_trait_id <- args$pheno
  df_hg38$region <- paste0('chr', unique(df_hg38$CHR),':', min(df_hg38$BP),'-', max(df_hg38$BP))
  message(paste('Old region:', args$region))
  message(paste('New region:', unique(df_hg38$region)))
  
  df_hg38$SNP <- paste0('chr', df_hg38$CHR, '_', df_hg38$BP, '_', df_hg38$A1, '_', df_hg38$A2)
  
  final_lbf <- inner_join(df_hg38, transposed_lbf, by=c('rsid', 'old_position', 'A1', 'A2'))
  final_lbf <- final_lbf[order(final_lbf$BP)]
  
  
  data_coloc5 <- final_lbf[, c('molecular_trait_id', 'region', 'SNP', 'CHR', 'BP', 'X1', 'X2',
                               'X3', 'X4', 'X5', 'X6', 'X7', 'X8', 'X9', 'X10')]
  
  colnames(data_coloc5) <- c('molecular_trait_id', 'region', 'variant', 'chromosome', 'position',
                             'lbf_variable1', 'lbf_variable2', 'lbf_variable3', 'lbf_variable4', 
                             'lbf_variable5', 'lbf_variable6', 'lbf_variable7', 'lbf_variable8', 
                             'lbf_variable9', 'lbf_variable10')
  
  write.table(data_coloc5, paste0(args$pheno, '_', unique(df_hg38$region), '_coloc5.tsv'), 
              row.names=F, sep='\t', quote=F)  
  message(paste('Coloc5 data is written.')) 
  
  data_coloc3 <- df_hg38[, c('molecular_trait_id', 'region', 'SNP', 'A1', 'A2', 'CHR', 'BP', 'maf',
                             'beta', 'se', 'an', 'ac', 'n', 'log10P', 'info')]
  
  colnames(data_coloc3) <- c('molecular_trait_id', 'region', 'variant', 'ref', 'alt', 'chromosome',
                             'position', 'maf', 'beta', 'se', 'an', 'ac', 'n', 'log10p', 'info')
  
  write.table(data_coloc3, paste0(args$pheno, '_', unique(df_hg38$region), '_coloc3.tsv'), 
              row.names=F, sep='\t', quote=F)  
  message(paste('Coloc3 data is written.'))

  df_rsid <- df$rsid
  variables <- cbind(df_rsid, res$variables)
  colnames(variables)[1] <- 'rsid'
  variables$old_position <- df$old_position
  variables$A1 <- df$A1
  variables$A2 <- df$A2
  
  pip <- as.data.frame(susie_obj$pip)
  colnames(pip) <- 'pip'
  pip$rsid <- rownames(pip)
  pip$old_position <- df$old_position
  pip$A1 <- df$A1
  pip$A2 <- df$A2
  
  variables <- merge(variables, pip, by=c('rsid', 'old_position', 'A1', 'A2'))
  
  variables_clpp <- variables %>% filter(cs > 0)
  
  cs <- res$cs
  if(is.null(cs)) {
    write.table(variables_clpp, paste0(args$pheno, '_NULL_clpp.tsv'),
                row.names=F, sep='\t')
  } else {
    variables_clpp <- merge(variables_clpp, cs, by='cs')
    variables_clpp_hg38 <- inner_join(df_hg38, variables_clpp, by=c('rsid', 'old_position', 'A1', 'A2'))
    
    message(paste('Number of variants in credible sets before liftover:', nrow(variables_clpp)))
    message(paste('Number of variants in credible sets after liftover:', nrow(variables_clpp_hg38)))
    
    variables_clpp_hg38 <- variables_clpp_hg38[order(variables_clpp_hg38$BP)]
    
    variables_clpp_hg38$z <- zsc(10^-(variables_clpp_hg38$log10P), variables_clpp_hg38$beta)
    variables_clpp_hg38$cs_index <- paste0('L', variables_clpp_hg38$cs)
    variables_clpp_hg38$cs_id <- paste0(variables_clpp_hg38$molecular_trait_id,
                                        '_', variables_clpp_hg38$region,
                                        '_', variables_clpp_hg38$cs_index)
    
    data_clpp <- variables_clpp_hg38[, c('molecular_trait_id', 'region', 'SNP', 'CHR', 'BP', 'A1', 'A2',
                                         'cs_id', 'cs_index', 'pip', 'z', 'cs_min_r2')]
    
    colnames(data_clpp) <- c('molecular_trait_id', 'region', 'variant', 'chromosome', 'position',
                             'ref', 'alt', 'cs_id', 'cs_index', 'pip', 'z', 'cs_min_r2')
    
    write.table(data_clpp, paste0(args$pheno, '_', unique(df_hg38$region), '_clpp.tsv'), 
                row.names=F, sep='\t', quote=F)
    
    message(paste('CLPP data is written.'))
  }

}

parser <- ArgumentParser()

parser$add_argument("--z", "-z", type = "character", required = TRUE)
parser$add_argument("--ld", type = "character")
parser$add_argument("--out", type = "character")
parser$add_argument("--log", type = "character")
parser$add_argument("--pheno", type = "character")
parser$add_argument("--region", type = "character")
parser$add_argument("--n-samples", "-n", type = "integer", required = TRUE)
parser$add_argument("--L", type = "integer", default = 10)
parser$add_argument("--min-cs-corr", default = 0.5, type = "double")
parser$add_argument("--low-purity-threshold", default = 0.5, type = "double")
parser$add_argument("--n-covariates", "-k", type = "integer")
parser$add_argument("--GRCh", type = "integer")


args <- parser$parse_args()

if (is.null(args$out)) {
  args$out <- "tmp"
}

if (is.null(args$log)) {
  args$log <- paste0(args$out, ".susie.log")
}

logfile <- file(args$log, open = "w")
sink(logfile, type = "output", split = TRUE)
sink(logfile, type = "message")
print(args)

print("Analysis started")
tryCatch(
  {
    main(args)
    print("Finished!")
  },
  error = function(e) {
    sink()
    message(as.character(e))
    sink(type = "message")
    stop(e)
  }
)
