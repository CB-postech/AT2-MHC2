library(future.apply)

# Function to normalize UMI counts across cells in single-cell RNA sequencing data
#
# This function performs random sampling of UMIs to normalize the total UMI count
# per cell to a target value. For cells with fewer UMIs than the target,
# it performs oversampling with replacement. For cells with more UMIs than
# the target, it performs downsampling without replacement.
#
# @param matrix A genes x cells matrix where each element represents UMI counts
# @param target_umi Integer specifying the desired number of UMIs per cell (default: 10000)
# @return A normalized matrix with exactly target_umi counts per cell
# normalized_matrix <- sampling_umi(raw_matrix, target_umi = 10000)

sampling_umi_optimized <- function(matrix, target_umi = 10000) {
  # Calculate total UMIs per cell
  col_sums <- colSums(matrix)
  
  for(i in colnames(matrix)) {
    # Calculate non-zero genes and their probabilities directly
    # Sample genes

    sampled_genes <- sample(
      which(matrix[, i] > 0),    # only sample expressed genes (exclude zero counts)
      size = target_umi,
      prob = matrix[which(matrix[, i] > 0), i] / col_sums[i], # sample probability determined by relative UMI count
      replace = TRUE
    )

    sampled_counts <- table(sampled_genes)
    matrix[as.numeric(names(sampled_counts)), i] <- as.vector(sampled_counts)
  }
  return(matrix)
}

# Function to perform chi-square test on a pseudo-bulk matrix
#
# This function takes a pseudo-bulk expression matrix and computes expected counts 
# based on given cell numbers. It then performs a chi-square test to determine 
# whether the observed expression significantly deviates from the expected distribution.
#
# @param observed_matrix A genes x cell-types matrix where each element represents UMI counts
# @param cell_numbers A vector specifying the number of cells in each cell type
# @return A list containing chi-square test results for each gene, including
#         statistic, p-value, observed counts, and expected counts.
# results <- perform_chisq_test(observed_matrix, cell_numbers)
perform_chisq_test <- function(observed_matrix, cell_numbers) {
  results <- list()
  
  for(gene in colnames(observed_matrix)) {
    observed <- observed_matrix[, gene]
    names(observed) <- rownames(observed_matrix)
    if (sum(observed) != 0) {
      expected <- sum(observed) * (cell_numbers / sum(cell_numbers))
      chisq_test <- chisq.test(x = observed, p = cell_numbers/sum(cell_numbers))
      
      n <- sum(observed)
      k <- length(observed)
      cramers_v <- sqrt(chisq_test$statistic / (n * (k - 1)))

      results[[gene]] <- list(
        statistic = chisq_test$statistic,
        p_value = chisq_test$p.value,
        observed = observed,
        expected = expected,
        cramers_v = cramers_v
      )
    }
    else {
      results[[gene]] <- list(
        statistic = NA,
        p_value = 1,
        cramers_v = NA,
        observed = observed,
        expected = observed
      )
    }

  }
  return(results)
}

# Function to create a summary matrix from chi-square test results
create_summary_matrix <- function(adjusted_results) {
  # Get all genes (keys from the list)
  genes <- names(adjusted_results)
  
  # Get condition names from the first gene's expected values
  conditions <- names(adjusted_results[[1]]$expected)
  
  # Create column names
  col_names <- c(
    "statistic",
    "p_value",
    "cramers_v",
    paste0(conditions, "_observed"),
    paste0(conditions, "_expected")
  )
  
  # Initialize matrix
  result_matrix <- matrix(
    nrow = length(genes),
    ncol = length(col_names),
    dimnames = list(genes, col_names)
  )
  
  # Fill matrix
  for (i in seq_along(genes)) {
    gene <- genes[i]
    result_matrix[i, "p_value"] <- adjusted_results[[gene]]$p_value
    result_matrix[i, "statistic"] <- adjusted_results[[gene]]$statistic
    result_matrix[i, "cramers_v"] <- adjusted_results[[gene]]$cramers_v
    result_matrix[i, paste0(conditions, "_observed")] <- adjusted_results[[gene]]$observed
    result_matrix[i, paste0(conditions, "_expected")] <- unname(adjusted_results[[gene]]$expected)
  }
  
  return(result_matrix)
}

basal_level_correction <- function(pseudobulk, target_umi_basal = 500) {
  for (i in colnames(pseudobulk)) {
    col = pseudobulk[, i]
    if (sum(col) >= target_umi_basal) {
      probs = col / sum(col)
      sampled_counts = table(sample(1:length(probs), size = target_umi_basal, prob = probs, replace = TRUE))

      # error handling
      # consider there are no counts in the cell in certain condition
      if (col[1] == 0) {
        pseudobulk[1, i] = 0
        pseudobulk[2, i] = target_umi_basal
      }
      else if (col[2] == 0) {
        pseudobulk[1, i] = target_umi_basal
        pseudobulk[2, i] = 0}
      else {
        pseudobulk[, i] = sampled_counts
      }
  }
  }
  return(pseudobulk)
}

# This function is utilty function to run chi-square test with UMI sampling
run_ChiSqaure_with_basal_exp_correction_mean <- function(so, condition, n_iter = 100, target_umi_depth = 10000, target_umi_basal = 500, num_cores = 1, use_batch = FALSE, batch_size = 10) {
    condition_number = so@meta.data[, condition] %>% table
    # > condition_number
    # floxed   dAT2 
    # 3578   4594 

    # 1. UMI sampling (depth correction)
    start_time <- Sys.time()
    all_results <- pbmcapply::pbmclapply(
        1:n_iter,
        function(i) sampling_umi_optimized(so@assays$RNA$counts, target_umi_depth),
        mc.cores = num_cores
    )
    
    mean_result <- Reduce('+', all_results) / n_iter
    rm(all_results)

    end_time = Sys.time()
    print('Sequencing depth correction time spent : ')
    print(end_time - start_time)

    # 2. building pseudobulk
    start_time <- Sys.time()
    pseudobulk_mean <- aggregate(t(mean_result), 
                            list(condition = so@meta.data[, condition]), 
                            sum)
    rownames(pseudobulk_mean) = pseudobulk_mean[[1]] 
    pseudobulk_mean = pseudobulk_mean[-1]
    end_time <- Sys.time()
    print('Pseudobulk building time spent : ')
    print(end_time - start_time)
    pseudobulk_mean <- round(pseudobulk_mean, 0)

    # 3. basal expression level correction
    start_time <- Sys.time()
    sampled_data_list <- pbmcapply::pbmclapply(
      1:n_iter,
      function(i) basal_level_correction(pseudobulk_mean, target_umi_basal),
      mc.cores = num_cores)

    arr <- array(unlist(sampled_data_list[1:5]), 
                dim = c(nrow(sampled_data_list[[1]]), ncol(sampled_data_list[[1]]), length(sampled_data_list)))
    median_matrix <- apply(arr, c(1, 2), median)
    rownames(median_matrix) <- rownames(pseudobulk_mean)
    colnames(median_matrix) <- colnames(pseudobulk_mean)

    end_time <- Sys.time()
    print('Basal expression depth correction time spent : ')
    print(end_time - start_time)

    # 4. chi-square test
    start_time <- Sys.time()

    chi_result <- perform_chisq_test(median_matrix, condition_number)
    chi_summary <- create_summary_matrix(chi_result)
    chi_summary <- data.frame(chi_summary)
    chi_summary$p_BH <- p.adjust(chi_summary$p_value, method = 'BH')
    chi_summary$adjust_p_values <- p.adjust(chi_summary$p_value, method = 'bonferroni')

    # 2.3. calculate differential standardized residuals (effect size)
    # only do if there is only two category in the condition
    if (so@meta.data[, condition] %>% levels %>% length == 2) { 
        
        control_name = levels(so@meta.data[, condition])[1]
        experiment_name = levels(so@meta.data[, condition])[2]

        chi_summary[!is.na(chi_summary$statistic), 'log2FC'] = unlist(log2(((pseudobulk_mean[experiment_name, ]) / (pseudobulk_mean[control_name, ]))[!is.na(chi_summary$statistic)]))
    }
    end_time <- Sys.time()
    print('Chi-square test time spent : ')
    print(end_time - start_time)
    return(list(mean_result = mean_result, chi_summary = chi_summary))
}
