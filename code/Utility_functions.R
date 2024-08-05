partial_spearman_rho <- function(formula, adj.method = "BH") {
  # This function calculates the partial Spearman correlation coefficients
  # between variables, controlling for confounding variables. The function requires a formula
  # indicating the dependent, independent, and confounding variables (e.g., Y ~ X | Z).
  # The input formula should specify three variables or dataframes:
  #   - The dependent variable (Y)
  #   - The independent variable (X)
  #   - The confounding variable(s) (Z)
  # The data for these variables must be available in the global environment.
  # The function returns a list containing:
  #   - rho: A matrix of partial Spearman correlation coefficients.
  #   - p.value: A matrix of associated p-values.
  #   - adjusted.p.value: A matrix of adjusted p-values (corrected for multiple testing).
  #   - CI: A data frame containing the upper and lower bounds of the confidence intervals for the correlations.
  # The 'adj.method' argument allows specifying the method for p-value adjustment, default is "BH" (Benjamini-Hochberg).

    # Extract the variable names from the formula and get them from the global environment
    data <- all.vars(formula) %>%
        purrr::set_names() %>%
        purrr::map_dfc(~ get(.x, envir = .GlobalEnv))

    # Extract variable names from formula
    formula_parts <- all.vars(formula)
    X_str <- formula_parts[1]
    Y_str <- formula_parts[2]
    Z_str <- formula_parts[3]

    # Extract variables from data
    X_data <- get(X_str, envir = .GlobalEnv)
    Y_data <- get(Y_str, envir = .GlobalEnv)
    Z_data <- get(Z_str, envir = .GlobalEnv)

    # Check if X and Y are vectors or have one column, and convert to dataframes if necessary
    X_data <- if (is.vector(X_data) || (is.matrix(X_data) && ncol(X_data) == 1)) as.data.frame(X_data) else X_data
    Y_data <- if (is.vector(Y_data) || (is.matrix(Y_data) && ncol(Y_data) == 1)) as.data.frame(Y_data) else Y_data

    # Create confounder list as a string
    confou.list <- paste(names(Z_data), collapse = "+")

    # Initialize matrices to store results
    dim_names <- list(colnames(X_data), colnames(Y_data))
    rho_values <- matrix(NA, ncol(X_data), ncol(Y_data), dimnames = dim_names)
    p_values <- matrix(NA, ncol(X_data), ncol(Y_data), dimnames = dim_names)
    CI_up <- matrix(NA, ncol(X_data), ncol(Y_data), dimnames = dim_names)
    CI_low <- matrix(NA, ncol(X_data), ncol(Y_data), dimnames = dim_names)

    # Compute partial Spearman correlation for each pair of variables
    for (i in 1:ncol(X_data)) {
        for (j in 1:ncol(Y_data)) {
            tmp.form <- formula(paste(names(X_data)[[i]], "|", names(Y_data)[[j]], "~", confou.list))
            rho_result <- PResiduals::partial_Spearman(tmp.form, data, fit.x = "lm", fit.y = "lm")
            rho_values[i, j] <- rho_result$TS$TB$ts
            p_values[i, j] <- rho_result$TS$TB$pval
            CI_up[i, j] <- rho_result$TS$TB$upper
            CI_low[i, j] <- rho_result$TS$TB$lower
        }
    }

    # Adjust p-values
    corrected_pvalues <- multtest::mt.rawp2adjp(as.vector(p_values), proc = adj.method)$adjp
    corrected_pvalues_matrix <- matrix(corrected_pvalues[, 2], ncol = ncol(Y_data), byrow = TRUE, dimnames = dim_names)

    # Prepare CI dataframe by merging results
    CI.df <- melt(CI_up, value.name = "CI_up") %>%
        mutate(pairs = paste0(Var1, "_", Var2)) %>%
        select(-c(Var1, Var2)) %>%
        merge(melt(CI_low, value.name = "CI_low") %>%
                  mutate(pairs = paste0(Var1, "_", Var2)) %>%
                  select(-c(Var1, Var2)), by = "pairs") %>%
        merge(melt(rho_values, value.name = "rho") %>%
                  mutate(pairs = paste0(Var1, "_", Var2)) %>%
                  select(-c(Var1, Var2)), by = "pairs") %>%
        merge(melt(p_values, value.name = "p.value") %>%
                  mutate(pairs = paste0(Var1, "_", Var2)) %>%
                  select(-c(Var1, Var2)), by = "pairs") %>%
        merge(melt(corrected_pvalues_matrix, value.name = "p.adj") %>%
                  mutate(pairs = paste0(Var1, "_", Var2)) %>%
                  select(-c(Var1, Var2)), by = "pairs")

    # Return the results
    return(list(rho = as.data.frame(rho_values),
                p.value = as.data.frame(p_values),
                adjusted.p.value = corrected_pvalues_matrix,
                CI = as.data.frame(CI.df)))
}
