#' getFittedResultByOneFactor 
#' @description
#' This function fits data to the following linear model: log10(y+1) ~ x
#' And return the slope data and p value for the slope
#' As it is a multiple test, the p value will be adjusted by BH method
#' @param count a character specifying the input count matrix, 
#' the first column should be y and the second column should be x
#' @param sample_sheet a character specifying the input sample sheet,
#' column names of sample_sheet must contain sample_col and conditon_col
#' @param sample_col a character specifying the column name in sample_sheet representing the sample name
#' @param condition_col a character specifying the column name in sample_sheet representing the testing conditions
#' @return return a dataframe speciyfing the beta and adjusted p value
#' @author Tsunghan Hsieh
#' @importFrom dplyr mutate left_join select
#' @importFrom assertthat assert_that
#' @export
#' @examples
#' \dontrun{
#' df <- getFittedResultByOneFactor(
#' count = count,
#' sample_sheet = sampleSheet_subset,
#' sample_col = "sample",
#' condition_col = "Weeks")
#' }

getFittedResultByOneFactor <- function(count, sample_sheet, sample_col, condition_col) {
  
  assertthat::assert_that(is.character(sample_col), 
                          msg = "The given argument 'sample_col' is not a string.\n")
  
  assertthat::assert_that(is.character(condition_col), 
                          msg = "The given argument 'condition_col' is not a string.\n")
  
  
  # Transform the dataframe ----------------------------------------------------
  ## Initiate a dataframe for results
  fit_result <- data.frame(matrix(ncol=3, nrow=0))
  colnames(fit_result) <- c("GENEID", "beta", "pValue")
  
  # Iterate --------------------------------------------------------------------
  ## Call a help function to do the fitting
  ## Apply the function to each row using lapply
  results_list <- lapply(1:nrow(count), .fit_row, count = count, sample_sheet = sample_sheet, sample_col = sample_col, condition_col = condition_col)
  
  ## Combine the list of data frames using do.call and rbind
  fit_result <- do.call(rbind, results_list)
  
  ## Convert fit_result to data.frame
  fit_result <- as.data.frame(fit_result)
  
  fit_result <- fit_result |>
    dplyr::mutate(adjPValue = p.adjust(pValue, method = "BH"))
  
  rownames(fit_result) <- fit_result$GENEID
  
  return(fit_result)
}

#' .fit_row
#' @description A helper function for fitting data into glm 
.fit_row <- function(i, count, sample_sheet, sample_col, condition_col) {
  df <- data.frame(sample = colnames(count), count = as.numeric(count[i,]))
  
  # Save the gene id
  GENEID <- rownames(count)[i]
  
  df <- df |>
    dplyr::mutate(log10count = log10(count + 1))
  
  # Subset sample sheet
  sample_sheet_subset <- sample_sheet[, c(sample_col, condition_col)]
  
  # Merge the data frames
  regression_df <- merge(sample_sheet_subset, df, by = sample_col)
  
  # Create the formula dynamically
  formula <- as.formula(paste("log10count ~", condition_col))
  
  # Run the glm
  glm_model <- glm(formula, data = regression_df, family = gaussian())
  
  # Extract results
  tmp_df <- data.frame(
    GENEID = GENEID,
    beta = coef(summary(glm_model))[, 1][condition_col],
    pValue = coef(summary(glm_model))[, 4][condition_col]
  )
  
  return(tmp_df)
}