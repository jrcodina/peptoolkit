#' Transform Amino Acid Appearance Probability into -1, 0, or 1
#'
#' This function transforms the counts of amino acids to a -1, 0, 1 matrix
#' based on a probability of appearance of each peptide in each position.
#'
#' @param x A data frame containing peptide sequences.
#' @param threshold The probability threshold to determine the transformation.
#' @param group A character string indicating which part of the data to consider. Either 'Best' or 'Worst'.
#' @param percentage The percentage of the data to consider, if group is specified.
#' @return A matrix with the same dimensions as the input where each cell
#' has been transformed to -1, 0, or 1 based on the probability threshold.
#' @examples
#' # Generate a mock data frame
#' peptide_data <- data.frame(Sequence = c("ACGT", "TGCA", "GATC", "CGAT"))
#'
#' # Apply the function to the mock data
#' appearance_to_binary(peptide_data, group = "Best", percentage = 0.5)
#' @export
#'
appearance_to_binary <- function(x, threshold = 1.65, group = "Best", percentage = 0.05) {
  # Check if the input data frame has the required column
  if (!"Sequence" %in% names(x)) {
    stop("Input data frame must contain a 'Sequence' column.")
  }

  total <- count_aa(x$Sequence)
  total_sample <- nrow(x)

  # Determine which part of the data to consider based on the 'group' parameter
  if (group == "Best") {
    sample_size = length(x$Sequence[1:(nrow(x)*percentage)])
    selected = count_aa(x$Sequence[1:(nrow(x)*percentage)])
  } else if (group == "Worst") {
    sample_size = length(x$Sequence[(nrow(x)-nrow(x)*percentage):nrow(x)])
    selected = count_aa(x$Sequence[(nrow(x)-nrow(x)*percentage):nrow(x)])
  } else {
    stop("Invalid group. Choose either 'Best' or 'Worst'.")
  }

  # Replace zeros in 'selected' with a small value to avoid division by zero
  selected[selected < 1] <- 1e-10

  # Calculate the probability
  probability_matrix <- ((selected/sample_size) - ((total-selected) / (total_sample-sample_size))) /
    (sqrt((total/total_sample)*(1-total/total_sample))* sqrt(1/sample_size + 1/(total_sample-sample_size)))

  # Replace Inf with -Inf to denote low probability
  probability_matrix[probability_matrix == Inf] <- -Inf

  # Transform the probabilities to -1, 0, 1 based on the threshold
  transformed <- ifelse(probability_matrix > threshold, 1, ifelse(probability_matrix < -threshold, -1, 0))

  return(transformed)
}
