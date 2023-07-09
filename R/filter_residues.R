#' Filter Peptides by Residue Counts
#'
#' This function counts the number of specified residues in each peptide sequence
#' and filters out the ones with more than the specified limit. It's defaults is
#' for filtering out small alliphatic residues.
#'
#' @param df A data frame containing peptide sequences.
#' @param sequence_col The name of the column that contains the sequences.
#' @param residues A character vector of residues to count.
#' @param max_residues The maximum number of allowed residues.
#' @return A filtered data frame.
#' @examples
#' # Generate a mock data frame
#' peptide_data <- data.frame(Sequence = c("AVILG", "VILGA", "ILGAV", "LGAVI"))
#' # Apply the function to the mock data
#' filter_residues(peptide_data, residues = c("A", "V", "I", "L", "G"), max_residues = 2)
#' @export
#'
filter_residues <- function(df, sequence_col = "Sequence", residues = c("A", "V", "I", "L", "G"), max_residues = 2) {
  # Ensure the sequence column exists
  if (!(sequence_col %in% names(df))) {
    stop(paste("Column", sequence_col, "not found in the data frame."))
  }

  # Count the number of specified residues in each sequence
  df$ResidueCount <- rowSums(sapply(residues, function(r) stringr::str_count(df[[sequence_col]], r)))

  # Filter out the sequences with more than the specified number of residues
  filtered_df <- dplyr::filter(df, df$ResidueCount <= max_residues)

  return(filtered_df)
}
