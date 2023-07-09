#' Count Amino Acids
#'
#' This function counts the occurrence of each of the 20 amino acids
#' at each of the first 'n' positions across a vector of peptide sequences.
#'
#' @param peptides A character vector of peptide sequences.
#' @param n The number of initial positions to consider in each peptide sequence.
#' @return A data frame with 'n' rows and 20 columns where each row represents a position
#' in the peptide sequence and each column represents an amino acid. Each cell
#' in the data frame contains the count of a particular amino acid at a particular position.
#' @examples
#' count_aa(c("ACDF", "BCDE", "ABCD"), n = 2)
#' @export
#'
count_aa <- function(peptides, n = 4) {
  # List of amino acids
  aminoAcids <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")

  # Create an empty data frame to store the counts
  aa_counts <- data.frame(matrix(nrow = n, ncol = 20))
  colnames(aa_counts) <- aminoAcids

  # Count the occurrence of each amino acid at each position
  for (w in 1:n) {
    position_aa <- substr(peptides, w, w)
    aa_counts[w,] <- table(factor(position_aa, levels = aminoAcids))
  }

  return(aa_counts)
}
