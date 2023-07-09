#' Increment Peptide Sequences
#'
#' This function generates new peptide sequences by adding each of the 20 amino acids
#' to each position of the input peptide or peptides.
#'
#' @param peptides A character vector of peptide sequences.
#' @param num_added The number of amino acids to be added to each position of the peptide.
#' @return A character vector of new peptide sequences.
#' @examples
#' increment(c("AC", "DE"))
#' increment("ACDE", num_added = 2)
#' @export
#'
increment <- function(peptides, num_added = 1){
  aminoAcids <- c('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')

  # Ensure the input is a character vector
  if (!is.character(peptides)) {
    stop("Input peptides must be a character vector.")
  }

  # Function to create new sequences from a single peptide
  create_sequences <- function(peptide, num_added) {
    original_length <- nchar(peptide)
    new_sequences <- c()
    for (i in 0:original_length) {
      a_combinations = expand.grid(rep(list(aminoAcids), num_added))
      for (j in 1:nrow(a_combinations)) {
        aa_comb = apply(a_combinations[j, , drop = FALSE], 1, paste, collapse = "")
        new_seq <- paste0(substr(peptide, 1, i), aa_comb, substr(peptide, i + 1, original_length))
        new_sequences <- c(new_sequences, new_seq)
      }
    }
    return(new_sequences)
  }

  # Generate new sequences from each peptide
  result <- unlist(lapply(peptides, create_sequences, num_added = num_added))
  return(result)
}
