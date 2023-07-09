#' Select Best vs Worst Peptides
#'
#' This function identifies the peptides from the function *appearance_to_binary*
#' that are 1 in one group and 0 or -1 in another group, and expands the grid
#' to all possible combinations.
#'
#' @param appearance_best A matrix with transformed counts for the 'best' group.
#' @param appearance_worst A matrix with transformed counts for the 'worst' group.
#' @return A data frame with combinations of 'best' peptides.
#' @examples
#' # Generate some mock data
#' appearance_best <- matrix(c(1, -1, 0, 1, -1), nrow = 5, ncol = 4)
#' appearance_worst <- matrix(c(-1, 1, 0, -1, 1), nrow = 5, ncol = 4)
#' # Call the function
#' select_best_vs_worst(appearance_best, appearance_worst)
#' @export
#'
select_best_vs_worst <- function(appearance_best, appearance_worst) {
  x <- as.data.frame(t(appearance_best))
  y <- as.data.frame(t(appearance_worst))

  poss <- list()
  for (j in 1:length(x)) {
    poss[[j]] <-  dplyr::setdiff(rownames(dplyr::filter(x, x[,j] == 1)),
                                 rownames(dplyr::filter(y, y[,j] == 1)))
  }

  best_pred <- expand.grid(poss)
  best_pred <- data.frame(Title = apply(best_pred, 1, paste, collapse = ""))

  return(best_pred)
}
