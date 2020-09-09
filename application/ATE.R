#' Compute Average Treatment Effect
#'
#' \code{ATE} computes average treatment effect of the matched subsamples by
#' a weighted average of the estimated treatment effects in each matched group.
#' The weight is the number of matched units.
#'
#' @param FLAME_object object returned by applying the FLAME algorithm
#'   (\code{\link{FLAME_bit}}, \code{\link{FLAME_PostgreSQL}}, or
#'   \code{\link{FLAME_SQLite}})
#' @examples
#' data(toy_data)
#' result <- FLAME::FLAME_bit(data = toy_data, holdout = toy_data)
#' FLAME::ATE(result)
#' @return average treatment effect (ATE) of the matched subsamples
#' @export

ATE <- function(FLAME_object) {
  size <- sapply(FLAME_object$MGs, length)
  return(sum(FLAME_object$CATE * size) / sum(size))
}