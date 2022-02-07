#'Print estimated Lifetime Attributable Risk for one person
#'
#'\code{print.LAR} is the basic function for printing class "LAR".
#'
#'@param x 'LAR', 'LAR_batch' or 'LAR_group' object.
#'@param digits the number of decimal points to print.
#'@param max.id the number of maximum of printing LAR results.
#'@param ... further arguments to be passed from or to other methods.
#'
#'@rdname print.LAR
#'@export
print.LAR <- function(x, digits=4, ...){
  result_lar <- x$LAR$total
  result_F_lar <- matrix(c(x$F_LAR$total, rep(x$BFR["total"], 3), (x$F_LAR$total + x$BFR["total"])),
                         nrow=3, ncol=3, byrow=TRUE, list(c("F.LAR", "BFR", "TFR"), c("Lower", "Mean", "Upper")))

  cat("LAR:", "\n")
  print(round(result_lar,digits=digits), ...)
  cat("\n")
  cat("Future LAR:", "\n")
  print(round(result_F_lar,digits=digits), ...)
  cat("---\n")
}

#'@rdname print.LAR
#'@export
print.LAR_batch <- function(x, digits=4, max.id=50, ...){
  id <- names(x)
  for(p in 1:ifelse(length(x)<=max.id, length(x), max.id)){
    cat("LAR result of", id[p], "\n\n")
    print(x[[p]], digits=digits, ...)
    cat("\n")
  }
  if(length(x)>max.id) cat("The results for", length(x)-max.id ,"people are omitted.")
}

#'@rdname print.LAR
#'@export
print.LAR_group <- function(x, digits=4, max.id=50, ...){
  id <- names(x)
  for(p in 1:ifelse(length(x)<=max.id, length(x), max.id)){
    cat("LAR result of", id[p], "\n\n")
    print(x[[p]], digits=digits, ...)
    cat("\n")
  }
  if(length(x)>max.id) cat("The results for", length(x)-max.id ,"groups are omitted.")
}
