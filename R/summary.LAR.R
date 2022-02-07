#'Summarize estimated Lifetime Attributable Risk for one person
#'
#'\code{summary.LAR} is the function for printing class "LAR".
#'@param object object of class 'LAR_batch' or LAR'.
#'@param digits the number of decimal points to print.
#'@param max.id the number of maximum of printing LAR results.
#'@param ... further arguments to be passed from or to other methods.
#'
#'@rdname summary.LAR
#'@export
summary.LAR <- function(object, digits=4, ...){
  result_lar <- matrix(0, nrow=length(object$LAR), ncol=5, dimnames=list(names(object$LAR), c("Lower", "Mean", "Upper", "LBR", "LFR")))
  result_F_lar <- matrix(0, nrow=length(object$F_LAR), ncol=5, dimnames=list(names(object$F_LAR), c("Lower", "Mean", "Upper", "BFR", "TFR")))


  for(i in 1:length(object$LAR)){
    result_lar[i,] <- c(object$LAR[[i]], object$LBR[i], object$LFR[i])
    result_F_lar[i,] <- c(object$F_LAR[[i]], object$BFR[i], object$TFR[i])
  }
  cat("Information:","\n")
  print(object$pinfo, row.names=FALSE)
  cat("\n")
  cat("LAR:", "\n")
  print(round(result_lar,digits=digits))
  cat("\n")
  cat("Future LAR:", "\n")
  print(round(result_F_lar,digits=digits))
  cat('\n')
  cat('Confidence Level:', paste0(round(object$ci, digits=3),"\n"))
  cat("Current Year:", object$current)
  cat('\n---\n')
}

#'@rdname summary.LAR
#'@export
summary.LAR_batch <- function(object, digits=4, max.id=50, ...){
  id <- names(object)
  for(p in 1:ifelse(length(object)<=max.id, length(object), max.id)){
    cat("summaries of LAR result : ", id[p], "\n\n")
    summary(object[[p]], digits=digits)
    cat("\n")
  }
  if(length(object)>max.id) cat("The results for", length(object)-max.id ,"people are omitted.")
}

#'@rdname summary.LAR
#'@export
summary.LAR_group <- function(object, digits=4, max.id=50, ...){
  id <- names(object)
  for(p in 1:ifelse(length(object)<=max.id, length(object), max.id)){
    cat("summaries of LAR result : Group", id[p], "\n\n")

    result_lar <- matrix(0, nrow=length(object[[p]]$LAR), ncol=5, dimnames=list(names(object[[p]]$LAR), c("Lower", "Mean", "Upper", "LBR", "LFR")))
    result_F_lar <- matrix(0, nrow=length(object[[p]]$F_LAR), ncol=5, dimnames=list(names(object[[p]]$F_LAR), c("Lower", "Mean", "Upper", "BFR", "TFR")))

    for(i in 1:length(object[[p]]$LAR)){
      result_lar[i,] <- c(object[[p]]$LAR[[i]], object[[p]]$LBR[i], object[[p]]$LFR[i])
      result_F_lar[i,] <- c(object[[p]]$F_LAR[[i]], object[[p]]$BFR[i], object[[p]]$TFR[i])
    }
    cat("Group Information:","\n")
    print(object[[p]]$pinfo, row.names=FALSE)
    cat("\n")
    cat("LAR:", "\n")
    print(round(result_lar,digits=digits))
    cat("\n")
    cat("Future LAR:", "\n")
    print(round(result_F_lar,digits=digits))
    cat('\n')
    cat('Confidence Level:', paste0(round(object[[p]]$ci, digits=3),"\n"))
    cat("Current Year:", object[[p]]$current)
    cat('\n---\n\n')
  }
  if(length(object)>max.id) cat("The results for", length(object)-max.id ,"groups are omitted.")
}
