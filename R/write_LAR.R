#'Write a LAR object
#'
#'Write 'LAR' object to CSV file
#'
#'@param x a 'LAR' object.
#'@param filename a string naming the file to save (.csv file)
#'@export
write_LAR <- function(x, filename) UseMethod('write_LAR')

#' @describeIn write_LAR  write an 'LAR' class object
#' @export
#' @importFrom utils write.csv
write_LAR.LAR <- function(x, filename){
  result_lar <- matrix(0, nrow=length(x$LAR), ncol=3, dimnames=list(names(x$LAR), c("Lower", "Mean", "Upper")))
  result_F_lar <- matrix(0, nrow=length(x$F_LAR), ncol=3, dimnames=list(names(x$F_LAR), c("F.Lower", "F.Mean", "F.Upper")))

  for(i in 1:length(x$LAR)){
    result_lar[i,] <- x$LAR[[i]]
    result_F_lar[i,] <- x$F_LAR[[i]]
  }
  write.csv(cbind(result_lar, result_F_lar, LBR=x$LBR, BFR=x$BFR, LFR=x$LFR, TFR=x$TFR), file=filename)
}


#' @describeIn write_LAR  write an 'LAR_batch' class object
#' @export
#' @importFrom utils write.csv
write_LAR.LAR_batch <- function(x, filename){
  site_list <- c('stomach', 'colon', 'liver', 'lung', 'breast', 'ovary', 'uterus', 'prostate', 'bladder', 'brain/cns', 'thyroid', 'remainder', 'oral', 'oesophagus', 'rectum', 'gallbladder', 'pancreas', 'kidney', 'leukemia')

  result_lar <- matrix(0, nrow=length(x), ncol=63, dimnames=list(NULL, paste0(rep(c("L_", "M_", "U_"), each=21), rep(c("total", "solid", site_list), 3))))
  result_F_lar <- matrix(0, nrow=length(x), ncol=63, dimnames=list(NULL, paste0(rep(c("F_L_", "F_M_", "F_U_"), each=21), rep(c("total", "solid", site_list), 3))))
  result_lbr <- matrix(0, nrow=length(x), ncol=21, dimnames=list(NULL, paste0("LBR_", c("total", 'solid', site_list))))
  result_bfr <- matrix(0, nrow=length(x), ncol=21, dimnames=list(NULL, paste0("BFR_", c("total", 'solid', site_list))))
  result_lfr <- matrix(0, nrow=length(x), ncol=21, dimnames=list(NULL, paste0("LFR_", c("total", 'solid', site_list))))
  result_tfr <- matrix(0, nrow=length(x), ncol=21, dimnames=list(NULL, paste0("TFR_", c("total", 'solid', site_list))))

  for(i in 1:length(x)){
    result_lar[i, paste0(rep(c("L_", "M_", "U_"), length(x[[i]]$LAR)), rep(names(x[[i]]$LAR), each=3))] <-
      unlist(lapply(x[[i]]$LAR, function(s) s[c("Lower", "Mean", "Upper")]))
    result_F_lar[i, paste0(rep(c("F_L_", "F_M_", "F_U_"), length(x[[i]]$F_LAR)), rep(names(x[[i]]$F_LAR), each=3))] <-
      unlist(lapply(x[[i]]$F_LAR, function(s) s[c("Lower", "Mean", "Upper")]))
    result_lbr[i, paste0("LBR_", names(x[[i]]$LBR))] <- x[[i]]$LBR
    result_bfr[i, paste0("BFR_", names(x[[i]]$BFR))] <- x[[i]]$BFR
    result_lfr[i, paste0("LFR_", names(x[[i]]$LFR))] <- x[[i]]$LFR
    result_tfr[i, paste0("TFR_", names(x[[i]]$TFR))] <- x[[i]]$TFR
  }
  write.csv(data.frame(PID=names(x),cbind(result_lar, result_F_lar, result_lbr, result_bfr, result_lfr, result_tfr)), file=filename, row.names=FALSE)
}


#' @describeIn write_LAR  write an 'LAR_group' class object
#' @export
#' @importFrom utils write.csv
write_LAR.LAR_group <- function(x, filename){
  site_list <- c('stomach', 'colon', 'liver', 'lung', 'breast', 'ovary', 'uterus', 'prostate', 'bladder', 'brain/cns', 'thyroid', 'remainder', 'oral', 'oesophagus', 'rectum', 'gallbladder', 'pancreas', 'kidney', 'leukemia')

  result_lar <- matrix(0, nrow=length(x), ncol=63, dimnames=list(NULL, paste0(rep(c("L_", "M_", "U_"), each=21), rep(c("total", "solid", site_list), 3))))
  result_F_lar <- matrix(0, nrow=length(x), ncol=63, dimnames=list(NULL, paste0(rep(c("F_L_", "F_M_", "F_U_"), each=21), rep(c("total", "solid", site_list), 3))))
  result_lbr <- matrix(0, nrow=length(x), ncol=21, dimnames=list(NULL, paste0("LBR_", c("total", 'solid', site_list))))
  result_bfr <- matrix(0, nrow=length(x), ncol=21, dimnames=list(NULL, paste0("BFR_", c("total", 'solid', site_list))))
  result_lfr <- matrix(0, nrow=length(x), ncol=21, dimnames=list(NULL, paste0("LFR_", c("total", 'solid', site_list))))
  result_tfr <- matrix(0, nrow=length(x), ncol=21, dimnames=list(NULL, paste0("TFR_", c("total", 'solid', site_list))))

  for(i in 1:length(x)){
    result_lar[i, paste0(rep(c("L_", "M_", "U_"), length(x[[i]]$LAR)), rep(names(x[[i]]$LAR), each=3))] <-
      unlist(lapply(x[[i]]$LAR, function(s) s[c("Lower", "Mean", "Upper")]))
    result_F_lar[i, paste0(rep(c("F_L_", "F_M_", "F_U_"), length(x[[i]]$F_LAR)), rep(names(x[[i]]$F_LAR), each=3))] <-
      unlist(lapply(x[[i]]$F_LAR, function(s) s[c("Lower", "Mean", "Upper")]))
    result_lbr[i, paste0("LBR_", names(x[[i]]$LBR))] <- x[[i]]$LBR
    result_bfr[i, paste0("BFR_", names(x[[i]]$BFR))] <- x[[i]]$BFR
    result_lfr[i, paste0("LFR_", names(x[[i]]$LFR))] <- x[[i]]$LFR
    result_tfr[i, paste0("TFR_", names(x[[i]]$TFR))] <- x[[i]]$TFR
  }
  write.csv(data.frame(GROUP=names(x),cbind(result_lar, result_F_lar, result_lbr, result_bfr, result_lfr, result_tfr)), file=filename, row.names=FALSE)
}

