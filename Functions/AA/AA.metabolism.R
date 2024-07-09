AA.metabolism <- function(ko3){
    aa <- c("ko00250","ko00260","ko00270","ko00290","ko00300",
            "ko00220","ko00330","ko00340","ko00350","ko00360","ko00380","ko00400")
    aam <- ko3[rownames(ko3) %in% aa,]
    aam <- aam[rowSums(aam[,1:(ncol(aam)-1)]) > 0,]
    aam <- aam[,c("Description",colnames(aam)[1:(ncol(aam)-1)])]
    return(aam)
}