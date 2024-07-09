Vitamin.metabolism <- function(ko3){
    aa <- c("ko00730","ko00740","ko00750","ko00780","ko00830")
    aam <- ko3[rownames(ko3) %in% aa,]
    aam <- aam[rowSums(aam[,1:(ncol(aam)-1)]) > 0,]
    aam <- aam[,c("Description",colnames(aam)[1:(ncol(aam)-1)])]
    aam$Description <- c("Vitamin B1 metabolism","Vitamin B2 metabolism",
                         "Vitamin B6 metabolism","Vitamin H metabolism",
                         "Vitamin A metabolism")
    return(aam)
}