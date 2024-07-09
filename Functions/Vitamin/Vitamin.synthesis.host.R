Vitamin.synthesis.host <- function(Gene){
    Vitamin.path <- list(VB1 = c("K01077","K01078","K23458","K14394","K06949","K22911"),
                    VB2 = "K00793",
                    VB6 = c("K07758","K13248","K05275","K18607","K00275","K23998"),
                    VH = c("K01012","K01435","K08351"),
                    VA = c("K13951","K13980","K13952","K13953","K00001","K00121","K15734",
                           "K13369","K11150","K11151","K11152","K11153","K11161","K11154",
                           "K11146","K11147","K11148","K11149"))
    temp <- c()
    for (i in 1:5) {
        d <- Vitamin.path[[i]]
        d <- ifelse(d %in% colnames(Gene),d,NA)
        d <- d[complete.cases(d)]
        a <- ifelse(length(d) == 0,0,1)
        temp <- c(temp,a)
    }
    result <- list(temp,Vitamin.path)
    return(result)
}