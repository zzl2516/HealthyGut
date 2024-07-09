Propionate.pathway <- function(ko){
    Propionate.ko <- c("K01899","K01900","K01902","K01903",
                       "K01847","K01848","K01849",
                       "K05606",
                       "K11264","K08426","K03416",
                       "K00626","K00627",
                       "K00232","K00248","K19745","K20143",
                       "K01026","K01905","K22224","K24012",
                       "K01895","K01908",
                       "K13788","K00625","K15024","K13923",
                       "K00925","K00932","K19697")
    
    Propionate.name <- c("LSC1","LSC2","sucD","sucC",
                         "MUT","mcmA1","mcmA2",
                         "epi",
                         "mmcD","G2A","mmcC5S",
                         "ACAT","DLAT",
                         "ACOX1","ACADS","acul","acrC",
                         "pct","acdA","acdB","acdAB",
                         "ACSS1_2","ACSS3",
                         "pta1","pta2","pta3","pduL",
                         "ackA","tdcD","pduW")
    
    Propionate <- c()
    for (i in 1:30) {
        if (Propionate.ko[i] %in% rownames(ko)) {
            C2 <- ko[Propionate.ko[i],]
        }else{
            C2 <- t(as.data.frame(rep(0,ncol(ko))))
            colnames(C2) <- colnames(ko)
        }
        Propionate <- rbind(Propionate,C2)
    }
    rownames(Propionate) <- Propionate.name
    Propionate <- as.data.frame(t(Propionate[,-ncol(Propionate)]))
    
    Propionate$acdAB <- (Propionate$acdA + Propionate$acdB)/2 + Propionate$acdAB
    Propionate <- Propionate[,-c(19,20)]
    Propionate$LSC1 <- (Propionate$LSC1 + Propionate$LSC2)/2
    Propionate <- Propionate[,-2]
    colnames(Propionate)[1] <- "LSC"
    Propionate$sucD <- (Propionate$sucD + Propionate$sucC)/2
    Propionate <- Propionate[,-3]
    colnames(Propionate)[2] <- "sucCD"
    Propionate$mcmA1 <- (Propionate$mcmA1 + Propionate$mcmA2)/2
    Propionate <- Propionate[,-5]
    colnames(Propionate)[4] <- "mcmA"
    Propionate$pta1 <- Propionate$pta1 + Propionate$pta2 + Propionate$pta3
    Propionate <- Propionate[,-c(20,21)]
    colnames(Propionate)[19] <- "PTA"
    Propionate <- as.data.frame(t(Propionate))
    Propionate$Gene <- rownames(Propionate)
    Propionate <- Propionate[,c("Gene",colnames(Propionate)[1:(ncol(Propionate)-1)])]
    Propionate <- Propionate[rowSums(Propionate[,2:ncol(Propionate)]) > 0,]
    return(Propionate)
}