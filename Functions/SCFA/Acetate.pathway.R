Acetate.pathway <- function(ko){
    Acetate.ko <- c("K00158",
                    "K17388","K15024","K00625","K04020",
                    "K00925","K01512",
                    "K01895","K01913",
                    "K01026","K18118","K01905","K22224","K24012",
                    "K01067",
                    "K00174","K00175",
                    "K00169","K00170","K00171","K00172","K00189","K03737",
                    "K00161","K00162","K00163","K00627",
                    "K00016","K00467",
                    "K00156",
                    "K00132","K18366","K04072","K04073","K04021",
                    "K00128","K00138","K00149","K14085",
                    "K00656")
    
    Acetate.name <- c("poxL",
                      "ROCK2","PTA","pta","eutD",
                      "ackA","acyP",
                      "ACSS1_2","ACN1",
                      "pct","aarC","acdA","acdB","acdAB",
                      "ACH1",
                      "korA","korB",
                      "porA","porB","porD","porC","porG","nifJ",
                      "pdhA","pdhB","aceE","aceF",
                      "LDH","L2MO",
                      "poxB",
                      "AcAldDH","bphJ","adhE","mhpF","eutE",
                      "ALDH","aldB","ALDH9A1","ALDH7A1",
                      "pflD")
    
    Acetate <- c()
    for (i in 1:40) {
        if (Acetate.ko[i] %in% rownames(ko)) {
            C2 <- ko[Acetate.ko[i],]
        }else{
            C2 <- t(as.data.frame(rep(0,ncol(ko))))
            colnames(C2) <- colnames(ko)
        }
        Acetate <- rbind(Acetate,C2)
    }
    rownames(Acetate) <- Acetate.name
    Acetate <- as.data.frame(t(Acetate[,-ncol(Acetate)]))
    
    Acetate$acdAB <- (Acetate$acdA + Acetate$acdB)/2 + Acetate$acdAB
    Acetate <- Acetate[,-c(12,13)]
    Acetate$PTA <- Acetate$PTA + Acetate$pta
    Acetate <- Acetate[,-4]
    Acetate$porA <- (Acetate$porA + Acetate$porB + Acetate$porC + Acetate$porD +
                         Acetate$porG)/5 + Acetate$nifJ
    Acetate <- Acetate[,-c(16:20)]
    colnames(Acetate)[15] <- "por"
    Acetate$pdhA <- (Acetate$pdhA + Acetate$pdhB + Acetate$aceE + Acetate$aceF)/4
    Acetate <- Acetate[,-c(17:19)]
    colnames(Acetate)[16] <- "pdh"
    Acetate$korA <- (Acetate$korA + Acetate$korB)/2
    Acetate <- Acetate[,-14]
    colnames(Acetate)[13] <- "kor"
    Acetate <- as.data.frame(t(Acetate))
    Acetate$Gene <- rownames(Acetate)
    Acetate <- Acetate[,c("Gene",colnames(Acetate)[1:(ncol(Acetate)-1)])]
    Acetate <- Acetate[rowSums(Acetate[,2:ncol(Acetate)]) > 0,]
    return(Acetate)
}