SCFA.synthesis <- function(ko){
    SCFA <- c("K00656",
              "K00925","K01512",
              "K01895","K01913",
              "K01026","K18118","K01905","K22224","K24012",
              "K01067",
              "K00128","K00138","K00149","K14085",
              "K00467",
              "K00156",
              "K01908",
              "K00932","K19697",
              "K01034","K01035","K19709",
              "K00929")
    SCFA.name <- c("pflD",
                   "ackA","acyP",
                   "ACSS1_2","ACN1",
                   "pct","aarC","acdA","acdB","acdAB",
                   "ACH1",
                   "ALDH","aldB","ALDH9A1","ALDH7A1",
                   "L2MO",
                   "poxB",
                   "ACSS3",
                   "tdcD","pduW",
                   "atoD","atoA","ydiF",
                   "buk")
    SCFA.type <- c("Formate",
                   "Acetate_Propionate","Acetate",
                   "Acetate_Propionate","Acetate",
                   "Acetate_Propionate","Acetate","Acetate_Propionate",
                   "Acetate",
                   "Acetate","Acetate","Acetate","Acetate",
                   "Acetate",
                   "Acetate",
                   "Propionate",
                   "Propionate","Propionate",
                   "Butyrate","Butyrate",
                   "Butyrate")
    SCFA.ko <- c()
    for (i in 1:24) {
        if (SCFA[i] %in% rownames(ko)) {
            C2 <- ko[SCFA[i],]
        }else{
            C2 <- rep(0,ncol(ko))
        }
        SCFA.ko <- rbind(SCFA.ko,C2)
    }
    rownames(SCFA.ko) <- SCFA.name
    SCFA.ko <- as.data.frame(t(SCFA.ko[,-ncol(SCFA.ko)]))
    
    SCFA.ko$acdAB <- (SCFA.ko$acdA + SCFA.ko$acdB)/2 + SCFA.ko$acdAB
    SCFA.ko <- SCFA.ko[,-c(8,9)]
    SCFA.ko$atoA <- (SCFA.ko$atoA + SCFA.ko$atoD)/2
    SCFA.ko <- SCFA.ko[,-19]
    colnames(SCFA.ko)[19] <- "atoAD"
    SCFA.ko <- as.data.frame(t(SCFA.ko))
    SCFA.ko$Gene <- rownames(SCFA.ko)
    SCFA.ko$Type <- SCFA.type
    SCFA.ko <- SCFA.ko[,c("Type","Gene",colnames(SCFA.ko)[1:(ncol(SCFA.ko)-2)])]
    SCFA.ko <- SCFA.ko[rowSums(SCFA.ko[,3:ncol(SCFA.ko)]) > 0,]
    return(SCFA.ko)
}