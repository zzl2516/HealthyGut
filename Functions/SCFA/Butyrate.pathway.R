Butyrate.pathway <- function(ko){
    Butyrate.ko <- c("K00626",
                     "K00074","K00022","K07516","K01825","K01782","K07514",
                     "K01692","K07515","K07511","K01715",
                     "K17829","K00209",
                     "K01034","K01035","K19709",
                     "K00634",
                     "K00929")
    
    Butyrate.name <- c("ACAT",
                       "HBDH","HADH","fadN","fadB","fadJ","EHHADH",
                       "echA","HADHA","ECHS1","crt",
                       "ccrA","fabV",
                       "atoD","atoA","ydiF",
                       "PTB",
                       "buk")
    
    Butyrate <- c()
    for (i in 1:18) {
        if (Butyrate.ko[i] %in% rownames(ko)) {
            C2 <- ko[Butyrate.ko[i],]
        }else{
            C2 <- t(as.data.frame(rep(0,ncol(ko))))
            colnames(C2) <- colnames(ko)
        }
        Butyrate <- rbind(Butyrate,C2)
    }
    rownames(Butyrate) <- Butyrate.name
    Butyrate <- as.data.frame(t(Butyrate[,-ncol(Butyrate)]))
    
    Butyrate$atoA <- (Butyrate$atoA + Butyrate$atoD)/2
    Butyrate <- Butyrate[,-14]
    colnames(Butyrate)[1] <- "atoAD"
    
    Butyrate <- as.data.frame(t(Butyrate))
    Butyrate$Gene <- rownames(Butyrate)
    Butyrate <- Butyrate[,c("Gene",colnames(Butyrate)[1:(ncol(Butyrate)-1)])]
    Butyrate <- Butyrate[rowSums(Butyrate[,2:ncol(Butyrate)]) > 0,]
    return(Butyrate)
}