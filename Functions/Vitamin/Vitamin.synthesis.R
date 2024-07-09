Vitamin.synthesis <- function(ko){
    Vitamin <- c("K01077","K01078","K23458","K14394","K06949","K22911",
                 "K00793",
                 "K07758","K13248","K05275","K18607","K00275","K23998",
                 "K01012","K01435","K08351",
                 "K13951","K13980","K13952","K13953","K00001","K00121","K15734",
                 "K13369","K11150","K11151","K11152","K11153","K11161","K11154",
                 "K11146","K11147","K11148","K11149")
    Vitamin.name <- c("phoA","PHO","PHO5_11_12","ACP1","rsgA","TH2",
                      "ribE",
                      "PDXP","PHOSPHO2","P4DH","pno","pdxH","PPOX",
                      "bioB","BTD","bisC",
                      "ADH1_7","ADH4","ADH6","adhP","adh","frmA","SDR16C5",
                      "HSD17B6","RDH8","RDH10","RDH11","RDH12","RDH13","RDH16",
                      "DHRS3","DHRS4","DHRS4L2","DHRS9")
    Vitamin.type <- c(rep("Vitamin B1",6),"Vitamin B2",rep("Vitamin B6",6),
                      rep("Vitamin H",3),rep("Vitamin A",18))
    Vitamin.ko <- c()
    for (i in 1:34) {
        if (Vitamin[i] %in% rownames(ko)) {
            C2 <- ko[Vitamin[i],]
        }else{
            C2 <- rep(0,ncol(ko))
        }
        Vitamin.ko <- rbind(Vitamin.ko,C2)
    }
    rownames(Vitamin.ko) <- Vitamin.name
    Vitamin.ko <- as.data.frame(t(Vitamin.ko[,-ncol(Vitamin.ko)]))
    Vitamin.ko <- as.data.frame(t(Vitamin.ko))
    Vitamin.ko$Gene <- rownames(Vitamin.ko)
    Vitamin.ko$Type <- Vitamin.type
    Vitamin.ko <- Vitamin.ko[,c("Type","Gene",colnames(Vitamin.ko)[1:(ncol(Vitamin.ko)-2)])]
    Vitamin.ko <- Vitamin.ko[rowSums(Vitamin.ko[,3:ncol(Vitamin.ko)]) > 0,]
    return(Vitamin.ko)
}