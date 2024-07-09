Immune.synthesis <- function(ko){
    Immune <- c("K00234","K00235","K00236","K00237","K25801","K00239","K00240",
                 "K00241","K00242","K18859","K18860","K00244","K00245","K00246",
                 "K00247","K00233","K25995","K25996","K01899","K01900","K01902",
                 "K01903","K18118","K00135","K08324","K00139","K17761","K15737",
                 "K05714","K01637","K14471","K14472","K10764","K01739",
                 "K00128","K14085","K00149","K11817","K22417","K01501","K01426",
                 "K21801","K11816",
                 "K04103","K11182","K00274",
                 "K01432","K14263","K07130","K00463","K00486","K00453",
                 "K01442","K22605","K15868","K15870","K15872","K15871","K15874",
                 "K07007",
                 "K18277","K21579","K20038","K20037","K22443","K22444")
    Immune.name <- c("SDHA","SDHB","SDHC","SDHD","SDH4","sdhA","sdhB","sdhC",
                      "sdhD","sdhD1","sdhD2","frdA","frdB","frdC","frdD","frdA1",
                      "frdB1","frdC1","LSC1","LSC2","sucD","sucC","aarC","gabD",
                      "sad","ALDH5A1","SSADH","csiD","mhpC","aceA","smtA","smtB",
                      "metZ","metB",
                      "ALDH","ALDH7A1","ALDH9A1","AAO1_2","AAO4","nitrilase",
                      "amiE","iaaH","YUCCA",
                      "ipdC","AOC1","MAO",
                      "AFMID","BNA7","kynB","IDO","KMO","TDO2",
                      "cbh","baiA","baiB","baiCD","baiE","baiF","baiI","baiN",
                      "tmm","grdH","cutC","cutD","cntA","cntB")
    Immune.type <- c(rep("Succinate",14),rep("Indole-3-acetic acid",9),
                     rep("Indole-3-aldehyde",3),rep("Kynurenine",6),
                     rep("Secondary bile acids",8),rep("Trimethylamine N-oxide",4))
    Immune.ko <- c()
    for (i in 1:66) {
        if (Immune[i] %in% rownames(ko)) {
            C2 <- ko[Immune[i],]
        }else{
            C2 <- as.data.frame(t(rep(0,ncol(ko))))
            colnames(C2) <- colnames(ko)
        }
        Immune.ko <- rbind(Immune.ko,C2)
    }
    rownames(Immune.ko) <- Immune.name
    Immune.ko <- as.data.frame(t(Immune.ko[,-ncol(Immune.ko)]))
    
    Immune.ko$cutC <- (Immune.ko$cutC + Immune.ko$cutD)/2
    Immune.ko$cntA <- (Immune.ko$cntA + Immune.ko$cntB)/2
    colnames(Immune.ko)[63] <- "cutCD"
    colnames(Immune.ko)[65] <- "cntAB"
    Immune.ko <- Immune.ko[,-c(64,66)]
    Immune.ko$SDHA <- (Immune.ko$SDHA + Immune.ko$SDHB + Immune.ko$SDHC + Immune.ko$SDHD + Immune.ko$SDH4)/5 + 
        (Immune.ko$sdhA + Immune.ko$sdhB + Immune.ko$sdhC + Immune.ko$sdhD + Immune.ko$sdhD1 + Immune.ko$sdhD2)/4 + 
        (Immune.ko$frdA + Immune.ko$frdA1 + Immune.ko$frdB + Immune.ko$frdB1 + Immune.ko$frdC + Immune.ko$frdC1 + Immune.ko$frdD)/4
    colnames(Immune.ko)[1] <- "SDH"
    Immune.ko <- Immune.ko[,-c(2:18)]
    Immune.ko$LSC1 <- (Immune.ko$LSC1 + Immune.ko$LSC2)/2
    colnames(Immune.ko)[2] <- "LSC"
    Immune.ko <- Immune.ko[,-3]
    Immune.ko$sucD <- (Immune.ko$sucD + Immune.ko$sucC)/2
    colnames(Immune.ko)[3] <- "sucCD"
    Immune.ko <- Immune.ko[,-4]
    Immune.ko$smtA<- (Immune.ko$smtA + Immune.ko$smtB)/2
    colnames(Immune.ko)[12] <- "smtAB"
    Immune.ko <- Immune.ko[,-13]
    
    Immune.ko <- as.data.frame(t(Immune.ko))
    Immune.ko$Gene <- rownames(Immune.ko)
    Immune.ko$Type <- Immune.type
    Immune.ko <- Immune.ko[,c("Type","Gene",colnames(Immune.ko)[1:(ncol(Immune.ko)-2)])]
    Immune.ko <- Immune.ko[rowSums(Immune.ko[,3:ncol(Immune.ko)]) > 0,]
    return(Immune.ko)
}