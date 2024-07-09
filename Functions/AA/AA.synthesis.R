AA.synthesis <- function(ko){
    AA <- c("K00259","K19244","K00814","K09758",
            "K01437","K01424","K13051","K01953",
            "K00260","K15371","K00261","K00262","K00264","K00265","K00266","K00294","K13821","K01425","K23265",
            "K00301","K00302","K00303","K00304","K00305","K00306","K00314","K21833","K21834","K25960","K25961","K00552","K00639","K20801",
            "K01079","K02203","K22305","K22528","K12235","K25316","K25317",
            "K01733",
            "K01755","K14681","K01478","K13240","K13241","K13242","K13247","K00613","K24954","K24955","K23233","K00296","K24952","K24953",
            "K00286","K01750","K19743","K19742","K01777","K01259","K11142",
            "K00013","K14152","K05604","K08660",
            "K04518","K14170","K01713","K05359",
            "K00220","K24018","K15226","K15227","K00500","K01688",
            "K01694","K01695","K01696","K06001","K00837",
            "K00830",
            "K05597",
            "K14454","K14455","K00811","K00812","K00813","K11358",
            "K00600",
            "K01620",
            "K00827","K14272",
            "K15849","K00817","K05821","K00270",
            "K01738","K13034","K17069","K10150","K00816","K01760","K14155","K01758","K17217",
            "K00815","K00838","K00832",
            "K00544","K00547","K00548","K24024","K00549","K08968","K01740","K08969","K23977",
            "K14260",
            "K00826","K00263","K00835",
            "K01586","K12526","K00290","K05831")
    AA.name <- c("ald","ala","ALT","asdA",
                 "aspA","ansB","ASRGL1","asnB",
                 "gudB","GDH2","gdhA","gdhA1","GLT1","gltB","gltD","P5CDH","putA","glsA","purQ",
                 "Sox","soxA","soxB","soxD","soxG","PIPOX","SARDH","dgcA","dgcB","etfA","etfB","GNMT","kbl","ItaA",
                 "serB","thrH","psp","ysaA","SRR","racX","alr",
                 "thrC",
                 "argH","argHA","arcA","NOS1","NOS2","NOS3","NOA1","GATM","ooxA","ooxB","ODH1","nos","noxA","noxB",
                 "proC","ocd","lhpl","lhpl1","prdF","pip","LAP3",
                 "hisD","HIS4","CNDP1","EIF2AK3",
                 "pheA2","pheA","pheC","ADT",
                 "tyrC","cpd","tyrAa","TyrAa1","phhA","TPL",
                 "TRP","trpA","trpB","trpB1","ISS1",
                 "AGXT",
                 "aspQ",
                 "GOT1","GOT2","ASP5","aspB","aspC","yhdR",
                 "glyA",
                 "ItaE",
                 "AGXT2","GGAT",
                 "PAT","hisC","ARO9","pdh",
                 "cysK","ATCYSC1","MET17","CysO","CCBL","metC","patB","CTH","mccB",
                 "TAT","ARO8","tyrB",
                 "BHMT","mmuM","metH","yitJ","metE","msrC","metY","mtnE","GTK",
                 "alaA",
                 "ilvE","LDH","avtA",
                 "lysA","lysAC","LYS1","lysK")
    AA.type <- c(rep("Alanine",4),
                 rep("Aspartate",4),
                 rep("Glutamate",10),
                 rep("Glycine",8),
                 rep("Serine",7),
                 "Threonine",
                 rep("Arginine",9),
                 rep("Proline",6),
                 rep("Histidine",4),
                 rep("Phenylalanine",4),
                 rep("Tyrosine",5),
                 rep("Tryptophan",2),
                 "Alanine_Glycine_Serine",
                 "Aspartate_Glutamate",
                 rep("Aspartate_Phenylalanine_Tyrosine_Cysteine",6),
                 "Glycine_Serine",
                 "Glycine_Threonine",
                 rep("Alanine_Glycine",2),
                 rep("Phenylalanine_Tyrosine",4),
                 rep("Cysteine",9),
                 rep("Phenylalanine_Tyrosine_Methionine",3),
                 rep("Methionine",9),
                 "Alanine_Valine",
                 rep("Valine_Isoleucine_Leucine",3),
                 rep("Lysine",4))
    AA.ko <- c()
    for (i in 1:126) {
        if (AA[i] %in% rownames(ko)) {
            C2 <- ko[AA[i],]
        }else{
            C2 <- rep(0,ncol(ko))
        }
        AA.ko <- rbind(AA.ko,C2)
    }
    rownames(AA.ko) <- AA.name
    AA.ko <- as.data.frame(t(AA.ko[,-ncol(AA.ko)]))
    
    AA.ko$gdhA <- AA.ko$gdhA + AA.ko$gdhA1
    AA.ko <- AA.ko[,-12]
    AA.ko$Sox <- AA.ko$Sox + (AA.ko$soxA + AA.ko$soxB + AA.ko$soxD + AA.ko$soxG)/4
    AA.ko <- AA.ko[,-c(20:23)]
    AA.ko$dgcA <- (AA.ko$dgcA + AA.ko$dgcB)/2
    colnames(AA.ko)[22] <- "dgcAB"
    AA.ko <- AA.ko[,-23]
    AA.ko$etfA <- (AA.ko$etfA + AA.ko$etfB)/2
    colnames(AA.ko)[23] <- "etfAB"
    AA.ko <- AA.ko[,-24]
    AA.ko$nos <- AA.ko$nos + (AA.ko$NOS1 + AA.ko$NOS2 + AA.ko$NOS3)/3
    AA.ko <- AA.ko[,-c(38:40)]
    AA.ko$ooxA <- (AA.ko$ooxA + AA.ko$ooxB)/2
    colnames(AA.ko)[40] <- "ooxAB"
    AA.ko <- AA.ko[,-41]
    AA.ko$noxA <- (AA.ko$noxA + AA.ko$noxB)/2
    colnames(AA.ko)[43] <- "noxAB"
    AA.ko <- AA.ko[,-44]
    AA.ko$lhpl <- AA.ko$lhpl + AA.ko$lhpl1
    AA.ko <- AA.ko[,-47]
    AA.ko$tyrAa <- AA.ko$tyrAa + AA.ko$TyrAa1
    AA.ko <- AA.ko[,-61]
    AA.ko$TRP <- AA.ko$TRP + (AA.ko$trpA + AA.ko$trpB + AA.ko$trpB)/2
    AA.ko <- AA.ko[,-c(64:66)]
    
    AA.ko <- as.data.frame(t(AA.ko))
    AA.ko$Gene <- rownames(AA.ko)
    AA.ko$Type <- AA.type
    AA.ko <- AA.ko[,c("Type","Gene",colnames(AA.ko)[1:(ncol(AA.ko)-2)])]
    AA.ko <- AA.ko[rowSums(AA.ko[,3:ncol(AA.ko)]) > 0,]
    return(AA.ko)
}