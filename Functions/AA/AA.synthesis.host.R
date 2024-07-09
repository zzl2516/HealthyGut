AA.synthesis.host <- function(Gene){
    AA.path <- list(Alanine = c("K00259","K19244","K00814","K09758"),
                    Aspartate = c("K01437","K01424","K013051","K01953"),
                    Glutamate = c("K00260","K15371","K00261","K00262","K00264","K00265","K00266","K00294","K13821","K01425","K23265"),
                    Glycine = c("K00301","K00302","K00303","K00304","K00305","K00306","K00314","K21833","K21834","K25960","K25961","K00552","K00639","K20801"),
                    Serine = c("K01079","K02203","K22305","K22528","K12235","K25316","K25317"),
                    Threonine = "K01733",
                    Arginine = c("K01755","K14681","K01478","K13240","K13241","K13242","K13247","K00613","K24954","K24955","K23233","K00296","K24952","K24953"),
                    Proline = c("K00286","K01750","K19743","K19742","K01777","K01259","K11142"),
                    Histidine = c("K00013","K14152","K05604","K08660"),
                    Phenylalanine = c("K04518","K14170","K01713","K05359"),
                    Tyrosine = c("K00220","K24018","K15226","K15227","K00500","K01688"),
                    Tryptophan = c("K01694","K01695","K01696","K06001","K00837"),
                    Alanine_Glycine_Serine = "K00830",
                    Aspartate_Glutamate = "K05597",
                    Aspartate_Phenylalanine_Tyrosine_Cysteine = c("K14454","K14455","K00811","K00812","K00813","K11358"),
                    Glycine_Serine = "K00600",
                    Glycine_Threonine = "K01620",
                    Alanine_Glycine = c("K00827","K14272"),
                    Phenylalanine_Tyrosine = c("K15849","K00817","K05821","K00270"),
                    Cysteine = c("K01738","K13034","K17069","K10150","K00816","K01760","K14155","K01758","K17217"),
                    Phenylalanine_Tyrosine_Methionine = c("K00815","K00838","K00832"),
                    Methionine = c("K00544","K00547","K00548","K24024","K00549","K08968","K01740","K08969","K23977"),
                    Alanine_Valine = "K14260",
                    Valine_Isoleucine_Leucine = c("K00826","K00263","K00835"),
                    Lysine = c("K01586","K12526","K00290","K05831"))
    temp <- c()
    for (i in 1:25) {
        d <- AA.path[[i]]
        d <- ifelse(d %in% colnames(Gene),d,NA)
        d <- d[complete.cases(d)]
        a <- ifelse(length(d) == 0,0,1)
        temp <- c(temp,a)
    }
    result <- list(temp,AA.path)
    return(result)
}