Immune.synthesis.host <- function(Gene){
    Immune.path <- list(Succinate = c("K00234","K00235","K00236","K00237","K25801","K00239","K00240",
                                      "K00241","K00242","K18859","K18860","K00244","K00245","K00246",
                                      "K00247","K00233","K25995","K25996","K01899","K01900","K01902",
                                      "K01903","K18118","K00135","K08324","K00139","K17761","K15737",
                                      "K05714","K01637","K14471","K14472","K10764","K01739"),
                        IAA = c("K00128","K14085","K00149","K11817","K22417","K01501","K01426",
                                "K21801","K11816"),
                        IAld = c("K04103","K11182","K00274"),
                        Kynurenine = c("K01432","K14263","K07130","K00463","K00486","K00453"),
                        SBA = c("K01442","K22605","K15868","K15870","K15872","K15871","K15874",
                                "K07007"),
                        TMAO = c("K18277","K21579","K20038","K20037","K22443","K22444"))
    temp <- c()
    for (i in 1:6) {
        d <- Immune.path[[i]]
        d <- ifelse(d %in% colnames(Gene),d,NA)
        d <- d[complete.cases(d)]
        a <- ifelse(length(d) == 0,0,1)
        temp <- c(temp,a)
    }
    result <- list(temp,Immune.path)
    return(result)
}