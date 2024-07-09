SCFA.host <- function(Gene){
    SCFA <- list(poxL = "K00158",ROCK2 = "K17388",PTA = c("K13788","K15024","K00625"),
                 eutD = "K04020",ackA = "K00925",acyP = "K01512",ACSS1_2 = "K01895",
                 ACN1 = "K01913", pct = "K01026",aarC = "K18118",
                 acdAB = c("K01905","K22224","K24012"),ACH1 = "K01067",
                 kor = c("K00174","K00175"),
                 por = c("K00169","K00170","K00171","K00172","K00189","K03737"),
                 pdh = c("K00161","K00162","K00163","K00627"),
                 LDH = "K00016",ASPH = "K00476",poxB = "K00156",
                 AcAldDH = "K00132",bphJ = "K18366",adhE = "K04072",
                 mhpF = "K04073",eutE = "K04021",ALDH = "K00128",aldB = "K00138",
                 ALDH9A1 = "K00149",ALDH7A1 = "K14085",pflD = "K00656",
                 LSC = c("K01899","K01900"),sucCD = c("K01902","K01903"),
                 MUT = "K01847",mcmA = c("K01848","K01849"),epi = "K05606",
                 mmcD = "K11264",G2A = "K08426",mmcC5S = "K03416",
                 ACAT = "K00626",DLAT = "K00627",ACOX1 = "K00232",ACADS = "K00248",
                 acul = "K19745",acrC = "K20143",ACSS3 = "K01908",
                 pduL = "K13923",tdcD = "K00932",pduW = "K19697",
                 HBDH = "K00074",HADH = "K00022",fadN = "K07516",fadB = "K01825",
                 fadJ = "K01782",EHHADH = "K07514",echA = "K01692",HADHA = "K07515",
                 ECHS1 = "K07511",crt = "K01715",ccrA =  "K17829",fabV = "K00209",
                 atoAD = c("K01034","K01035"), ydiF = "K19709",
                 PTB = "K00634",buk = "K00929",L2MO = "K00467")
    temp <- c()
    for (i in 1:63) {
        d <- SCFA[[i]]
        d <- ifelse(d %in% colnames(Gene),d,NA)
        d <- d[complete.cases(d)]
        a <- ifelse(length(d) == 0,0,1)
        temp <- c(temp,a)
    }
    result <- list(temp,SCFA)
    return(result)
}