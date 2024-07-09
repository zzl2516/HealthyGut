abun.heatmap.s2 <- function(ARG_type,group,Sample_numb){
    ARG_type2 <- ARG_type[which(rowSums(ARG_type[,3:ncol(ARG_type)])>0),]
    ARG_type3 <- melt(ARG_type2)
    colnames(ARG_type2) <- c("Type","Gene","variable","value")
    
    colnames(ARG_type3)[4] <- "Value"
    ARG_type3$Value <- log10(ARG_type3$Value)
    a <- floor(min(ARG_type3$Value[ARG_type3$Value != -Inf]))
    b <- floor(max(ARG_type3$Value[ARG_type3$Value != -Inf]))
    ARG_type3$Value <- ARG_type3$Value - a
    ARG_type3[ARG_type3==-Inf] <- 0
    if (b > 0) {
        i <- seq(a+1,b+1,1) 
    }else{
        i <- seq(a+1,b,1) 
    }
    
    if (b > 0) {
        breaks = c(0:(floor(max(ARG_type3$Value))+1)) 
    }else{
        breaks = c(0:floor(max(ARG_type3$Value))) 
    }
    
    Gene.type <- data.frame(Gene = c("ald","ala","ALT","asdA",
                                     "aspA","ansB","ASRGL1","asnB",
                                     "gudB","GDH2","gdhA","GLT1","gltB","gltD","P5CDH","putA","glsA","purQ",
                                     "Sox","PIPOX","SARDH","dgcAB","etfAB","GNMT","kbl","ItaA",
                                     "serB","thrH","psp","ysaA","SRR","racX","alr",
                                     "thrC",
                                     "argH","argHA","arcA","NOA1","GATM","ooxAB","ODH1","nos","noxAB",
                                     "proC","ocd","lhpl","prdF","pip","LAP3",
                                     "hisD","HIS4","CNDP1","EIF2AK3",
                                     "pheA2","pheA","pheC","ADT",
                                     "tyrC","cpd","tyrAa","phhA","TPL",
                                     "TRP","ISS1",
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
                                     "lysA","lysAC","LYS1","lysK"),
                            Type = c(rep("Alanine",4),
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
                                     rep("Valine",4)))
    
    ARG_type3 <- merge(ARG_type3,Gene.type)
    
    ARG_type_abun <- ggplot(ARG_type3,aes(variable,Gene)) +
        geom_tile(aes(fill = Value),colour = "white") +
        scale_fill_gradientn(name = "Abundance",
                             colours = colorRampPalette(rev(brewer.pal(n = 7,name = "RdYlBu")))(100),
                             breaks = breaks,
                             labels = c("N.D.",10^i)) +
        facet_grid(Type~.,scales = "free",space = "free") + 
        theme_bw()+
        theme(panel.grid=element_blank(),
              axis.ticks.length = unit(0.4,"lines"), 
              axis.ticks = element_line(color='black'),
              axis.line = element_line(colour = "black"), 
              axis.title.x=element_blank(),
              axis.title.y = element_blank(),
              axis.text.y=element_text(colour='black',size=8),
              axis.text.x=element_text(colour = "black",size = 8,
                                       angle = 45,hjust = 1,vjust = 1),
              legend.title = element_text(color = "black",face = "bold",size = 12),
              legend.text = element_text(color = "black",size = 12,face = "bold"),
              strip.text.y = element_text(angle = 0,color = "black",size = 10,
                                          face = "bold")) +
        guides(fill = guide_colorbar(barheight = 5))
    result <- ARG_type_abun
    return(result)
}