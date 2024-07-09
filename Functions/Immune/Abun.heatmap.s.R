abun.heatmap.s <- function(ARG_type,group,Sample_numb){
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
    
    Gene.type <- data.frame(Gene = c("SDH","LSC","sucCD","aarC","gabD","sad",
                                     "ALDH5A1","SSADH","csiD","mhpC","aceA",
                                     "smtAB","metZ","metB","ALDH","ALDH7A1",
                                     "ALDH9A1","AAO1_2","AAO4","nitrilase",
                                     "amiE","iaaH","YUCCA","ipdC","AOC1","MAO",
                                     "AFMID","BNA7","kynB","IDO","KMO","TDO2",
                                     "cbh","baiA","baiB","baiCD","baiE","baiF",
                                     "baiI","baiN", "tmm","grdH","cutCD","cntAB"),
                            Type = c(rep("Succinate",14),rep("Indole-3-acetic acid",9),
                                     rep("Indole-3-aldehyde",3),rep("Kynurenine",6),
                                     rep("Secondary bile acids",8),rep("Trimethylamine N-oxide",4)))
    
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