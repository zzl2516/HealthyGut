abun.heatmap.g <- function(ARG_type,Group_numb){
    ARG_type2 <- ARG_type[which(rowSums(ARG_type[,3:ncol(ARG_type)])>0),]
    ARG_type2 <- melt(ARG_type2)
    colnames(ARG_type2) <- c("Type","Gene","variable","value")
    colnames(group) <- c("variable","Group")
    ARG_type3 <- merge(ARG_type2,group)
    
    ARG_type4 <- spread(ARG_type3,Gene,value)
    ARG_type4 <- ARG_type4[,3:ncol(ARG_type4)]
    type.sample <- colnames(ARG_type4[3:ncol(ARG_type4)])
    
    ARG_type3 <- ARG_type3 %>%
        group_by(Gene,Group) %>%
        summarise(Value = mean(value))
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
    
    Gene.type <- data.frame(Gene = c("pflD","ackA","acyP","ACSS1_2","ACN1",
                                      "pct","aarC","acdAB","ACH1","ALDH","aldB",
                                      "ALDH9A1","ALDH7A1","L2MO","poxB","ACSS3",
                                      "tdcD","pduW","atoAD","ydiF","buk"),
                            Type = c("Formate","Acetate_Propionate","Acetate",
                                     "Acetate_Propionate","Acetate",
                                     "Acetate_Propionate","Acetate",
                                     "Acetate_Propionate","Acetate",
                                     "Acetate","Acetate","Acetate","Acetate",
                                     "Acetate","Acetate","Propionate",
                                     "Propionate","Propionate","Butyrate",
                                     "Butyrate","Butyrate"))
    Gene.type$Type <- factor(Gene.type$Type,levels = c("Formate","Acetate",
                                                       "Propionate","Acetate_Propionate",
                                                       "Butyrate"))
    ARG_type3 <- merge(ARG_type3,Gene.type)
    
    ARG_type5 <- spread(ARG_type3,Group,Value)
    for (j in 4:ncol(ARG_type5)) {
        ARG_type5[,j] <- ARG_type5[,j]/ARG_type5[,3]
    }
    ARG_type5 <- melt(ARG_type5)
    colnames(ARG_type5)[4] <- "value2"
    ARG_type5[ARG_type5$variable == levels(group$Group)[1],"value2"] <- NA
    colnames(ARG_type5)[3] <- "Group"
    ARG_type3 <- merge(ARG_type3,ARG_type5)
    p <- ggplot(ARG_type3,aes(Group,Gene)) +
        geom_tile(aes(fill = Value,colour = "white")) +
        scale_fill_gradientn(name = "Abundance",
                             colours = colorRampPalette(brewer.pal(n = 9,name = "YlGn"))(100),
                             breaks = breaks,
                             labels = c("N.D.",10^i)) +
        geom_text(aes(label = round(value2,2)),
                  color = ifelse(ARG_type3$Value > max(ARG_type3$Value/3),"white","black")) +
        facet_grid(Type~.,scales = "free",space = "free") + 
        theme_bw()+
        theme(panel.grid=element_blank(),
              axis.ticks.length = unit(0.4,"lines"), 
              axis.ticks = element_line(color='black'),
              axis.line = element_line(colour = "black"), 
              axis.title.x=element_blank(),
              axis.title.y = element_blank(),
              axis.text.y=element_text(colour='black',size=10,face = "bold.italic"),
              axis.text.x=element_text(colour = "black",size = 14,face = "bold",
                                       angle = 45,hjust = 1,vjust = 1),
              legend.title = element_text(color = "black",face = "bold",size = 10),
              legend.text = element_text(color = "black",size = 10,face = "bold"),
              strip.text.y = element_text(angle = 0,color = "black",size = 10,
                                          face = "bold")) +
        guides(fill = guide_colorbar(barheight = 10),colour = "none")

    return(p)
}