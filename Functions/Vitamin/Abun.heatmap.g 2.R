abun.heatmap.g2 <- function(ARG_type,group,Group_numb){
    ARG_type2 <- ARG_type[which(rowSums(ARG_type[,3:ncol(ARG_type)])>0),]
    ARG_type2 <- melt(ARG_type2)
    colnames(ARG_type2) <- c("Type","Gene","variable","value")
    colnames(group) <- c("variable","Group")
    ARG_type3 <- merge(ARG_type2,group)
    
    ARG_type4 <- spread(ARG_type3,Gene,value)
    ARG_type4 <- ARG_type4[,3:ncol(ARG_type4)]
    type.sample <- colnames(ARG_type4[2:ncol(ARG_type4)])
    p.value <- c()
    if (length(levels(ARG_type4$Group)) == 2) {
        for (i in type.sample) {
            fit1 <- wilcox.test(as.formula(sprintf("`%s` ~ Group",i)),
                           data = ARG_type4)
            p.value <- c(p.value,fit1$p.value)
        }
    }else{
        for (i in type.sample) {
            fit1 <- aov(as.formula(sprintf("`%s` ~ Group",i)),
                        data = ARG_type4)
            p.value <- c(p.value,summary(fit1)[[1]][["Pr(>F)"]][[1]])
        }
    }
    
    p.value <- as.data.frame(cbind(type.sample,p.value))
    colnames(p.value) <- c("Gene","P")
    p.value$P <- as.numeric(p.value$P)
    p.value$value <- ifelse(p.value$P < 0.05,1,NA)
    p.value$value2 <- rep(NA,nrow(p.value))
    
    ARG_type3 <- ARG_type3 %>%
        group_by(Gene,Group) %>%
        summarise(Value = mean(value))
    ARG_type6 <- ARG_type3
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
    
    Gene.type <- data.frame(Gene = c("phoA","PHO","PHO5_11_12","ACP1","rsgA","TH2",
                                     "ribE",
                                     "PDXP","PHOSPHO2","P4DH","pno","pdxH","PPOX",
                                     "bioB","BTD","bisC",
                                     "ADH1_7","ADH4","ADH6","adhP","adh","frmA","SDR16C5",
                                     "HSD17B6","RDH8","RDH10","RDH11","RDH12","RDH13","RDH16",
                                     "DHRS3","DHRS4","DHRS4L2","DHRS9"),
                            Type = c(rep("Vitamin B1",6),"Vitamin B2",rep("Vitamin B6",6),
                                     rep("Vitamin H",3),rep("Vitamin A",18)))
    
    ARG_type3 <- merge(ARG_type3,Gene.type)
    ARG_type6 <- merge(ARG_type6,Gene.type)
    p.value <- merge(p.value,Gene.type)
    
    ARG_type_abun <- ggplot(ARG_type3,aes(Group,Gene)) +
        geom_tile(aes(fill = Value),colour = "white") +
        geom_point(data = p.value,aes(x = Group_numb + 1,y = Gene,size = value),
                   color = "#B2182B",show.legend = FALSE) +
        geom_point(data = p.value,aes(x = Group_numb + 1.5,y = Gene),color = "white",
                   show.legend = FALSE) +
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
              axis.text.y=element_text(colour='black',size=10,face = "bold.italic"),
              axis.text.x=element_text(colour = "black",size = 12,face = "bold",
                                       angle = 45,hjust = 1,vjust = 1),
              legend.title = element_text(color = "black",face = "bold",size = 12),
              legend.text = element_text(color = "black",size = 12,face = "bold"),
              strip.text.y = element_text(angle = 0,color = "black",size = 10,
                                          face = "bold")) +
        guides(fill = guide_colorbar(barheight = 10))
    
    ARG_type5 <- spread(ARG_type6,Group,Value)
    for (j in 4:ncol(ARG_type5)) {
        ARG_type5[,j] <- ARG_type5[,j]/ARG_type5[,3]
    }
    ARG_type5 <- melt(ARG_type5)
    colnames(ARG_type5)[4] <- "value2"
    ARG_type5[ARG_type5$variable == levels(group$Group)[1],"value2"] <- 1
    colnames(ARG_type5)[3] <- "Group"
    ARG_type3 <- merge(ARG_type3,ARG_type5)
    p <- ggplot(ARG_type3,aes(Group,Gene)) +
        geom_tile(aes(fill = Value,colour = "white")) +
        geom_point(data = p.value,aes(x = Group_numb + 1,y = Gene,size = value),
                   color = "#B2182B",show.legend = FALSE) +
        geom_point(data = p.value,aes(x = Group_numb + 1.5,y = Gene),color = "white",
                   show.legend = FALSE) +
        scale_fill_gradientn(name = "Abundance",
                             colours = colorRampPalette(brewer.pal(n = 9,name = "YlGn"))(100),
                             breaks = breaks,
                             labels = c("N.D.",10^i)) +
        geom_text(aes(label = sprintf("%0.2f",round(value2,digits = 2))),
                  color = ifelse(ARG_type3$Value > max(ARG_type3$Value/2),"white","black")) +
        facet_grid(Type~.,scales = "free",space = "free") + 
        theme_bw()+
        theme(panel.grid=element_blank(),
              axis.ticks.length = unit(0.4,"lines"), 
              axis.ticks = element_line(color='black'),
              axis.line = element_line(colour = "black"), 
              axis.title.x=element_blank(),
              axis.title.y = element_blank(),
              axis.text.y=element_text(colour='black',size=10,face = "bold.italic"),
              axis.text.x=element_text(colour = "black",size = 12,face = "bold",
                                       angle = 45,hjust = 1,vjust = 1),
              legend.title = element_text(color = "black",face = "bold",size = 10),
              legend.text = element_text(color = "black",size = 10,face = "bold"),
              strip.text.y = element_text(angle = 0,color = "black",size = 10,
                                          face = "bold")) +
        guides(fill = guide_colorbar(barheight = 10),colour = "none")
    
    x <- c("a","b")
    y <- c("a","a")
    test.b <- c()
    bb <- ARG_type4
    aa <- c()
    if (length(levels(ARG_type4$Group)) == 2) {
        for (j in type.sample) {
            fit1 <- wilcox.test(as.formula(sprintf("`%s` ~ Group",j)),
                           data = bb)
            if (fit1$p.value == "NaN") {
                fit1$p.value = 1
            }
            test <- if(fit1$p.value < 0.05){
                x
            }else{
                y
            }
            test.b <- cbind(test.b,test)
            rownames(test.b) <- levels(ARG_type4$Group)
            aa <- c(aa,fit1$p.value)
        }
    }else{
        for (j in type.sample) {
            fit1 <- aov(as.formula(sprintf("`%s` ~ Group",j)),
                        data = bb)
            tuk1<-glht(fit1,linfct=mcp(Group="Tukey"))
            res1 <- cld(tuk1,alpah=0.05)
            test.b <- cbind(test.b,res1$mcletters$Letters)
            aa <- rbind(aa,summary(tuk1)$test$pvalue)
        }
    }
    colnames(test.b) <- type.sample
    test.b <- melt(test.b)
    colnames(test.b) <- c("Group","Gene","value")
    test.b <- merge(test.b,Gene.type)
    
    ARG_type_abun2 <- ggplot(ARG_type3,aes(Group,Gene)) +
        geom_tile(aes(fill = Value),colour = "white") +
        geom_text(data = test.b,aes(Group,Gene,label = value),
                  color = "black",fontface = "bold") +
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
              axis.text.y=element_text(colour='black',size=10,face = "bold.italic"),
              axis.text.x=element_text(colour = "black",size = 12,face = "bold",
                                       angle = 45,hjust = 1,vjust = 1),
              legend.title = element_text(color = "black",face = "bold",size = 12),
              legend.text = element_text(color = "black",size = 12,face = "bold"),
              strip.text.y = element_text(angle = 0,color = "black",size = 10,
                                          face = "bold")) +
        guides(fill = guide_colorbar(barheight = 10))
    
    if (length(levels(group$Group)) == 2) {
        test.result <- data.frame(Index = type.sample,
                                  pvalue = aa)
    }else{
        colnames(aa) <- names(summary(tuk1)$test$coefficients)
        test.result <- as.data.frame(aa)
        test.result$Index <- type.sample
        test.result <- test.result[,c("Index",colnames(test.result)[1:(ncol(test.result)-1)])]
    }
    
    result <- list(p.value,ARG_type_abun,p,ARG_type_abun2,test.result)
    return(result)
}