# Pipeline for AA synthesis analysis based on metagenome
# All codes are based on R 4.0.2
library(ggplot2)
library(ggimage)
library(reshape2)
library(RColorBrewer)
library(multcomp)
library(dplyr)
library(tidyverse)
library(viridis)
library(cowplot)
library(grid)
library(tidyr)
library(pheatmap)
library(vegan)
library(ape)
library(ggrepel)
library(FactoMineR)

## Input
source("Functions/Dataprocess/Loading.R")
cbbPalette <- Loading_Color()
group <- Loading_Group("Input/group.txt")
Group_numb <- length(unique(group[,2]))
Sample_numb <- length(unique(group[,1]))
ko <- Loading_Table("Input/kegg.profile.entry.xls")
Gene <- Loading_Gene("Input/uniqGeneSet.kegg.category.xls")
tax <- Loading_tax("Input/uniqGeneSet.NR.anno.xls")
abundance <- Loading_abun("Input/profile.txt")
path <- Loading_Table("Input/kegg.profile.pathway.xls")
abundance <- abundance[,c("V1",colnames(abundance)[colnames(abundance) %in% group$ID])]
ko <- ko[,c(colnames(ko)[colnames(ko) %in% group$ID],"Description")]
path <- path[,c(colnames(path)[colnames(path) %in% group$ID],"Description")]

dir.create("Results")

# Amino acid metabolism
dir.create("Results/AA")

## Gene identification
dir.create("Results/AA/Gene")

## Pathways of amino acids
dir.create("Results/AA/Gene/AA_metabolism")
source("Functions/AA/AA.metabolism.R")
AA.path <- AA.metabolism(path)
write.table(AA.path,"Results/AA/Gene/AA_metabolism/AA_metabolism_pathway.txt",
            sep = "\t",row.names = FALSE)

### Group heatmap for amino acid metabolism pathways
source("Functions/AA/Abun.heatmap.g.R")
result <- abun.heatmap.g(AA.path,group,Group_numb)
pdf(file = "Results/AA/Gene/AA_metabolism/AA_metabolism_pathway_group.pdf",
    width = 4.1 + 0.4*(Group_numb + 1),
    height = max(str_length(levels(group$Group)))/7 + nrow(result[[1]])*0.2)
result[[2]]
dev.off()

pdf(file = "Results/AA/Gene/AA_metabolism/AA_metabolism_pathway_FC.pdf",
    width = 4.85 + 0.5*(Group_numb + 1),
    height = max(str_length(levels(group$Group)))/7 + nrow(result[[1]])*0.2)
result[[3]]
dev.off()

pdf(file = "Results/AA/Gene/AA_metabolism/AA_metabolism_pathway_test.pdf",
    width = 4.1 + 0.4*Group_numb,
    height = max(str_length(levels(group$Group)))/7 + nrow(result[[1]])*0.2)
result[[4]]
dev.off()

write.table(result[[5]],
            "Results/AA/Gene/AA_metabolism/AA_metabolism_pathway_test.txt",
            sep = "\t",row.names = FALSE)

### Sample heatmap for amino acid metabolism pathways
source("Functions/AA/Abun.heatmap.s.R")
result <- abun.heatmap.s(AA.path,group,Sample_numb)
pdf(file = "Results/AA/Gene/AA_metabolism/AA_metabolism_pathway_sample.pdf",
    width = 4.1 + 0.3*Sample_numb,
    height = max(str_length(group$ID))/7 + nrow(AA.path)*0.2)
result
dev.off()

## Gene abundance in amino acid synthesis
dir.create("Results/AA/Gene/AA_synthesis")
source("Functions/AA/AA.synthesis.R")
source("Functions/AA/Abun.heatmap.g 2.R")
source("Functions/AA/Abun.heatmap.s 2.R")
AA.ko <- AA.synthesis(ko)
write.table(AA.ko,"Results/AA/Gene/AA_synthesis/AA_synthesis_gene.txt",
            sep = "\t",row.names = FALSE)

result <- abun.heatmap.g2(AA.ko,group,Group_numb)
pdf(file = "Results/AA/Gene/AA_synthesis/AA_synthesis_group.pdf",
    width = 5 + 0.4*(Group_numb + 1),
    height = max(str_length(levels(group$Group)))/7 + nrow(result[[1]])*0.2)
result[[2]]
dev.off()

pdf(file = "Results/AA/Gene/AA_synthesis/AA_synthesis_FC.pdf",
    width = 5 + 0.5*(Group_numb + 1),
    height = max(str_length(levels(group$Group)))/7 + nrow(result[[1]])*0.2)
result[[3]]
dev.off()

pdf(file = "Results/AA/Gene/AA_synthesis/AA_synthesis_test.pdf",
    width = 5 + 0.4*Group_numb,
    height = max(str_length(levels(group$Group)))/7 + nrow(result[[1]])*0.2)
result[[4]]
dev.off()

write.table(result[[5]],
            "Results/AA/Gene/AA_synthesis/AA_synthesis_test.txt",
            sep = "\t",row.names = FALSE)

result <- abun.heatmap.s2(AA.ko,group,Sample_numb)
pdf(file = "Results/AA/Gene/AA_synthesis/AA_synthesis_sample.pdf",
    width = 5.1 + 0.2*Sample_numb,
    height = max(str_length(group$ID))/7 + nrow(AA.ko)*0.2)
result
dev.off()

### Host_relative stat
tax.name <- c("phylum","class","order","family","genus","species")
source("Functions/AA/AA.synthesis.host.R")
AA.h <- AA.synthesis.host(Gene)

#### relative_Group
dir.create("Results/AA/Host_relative_Group")
source("Functions/AA/Host.realtive.group.R")
for (i in 1:25) {
    if (names(AA.h[[2]])[i] %in% AA.ko$Type & AA.h[[1]][[i]] == 1) {
        dir.create(paste("Results/AA/Host_relative_Group/",names(AA.h[[2]])[i],sep = ""))
        result <- Host.relative.group(Gene,tax,abundance,group,AA.h[[2]][[i]])
        title <- names(AA.h[[2]][i])
        for (j in 1:6) {
            write.csv(result[[j]],
                      paste("Results/AA/Host_relative_Group/",names(AA.h[[2]])[i],
                            "/",names(AA.h[[2]])[i],"_",tax.name[j],".csv",sep = ""))
            aa <- ifelse(nrow(result[[j]]) > 5,
                         max(str_length(rownames(result[[j]])[1:5])),
                         max(str_length(rownames(result[[j]]))))
            if (nrow(result[[j]]) > 6) {
                result[[j]] <- result[[j]][1:5,]
            }
            pdf(paste("Results/AA/Host_relative_Group/",names(AA.h[[2]])[i],
                     "/",names(AA.h[[2]])[i],"_",tax.name[j],".pdf",sep = ""),
                width = ifelse(aa > 40,aa*0.1 + 0.5*Group_numb,3.2 + 0.5*Group_numb),
                height = 3.5)
            pheatmap(result[[j]],fontsize = 10,cluster_rows = FALSE,fontface = "bold",cluster_cols = FALSE,
                     cellwidth = 30,cellheight = 20,legend = FALSE,breaks = c(seq(0,1,by = 0.01)),
                     color = colorRampPalette(c("white","Red"))(100),border_color = "black")
            grid.text(title,hjust = 0.5,x = 0.5,y = 0.92,gp = gpar(font = 2,size = 1.2))
            dev.off()
        }
    }
}

#### Legend
df <- data.frame(A = c(0,0.1,0.2,0.3,0.4,0.5),
                 B = c(0.5,0.6,0.7,0.8,0.9,1))
df <- t(df)
df1 <- data.frame(A = c("0%","10%","20%","30%","40%","50%"),
                  B = c("50%","60%","70%","80%","90%","100%"))
df1 <- t(df1)
bk <- c(seq(0,1,by = 0.01))
pheatmap(df,fontsize = 30,cluster_rows = FALSE,fontface = "bold",cluster_cols = FALSE,
         cellwidth = 80,cellheight = 50,legend = FALSE,breaks = bk,show_rownames = FALSE,
         color = colorRampPalette(c("white","Red"))(100),show_colnames = FALSE,
         display_numbers = df1,number_color = "black",border_color = "black",
         filename = "Results/AA/Host_relative_Group/Legend_relative_Group.pdf",width = 7,height = 2)

### AA host stat absolute_Group
dir.create("Results/AA/Host_absolute_Group")
source("Functions/AA/Host.absolute.group.R")
for (i in 1:25) {
    if (names(AA.h[[2]])[i] %in% AA.ko$Type & AA.h[[1]][[i]] == 1) {
        dir.create(paste("Results/AA/Host_absolute_Group/",names(AA.h[[2]])[i],sep = ""))
        result <- Host.absolute.group(Gene,tax,abundance,group,AA.h[[2]][[i]])
        title <- names(AA.h[[2]][i])
        for (j in 1:6) {
            write.csv(result[[j]],
                      paste("Results/AA/Host_absolute_Group/",names(AA.h[[2]])[i],
                            "/",names(AA.h[[2]])[i],"_",tax.name[j],".csv",sep = ""))
            aa <- ifelse(nrow(result[[j]]) > 5,
                         max(str_length(rownames(result[[j]])[1:5])),
                         max(str_length(rownames(result[[j]]))))
            if (nrow(result[[j]]) > 6) {
                result[[j]] <- result[[j]][1:5,]
            }
            pdf(paste("Results/AA/Host_absolute_Group/",names(AA.h[[2]])[i],
                      "/",names(AA.h[[2]])[i],"_",tax.name[j],".pdf",sep = ""),
                width = ifelse(aa > 40,aa*0.1 + 0.5*Group_numb,3.2 + 0.5*Group_numb),
                height = 3.5)
            pheatmap(result[[j]],fontsize = 10,cluster_rows = FALSE,fontface = "bold",cluster_cols = FALSE,
                     cellwidth = 30,cellheight = 20,
                     color = colorRampPalette(c("white","Red"))(100),border_color = "black")
            grid.text(title,hjust = 0.5,x = 0.5,y = 0.92,gp = gpar(font = 2,size = 1.2))
            dev.off()
        }
    }
}

#### relative_Ungroup
dir.create("Results/AA/Host_relative_Ungroup")
source("Functions/AA/Host.relative.R")
for (i in 1:25) {
    if (names(AA.h[[2]])[i] %in% AA.ko$Type & AA.h[[1]][[i]] == 1) {
        dir.create(paste("Results/AA/Host_relative_Ungroup/",names(AA.h[[2]])[i],sep = ""))
        result <- Host.relative(Gene,tax,abundance,AA.h[[2]][[i]])
        title <- names(AA.h[[2]][i])
        for (j in 1:6) {
            write.csv(result[[j]],
                      paste("Results/AA/Host_relative_Ungroup/",names(AA.h[[2]])[i],
                            "/",names(AA.h[[2]])[i],"_",tax.name[j],".csv",sep = ""))
            aa <- ifelse(nrow(result[[j]]) > 5,
                         max(str_length(rownames(result[[j]])[1:5])),
                         max(str_length(rownames(result[[j]]))))
            if (nrow(result[[j]]) > 6) {
                result[[j]] <- result[[j]][1:5,]
            }
            pdf(paste("Results/AA/Host_relative_Ungroup/",names(AA.h[[2]])[i],
                      "/",names(AA.h[[2]])[i],"_",tax.name[j],".pdf",sep = ""),
                width = ifelse(aa > 40,aa*0.1 + 0.6*Sample_numb,3.2 + 0.6*Sample_numb),
                height = 3.5)
            pheatmap(result[[j]],fontsize = 10,cluster_rows = FALSE,fontface = "bold",cluster_cols = FALSE,
                     cellwidth = 30,cellheight = 20,legend = FALSE,breaks = c(seq(0,1,by = 0.01)),
                     color = colorRampPalette(c("white","Red"))(100),border_color = "black")
            grid.text(title,hjust = 0.5,x = 0.5,y = 0.92,gp = gpar(font = 2,size = 1.2))
            dev.off()
        }
    }
}

#### Legend
df <- data.frame(A = c(0,0.1,0.2,0.3,0.4,0.5),
                 B = c(0.5,0.6,0.7,0.8,0.9,1))
df <- t(df)
df1 <- data.frame(A = c("0%","10%","20%","30%","40%","50%"),
                  B = c("50%","60%","70%","80%","90%","100%"))
df1 <- t(df1)
bk <- c(seq(0,1,by = 0.01))
pheatmap(df,fontsize = 30,cluster_rows = FALSE,fontface = "bold",cluster_cols = FALSE,
         cellwidth = 80,cellheight = 50,legend = FALSE,breaks = bk,show_rownames = FALSE,
         color = colorRampPalette(c("white","Red"))(100),show_colnames = FALSE,
         display_numbers = df1,number_color = "black",border_color = "black",
         filename = "Results/AA/Host_relative_Ungroup/Legend_relative_Ungroup.pdf",width = 7,height = 2)

### AA host stat absolute_Ungroup
dir.create("Results/AA/Host_absolute_Ungroup")
source("Functions/AA/Host.absolute.R")
for (i in 1:25) {
    if (names(AA.h[[2]])[i] %in% AA.ko$Type & AA.h[[1]][[i]] == 1) {
        dir.create(paste("Results/AA/Host_absolute_Ungroup/",names(AA.h[[2]])[i],sep = ""))
        result <- Host.absolute(Gene,tax,abundance,AA.h[[2]][[i]])
        title <- names(AA.h[[2]][i])
        for (j in 1:6) {
            write.csv(result[[j]],
                      paste("Results/AA/Host_absolute_Ungroup/",names(AA.h[[2]])[i],
                            "/",names(AA.h[[2]])[i],"_",tax.name[j],".csv",sep = ""))
            aa <- ifelse(nrow(result[[j]]) > 5,
                         max(str_length(rownames(result[[j]])[1:5])),
                         max(str_length(rownames(result[[j]]))))
            if (nrow(result[[j]]) > 6) {
                result[[j]] <- result[[j]][1:5,]
            }
            pdf(paste("Results/AA/Host_absolute_Ungroup/",names(AA.h[[2]])[i],
                      "/",names(AA.h[[2]])[i],"_",tax.name[j],".pdf",sep = ""),
                width = ifelse(aa > 40,aa*0.1 + 0.6*Sample_numb,3.5 + 0.6*Sample_numb),
                height = 3.5)
            pheatmap(result[[j]],fontsize = 10,cluster_rows = FALSE,fontface = "bold",cluster_cols = FALSE,
                     cellwidth = 30,cellheight = 20,
                     color = colorRampPalette(c("white","Red"))(100),border_color = "black")
            grid.text(title,hjust = 0.5,x = 0.5,y = 0.92,gp = gpar(font = 2,size = 1.2))
            dev.off()
        }
    }
}

## SCFA
dir.create("Results/SCFA")

## Gene identification
dir.create("Results/SCFA/Gene")

## Gene of SCFA synthesis key gene
dir.create("Results/SCFA/Gene/Synthesis")
source("Functions/SCFA/SCFA.synthesis.R")
SCFA.ko <- SCFA.synthesis(ko)
write.table(SCFA.ko,"Results/SCFA/Gene/Synthesis/SCFA_key_synthesis_gene.txt",
            sep = "\t",row.names = FALSE)

### Group heatmap for SCFA synthesis key gene
source("Functions/SCFA/Abun.heatmap.g.R")
result <- abun.heatmap.g(SCFA.ko,group,Group_numb)
pdf(file = "Results/SCFA/Gene/Synthesis/SCFA_key_synthesis_gene_group.pdf",
    width = 3.4 + 0.4*(Group_numb + 1),
    height = max(str_length(levels(group$Group)))/7 + nrow(result[[1]])*0.22)
result[[2]]
dev.off()

pdf(file = "Results/SCFA/Gene/Synthesis/SCFA_key_synthesis_gene_FC.pdf",
    width = 3.4 + 0.4*(Group_numb + 1),
    height = max(str_length(levels(group$Group)))/7 + nrow(result[[1]])*0.22)
result[[3]]
dev.off()

pdf(file = "Results/SCFA/Gene/Synthesis/SCFA_key_synthesis_gene_test.pdf",
    width = 3.4 + 0.4*Group_numb,
    height = max(str_length(levels(group$Group)))/7 + nrow(result[[1]])*0.22)
result[[4]]
dev.off()

write.table(result[[5]],
            "Results/SCFA/Gene/Synthesis/SCFA_key_synthesis_gene_test.txt",
            sep = "\t",row.names = FALSE)

### Sample heatmap for SCFA synthesis key gene
source("Functions/SCFA/Abun.heatmap.s.R")
result <- abun.heatmap.s(SCFA.ko,group,Sample_numb)
pdf(file = "Results/SCFA/Gene/Synthesis/SCFA_key_synthesis_gene_sample.pdf",
    width = 3.4 + 0.3*Sample_numb,
    height = max(str_length(group$ID))/7 + nrow(SCFA.ko)*0.2)
result
dev.off()

## Gene abundance in synthesis pathway
### Acetate
dir.create("Results/SCFA/Gene/Acetate")
source("Functions/SCFA/Acetate.pathway.R")
Acetate.ko <- Acetate.pathway(ko)
write.table(Acetate.ko,"Results/SCFA/Gene/Acetate/Acetate_pathway.txt",
            sep = "\t",row.names = FALSE)

source("Functions/SCFA/Fold change.R")
result <- fold.change(Acetate.ko,group)
write.table(result[[1]],"Results/SCFA/Gene/Acetate/Gene_fold_change.txt",sep = "\t",
            row.names = FALSE)
pdf(file = "Results/SCFA/Gene/Acetate/Gene_fold_change.pdf",
    width = 0.5*(Group_numb + 1) + 0.8,
    height = max(str_length(levels(group$Group)))/7 + nrow(Acetate.ko)*0.2)
result[[3]]
dev.off()

Acetate.img <- "Functions/SCFA/Acetate.pdf"
A1 <- result[[2]] + ylim(Acetate.ko[Acetate.ko$Gene %in% c("pflD","kor","por","pdh"),1])
A2 <- result[[2]] + ylim(Acetate.ko[Acetate.ko$Gene %in% c("poxL"),1])
A3 <- result[[2]] + ylim(Acetate.ko[Acetate.ko$Gene %in% c("ROCK2","PTA","eutD"),1])
A4 <- result[[2]] + ylim(Acetate.ko[Acetate.ko$Gene %in% c("ackA","acyP"),1])
A5 <- result[[2]] + ylim(Acetate.ko[Acetate.ko$Gene %in% c("pct","aarC","acdAB","ACH1","ACSS1_2","ACN1"),1])
A6 <- result[[2]] + ylim(Acetate.ko[Acetate.ko$Gene %in% c("AcAldDH","bphJ","adhE","mhpF","eutE"),1])
A7 <- result[[2]] + ylim(Acetate.ko[Acetate.ko$Gene %in% c("ALDH","aldB","ALDH9A1","ALDH7A1"),1])
A8 <- result[[2]] + ylim(Acetate.ko[Acetate.ko$Gene %in% c("LDH"),1])
A9 <- result[[2]] + ylim(Acetate.ko[Acetate.ko$Gene %in% c("L2MO"),1])
A10 <- result[[2]] + ylim(Acetate.ko[Acetate.ko$Gene %in% c("poxB"),1])

pdf("Results/SCFA/Gene/Acetate/Synthesis_pathway.pdf",width = 13,height = 7)
gg <- ggplot()
ggbackground(gg,Acetate.img)
if (nrow(Acetate.ko[Acetate.ko$Gene %in% c("pflD","kor","por","pdh"),]) > 0) {
    print(A1,vp = viewport(width = 0.025*(Group_numb + 2),
                           height = 0.025*(nrow(Acetate.ko[Acetate.ko$Gene %in% c("pflD","kor","por","pdh"),]) + 1),
                           x = 0.274,y = 0.53,just = "top"))
}
if (nrow(Acetate.ko[Acetate.ko$Gene %in% c("poxL"),]) > 0) {
    print(A2,vp = viewport(width = 0.025*(Group_numb + 2),
                           height = 0.025*(nrow(Acetate.ko[Acetate.ko$Gene %in% c("poxL"),]) + 1),
                           x = 0.26,y = 0.62,just = "right"))
}
if (nrow(Acetate.ko[Acetate.ko$Gene %in% c("ROCK2","PTA","eutD"),]) > 0) {
    print(A3,vp = viewport(width = 0.025*(Group_numb + 2),
                           height = 0.025*(nrow(Acetate.ko[Acetate.ko$Gene %in% c("ROCK2","PTA","eutD"),]) + 1),
                           x = 0.29,y = 0.665,just = "left"))
}
if (nrow(Acetate.ko[Acetate.ko$Gene %in% c("ackA","acyP"),]) > 0) {
    print(A4,vp = viewport(width = 0.025*(Group_numb + 2),
                           height = 0.025*(nrow(Acetate.ko[Acetate.ko$Gene %in% c("ackA","acyP"),]) + 1),
                           x = 0.805,y = 0.833,just = "right"))
}
if (nrow(Acetate.ko[Acetate.ko$Gene %in% c("pct","aarC","acdAB","ACH1","ACSS1_2","ACN1"),]) > 0) {
    print(A5,vp = viewport(width = 0.028*(Group_numb + 2),
                           height = 0.025*(nrow(Acetate.ko[Acetate.ko$Gene %in% c("pct","aarC","acdAB","ACH1","ACSS1_2","ACN1"),]) + 1),
                           x = 0.5745,y = 0.5255))
}
if (nrow(Acetate.ko[Acetate.ko$Gene %in% c("AcAldDH","bphJ","adhE","mhpF","eutE"),]) > 0) {
    print(A6,vp = viewport(width = 0.028*(Group_numb + 2),
                           height = 0.025*(nrow(Acetate.ko[Acetate.ko$Gene %in% c("AcAldDH","bphJ","adhE","mhpF","eutE"),]) + 1),
                           x = 0.48,y = 0.38,just = "right"))
}
if (nrow(Acetate.ko[Acetate.ko$Gene %in% c("ALDH","aldB","ALDH9A1","ALDH7A1"),]) > 0) {
    print(A7,vp = viewport(width = 0.025*(Group_numb + 2),
                           height = 0.025*(nrow(Acetate.ko[Acetate.ko$Gene %in% c("ALDH","aldB","ALDH9A1","ALDH7A1"),]) + 1),
                           x = 0.675,y = 0.395,just = "left"))
}
if (nrow(Acetate.ko[Acetate.ko$Gene %in% c("LDH"),]) > 0) {
    print(A8,vp = viewport(width = 0.025*(Group_numb + 2),
                           height = 0.025*(nrow(Acetate.ko[Acetate.ko$Gene %in% c("LDH"),]) + 1),
                           x = 0.446,y = 0.65,just = "left"))
}
if (nrow(Acetate.ko[Acetate.ko$Gene %in% c("L2MO"),]) > 0) {
    print(A9,vp = viewport(width = 0.025*(Group_numb + 2),
                           height = 0.025*(nrow(Acetate.ko[Acetate.ko$Gene %in% c("L2MO"),]) + 1),
                           x = 0.712,y = 0.65,just = "right"))
}
if (nrow(Acetate.ko[Acetate.ko$Gene %in% c("poxB"),]) > 0) {
    print(A10,vp = viewport(width = 0.025*(Group_numb + 2),
                            height = 0.025*(nrow(Acetate.ko[Acetate.ko$Gene %in% c("poxB"),]) + 1),
                            x = 0.25,y = 0.3,just = "right"))
}
dev.off()

### Propionate
dir.create("Results/SCFA/Gene/Propionate")
source("Functions/SCFA/Propionate.pathway.R")
Propionate.ko <- Propionate.pathway(ko)
write.table(Propionate.ko,"Results/SCFA/Gene/Propionate/Propionate_pathway.txt",
            sep = "\t",row.names = FALSE)

source("Functions/SCFA/Fold change.R")
result <- fold.change(Propionate.ko,group)
write.table(result[[1]],"Results/SCFA/Gene/Propionate/Gene_fold_change.txt",sep = "\t",
            row.names = FALSE)
pdf(file = "Results/SCFA/Gene/Propionate/Gene_fold_change.pdf",
    width = 0.5*(Group_numb + 1) + 0.8,
    height = max(str_length(levels(group$Group)))/7 + nrow(Propionate.ko)*0.2)
result[[3]]
dev.off()

Propionate.img <- "Functions/SCFA/Propionate.pdf"
P1 <- result[[2]] + ylim(Propionate.ko[Propionate.ko$Gene %in% c("LCS","sucCD"),1])
P2 <- result[[2]] + ylim(Propionate.ko[Propionate.ko$Gene %in% c("MUT","mcmA"),1])
P3 <- result[[2]] + ylim(Propionate.ko[Propionate.ko$Gene %in% c("epi"),1])
P4 <- result[[2]] + ylim(Propionate.ko[Propionate.ko$Gene %in% c("mmcD","G2A","mmcC5S"),1])
P5 <- result[[2]] + ylim(Propionate.ko[Propionate.ko$Gene %in% c("pct"),1])
P6 <- result[[2]] + ylim(Propionate.ko[Propionate.ko$Gene %in% c("ACAT","DLAT"),1])
P7 <- result[[2]] + ylim(Propionate.ko[Propionate.ko$Gene %in% c("ACOX1","ACADS","acul","acrC"),1])
P8 <- result[[2]] + ylim(Propionate.ko[Propionate.ko$Gene %in% c("pct","acdAB","ACSS1_2","ACSS3"),1])
P9 <- result[[2]] + ylim(Propionate.ko[Propionate.ko$Gene %in% c("PTA","pduL"),1])
P10 <- result[[2]] + ylim(Propionate.ko[Propionate.ko$Gene %in% c("ackA","tdcD","pduW"),1])

pdf("Results/SCFA/Gene/Propionate/Synthesis_pathway.pdf",width = 13,height = 7)
gg <- ggplot()
ggbackground(gg,Propionate.img)
if (nrow(Propionate.ko[Propionate.ko$Gene %in% c("LCS","sucCD"),])  > 0) {
    print(P1,vp = viewport(width = 0.028*(Group_numb + 2),
                           height = 0.028*(nrow(Propionate.ko[Propionate.ko$Gene %in% c("LCS","sucCD"),]) + 1),
                           x = 0.3,y = 0.274,just = "right"))
}
if (nrow(Propionate.ko[Propionate.ko$Gene %in% c("MUT","mcmA"),])  > 0) {
    print(P2,vp = viewport(width = 0.028*(Group_numb + 2),
                           height = 0.028*(nrow(Propionate.ko[Propionate.ko$Gene %in% c("MUT","mcmA"),]) + 1),
                           x = 0.3,y = 0.474,just = "right"))
}
if (nrow(Propionate.ko[Propionate.ko$Gene %in% c("epi"),]) > 0) {
    print(P3,vp = viewport(width = 0.025*(Group_numb + 2),
                           height = 0.028*(nrow(Propionate.ko[Propionate.ko$Gene %in% c("epi"),]) + 1),
                           x = 0.3,y = 0.676,just = "right"))
}
if (nrow(Propionate.ko[Propionate.ko$Gene %in% c("mmcD","G2A","mmcC5S"),]) > 0) {
    print(P4,vp = viewport(width = 0.03*(Group_numb + 2),
                           height = 0.028*(nrow(Propionate.ko[Propionate.ko$Gene %in% c("mmcD","G2A","mmcC5S"),]) + 1),
                           x = 0.41,y = 0.834,just = "bottom"))
}
if (nrow(Propionate.ko[Propionate.ko$Gene %in% c("pct"),]) > 0) {
    print(P5,vp = viewport(width = 0.025*(Group_numb + 2),
                           height = 0.028*(nrow(Propionate.ko[Propionate.ko$Gene %in% c("pct"),]) + 1),
                           x = 0.52,y = 0.274,just = "left"))
}
if (nrow(Propionate.ko[Propionate.ko$Gene %in% c("ACAT","DLAT"),]) > 0) {
    print(P6,vp = viewport(width = 0.028*(Group_numb + 2),
                           height = 0.028*(nrow(Propionate.ko[Propionate.ko$Gene %in% c("ACAT","DLAT"),]) + 1),
                           x = 0.52,y = 0.474,just = "left"))
}
if (nrow(Propionate.ko[Propionate.ko$Gene %in% c("ACOX1","ACADS","acul","acrC"),]) > 0) {
    print(P7,vp = viewport(width = 0.028*(Group_numb + 2),
                           height = 0.025*(nrow(Propionate.ko[Propionate.ko$Gene %in% c("ACOX1","ACADS","acul","acrC"),]) + 1),
                           x = 0.51,y = 0.676,just = "right"))
}
if (nrow(Propionate.ko[Propionate.ko$Gene %in% c("pct","acdAB","ACSS1_2","ACSS3"),]) > 0) {
    print(P8,vp = viewport(width = 0.03*(Group_numb + 2),
                           height = 0.025*(nrow(Propionate.ko[Propionate.ko$Gene %in% c("pct","acdAB","ACSS1_2","ACSS3"),]) + 1),
                           x = 0.6,y = 0.834,just = "bottom"))
}
if (nrow(Propionate.ko[Propionate.ko$Gene %in% c("PTA","pduL"),]) > 0) {
    print(P9,vp = viewport(width = 0.028*(Group_numb + 2),
                           height = 0.025*(nrow(Propionate.ko[Propionate.ko$Gene %in% c("PTA","pduL"),]) + 1),
                           x = 0.585,y = 0.676,just = "left"))
}
if (nrow(Propionate.ko[Propionate.ko$Gene %in% c("ackA","tdcD","pduW"),]) > 0) {
    print(P10,vp = viewport(width = 0.028*(Group_numb + 2),
                            height = 0.025*(nrow(Propionate.ko[Propionate.ko$Gene %in% c("ackA","tdcD","pduW"),]) + 1),
                            x = 0.715,y = 0.676,just = "left"))
}
dev.off()

### Butyrate
dir.create("Results/SCFA/Gene/Butyrate")
source("Functions/SCFA/Butyrate.pathway.R")
Butyrate.ko <- Butyrate.pathway(ko)
write.table(Butyrate.ko,"Results/SCFA/Gene/Butyrate/Butyrate_pathway.txt",
            sep = "\t",row.names = FALSE)

source("Functions/SCFA/Fold change.R")
result <- fold.change(Butyrate.ko,group)
write.table(result[[1]],"Results/SCFA/Gene/Butyrate/Gene_fold_change.txt",sep = "\t",
            row.names = FALSE)
pdf(file = "Results/SCFA/Gene/Butyrate/Gene_fold_change.pdf",
    width = 0.5*(Group_numb + 1) + 0.8,
    height = max(str_length(levels(group$Group)))/7 + nrow(Butyrate.ko)*0.2)
result[[3]]
dev.off()

Butyrate.img <- "Functions/SCFA/Butyrate.pdf"
B1 <- result[[2]] + ylim(Butyrate.ko[Butyrate.ko$Gene %in% c("ACAT"),1])
B2 <- result[[2]] + ylim(Butyrate.ko[Butyrate.ko$Gene %in% c("fadB","fadJ","EHHADH"),1])
B3 <- result[[2]] + ylim(Butyrate.ko[Butyrate.ko$Gene %in% c("HBDH","HADH","fadN"),1])
B4 <- result[[2]] + ylim(Butyrate.ko[Butyrate.ko$Gene %in% c("echA","HADHA","ECHS1","crt"),1])
B5 <- result[[2]] + ylim(Butyrate.ko[Butyrate.ko$Gene %in% c("ccrA","fabV"),1])
B6 <- result[[2]] + ylim(Butyrate.ko[Butyrate.ko$Gene %in% c("atoAD","ydiF"),1])
B7 <- result[[2]] + ylim(Butyrate.ko[Butyrate.ko$Gene %in% c("PTB"),1])
B8 <- result[[2]] + ylim(Butyrate.ko[Butyrate.ko$Gene %in% c("buk"),1])

pdf("Results/SCFA/Gene/Butyrate/Synthesis_pathway.pdf",width = 13,height = 7)
gg <- ggplot()
ggbackground(gg,Butyrate.img)
if (nrow(Butyrate.ko[Butyrate.ko$Gene %in% c("ACAT"),]) > 0) {
    print(B1,vp = viewport(width = 0.025*(Group_numb + 2),
                           height = 0.025*(nrow(Butyrate.ko[Butyrate.ko$Gene %in% c("ACAT"),]) + 1),
                           x = 0.41,y = 0.92,just = "bottom"))
}
if (nrow(Butyrate.ko[Butyrate.ko$Gene %in% c("fadB","fadJ","EHHADH"),]) > 0) {
    print(B2,vp = viewport(width = 0.025*(Group_numb + 2),
                           height = 0.025*(nrow(Butyrate.ko[Butyrate.ko$Gene %in% c("fadB","fadJ","EHHADH"),]) + 1),
                           x = 0.55,y = 0.71,just = "left"))
}
if (nrow(Butyrate.ko[Butyrate.ko$Gene %in% c("HBDH","HADH","fadN"),]) > 0) {
    print(B3,vp = viewport(width = 0.028*(Group_numb + 2),
                           height = 0.025*(nrow(Butyrate.ko[Butyrate.ko$Gene %in% c("HBDH","HADH","fadN"),]) + 1),
                           x = 0.49,y = 0.78,just = "right"))
}
if (nrow(Butyrate.ko[Butyrate.ko$Gene %in% c("echA","HADHA","ECHS1","crt"),]) > 0) {
    print(B4,vp = viewport(width = 0.03*(Group_numb + 2),
                           height = 0.025*(nrow(Butyrate.ko[Butyrate.ko$Gene %in% c("echA","HADHA","ECHS1","crt"),]) + 1),
                           x = 0.45,y = 0.6,just = "right"))
}
if (nrow(Butyrate.ko[Butyrate.ko$Gene %in% c("ccrA","fabV"),]) > 0) {
    print(B5,vp = viewport(width = 0.028*(Group_numb + 2),
                           height = 0.025*(nrow(Butyrate.ko[Butyrate.ko$Gene %in% c("ccrA","fabV"),]) + 1),
                           x = 0.55,y = 0.45,just = "left"))
}
if (nrow(Butyrate.ko[Butyrate.ko$Gene %in% c("atoAD","ydiF"),]) > 0) {
    print(B6,vp = viewport(width = 0.028*(Group_numb + 2),
                           height = 0.025*(nrow(Butyrate.ko[Butyrate.ko$Gene %in% c("atoAD","ydiF"),]) + 1),
                           x = 0.55,y = 0.23,just = "left"))
}
if (nrow(Butyrate.ko[Butyrate.ko$Gene %in% c("PTB"),]) > 0) {
    print(B7,vp = viewport(width = 0.025*(Group_numb + 2),
                           height = 0.025*(nrow(Butyrate.ko[Butyrate.ko$Gene %in% c("PTB"),]) + 1),
                           x = 0.45,y = 0.38,just = "right"))
}
if (nrow(Butyrate.ko[Butyrate.ko$Gene %in% c("buk"),]) > 0) {
    print(B8,vp = viewport(width = 0.025*(Group_numb + 2),
                           height = 0.025*(nrow(Butyrate.ko[Butyrate.ko$Gene %in% c("buk"),]) + 1),
                           x = 0.45,y = 0.08,just = "right"))
}
dev.off()

### Host_relative stat
tax.name <- c("phylum","class","order","family","genus","species")
source("Functions/SCFA/SCFA.host.R")
SCFA.h <- SCFA.host(Gene)
SCFA.abundance <- as.data.frame(rbind(SCFA.ko[,-1],Acetate.ko,Propionate.ko,Butyrate.ko))
SCFA.abundance <- aggregate(SCFA.abundance[,2:ncol(SCFA.abundance)],
                            list(SCFA.abundance$Gene),sum)
#### relative_Group
dir.create("Results/SCFA/Host_relative_Group")
source("Functions/SCFA/Host.realtive.group.R")
for (i in 1:63) {
    if (names(SCFA.h[[2]])[i] %in% SCFA.abundance$Group.1 & SCFA.h[[1]][[i]] == 1) {
        dir.create(paste("Results/SCFA/Host_relative_Group/",names(SCFA.h[[2]])[i],sep = ""))
        result <- Host.relative.group(Gene,tax,abundance,group,SCFA.h[[2]][[i]])
        title <- names(SCFA.h[[2]][i])
        for (j in 1:6) {
            write.csv(result[[j]],
                      paste("Results/SCFA/Host_relative_Group/",names(SCFA.h[[2]])[i],
                            "/",names(SCFA.h[[2]])[i],"_",tax.name[j],".csv",sep = ""))
            aa <- ifelse(nrow(result[[j]]) > 5,
                         max(str_length(rownames(result[[j]])[1:5])),
                         max(str_length(rownames(result[[j]]))))
            if (nrow(result[[j]]) > 6) {
                result[[j]] <- result[[j]][1:5,]
            }
            pdf(paste("Results/SCFA/Host_relative_Group/",names(SCFA.h[[2]])[i],
                      "/",names(SCFA.h[[2]])[i],"_",tax.name[j],".pdf",sep = ""),
                width = ifelse(aa > 40,aa*0.1 + 0.5*Group_numb,3.2 + 0.5*Group_numb),
                height = 3.5)
            pheatmap(result[[j]],fontsize = 10,cluster_rows = FALSE,fontface = "bold",cluster_cols = FALSE,
                     cellwidth = 30,cellheight = 20,legend = FALSE,breaks = c(seq(0,1,by = 0.01)),
                     color = colorRampPalette(c("white","Red"))(100),border_color = "black")
            grid.text(title,hjust = 0.5,x = 0.5,y = 0.92,gp = gpar(font = 2,size = 1.2))
            dev.off()
        }
    }
}

#### Legend
df <- data.frame(A = c(0,0.1,0.2,0.3,0.4,0.5),
                 B = c(0.5,0.6,0.7,0.8,0.9,1))
df <- t(df)
df1 <- data.frame(A = c("0%","10%","20%","30%","40%","50%"),
                  B = c("50%","60%","70%","80%","90%","100%"))
df1 <- t(df1)
bk <- c(seq(0,1,by = 0.01))
pheatmap(df,fontsize = 30,cluster_rows = FALSE,fontface = "bold",cluster_cols = FALSE,
         cellwidth = 80,cellheight = 50,legend = FALSE,breaks = bk,show_rownames = FALSE,
         color = colorRampPalette(c("white","Red"))(100),show_colnames = FALSE,
         display_numbers = df1,number_color = "black",border_color = "black",
         filename = "Results/SCFA/Host_relative_Group/Legend_relative_Group.pdf",width = 7,height = 2)

### SCFA host stat absolute_Group
dir.create("Results/SCFA/Host_absolute_Group")
source("Functions/SCFA/Host.absolute.group.R")
for (i in 1:63) {
    if (names(SCFA.h[[2]])[i] %in% SCFA.abundance$Group.1 & SCFA.h[[1]][[i]] == 1) {
        dir.create(paste("Results/SCFA/Host_absolute_Group/",names(SCFA.h[[2]])[i],sep = ""))
        result <- Host.absolute.group(Gene,tax,abundance,group,SCFA.h[[2]][[i]])
        title <- names(SCFA.h[[2]][i])
        for (j in 1:6) {
            write.csv(result[[j]],
                      paste("Results/SCFA/Host_absolute_Group/",names(SCFA.h[[2]])[i],
                            "/",names(SCFA.h[[2]])[i],"_",tax.name[j],".csv",sep = ""))
            aa <- ifelse(nrow(result[[j]]) > 5,
                         max(str_length(rownames(result[[j]])[1:5])),
                         max(str_length(rownames(result[[j]]))))
            if (nrow(result[[j]]) > 6) {
                result[[j]] <- result[[j]][1:5,]
            }
            pdf(paste("Results/SCFA/Host_absolute_Group/",names(SCFA.h[[2]])[i],
                      "/",names(SCFA.h[[2]])[i],"_",tax.name[j],".pdf",sep = ""),
                width = ifelse(aa > 40,aa*0.1 + 0.5*Group_numb,3.2 + 0.5*Group_numb),
                height = 3.5)
            pheatmap(result[[j]],fontsize = 10,cluster_rows = FALSE,fontface = "bold",cluster_cols = FALSE,
                     cellwidth = 30,cellheight = 20,
                     color = colorRampPalette(c("white","Red"))(100),border_color = "black")
            grid.text(title,hjust = 0.5,x = 0.5,y = 0.92,gp = gpar(font = 2,size = 1.2))
            dev.off()
        }
    }
}

#### relative_Ungroup
dir.create("Results/SCFA/Host_relative_Ungroup")
source("Functions/SCFA/Host.relative.R")
for (i in 1:63) {
    if (names(SCFA.h[[2]])[i] %in% SCFA.abundance$Group.1 & SCFA.h[[1]][[i]] == 1) {
        dir.create(paste("Results/SCFA/Host_relative_Ungroup/",names(SCFA.h[[2]])[i],sep = ""))
        result <- Host.relative(Gene,tax,abundance,SCFA.h[[2]][[i]])
        title <- names(SCFA.h[[2]][i])
        for (j in 1:6) {
            write.csv(result[[j]],
                      paste("Results/SCFA/Host_relative_Ungroup/",names(SCFA.h[[2]])[i],
                            "/",names(SCFA.h[[2]])[i],"_",tax.name[j],".csv",sep = ""))
            aa <- ifelse(nrow(result[[j]]) > 5,
                         max(str_length(rownames(result[[j]])[1:5])),
                         max(str_length(rownames(result[[j]]))))
            if (nrow(result[[j]]) > 6) {
                result[[j]] <- result[[j]][1:5,]
            }
            pdf(paste("Results/SCFA/Host_relative_Ungroup/",names(SCFA.h[[2]])[i],
                      "/",names(SCFA.h[[2]])[i],"_",tax.name[j],".pdf",sep = ""),
                width = ifelse(aa > 40,aa*0.1 + 0.6*Sample_numb,3.2 + 0.6*Sample_numb),
                height = 3.5)
            pheatmap(result[[j]],fontsize = 10,cluster_rows = FALSE,fontface = "bold",cluster_cols = FALSE,
                     cellwidth = 30,cellheight = 20,legend = FALSE,breaks = c(seq(0,1,by = 0.01)),
                     color = colorRampPalette(c("white","Red"))(100),border_color = "black")
            grid.text(title,hjust = 0.5,x = 0.5,y = 0.92,gp = gpar(font = 2,size = 1.2))
            dev.off()
        }
    }
}

#### Legend
df <- data.frame(A = c(0,0.1,0.2,0.3,0.4,0.5),
                 B = c(0.5,0.6,0.7,0.8,0.9,1))
df <- t(df)
df1 <- data.frame(A = c("0%","10%","20%","30%","40%","50%"),
                  B = c("50%","60%","70%","80%","90%","100%"))
df1 <- t(df1)
bk <- c(seq(0,1,by = 0.01))
pheatmap(df,fontsize = 30,cluster_rows = FALSE,fontface = "bold",cluster_cols = FALSE,
         cellwidth = 80,cellheight = 50,legend = FALSE,breaks = bk,show_rownames = FALSE,
         color = colorRampPalette(c("white","Red"))(100),show_colnames = FALSE,
         display_numbers = df1,number_color = "black",border_color = "black",
         filename = "Results/SCFA/Host_relative_Ungroup/Legend_relative_Ungroup.pdf",width = 7,height = 2)

### SCFA host stat absolute_Ungroup
dir.create("Results/SCFA/Host_absolute_Ungroup")
source("Functions/SCFA/Host.absolute.R")
for (i in 1:63) {
    if (names(SCFA.h[[2]])[i] %in% SCFA.abundance$Group.1 & SCFA.h[[1]][[i]] == 1) {
        dir.create(paste("Results/SCFA/Host_absolute_Ungroup/",names(SCFA.h[[2]])[i],sep = ""))
        result <- Host.absolute(Gene,tax,abundance,SCFA.h[[2]][[i]])
        title <- names(SCFA.h[[2]][i])
        for (j in 1:6) {
            write.csv(result[[j]],
                      paste("Results/SCFA/Host_absolute_Ungroup/",names(SCFA.h[[2]])[i],
                            "/",names(SCFA.h[[2]])[i],"_",tax.name[j],".csv",sep = ""))
            aa <- ifelse(nrow(result[[j]]) > 5,
                         max(str_length(rownames(result[[j]])[1:5])),
                         max(str_length(rownames(result[[j]]))))
            if (nrow(result[[j]]) > 6) {
                result[[j]] <- result[[j]][1:5,]
            }
            pdf(paste("Results/SCFA/Host_absolute_Ungroup/",names(SCFA.h[[2]])[i],
                      "/",names(SCFA.h[[2]])[i],"_",tax.name[j],".pdf",sep = ""),
                width = ifelse(aa > 40,aa*0.1 + 0.6*Sample_numb,3.5 + 0.6*Sample_numb),
                height = 3.5)
            pheatmap(result[[j]],fontsize = 10,cluster_rows = FALSE,fontface = "bold",cluster_cols = FALSE,
                     cellwidth = 30,cellheight = 20,
                     color = colorRampPalette(c("white","Red"))(100),border_color = "black")
            grid.text(title,hjust = 0.5,x = 0.5,y = 0.92,gp = gpar(font = 2,size = 1.2))
            dev.off()
        }
    }
}

# Vitamin metabolism
dir.create("Results/Vitamin")

## Gene identification
dir.create("Results/Vitamin/Gene")

## Pathways of amino acids
dir.create("Results/Vitamin/Gene/Vitamin_metabolism")
source("Functions/Vitamin/Vitamin.metabolism.R")
Vitamin.path <- Vitamin.metabolism(path)
write.table(Vitamin.path,"Results/Vitamin/Gene/Vitamin_metabolism/Vitamin_metabolism_pathway.txt",
            sep = "\t",row.names = FALSE)

### Group heatmap for vitamin metabolism pathways
source("Functions/Vitamin/Abun.heatmap.g.R")
result <- abun.heatmap.g(Vitamin.path,group,Group_numb)
pdf(file = "Results/Vitamin/Gene/Vitamin_metabolism/Vitamin_metabolism_pathway_group.pdf",
    width = 3.55 + 0.4*(Group_numb + 1),
    height = max(str_length(levels(group$Group)))/7 + nrow(result[[1]])*0.3)
result[[2]]
dev.off()

pdf(file = "Results/Vitamin/Gene/Vitamin_metabolism/Vitamin_metabolism_pathway_FC.pdf",
    width = 3.55 + 0.5*(Group_numb + 1),
    height = max(str_length(levels(group$Group)))/7 + nrow(result[[1]])*0.3)
result[[3]]
dev.off()

pdf(file = "Results/Vitamin/Gene/Vitamin_metabolism/Vitamin_metabolism_pathway_test.pdf",
    width = 3.55 + 0.4*Group_numb,
    height = max(str_length(levels(group$Group)))/7 + nrow(result[[1]])*0.3)
result[[4]]
dev.off()

write.table(result[[5]],
            "Results/Vitamin/Gene/Vitamin_metabolism/Vitamin_metabolism_pathway_test.txt",
            sep = "\t",row.names = FALSE)

### Sample heatmap for vitamin metabolism pathways
source("Functions/Vitamin/Abun.heatmap.s.R")
result <- abun.heatmap.s(Vitamin.path,group,Sample_numb)
pdf(file = "Results/Vitamin/Gene/Vitamin_metabolism/Vitamin_metabolism_pathway_sample.pdf",
    width = 3.55 + 0.3*Sample_numb,
    height = max(str_length(group$ID))/7 + nrow(Vitamin.path)*0.3)
result
dev.off()

## Gene abundance in vitamin synthesis
dir.create("Results/Vitamin/Gene/Vitamin_synthesis")
source("Functions/Vitamin/Vitamin.synthesis.R")
source("Functions/Vitamin/Abun.heatmap.g 2.R")
source("Functions/Vitamin/Abun.heatmap.s 2.R")
Vitamin.ko <- Vitamin.synthesis(ko)
write.table(Vitamin.ko,"Results/Vitamin/Gene/Vitamin_synthesis/Vitamin_synthesis_gene.txt",
            sep = "\t",row.names = FALSE)

result <- abun.heatmap.g2(Vitamin.ko,group,Group_numb)
pdf(file = "Results/Vitamin/Gene/Vitamin_synthesis/Vitamin_synthesis_group.pdf",
    width = 3 + 0.4*(Group_numb + 1),
    height = max(str_length(levels(group$Group)))/7 + nrow(result[[1]])*0.25)
result[[2]]
dev.off()

pdf(file = "Results/Vitamin/Gene/Vitamin_synthesis/Vitamin_synthesis_FC.pdf",
    width = 3 + 0.45*(Group_numb + 1),
    height = max(str_length(levels(group$Group)))/7 + nrow(result[[1]])*0.25)
result[[3]]
dev.off()

pdf(file = "Results/Vitamin/Gene/Vitamin_synthesis/Vitamin_synthesis_test.pdf",
    width = 3 + 0.4*Group_numb,
    height = max(str_length(levels(group$Group)))/7 + nrow(result[[1]])*0.25)
result[[4]]
dev.off()

write.table(result[[5]],
            "Results/Vitamin/Gene/Vitamin_synthesis/Vitamin_synthesis_test.txt",
            sep = "\t",row.names = FALSE)

result <- abun.heatmap.s2(Vitamin.ko,group,Sample_numb)
pdf(file = "Results/Vitamin/Gene/Vitamin_synthesis/Vitamin_synthesis_sample.pdf",
    width = 3 + 0.2*Sample_numb,
    height = max(str_length(group$ID))/7 + nrow(Vitamin.ko)*0.2)
result
dev.off()

### Host_relative stat
tax.name <- c("phylum","class","order","family","genus","species")
source("Functions/Vitamin/Vitamin.synthesis.host.R")
Vitamin.h <- Vitamin.synthesis.host(Gene)
Vitamin.ko$Type <- gsub("itamin ","",Vitamin.ko$Type)

#### relative_Group
dir.create("Results/Vitamin/Host_relative_Group")
source("Functions/Vitamin/Host.realtive.group.R")
for (i in 1:5) {
    if (names(Vitamin.h[[2]])[i] %in% Vitamin.ko$Type & Vitamin.h[[1]][[i]] == 1) {
        dir.create(paste("Results/Vitamin/Host_relative_Group/",names(Vitamin.h[[2]])[i],sep = ""))
        result <- Host.relative.group(Gene,tax,abundance,group,Vitamin.h[[2]][[i]])
        title <- names(Vitamin.h[[2]][i])
        for (j in 1:6) {
            write.csv(result[[j]],
                      paste("Results/Vitamin/Host_relative_Group/",names(Vitamin.h[[2]])[i],
                            "/",names(Vitamin.h[[2]])[i],"_",tax.name[j],".csv",sep = ""))
            Vitamin <- ifelse(nrow(result[[j]]) > 5,
                              max(str_length(rownames(result[[j]])[1:5])),
                              max(str_length(rownames(result[[j]]))))
            if (nrow(result[[j]]) > 6) {
                result[[j]] <- result[[j]][1:5,]
            }
            pdf(paste("Results/Vitamin/Host_relative_Group/",names(Vitamin.h[[2]])[i],
                      "/",names(Vitamin.h[[2]])[i],"_",tax.name[j],".pdf",sep = ""),
                width = ifelse(Vitamin > 40,Vitamin*0.1 + 0.5*Group_numb,3.2 + 0.5*Group_numb),
                height = 3.5)
            pheatmap(result[[j]],fontsize = 10,cluster_rows = FALSE,fontface = "bold",cluster_cols = FALSE,
                     cellwidth = 30,cellheight = 20,legend = FALSE,breaks = c(seq(0,1,by = 0.01)),
                     color = colorRampPalette(c("white","Red"))(100),border_color = "black")
            grid.text(title,hjust = 0.5,x = 0.5,y = 0.92,gp = gpar(font = 2,size = 1.2))
            dev.off()
        }
    }
}

#### Legend
df <- data.frame(A = c(0,0.1,0.2,0.3,0.4,0.5),
                 B = c(0.5,0.6,0.7,0.8,0.9,1))
df <- t(df)
df1 <- data.frame(A = c("0%","10%","20%","30%","40%","50%"),
                  B = c("50%","60%","70%","80%","90%","100%"))
df1 <- t(df1)
bk <- c(seq(0,1,by = 0.01))
pheatmap(df,fontsize = 30,cluster_rows = FALSE,fontface = "bold",cluster_cols = FALSE,
         cellwidth = 80,cellheight = 50,legend = FALSE,breaks = bk,show_rownames = FALSE,
         color = colorRampPalette(c("white","Red"))(100),show_colnames = FALSE,
         display_numbers = df1,number_color = "black",border_color = "black",
         filename = "Results/Vitamin/Host_relative_Group/Legend_relative_Group.pdf",width = 7,height = 2)

### Vitamin host stat absolute_Group
dir.create("Results/Vitamin/Host_absolute_Group")
source("Functions/Vitamin/Host.absolute.group.R")
for (i in 1:5) {
    if (names(Vitamin.h[[2]])[i] %in% Vitamin.ko$Type & Vitamin.h[[1]][[i]] == 1) {
        dir.create(paste("Results/Vitamin/Host_absolute_Group/",names(Vitamin.h[[2]])[i],sep = ""))
        result <- Host.absolute.group(Gene,tax,abundance,group,Vitamin.h[[2]][[i]])
        title <- names(Vitamin.h[[2]][i])
        for (j in 1:6) {
            write.csv(result[[j]],
                      paste("Results/Vitamin/Host_absolute_Group/",names(Vitamin.h[[2]])[i],
                            "/",names(Vitamin.h[[2]])[i],"_",tax.name[j],".csv",sep = ""))
            Vitamin <- ifelse(nrow(result[[j]]) > 5,
                              max(str_length(rownames(result[[j]])[1:5])),
                              max(str_length(rownames(result[[j]]))))
            if (nrow(result[[j]]) > 6) {
                result[[j]] <- result[[j]][1:5,]
            }
            pdf(paste("Results/Vitamin/Host_absolute_Group/",names(Vitamin.h[[2]])[i],
                      "/",names(Vitamin.h[[2]])[i],"_",tax.name[j],".pdf",sep = ""),
                width = ifelse(Vitamin > 40,Vitamin*0.1 + 0.5*Group_numb,3.2 + 0.5*Group_numb),
                height = 3.5)
            pheatmap(result[[j]],fontsize = 10,cluster_rows = FALSE,fontface = "bold",cluster_cols = FALSE,
                     cellwidth = 30,cellheight = 20,
                     color = colorRampPalette(c("white","Red"))(100),border_color = "black")
            grid.text(title,hjust = 0.5,x = 0.5,y = 0.92,gp = gpar(font = 2,size = 1.2))
            dev.off()
        }
    }
}

#### relative_Ungroup
dir.create("Results/Vitamin/Host_relative_Ungroup")
source("Functions/Vitamin/Host.relative.R")
for (i in 1:5) {
    if (names(Vitamin.h[[2]])[i] %in% Vitamin.ko$Type & Vitamin.h[[1]][[i]] == 1) {
        dir.create(paste("Results/Vitamin/Host_relative_Ungroup/",names(Vitamin.h[[2]])[i],sep = ""))
        result <- Host.relative(Gene,tax,abundance,Vitamin.h[[2]][[i]])
        title <- names(Vitamin.h[[2]][i])
        for (j in 1:6) {
            write.csv(result[[j]],
                      paste("Results/Vitamin/Host_relative_Ungroup/",names(Vitamin.h[[2]])[i],
                            "/",names(Vitamin.h[[2]])[i],"_",tax.name[j],".csv",sep = ""))
            Vitamin <- ifelse(nrow(result[[j]]) > 5,
                              max(str_length(rownames(result[[j]])[1:5])),
                              max(str_length(rownames(result[[j]]))))
            if (nrow(result[[j]]) > 6) {
                result[[j]] <- result[[j]][1:5,]
            }
            pdf(paste("Results/Vitamin/Host_relative_Ungroup/",names(Vitamin.h[[2]])[i],
                      "/",names(Vitamin.h[[2]])[i],"_",tax.name[j],".pdf",sep = ""),
                width = ifelse(Vitamin > 40,Vitamin*0.1 + 0.6*Sample_numb,3.2 + 0.6*Sample_numb),
                height = 3.5)
            pheatmap(result[[j]],fontsize = 10,cluster_rows = FALSE,fontface = "bold",cluster_cols = FALSE,
                     cellwidth = 30,cellheight = 20,legend = FALSE,breaks = c(seq(0,1,by = 0.01)),
                     color = colorRampPalette(c("white","Red"))(100),border_color = "black")
            grid.text(title,hjust = 0.5,x = 0.5,y = 0.92,gp = gpar(font = 2,size = 1.2))
            dev.off()
        }
    }
}

#### Legend
df <- data.frame(A = c(0,0.1,0.2,0.3,0.4,0.5),
                 B = c(0.5,0.6,0.7,0.8,0.9,1))
df <- t(df)
df1 <- data.frame(A = c("0%","10%","20%","30%","40%","50%"),
                  B = c("50%","60%","70%","80%","90%","100%"))
df1 <- t(df1)
bk <- c(seq(0,1,by = 0.01))
pheatmap(df,fontsize = 30,cluster_rows = FALSE,fontface = "bold",cluster_cols = FALSE,
         cellwidth = 80,cellheight = 50,legend = FALSE,breaks = bk,show_rownames = FALSE,
         color = colorRampPalette(c("white","Red"))(100),show_colnames = FALSE,
         display_numbers = df1,number_color = "black",border_color = "black",
         filename = "Results/Vitamin/Host_relative_Ungroup/Legend_relative_Ungroup.pdf",width = 7,height = 2)

### Vitamin host stat absolute_Ungroup
dir.create("Results/Vitamin/Host_absolute_Ungroup")
source("Functions/Vitamin/Host.absolute.R")
for (i in 1:5) {
    if (names(Vitamin.h[[2]])[i] %in% Vitamin.ko$Type & Vitamin.h[[1]][[i]] == 1) {
        dir.create(paste("Results/Vitamin/Host_absolute_Ungroup/",names(Vitamin.h[[2]])[i],sep = ""))
        result <- Host.absolute(Gene,tax,abundance,Vitamin.h[[2]][[i]])
        title <- names(Vitamin.h[[2]][i])
        for (j in 1:6) {
            write.csv(result[[j]],
                      paste("Results/Vitamin/Host_absolute_Ungroup/",names(Vitamin.h[[2]])[i],
                            "/",names(Vitamin.h[[2]])[i],"_",tax.name[j],".csv",sep = ""))
            Vitamin <- ifelse(nrow(result[[j]]) > 5,
                              max(str_length(rownames(result[[j]])[1:5])),
                              max(str_length(rownames(result[[j]]))))
            if (nrow(result[[j]]) > 6) {
                result[[j]] <- result[[j]][1:5,]
            }
            pdf(paste("Results/Vitamin/Host_absolute_Ungroup/",names(Vitamin.h[[2]])[i],
                      "/",names(Vitamin.h[[2]])[i],"_",tax.name[j],".pdf",sep = ""),
                width = ifelse(Vitamin > 40,Vitamin*0.1 + 0.6*Sample_numb,3.5 + 0.6*Sample_numb),
                height = 3.5)
            pheatmap(result[[j]],fontsize = 10,cluster_rows = FALSE,fontface = "bold",cluster_cols = FALSE,
                     cellwidth = 30,cellheight = 20,
                     color = colorRampPalette(c("white","Red"))(100),border_color = "black")
            grid.text(title,hjust = 0.5,x = 0.5,y = 0.92,gp = gpar(font = 2,size = 1.2))
            dev.off()
        }
    }
}

# Immune
dir.create("Results/Innate_immune")

## Gene identification
dir.create("Results/Innate_immune/Gene")

## Gene abundance in Immune synthesis
dir.create("Results/Innate_immune/Gene/Immune_factor_synthesis")
source("Functions/Immune/Immune.synthesis.R")
source("Functions/Immune/Abun.heatmap.g.R")
source("Functions/Immune/Abun.heatmap.s.R")
Immune.ko <- Immune.synthesis(ko)
write.table(Immune.ko,
            "Results/Innate_immune/Gene/Immune_factor_synthesis/Immune_synthesis_gene.txt",
            sep = "\t",row.names = FALSE)

result <- abun.heatmap.g(Immune.ko,group,Group_numb)
pdf(file = "Results/Innate_immune/Gene/Immune_factor_synthesis/Immune_synthesis_group.pdf",
    width = 3.55 + 0.4*(Group_numb + 1),
    height = max(str_length(levels(group$Group)))/7 + nrow(result[[1]])*0.25)
result[[2]]
dev.off()

pdf(file = "Results/Innate_immune/Gene/Immune_factor_synthesis/Immune_synthesis_FC.pdf",
    width = 3.55 + 0.45*(Group_numb + 1),
    height = max(str_length(levels(group$Group)))/7 + nrow(result[[1]])*0.25)
result[[3]]
dev.off()

pdf(file = "Results/Innate_immune/Gene/Immune_factor_synthesis/Immune_synthesis_test.pdf",
    width = 3.55 + 0.4*Group_numb,
    height = max(str_length(levels(group$Group)))/7 + nrow(result[[1]])*0.25)
result[[4]]
dev.off()

write.table(result[[5]],
            "Results/Innate_immune/Gene/Immune_factor_synthesis/Immune_synthesis_test.txt",
            sep = "\t",row.names = FALSE)

result <- abun.heatmap.s(Immune.ko,group,Sample_numb)
pdf(file = "Results/Innate_immune/Gene/Immune_factor_synthesis/Immune_synthesis_sample.pdf",
    width = 3.85 + 0.2*Sample_numb,
    height = max(str_length(group$ID))/7 + nrow(Immune.ko)*0.2)
result
dev.off()

### Host_relative stat
tax.name <- c("phylum","class","order","family","genus","species")
source("Functions/Immune/Immune.synthesis.host.R")
Immune.h <- Immune.synthesis.host(Gene)
names(Immune.h[[2]]) <- c("Succinate","Indole-3-acetic acid","Indole-3-aldehyde",
                          "Kynurenine","Secondary bile acids","Trimethylamine N-oxide")

#### relative_Group
dir.create("Results/Innate_immune/Host_relative_Group")
source("Functions/Immune/Host.realtive.group.R")
for (i in 1:5) {
    if (names(Immune.h[[2]])[i] %in% Immune.ko$Type & Immune.h[[1]][[i]] == 1) {
        dir.create(paste("Results/Innate_immune/Host_relative_Group/",names(Immune.h[[2]])[i],sep = ""))
        result <- Host.relative.group(Gene,tax,abundance,group,Immune.h[[2]][[i]])
        title <- names(Immune.h[[2]][i])
        for (j in 1:6) {
            write.csv(result[[j]],
                      paste("Results/Innate_immune/Host_relative_Group/",names(Immune.h[[2]])[i],
                            "/",names(Immune.h[[2]])[i],"_",tax.name[j],".csv",sep = ""))
            Immune <- ifelse(nrow(result[[j]]) > 5,
                             max(str_length(rownames(result[[j]])[1:5])),
                             max(str_length(rownames(result[[j]]))))
            if (nrow(result[[j]]) > 6) {
                result[[j]] <- result[[j]][1:5,]
            }
            pdf(paste("Results/Innate_immune/Host_relative_Group/",names(Immune.h[[2]])[i],
                      "/",names(Immune.h[[2]])[i],"_",tax.name[j],".pdf",sep = ""),
                width = ifelse(Immune > 40,Immune*0.1 + 0.5*Group_numb,3.2 + 0.5*Group_numb),
                height = 3.5)
            pheatmap(result[[j]],fontsize = 10,cluster_rows = FALSE,fontface = "bold",cluster_cols = FALSE,
                     cellwidth = 30,cellheight = 20,legend = FALSE,breaks = c(seq(0,1,by = 0.01)),
                     color = colorRampPalette(c("white","Red"))(100),border_color = "black")
            grid.text(title,hjust = 0.5,x = 0.5,y = 0.92,gp = gpar(font = 2,size = 1.2))
            dev.off()
        }
    }
}

#### Legend
df <- data.frame(A = c(0,0.1,0.2,0.3,0.4,0.5),
                 B = c(0.5,0.6,0.7,0.8,0.9,1))
df <- t(df)
df1 <- data.frame(A = c("0%","10%","20%","30%","40%","50%"),
                  B = c("50%","60%","70%","80%","90%","100%"))
df1 <- t(df1)
bk <- c(seq(0,1,by = 0.01))
pheatmap(df,fontsize = 30,cluster_rows = FALSE,fontface = "bold",cluster_cols = FALSE,
         cellwidth = 80,cellheight = 50,legend = FALSE,breaks = bk,show_rownames = FALSE,
         color = colorRampPalette(c("white","Red"))(100),show_colnames = FALSE,
         display_numbers = df1,number_color = "black",border_color = "black",
         filename = "Results/Innate_immune/Host_relative_Group/Legend_relative_Group.pdf",width = 7,height = 2)

### Immune host stat absolute_Group
dir.create("Results/Innate_immune/Host_absolute_Group")
source("Functions/Immune/Host.absolute.group.R")
for (i in 1:5) {
    if (names(Immune.h[[2]])[i] %in% Immune.ko$Type & Immune.h[[1]][[i]] == 1) {
        dir.create(paste("Results/Innate_immune/Host_absolute_Group/",names(Immune.h[[2]])[i],sep = ""))
        result <- Host.absolute.group(Gene,tax,abundance,group,Immune.h[[2]][[i]])
        title <- names(Immune.h[[2]][i])
        for (j in 1:6) {
            write.csv(result[[j]],
                      paste("Results/Innate_immune/Host_absolute_Group/",names(Immune.h[[2]])[i],
                            "/",names(Immune.h[[2]])[i],"_",tax.name[j],".csv",sep = ""))
            Immune <- ifelse(nrow(result[[j]]) > 5,
                             max(str_length(rownames(result[[j]])[1:5])),
                             max(str_length(rownames(result[[j]]))))
            if (nrow(result[[j]]) > 6) {
                result[[j]] <- result[[j]][1:5,]
            }
            pdf(paste("Results/Innate_immune/Host_absolute_Group/",names(Immune.h[[2]])[i],
                      "/",names(Immune.h[[2]])[i],"_",tax.name[j],".pdf",sep = ""),
                width = ifelse(Immune > 40,Immune*0.1 + 0.5*Group_numb,3.2 + 0.5*Group_numb),
                height = 3.5)
            pheatmap(result[[j]],fontsize = 10,cluster_rows = FALSE,fontface = "bold",cluster_cols = FALSE,
                     cellwidth = 30,cellheight = 20,
                     color = colorRampPalette(c("white","Red"))(100),border_color = "black")
            grid.text(title,hjust = 0.5,x = 0.5,y = 0.92,gp = gpar(font = 2,size = 1.2))
            dev.off()
        }
    }
}

#### relative_Ungroup
dir.create("Results/Innate_immune/Host_relative_Ungroup")
source("Functions/Immune/Host.relative.R")
for (i in 1:5) {
    if (names(Immune.h[[2]])[i] %in% Immune.ko$Type & Immune.h[[1]][[i]] == 1) {
        dir.create(paste("Results/Innate_immune/Host_relative_Ungroup/",names(Immune.h[[2]])[i],sep = ""))
        result <- Host.relative(Gene,tax,abundance,Immune.h[[2]][[i]])
        title <- names(Immune.h[[2]][i])
        for (j in 1:6) {
            write.csv(result[[j]],
                      paste("Results/Innate_immune/Host_relative_Ungroup/",names(Immune.h[[2]])[i],
                            "/",names(Immune.h[[2]])[i],"_",tax.name[j],".csv",sep = ""))
            Immune <- ifelse(nrow(result[[j]]) > 5,
                             max(str_length(rownames(result[[j]])[1:5])),
                             max(str_length(rownames(result[[j]]))))
            if (nrow(result[[j]]) > 6) {
                result[[j]] <- result[[j]][1:5,]
            }
            pdf(paste("Results/Innate_immune/Host_relative_Ungroup/",names(Immune.h[[2]])[i],
                      "/",names(Immune.h[[2]])[i],"_",tax.name[j],".pdf",sep = ""),
                width = ifelse(Immune > 40,Immune*0.1 + 0.6*Sample_numb,3.2 + 0.6*Sample_numb),
                height = 3.5)
            pheatmap(result[[j]],fontsize = 10,cluster_rows = FALSE,fontface = "bold",cluster_cols = FALSE,
                     cellwidth = 30,cellheight = 20,legend = FALSE,breaks = c(seq(0,1,by = 0.01)),
                     color = colorRampPalette(c("white","Red"))(100),border_color = "black")
            grid.text(title,hjust = 0.5,x = 0.5,y = 0.92,gp = gpar(font = 2,size = 1.2))
            dev.off()
        }
    }
}

#### Legend
df <- data.frame(A = c(0,0.1,0.2,0.3,0.4,0.5),
                 B = c(0.5,0.6,0.7,0.8,0.9,1))
df <- t(df)
df1 <- data.frame(A = c("0%","10%","20%","30%","40%","50%"),
                  B = c("50%","60%","70%","80%","90%","100%"))
df1 <- t(df1)
bk <- c(seq(0,1,by = 0.01))
pheatmap(df,fontsize = 30,cluster_rows = FALSE,fontface = "bold",cluster_cols = FALSE,
         cellwidth = 80,cellheight = 50,legend = FALSE,breaks = bk,show_rownames = FALSE,
         color = colorRampPalette(c("white","Red"))(100),show_colnames = FALSE,
         display_numbers = df1,number_color = "black",border_color = "black",
         filename = "Results/Innate_immune/Host_relative_Ungroup/Legend_relative_Ungroup.pdf",width = 7,height = 2)

### Immune host stat absolute_Ungroup
dir.create("Results/Innate_immune/Host_absolute_Ungroup")
source("Functions/Immune/Host.absolute.R")
for (i in 1:5) {
    if (names(Immune.h[[2]])[i] %in% Immune.ko$Type & Immune.h[[1]][[i]] == 1) {
        dir.create(paste("Results/Innate_immune/Host_absolute_Ungroup/",names(Immune.h[[2]])[i],sep = ""))
        result <- Host.absolute(Gene,tax,abundance,Immune.h[[2]][[i]])
        title <- names(Immune.h[[2]][i])
        for (j in 1:6) {
            write.csv(result[[j]],
                      paste("Results/Innate_immune/Host_absolute_Ungroup/",names(Immune.h[[2]])[i],
                            "/",names(Immune.h[[2]])[i],"_",tax.name[j],".csv",sep = ""))
            Immune <- ifelse(nrow(result[[j]]) > 5,
                             max(str_length(rownames(result[[j]])[1:5])),
                             max(str_length(rownames(result[[j]]))))
            if (nrow(result[[j]]) > 6) {
                result[[j]] <- result[[j]][1:5,]
            }
            pdf(paste("Results/Innate_immune/Host_absolute_Ungroup/",names(Immune.h[[2]])[i],
                      "/",names(Immune.h[[2]])[i],"_",tax.name[j],".pdf",sep = ""),
                width = ifelse(Immune > 40,Immune*0.1 + 0.6*Sample_numb,3.5 + 0.6*Sample_numb),
                height = 3.5)
            pheatmap(result[[j]],fontsize = 10,cluster_rows = FALSE,fontface = "bold",cluster_cols = FALSE,
                     cellwidth = 30,cellheight = 20,
                     color = colorRampPalette(c("white","Red"))(100),border_color = "black")
            grid.text(title,hjust = 0.5,x = 0.5,y = 0.92,gp = gpar(font = 2,size = 1.2))
            dev.off()
        }
    }
}
