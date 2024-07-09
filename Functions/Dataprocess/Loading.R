
Loading_Color <- function() {
    cbbPalette <- c("#B2182B","#56B4E9","#E69F00","#009E73","#F0E442","#0072B2",
                "#D55E00","#CC79A7","#CC6666","#9999CC","#66CC99","#999999",
                "#ADD1E5")
    return(cbbPalette)
}


Loading_Group <- function(file_dir) {
    group <- read.table(file = file_dir,header = TRUE,sep = "\t")
    colnames(group) <- c("ID","Group")
    group$Group <- factor(group$Group,levels = unique(group$Group))
    return(group)
}

Loading_Table <- function(table_dir) {
    bac <- read.table(file = table_dir,header = FALSE,row.names = 1,sep = "\t",quote = "")
    colnames(bac) <- bac[1,]
    bac <- bac[-1,]
    for (i in 1:(ncol(bac)-1)) {
        bac[,i] <- as.numeric(bac[,i])
    }
    return(bac)
}

Loading_Gene <- function(table_dir){
    Gene <- read.csv(table_dir,header = TRUE,sep = "\t",quote = "")
    Gene <- Gene[!duplicated(Gene[,1:2]),1:2]
    Gene <- Gene %>%
        group_by(Entry) %>%
        mutate(index = row_number()) %>%
        pivot_wider(names_from = Entry,
                    values_from = GeneID) %>%
        select(-index)
    return(Gene)
}

Loading_tax <- function(table_dir){
    tax <- read.csv(table_dir,header = TRUE,sep = "\t")
    tax <- tax[,c(1,7)]
    colnames(tax) <- c("V1","V2")
    tax$V2 <- gsub("d__","k__",tax$V2)
    return(tax)
}

Loading_abun <- function(table_dir){
    abundance <- read.table(table_dir,header = FALSE,row.names = 1,sep = "\t")
    colnames(abundance) <- abundance[1,]
    abundance <- abundance[-1,]
    for (i in 1:ncol(abundance)) {
        abundance[,i] <- as.numeric(abundance[,i])
    }
    abundance$V1 <- rownames(abundance)
    abundance <- abundance[,c("V1",colnames(abundance)[1:(ncol(abundance)-1)])]
    return(abundance)
}



