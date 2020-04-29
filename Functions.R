####
####
#### A few useful functions used in the analysis



##------------------------------------------------------------------------##
##------------------------------------------------------------------------##

## MA plots between two columns of a matrix. Also calculates the proportion of genes that are >2-fold against the logCount.
MA_BPlot <- function(data, col1, col2) {
    arguments <- as.list(match.call())
    y = eval(as.name(arguments$col2), data)
    x = eval(as.name(arguments$col1), data)
    M = log2(x/y)
    A = 0.5*log2(x*y);
    res = data.frame(row.names(data),x=M, y=A, row.names = 1)
    res$bins <- cut(x=res$y, breaks=seq(from=0, to=20, by = 0.5), labels=as.character(seq(0.5,20,0.5)))
    res$count = 1
    bad.res = subset(res, x >= 1 | x <= -1)
    bad.res = subset(bad.res, x != Inf & x != (-Inf))
    badBygroup = tapply(bad.res$count, bad.res$bins, sum)
    allBygroup = tapply(res$count, res$bins, sum)
    off.2fold = badBygroup/allBygroup
    off.2fold = data.frame(off.2fold)
    off.2fold <- cbind(log2Bin = row.names(off.2fold), off.2fold)
    rownames(off.2fold) <- NULL
    off.2fold$log2Bin <- factor(off.2fold$log2Bin, levels=as.character(seq(0.5,20,0.5)))
    b.plot <- ggplot(off.2fold, aes(log2Bin, off.2fold)) + geom_bar(stat = "identity", fill="#F8766D", colour = "#00BFC4") + labs(y="prop. >2-fold", x="log2 read count bin") + theme(legend.position="none") + labs(title = paste(col1, "vs", col2, sep = " "), face = "bold")
    res$FC = NA
    res$FC[res$x > 1 ] = "above" 
    res$FC[res$x < -1 ] = "above"
    ma.plot = ggplot(res, aes(y, x, colour = FC)) + geom_point() + labs(y = "M [log2(x/y)]", x = "A [0.5xlog2(xy)]")+ scale_x_continuous(limits=c(0,20)) + theme(legend.position="none")
    #ma.plot = qplot(y, x, data = res, colour = FC, ylab = "M [log2(x/y)]", xlab = "A [0.5xlog2(xy)]") + scale_x_continuous(limits=c(0,20)) + theme(legend.position="none")
    plots = plot_grid(b.plot, ma.plot, ncol = 1)
    return(plots)
}

##------------------------------------------------------------------------##
##------------------------------------------------------------------------##

## cluster object (can plot with cummerbund's csClusterPlot)
csCluster_c<-function(object, k, logMode=T, method='none', pseudocount=1,...){
    require(cluster)
    m<-as.data.frame(object)
    m<-m[rowSums(m)>0,]
    if(logMode){
        m<-log10(m+pseudocount)
    }
    
    if(!is.function(method)){
        method = function(mat){JSdist(makeprobs(t(m)))}	
    }		
    n<-method(m)
    clusters<-pam(n,k, ...)
    #clsuters<-pamk(n,krange=2:20)
    class(clusters)<-"list"
    clusters$fpkm<-m
    clusters
}

##------------------------------------------------------------------------##
##------------------------------------------------------------------------##

geneBoxPlot <- function (tpmTable, gene, show_reps = F) {
        if (grepl("FBgn", gene)) {
            description <- subset(snapshots, FBgn_ID == gene)$GeneName
            geneName <- subset(tpmTable, gene_id == gene)$gene_symbol
            cg_id <- subset(tpmTable, gene_id == gene)$annotation_ID
            
            p <- ggplot(subset(tpmTable, gene_id == gene), aes(Male, TPM, fill = Status)) +
                geom_boxplot(outlier.size = 0, width = 0.3) +
                geom_point(pch = 21, position = position_jitterdodge(), aes(fill = Status)) + 
                facet_grid(.~Female, scales = "free_x", space = "free_x") +
                labs(title = paste(gene, " (", geneName, ")", sep = ""), subtitle = paste(cg_id, ": ", description, sep = "")) + 
                theme_bw() + 
                theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
                    axis.text.y = element_text(size = 14),
                    strip.text = element_text(size = 14, face = "bold")) +
                scale_fill_manual(values = c("#3f5a2a","#ffb200"))
        } else {
            fbgn_id <- subset(tpmTable, gene_symbol == gene)$gene_id
            description <- subset(snapshots, GeneSymbol == gene)$GeneName
            cg_id <- subset(tpmTable, gene_symbol == gene)$annotation_ID
            
            p <- ggplot(subset(tpmTable, gene_symbol == gene), aes(Male, TPM, fill = Status)) +
                geom_boxplot(outlier.size = 0, width = 0.3) +
                geom_point(pch = 21, position = position_jitterdodge(), aes(fill = Status)) + 
                facet_grid(.~Female, scales = "free_x", space = "free_x") +
                labs(title = paste(gene, " (", fbgn_id, ")", sep = ""), subtitle = paste(cg_id, ": ", description, sep = "")) + 
                theme_bw() + 
                theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
                    axis.text.y = element_text(size = 14),
                    strip.text = element_text(size = 14, face = "bold")) +
                scale_fill_manual(values = c("#3f5a2a","#ffb200"))
        }
        if (show_reps){
            p <- p + geom_text_repel(aes(label=replicate), force = 10)
        } 
    return(p)
}

##------------------------------------------------------------------------##
##------------------------------------------------------------------------##

geneBoxPlot_mon <- function (tpmTable, gene, show_reps = F) {

        if (grepl("FBgn", gene)) {
            description <- subset(snapshots, FBgn_ID == gene)$GeneName
            geneName <- subset(tpmTable, gene_id == gene)$gene_symbol
            cg_id <- subset(tpmTable, gene_id == gene)$annotation_ID
            
            p <- ggplot(subset(tpmTable, gene_id == gene), aes(Male, TPM, fill = Status, colour = Status)) +
                geom_boxplot(outlier.size = 0, width = 0.3, lwd = 0.1) +
                geom_point(pch = 21, position = position_jitterdodge(), aes(fill = Status)) + 
                facet_grid(.~Female, scales = "free_x", space = "free_x") +
                labs(title = paste(gene, " (", geneName, ")", sep = ""), subtitle = paste(cg_id, ": ", description, sep = "")) + 
                theme_bw() + 
                theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
                    axis.text.y = element_text(size = 14),
                    strip.text = element_text(size = 14, face = "bold")) +
                scale_fill_manual(values = c("#028fc2","#ffb200")) +
                scale_colour_manual(values = c("#d3004a", "#3f5a2a"))
        } else {
            fbgn_id <- subset(tpmTable, gene_symbol == gene)$gene_id
            description <- subset(snapshots, GeneSymbol == gene)$GeneName
            cg_id <- subset(tpmTable, gene_symbol == gene)$annotation_ID
            
            p <- ggplot(subset(tpmTable, gene_symbol == gene), aes(Male, TPM, fill = Status, colour = Status)) +
                geom_boxplot(outlier.size = 0, width = 0.3) +
                geom_point(pch = 21, position = position_jitterdodge(), aes(fill = Status)) + 
                facet_grid(.~Female, scales = "free_x", space = "free_x") +
                labs(title = paste(gene, " (", fbgn_id, ")", sep = ""), subtitle = paste(cg_id, ": ", description, sep = "")) + 
                theme_bw() + 
                theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
                    axis.text.y = element_text(size = 14),
                    strip.text = element_text(size = 14, face = "bold")) +
                scale_fill_manual(values = c("#028fc2","#ffb200")) +
                scale_colour_manual(values = c("#d3004a", "#3f5a2a"))
        }
        if (show_reps){
            p <- p + geom_text_repel(aes(label=replicate, colour = Status), force = 10)
        } 
    return(p)
}

##------------------------------------------------------------------------##
##------------------------------------------------------------------------##

geneBarPlot_mon <- function (tpmTable, gene, show_reps = F, show_rep_names = F) {

        if (grepl("FBgn", gene)) {
            description <- subset(snapshots, FBgn_ID == gene)$GeneName
            geneName <- subset(tpmTable, gene_id == gene)$gene_symbol
            cg_id <- subset(tpmTable, gene_id == gene)$annotation_ID
            
            tmpDF = subset(tpmTable, gene_id == gene)
            tmpDF.se = summarySE(tmpDF, measurevar = "TPM", groupvars = c("gene_id", "gene_symbol", "annotation_ID", "Sample", "Female", "Male"))

            p <- ggplot() +
                geom_bar(data = tmpDF.se, mapping = aes(Female, TPM, colour = Male, fill = Male), stat = "identity", position = position_dodge(.6), size = 0.7, width = 0.5) +
                geom_errorbar(data = tmpDF.se, aes(Female, colour = Male, ymin=TPM-se, ymax=TPM+se), width=.2, position=position_dodge(.6), size = 1) + 
                labs(title = paste(gene, " (", geneName, ")", sep = ""), subtitle = paste(cg_id, ": ", description, sep = "")) + 
                theme(axis.text.x = element_text(size = 13),
                    axis.text.y = element_text(size = 14)) +
                scale_colour_manual(values = c("#6d446e","#7cdc66","#f7006c")) +
                scale_fill_manual(values = c("#01c5a8","#989000","#a7ccde"))
                
        } else {
            fbgn_id <- subset(tpmTable, gene_symbol == gene)$gene_id
            description <- subset(snapshots, GeneSymbol == gene)$GeneName
            cg_id <- subset(tpmTable, gene_symbol == gene)$annotation_ID
            
            tmpDF = subset(tpmTable, gene_symbol == gene)
            tmpDF.se = summarySE(tmpDF, measurevar = "TPM", groupvars = c("gene_id", "gene_symbol", "annotation_ID", "Sample", "Female", "Male"))

            p <- ggplot() +
                geom_bar(data = tmpDF.se, mapping = aes(Female, TPM, colour = Male, fill = Male), stat = "identity", position = position_dodge(.6), size = 0.7, width = 0.5) +
                geom_errorbar(data = tmpDF.se, aes(Female, colour = Male, ymin=TPM-se, ymax=TPM+se), width=.2, position=position_dodge(.6), size = 1) + 
                labs(title = paste(gene, " (", fbgn_id, ")", sep = ""), subtitle = paste(cg_id, ": ", description, sep = "")) + 
                theme(axis.text.x = element_text(size = 13),
                    axis.text.y = element_text(size = 14)) +
                scale_colour_manual(values = c("#6d446e","#7cdc66","#f7006c")) +
                scale_fill_manual(values = c("#01c5a8","#989000","#a7ccde"))
        }
        if (show_reps){
            p <- p + geom_point(data = tmpDF, mapping = aes(Female, TPM, colour = Male), position = position_jitterdodge(jitter.width = 0,dodge.width = 0.6), size = 0.8)

            if (show_rep_names) {
                p <- p + geom_text_repel(data = tmpDF, aes(label=replicate, colour = Male), force = 1)
            }
        } 
    return(p)
}
##------------------------------------------------------------------------##
##------------------------------------------------------------------------##

geneBoxPlot_mon_byFemale <- function (tpmTable, gene, female, show_reps = F) 
{
    if (grepl("FBgn", gene)) {
        description <- subset(snapshots, FBgn_ID == gene)$GeneName
        geneName <- subset(tpmTable, gene_id == gene)$gene_symbol
        cg_id <- subset(tpmTable, gene_id == gene)$annotation_ID
        p <- ggplot(subset(tpmTable, gene_id == gene & Female == female), aes(Male, 
            TPM, fill = Status, colour = Status)) + geom_boxplot(outlier.size = 0, 
            width = 0.3, lwd = 0.1) + geom_point(pch = 21, position = position_jitterdodge(), 
            aes(fill = Status)) + labs(title = paste(gene, " (", 
            geneName, ")", sep = ""), subtitle = paste(cg_id, 
            ": ", description, sep = "")) + theme_bw() + theme(axis.text.x = element_text(angle = 45, 
            hjust = 1, size = 13), axis.text.y = element_text(size = 14), 
            strip.text = element_text(size = 14, face = "bold")) + 
            scale_fill_manual(values = c("#028fc2", "#ffb200")) + 
            scale_colour_manual(values = c("#d3004a", "#3f5a2a"))
    }
    else {
        fbgn_id <- subset(tpmTable, gene_symbol == gene)$gene_id
        description <- subset(snapshots, GeneSymbol == gene)$GeneName
        cg_id <- subset(tpmTable, gene_symbol == gene)$annotation_ID
        p <- ggplot(subset(tpmTable, gene_symbol == gene & Female == female), aes(Male, 
            TPM, fill = Status, colour = Status)) + geom_boxplot(outlier.size = 0, 
            width = 0.3) + geom_point(pch = 21, position = position_jitterdodge(), 
            aes(fill = Status)) + labs(title = paste(gene, " (", 
            fbgn_id, ")", sep = ""), subtitle = paste(cg_id, 
            ": ", description, sep = "")) + theme_bw() + theme(axis.text.x = element_text(angle = 45, 
            hjust = 1, size = 13), axis.text.y = element_text(size = 14), 
            strip.text = element_text(size = 14, face = "bold")) + 
            scale_fill_manual(values = c("#028fc2", "#ffb200")) + 
            scale_colour_manual(values = c("#d3004a", "#3f5a2a"))
    }
    if (show_reps) {
        p <- p + geom_text_repel(aes(label = replicate, colour = Status), 
            force = 10)
    }
    return(p)
}

##------------------------------------------------------------------------##
##------------------------------------------------------------------------##

geneBoxPlot_fa2 <- function (tpmTable, gene) {

        if (grepl("FBgn", gene)) {
            description <- subset(snapshots, FBgn_ID == gene)$GeneName
            geneName <- subset(tpmTable, ref_gene_id == gene)$gene_name
            detail <- subset(snapshots, FBgn_ID == gene)$gene_snapshot_text
            
            p <- ggplot(subset(tpmTable, ref_gene_id == gene), aes(Library_Name, TPM, fill = Sex, colour = Sex)) + 
                geom_boxplot(lwd = 0.3) + 
#                 geom_jitter() +
                facet_grid(.~dev_stage, scale = "free_x", space = "free_x") + 
                labs(title = paste(gene, " (", geneName, "): ", description, sep = ""), subtitle = paste(str_wrap(detail, width = 110), sep = "")) + 
                # theme_monokai_full() +
                theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                scale_fill_brewer(palette="Dark2") +
                scale_colour_brewer(palette="Set1")
            } else {
            description <- subset(snapshots, GeneSymbol == gene)$GeneName
            geneName <- subset(tpmTable, gene_name == gene)$ref_gene_id
            detail <- subset(snapshots, GeneSymbol == gene)$gene_snapshot_text
            
            p <- ggplot(subset(tpmTable, gene_name == gene), aes(Library_Name, TPM, fill = Sex, colour = Sex)) + 
                geom_boxplot(lwd = 0.3) + 
#                 geom_jitter() +
                facet_grid(.~dev_stage, scale = "free_x", space = "free_x") + 
                labs(title = paste(gene, " (", geneName, "): ", description, sep = ""), subtitle = paste(str_wrap(detail, width = 110), sep = "")) + 
                # theme_monokai_full() +
                theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                scale_fill_brewer(palette="Dark2") +
                scale_colour_brewer(palette="Set1")
        }
    return(p)
}

##------------------------------------------------------------------------##
##------------------------------------------------------------------------##

heatmap_fa2 <- function(matrix, gene_list, fly_atlas = F, title) {

if (fly_atlas) {
    ## make an annotation bar for "sex"
    col_annot = unique(select(sample.info, Sample_Name, Sex, dev_stage))
    rownames(col_annot) = col_annot$Sample_Name
    rownames(col_annot) = gsub("_", " ", rownames(col_annot))
    col_annot = subset(col_annot, select = c("Sex", "dev_stage"))

    # # set colors
    mat_colors = list(Sex = brewer.pal(3, "Set1"), dev_stage = brewer.pal(2, "Set2"))
    names(mat_colors$Sex) = unique(sample.info$Sex)
    names(mat_colors$dev_stage) = unique(sample.info$dev_stage)
} else {
    col_annot = unique(select(sampleInfo, Replicate, Status, Female, Male, Handler))
    rownames(col_annot) = col_annot$Replicate
    rownames(col_annot) = gsub("_", " ", rownames(col_annot))
    col_annot = subset(col_annot, select = c("Status", "Female", "Male", "Handler"))

    # # set colors
    # mat_colors <- list(group = brewer.pal(3, "Set1"))
    # names(mat_colors$group) <- unique(col_groups)
}
    

## process tmm matrix
data = subset(matrix, rownames(matrix) %in% gene_list)
data = log2(data+1)
data = as.data.frame(t(scale(t(data), scale=F)))
data[data < -2] = -2
data[data > 2] = 2
init_cols = colnames(data)
data$FBgn_ID = rownames(data)
data = merge(data, FBgn_to_symbol, by.x = "FBgn_ID", by.y = "primary_FBgn", all.x = T)
rownames(data) = data$gene_symbol
data = subset(data, select = init_cols)
colnames(data) = gsub("_", " ", colnames(data))

# plotter
p <- pheatmap(
  mat               = data,
  main              = title,
  color             = inferno(100),
  border_color      = NA,
  show_colnames     = TRUE,
  show_rownames     = TRUE,
  annotation_col    = col_annot,
  drop_levels       = TRUE,
#   cluster_col    = FALSE,
  annotation_names_row = F,
  fontsize          = 8    
)
    return(p)
    }


heatmap_mean <- function(tpmTable, gene_list, title) {

    col_annot = unique(select(sampleInfo, Sample, Status, Female, Male))
    rownames(col_annot) = col_annot$Sample
    col_annot = subset(col_annot, select = c("Status", "Female", "Male"))

    tmpMat<-cast(subset(tpmTable, gene_id %in% gene_list), gene_id~Sample, value ="TPM", fun.aggregate = mean)
    data <- tmpMat[,-1]
    rownames(data) <- tmpMat[,1]
    
## process tmm matrix
data = log2(data+1)
data = as.data.frame(t(scale(t(data), scale=F)))
data[data < -2] = -2
data[data > 2] = 2
init_cols = colnames(data)
data$FBgn_ID = rownames(data)
data = merge(data, FBgn_to_symbol, by.x = "FBgn_ID", by.y = "primary_FBgn", all.x = T)
rownames(data) = data$gene_symbol
data = subset(data, select = init_cols)
colnames(data) = gsub("_", " ", colnames(data))

# plotter
p <- pheatmap(
  mat               = data,
  main              = title,
  color             = inferno(100),
  border_color      = NA,
  show_colnames     = TRUE,
  show_rownames     = TRUE,
  annotation_col    = col_annot,
  drop_levels       = TRUE,
#   cluster_col    = FALSE,
  annotation_names_row = F,
  fontsize          = 8    
)
    return(p)
    }


geneBarPlot <- function (tpmTable, gene, show_reps = F, show_rep_names = F) 
{
    if (grepl("FBgn", gene)) {
        description <- subset(snapshots, FBgn_ID == gene)$GeneName
        geneName <- subset(tpmTable, gene_id == gene)$gene_symbol
        cg_id <- subset(tpmTable, gene_id == gene)$annotation_ID
        tmpDF = subset(tpmTable, gene_id == gene)
        tmpDF.se = summarySE(tmpDF, measurevar = "TPM", groupvars = c("gene_id", 
            "gene_symbol", "annotation_ID", "Sample", "Female", 
            "Male"))
        p <- ggplot() + geom_bar(data = tmpDF.se, mapping = aes(Female, 
            TPM, fill = Male), stat = "identity", 
            position = position_dodge(0.6), size = 0.7, width = 0.5) + 
            geom_errorbar(data = tmpDF.se, aes(Female, 
                ymin = TPM - se, ymax = TPM + se, colour = Male), width = 0.2, 
                position = position_dodge(0.6), size = 1) + labs(title = paste(gene, 
            " (", geneName, ")", sep = ""), subtitle = paste(cg_id, 
            ": ", description, sep = "")) + facet_grid(.~Female, space = "free_x", scale = "free_x") + theme(axis.text.x = element_text(size = 13), 
            axis.text.y = element_text(size = 14)) + scale_fill_manual(values = c("#003440", 
            "#003440","#003440"))
    }
    else {
        fbgn_id <- subset(tpmTable, gene_symbol == gene)$gene_id
        description <- subset(snapshots, GeneSymbol == gene)$GeneName
        cg_id <- subset(tpmTable, gene_symbol == gene)$annotation_ID
        tmpDF = subset(tpmTable, gene_symbol == gene)
        tmpDF.se = summarySE(tmpDF, measurevar = "TPM", groupvars = c("gene_id", 
            "gene_symbol", "annotation_ID", "Sample", "Female", 
            "Male"))
        p <- ggplot() + geom_bar(data = tmpDF.se, mapping = aes(Female, 
            TPM, fill = Male), stat = "identity", 
            position = position_dodge(0.6), size = 0.7, width = 0.5) + 
            geom_errorbar(data = tmpDF.se, aes(Female, 
                ymin = TPM - se, ymax = TPM + se, colour = Male), width = 0.2, 
                position = position_dodge(0.6), size = 1) + labs(title = paste(gene, 
            " (", fbgn_id, ")", sep = ""), subtitle = paste(cg_id, 
            ": ", description, sep = "")) + facet_grid(.~Female, space = "free_x", scale = "free_x") + theme(axis.text.x = element_text(size = 13), 
            axis.text.y = element_text(size = 14)) + scale_fill_manual(values = c("gray", 
            "#58b8f0","#ccae50")) + scale_colour_manual(values = c("#003440", 
            "#003440","#003440"))
    }
    if (show_reps) {
        p <- p + geom_point(data = tmpDF, mapping = aes(Female, 
            TPM, colour = Male), position = position_jitterdodge(jitter.width = 0.1, 
            dodge.width = 0.6), size = 0.8)
        if (show_rep_names) {
            p <- p + geom_text_repel(data = tmpDF, aes(label = replicate, 
                colour = Male), force = 1)
        }
    }
    return(p)
}


geneBarPlot_byFemale <- function (tpmTable, gene, show_reps = F, show_rep_names = F, female = "axenic") 
{
    if (grepl("FBgn", gene)) {
        description <- subset(snapshots, FBgn_ID == gene)$GeneName
        geneName <- subset(tpmTable, gene_id == gene)$gene_symbol
        cg_id <- subset(tpmTable, gene_id == gene)$annotation_ID
        tmpDF = subset(tpmTable, gene_id == gene)
        tmpDF.se = summarySE(tmpDF, measurevar = "TPM", groupvars = c("gene_id", 
            "gene_symbol", "annotation_ID", "Sample", "Female", 
            "Male"))
        p <- ggplot() + geom_bar(data = filter(tmpDF.se, Female == female), mapping = aes(Female, 
            TPM, fill = Male), stat = "identity", 
            position = position_dodge(0.6), size = 0.7, width = 0.5) + 
            geom_errorbar(data = filter(tmpDF.se, Female == female), aes(Female, 
                ymin = TPM - se, ymax = TPM + se, colour = Male), width = 0.2, 
                position = position_dodge(0.6), size = 1) + labs(title = paste(gene, 
            " (", geneName, ")", sep = ""), subtitle = paste(cg_id, 
            ": ", description, sep = "")) + theme(axis.text.x = element_text(size = 13), 
            axis.text.y = element_text(size = 14)) + scale_fill_manual(values = c("#003440", 
            "#003440","#003440"))
    }
    else {
        fbgn_id <- subset(tpmTable, gene_symbol == gene)$gene_id
        description <- subset(snapshots, GeneSymbol == gene)$GeneName
        cg_id <- subset(tpmTable, gene_symbol == gene)$annotation_ID
        tmpDF = subset(tpmTable, gene_symbol == gene)
        tmpDF.se = summarySE(tmpDF, measurevar = "TPM", groupvars = c("gene_id", 
            "gene_symbol", "annotation_ID", "Sample", "Female", 
            "Male"))
        p <- ggplot() + geom_bar(data = filter(tmpDF.se, Female == female), mapping = aes(Female, 
            TPM, fill = Male), stat = "identity", 
            position = position_dodge(0.6), size = 0.7, width = 0.5) + 
            geom_errorbar(data = filter(tmpDF.se, Female == female), aes(Female, 
                ymin = TPM - se, ymax = TPM + se, colour = Male), width = 0.2, 
                position = position_dodge(0.6), size = 1) + labs(title = paste(gene, 
            " (", fbgn_id, ")", sep = ""), subtitle = paste(cg_id, 
            ": ", description, sep = "")) + theme(axis.text.x = element_text(size = 13), 
            axis.text.y = element_text(size = 14)) + scale_fill_manual(values = c("gray", 
            "#58b8f0","#ccae50")) + scale_colour_manual(values = c("#003440", 
            "#003440","#003440"))
    }
    if (show_reps) {
        p <- p + geom_point(data = filter(tmpDF, Female == female), mapping = aes(Female, 
            TPM, colour = Male), position = position_jitterdodge(jitter.width = 0.1, 
            dodge.width = 0.6), size = 0.8)
        if (show_rep_names) {
            p <- p + geom_text_repel(data = filter(tmpDF, Female == female), aes(label = replicate, 
                colour = Male), force = 1)
        }
    }
    return(p)
}


 geneBarPlot_mon_byFemale <- function (tpmTable, gene, show_reps = F, show_rep_names = F, female = "axenic") 
{
    if (grepl("FBgn", gene)) {
        description <- subset(snapshots, FBgn_ID == gene)$GeneName
        geneName <- subset(tpmTable, gene_id == gene)$gene_symbol
        cg_id <- subset(tpmTable, gene_id == gene)$annotation_ID
        tmpDF = subset(tpmTable, gene_id == gene)
        tmpDF.se = summarySE(tmpDF, measurevar = "TPM", groupvars = c("gene_id", 
            "gene_symbol", "annotation_ID", "Sample", "Female", 
            "Male"))
        p <- ggplot() + geom_bar(data = filter(tmpDF.se, Female == female), mapping = aes(Female, 
            TPM, fill = Male), stat = "identity", 
            position = position_dodge(0.6), size = 0.7, width = 0.5) + 
            geom_errorbar(data = filter(tmpDF.se, Female == female), aes(Female, 
                ymin = TPM - se, ymax = TPM + se), width = 0.2, 
                position = position_dodge(0.6), size = 1) + labs(title = paste(gene, 
            " (", geneName, ")", sep = ""), subtitle = paste(cg_id, 
            ": ", description, sep = "")) + theme(axis.text.x = element_text(size = 13), 
            axis.text.y = element_text(size = 14)) + scale_fill_manual(values = c("gray", 
            "#58b8f0","#ccae50"))
    }
    else {
        fbgn_id <- subset(tpmTable, gene_symbol == gene)$gene_id
        description <- subset(snapshots, GeneSymbol == gene)$GeneName
        cg_id <- subset(tpmTable, gene_symbol == gene)$annotation_ID
        tmpDF = subset(tpmTable, gene_symbol == gene)
        tmpDF.se = summarySE(tmpDF, measurevar = "TPM", groupvars = c("gene_id", 
            "gene_symbol", "annotation_ID", "Sample", "Female", 
            "Male"))
        p <- ggplot() + geom_bar(data = filter(tmpDF.se, Female == female), mapping = aes(Female, 
            TPM, fill = Male), stat = "identity", 
            position = position_dodge(0.6), size = 0.7, width = 0.5) + 
            geom_errorbar(data = filter(tmpDF.se, Female == female), aes(Female, 
                ymin = TPM - se, ymax = TPM + se), width = 0.2, 
                position = position_dodge(0.6), size = 1) + labs(title = paste(gene, 
            " (", fbgn_id, ")", sep = ""), subtitle = paste(cg_id, 
            ": ", description, sep = "")) + theme(axis.text.x = element_text(size = 13), 
            axis.text.y = element_text(size = 14)) + scale_fill_manual(values = c("gray", 
            "#58b8f0","#ccae50"))
    }
    if (show_reps) {
        p <- p + geom_point(data = filter(tmpDF, Female == female), mapping = aes(Female, 
            TPM, colour = Male), position = position_jitterdodge(jitter.width = 0, 
            dodge.width = 0.6), size = 0.8)
        if (show_rep_names) {
            p <- p + geom_text_repel(data = filter(tmpDF, Female == female), aes(label = replicate, 
                colour = Male), force = 1)
        }
    }
    return(p)
}

##------------------------------------------------------------------------##
##------------------------------------------------------------------------##

 edge.DE <- function(fit, contrast, LRT = F) {

    if(LRT){
        lrt <- glmLRT(fit, contrast = contrast)
        } else {
            lrt <- glmQLFTest(fit, contrast = contrast)
        }
    
    lrt.tTags <- topTags(lrt, n = NULL)
    lrt.table <- lrt.tTags$table
    lrt.table$sig = ifelse(lrt.table$FDR < 0.05 & (lrt.table$logFC > 1 | lrt.table$logFC < -1), "yes", "no")
    lrt.table$direction = ifelse(lrt.table$sig == "yes" & lrt.table$logFC > 1, "Up", ifelse(lrt.table$sig == "yes" &    lrt.table$logFC < -1, "Down", "n.s."))
    lrt.table$gene = rownames(lrt.table)
    lrt.table = merge(lrt.table, FBgn_to_symbol, by.x = "gene", by.y = "primary_FBgn", all.x = T)
    return(lrt.table)
}

##------------------------------------------------------------------------##
##------------------------------------------------------------------------##

c.mds.gg <- function(dge){
    mdsObj <- plotMDS(dge, plot = F)$cmdscale.out
    mdsObj <- as.data.frame(as.matrix(mdsObj))
    mdsObj$replicate <- rownames(mdsObj)
    colnames(mdsObj) = c("dim1", "dim2", "replicate")
    mdsObj = merge(mdsObj, sampleInfo, by.x = "replicate", by.y = "Replicate")
    mdsObj$replicate_num = gsub(".*_", "", mdsObj$replicate)

    p <- ggscatter(mdsObj, 
              x = "dim1", 
              y = "dim2",
              color = "Sample",
              shape = "Female",
              size = 3.5,
              alpha = 0.8, 
              ggtheme = theme_bw(),
              repel = "Time",) + 
                theme(axis.text = element_text(size = 10), legend.title = element_blank(), axis.title = element_text(size = 12), legend.text = element_text(size = 12)) +
                geom_text_repel(aes(label=replicate_num), 
                    force = 10, 
                    box.padding = unit(0.35, "lines"), 
                    point.padding = unit(0.5, "lines"), 
                    fontface = "bold", 
                    size = 3,
                    colour = "grey") +
                labs ( x = "Dimension 1", y = "Dimension 2")
    return(p)
}

##------------------------------------------------------------------------##
##------------------------------------------------------------------------##

diag.plot <- function(edge.DE.obj1, edge.DE.obj2){
    df = merge(edge.DE.obj1, edge.DE.obj2, by.x = "gene", by.y = "gene", all.x = T)

    x_axis <- deparse(substitute(edge.DE.obj1))
    x_axis <- gsub("lrt.", "logFC ", x_axis)
    x_axis <- gsub(".table", "", x_axis)
    y_axis <- deparse(substitute(edge.DE.obj2))
    y_axis <- gsub("lrt.", "logFC ", y_axis)
    y_axis <- gsub(".table", "", y_axis)

    p <- ggplot() + 
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.8, colour = "grey") + 
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.8, colour = "grey") +
    geom_abline(slope = 1, linetype = "dashed", alpha = 0.8, colour = "grey") +
    geom_point(data = subset(df, sig.x == "no" & sig.y == "no"),
               aes(logFC.x, logFC.y), 
               colour = "gray",
               alpha = 0.4) +
    geom_point(data = subset(df, sig.x == "yes" & sig.y == "no"),
               aes(logFC.x, logFC.y, colour = "#007954"), 
               alpha = 0.8,
               size = 3) +
    geom_point(data = subset(df, sig.x == "no" & sig.y == "yes"),
               aes(logFC.x, logFC.y, colour = "#c9578c"), 
               alpha = 0.8,
               size = 3) +
    geom_point(data = subset(df, sig.x == "yes" & sig.y == "yes"),
               aes(logFC.x, logFC.y, colour = "orange"), 
               alpha = 0.8,
               size = 3) +
    theme_bw() +
    theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14)) +
    scale_colour_manual(name = 'Status', values =c('#007954'='#007954','#c9578c'='#c9578c','orange'='orange'), labels = c('sig in X','sig in Y','normal')) +
    labs(x = x_axis, y = y_axis)
    return(p)
}

##------------------------------------------------------------------------##
##------------------------------------------------------------------------##

edge.DE.gg <- function(edge.DE.obj) {
    p <- ggplot(edge.DE.obj, aes(logFC, -log10(PValue), colour = sig)) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.8, colour = "grey") +
    geom_point(alpha = 0.5) +
    theme_bw()
    return(p)
}

##------------------------------------------------------------------------##
##------------------------------------------------------------------------##

edge.DE.gg.size <- function (edge.DE.obj) 
{
    p <- ggplot(edge.DE.obj, aes(logFC, -log10(PValue), shape = sig, colour = DE.cat)) + 
        geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.8, 
            colour = "grey") + geom_point(alpha = 0.5) + theme_bw()
    return(p)
}

##------------------------------------------------------------------------##
##------------------------------------------------------------------------##

## extract gene IDs based on GO term and enriched gene list:
extract_GO_genes = function(go_term, gene_set){
    rownames(subset(GO_info, row.names(GO_info) %in% gene_set & grepl(go_term, GO_info$term)))
}

##------------------------------------------------------------------------##
##------------------------------------------------------------------------##

## extract GO term descriptions for GOseq:
get_GO_term_descr =  function(x) {
    d = 'none';
    go_info = GOTERM[[x]];
    if (length(go_info) >0) { d = paste(Ontology(go_info), Term(go_info), sep=' ');}
    return(d);
}

##------------------------------------------------------------------------##
##------------------------------------------------------------------------##

corr_eqn <- function(x,y, digits = 2) {
  corr_coef <- round(cor(x, y, use = "pairwise.complete.obs"), digits = digits)
  paste("italic(R)^2 == ", corr_coef)
}

##------------------------------------------------------------------------##
##------------------------------------------------------------------------##

## A function to calculate the tissue specificity index (based on CummerBund's S function)
calcSpecificity<-function(matrix,logMode=T,pseudocount=1,relative=FALSE){
    tpms<-matrix
    if(logMode){
        tpms<-log10(tpms+pseudocount)
    }
    tpms<-t(makeprobs(t(tpms)))
    d<-diag(ncol(tpms))
    res<-apply(d,MARGIN=1,function(q){
        JSdistFromP(tpms,q)
    })
    colnames(res)<-paste(colnames(tpms))
    
    if(relative){
        res<-res/max(res)
    }
    1-res
}

##------------------------------------------------------------------------##
##------------------------------------------------------------------------##

## estimate number of cluster for K-means
findK<-function(object, k.range=c(2:20), logMode=T, pseudocount=1,...){
    require(cluster)
    m<-as.data.frame(object)
    m<-m[rowSums(m)>0,]
    if(logMode){
        m<-log10(m+pseudocount)
    }
    n<-JSdist(makeprobs(t(m)))
    myWidths<-c()
    for (k in k.range){
        #print(k)
        myWidths<-c(myWidths,pam(n,k,...)$silinfo$avg.width)
    }
    plot(k.range,myWidths)
}

##------------------------------------------------------------------------##
##------------------------------------------------------------------------##

## Modifications of functions to compare groups of lists 
## (from http://stackoverflow.com/questions/23559371/how-to-get-the-list-of-items-in-venn-diagram-in-r)
Intersect <- function (x) {  
    # Multiple set version of intersect
    # x is a list
    if (length(x) == 1) {
        unlist(x)
    } else if (length(x) == 2) {
        intersect(x[[1]], x[[2]])
    } else if (length(x) > 2){
        intersect(x[[1]], Intersect(x[-1]))
    }
}
#
Union <- function (x) {  
    # Multiple set version of union
    # x is a list
    if (length(x) == 1) {
        unlist(x)
    } else if (length(x) == 2) {
        union(x[[1]], x[[2]])
    } else if (length(x) > 2) {
        union(x[[1]], Union(x[-1]))
    }
}
#
Setdiff <- function (x, y) {
    # Remove the union of the y's from the common x's. 
    # x and y are lists of characters.
    xx <- Intersect(x)
    yy <- Union(y)
    setdiff(xx, yy)
}

##------------------------------------------------------------------------##
##------------------------------------------------------------------------##

## Ouput the color IDs used by ggplot
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}

##------------------------------------------------------------------------##
##------------------------------------------------------------------------##

# TPM function
tpm <- function(counts, lengths) {
    rate <- counts / lengths
    rate / sum(rate) * 1e6
}

##------------------------------------------------------------------------##
##------------------------------------------------------------------------##

## Miscellaneous operators
'%!in%' <- function(x,y)!('%in%'(x,y))

##------------------------------------------------------------------------##
##------------------------------------------------------------------------##

theme_monokai_full <- function(base_size = 14, base_family = ""){
  color.background = "#232323"
  color.grid.major = "#232323"
  color.text = "#ffffff"
  color.axis = "#ffffff"

  theme_bw(base_size=base_size) +
    theme(

      panel.background=element_rect(fill=color.background, color=NA),
      plot.background=element_rect(fill=color.background, color=color.background),
      panel.border=element_rect(color=color.axis),

      panel.grid.major=element_line(color="grey95",size=.1, linetype=3),
      panel.grid.minor=element_blank(),
#       axis.line.x=element_line(color=color.grid.major, size=1),
#       axis.line.y=element_line(color=color.grid.major, size=1),
#       axis.ticks=element_line(color=NA),

      legend.background = element_rect(fill=color.background),
      legend.key = element_rect(fill=color.background, color=NA),
      legend.text = element_text(size=rel(.8),color=color.text),#color.axis.title),
      legend.title = element_text(color=color.text),

      plot.title=element_text(color=color.text, size=rel(1.2)),
      plot.subtitle=element_text(color=color.text, size=rel(0.8)),
      axis.text.x=element_text(size=rel(.95),color=color.text),
      axis.text.y=element_text(size=rel(.95),color=color.text),
      axis.title.x=element_text(size=rel(1),color=color.text, vjust=0),
      axis.title.y=element_text(size=rel(1),color=color.text, vjust=1.25),
      strip.background = element_rect(fill = color.background, colour = color.axis),
      strip.text = element_text(colour = color.text)
    )

}

###--------------------------------#####
theme_black_full <- function(base_size = 14, base_family = ""){
  color.background = "#000000"
  color.grid.major = "#000000"
  color.text = "#ffffff"
  color.axis = "#ffffff"

  theme_bw(base_size=base_size) +
    theme(

      panel.background=element_rect(fill=color.background, color=NA),
      plot.background=element_rect(fill=color.background, color=color.background),
      panel.border=element_rect(color=color.axis),

      panel.grid.major=element_line(color="grey95",size=.1, linetype=3),
      panel.grid.minor=element_blank(),
#       axis.line.x=element_line(color=color.grid.major, size=1),
#       axis.line.y=element_line(color=color.grid.major, size=1),
#       axis.ticks=element_line(color=NA),

      legend.background = element_rect(fill=color.background),
      legend.key = element_rect(fill=color.background, color=NA),
      legend.text = element_text(size=rel(.8),color=color.text),#color.axis.title),
      legend.title = element_text(color=color.text),

      plot.title=element_text(color=color.text, size=rel(1.2)),
      plot.subtitle=element_text(color=color.text, size=rel(0.8)),
      axis.text.x=element_text(size=rel(.95),color=color.text),
      axis.text.y=element_text(size=rel(.95),color=color.text),
      axis.title.x=element_text(size=rel(1),color=color.text, vjust=0),
      axis.title.y=element_text(size=rel(1),color=color.text, vjust=1.25),
      strip.background = element_rect(fill = color.background, colour = color.axis),
      strip.text = element_text(colour = color.text)
    )

}


plot.qq <- function(vec, title.str="qqplot", hit.idx=NULL) {
  col.vec <- rep(1, length(vec))
  if (length(hit.idx) > 0) {
    col.vec[1:length(hit.idx)] <- 2
    col.vec <- rev(col.vec)
  }	
  qqplot(-log10(ppoints(length(vec))), -log10(vec), ylab="-log10(observed)", xlab="-log10(expected)", main=title.str, col=col.vec)
  abline(a=0, b=1)
}

ggplot.qq<- function(vec, hit.idx=NULL) {
    col.vec <- rep(1, length(vec))
  if (length(hit.idx) > 0) {
    col.vec[1:length(hit.idx)] <- 2
    col.vec <- rev(col.vec)
  } 
    matrix <- data.frame(expected = -log10(ppoints(length(vec))), observed = -log10(vec), sigLabel = rev(col.vec))
    matrix$sig = ifelse(matrix$sigLabel == "1", "no", "yes")
    p <- ggplot(matrix, aes(expected, observed, colour = sig)) + 
        geom_point(size = 2, alpha = 0.75) + 
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "white") + 
        theme(legend.position = "none") +
        scale_colour_manual(values = c("gray", "red"))
    return(p)
}