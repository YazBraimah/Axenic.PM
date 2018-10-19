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

geneBoxPlot <- function (gene, show_reps = F) {
        if (grepl("FBgn", gene)) {
            description <- subset(snapshots, FBgn_ID == gene)$GeneName
            geneName <- subset(tpm.table, gene_id == gene)$gene_symbol
            cg_id <- subset(tpm.table, gene_id == gene)$annotation_ID
            
            p <- ggplot(subset(tpm.table, gene_id == gene), aes(Male, TPM, fill = Status)) +
                geom_boxplot(outlier.size = 0, width = 0.3, alpha = 0.5) +
                geom_point(pch = 21, position = position_jitterdodge(), aes(fill = Status)) + 
                facet_grid(.~Female, scales = "free_x", space = "free_x") +
                labs(title = paste(gene, " (", geneName, ")", sep = ""), subtitle = paste(cg_id, ": ", description, sep = "")) + 
                theme_bw() + 
                theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
                    axis.text.y = element_text(size = 14),
                    strip.text = element_text(size = 14, face = "bold")) +
                scale_fill_manual(values = c("#3f5a2a","#ffb200"))
        } else {
            fbgn_id <- subset(tpm.table, gene_symbol == gene)$gene_id
            description <- subset(snapshots, GeneSymbol == gene)$GeneName
            cg_id <- subset(tpm.table, gene_symbol == gene)$annotation_ID
            
            p <- ggplot(subset(tpm.table, gene_symbol == gene), aes(Male, TPM, fill = Status)) +
                geom_boxplot(outlier.size = 0, width = 0.3, alpha = 0.5) +
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

 edge.DE <- function(fit, contrast) {
    lrt <- glmLRT(fit, contrast = contrast)
    lrt.tTags <- topTags(lrt, n = NULL)
    lrt.table <- lrt.tTags$table
    lrt.table$sig = ifelse(lrt.table$FDR < 0.01 & (lrt.table$logFC > 1 | lrt.table$logFC < -1), "yes", "no")
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
                    size = 3) +
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
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) + 
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    geom_abline(slope = 1, linetype = "dashed", alpha = 0.5) +
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
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    geom_point(alpha = 0.5) +
    theme_bw()
    return(p)
}

##------------------------------------------------------------------------##
##------------------------------------------------------------------------##

## extract gene IDs based on GO term and enriched gene list:
extract_GO_genes = function(go_term, gene_set, trinity = FALSE){
    if (trinity){
        rownames(subset(GOinfo_pasa, row.names(GOinfo_pasa) %in% gene_set & grepl(go_term, GOinfo_pasa$V2)))
    } else {
    rownames(subset(GOinfo_annotated, row.names(GOinfo_annotated) %in% gene_set & grepl(go_term, GOinfo_annotated$V2)))
    }
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
