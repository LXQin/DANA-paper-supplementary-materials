### Include Packages -----------------------------------------------------------

# Plots
library(ggplot2)  # plots
library(gridExtra)  # Layout for ggplot2
library(ggnewscale)  # two scales for one plot
library(corrplot)  # correlation plots
library(stargazer)  # ascii table
library(plotly)  # interactive plots
library(ggrepel)  # repelling labels in scatter plots
# Precision Estimation
library(glmnet)  # lasso estimator
library(huge)  # precision estimation 
library(Rcpp)  # Requirement for FastGGM
library(FastGGM)  # preciison estimation method
# Normalization Methods
library(edgeR) 
library(DESeq)
library(PoissonSeq)
library(sva)
library(RUVSeq)
library(vsn)  # variance stabilization (array data)
# Statistics
library(DescTools) # concordance correlation coeff (CCC)
library(ffpe)  # used for CATplot function


### Auxiliary Functions --------------------------------------------------------




#' Function for building a coordinates data frame for miRNA genome locations
#' based on miRNA data from miRBase. MiRNAs with ambiguous
#' location on the genome are excluded from the data frame
#' @param file Path to the miRBase database file. E.g. "hsa.gff3"
#' @return Coordinate data frame that relates each miRNA to its location on the 
#' genome.
mirbase.coords <- function(file){
  # read the coordinate file
  coordData <- read.table(file, header=FALSE, stringsAsFactors = FALSE, comment.char="#")
  # keep only essential (needed) data columns
  coordData <- coordData[, c(3, 1, 4, 5, 9)]
  names(coordData) <- c("type", "chr", "start", "end", "info")
  # The info column contains ID, Alias, and name of the miRNA
  info <- coordData[["info"]]  # for simplicity
  
  # extract information from the info column
  idRead <- character(length(info))
  aliasRead <- character(length(info))
  nameRead <- character(length(info))
  # read all information from the info string
  for (i in 1:length(info)) {
    # information form the info string
    line <- read.csv(text=info[i], sep=";", stringsAsFactors = FALSE, 
                     header=FALSE)
    # ignore the identifier part of the information (e.g. "ID=")
    idRead[i] <- substr(line[1], 4, 512)
    aliasRead[i] <- substr(line[2], 7, 512)
    nameRead[i] <- tolower(substr(line[3], 6, 512))
  }
  # add information to coordData and remove old info column
  coords <- cbind(coordData[1:4], id=idRead, name=nameRead)
  # Set the names of the entries of the data frame to the unique accession ids
  rownames(coords) <- coords$id
  rm(coordData, line, aliasRead, idRead, info, nameRead, i) # environment clean-up
  coords$id <- as.character(coords$id)
  coords$name <- tolower(as.character(coords$name))
  coords$type <- tolower(as.character(coords$type))
  
  # find all duplicates
  duplicated.genes <- unique(coords$name[which(duplicated(coords$name))])
  # l <- lapply(duplicated.genes, function(a) coords[which(coords$name==a), ])
  keep <- rep(TRUE, dim(coords)[1])
  for(dg in duplicated.genes){
    # if the duplication is caused by a pair of miRNA and miRNA_primary_transcript
    # then keep both entries, otherwiese remove from the data frame
    pos <- which(coords$name==dg)
    if(length(pos) == 2 && (coords$type[pos[1]] != coords$type[pos[2]])) {
      # keep mirna location
      keep[pos[coords$type[pos]!="mirna"]] <- FALSE
    } else {
      # multiples elements of the same miRNA or miRNA_primary_transcript => Remove
      keep[pos] <- FALSE
    }
  }
  coords <- coords[keep, ]
  rownames(coords) <- coords$name
  return(coords)
}



#' Computes partial correlation matrix from a given precision matrix
#' @param P A square precision matrix
#' @return Partial correlation matrix computed from P
prec2corr <- function(P) {
  
  # Test if matrix is square
  if(dim(P)[1] != dim(P)[2]) {
    stop("Error: Function prec2cor expected a square matrix.")
  }
  
  n <- dim(P)[1]  # dimension of matrix
  C <-matrix(1, nrow=n, ncol=n)
  
  # compute the partial correlations
  for (i in 1:n) {
    for (j in 1:n) {
      C[i,j] <- -P[i,j] / sqrt(P[i,i] * P[j,j])
    }
  }
  
  # catch any correlations >1 or <-1 and cap the values
  C[C > 1] <- 1
  C[C < (-1)] <- -1
  
  rownames(C) <- rownames(P)
  colnames(C) <- colnames(P)
  
  return(C)
}



#' Sorts a given set of genes by chromosomes and by gene location on 
#' each chromosome. Genes that are not found in coords are appended at the end.
#' @param genes List of gene names
#' @param coords Coordinates data.frame that specifies the chromosome and 
#' base-pair start location of each gene.
#' @return List of sorted genes.
sort.on.chromosome <- function(genes, coords) {
  ## Sort by chromosome
  chromosomes <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", 
                   "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", 
                   "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", 
                   "chr21","chr22", "chrX", "chrY")
  genes.sorted <- genes[order(match(coords[genes, "chr"], chromosomes))]
  
  ## Sort by base-pair location on each chromosome
  for(chr in chromosomes) {
    idx <- which(coords[genes.sorted, "chr"] == chr)
    genes.sorted[idx] <- genes.sorted[idx][order(coords[genes.sorted[idx], "start"])]
  }
  
  return(genes.sorted)
}



#' computes the graph corresponding to a given precision or correlation matrix
#' @param X precision/correlation matrix
#' @param sym Symmetrize the output graphs, with the "and" or "or" rule. 
#' Default is "or".
#' @returns graph connectivity/adjacency matrix
to.graph <- function(X, sym="or") {
  res <- NULL
  if(sym=="or") {
    res <- 1*(X | (t(X)))
  } else if(sym == "and") {
    res <- 1*(X & (t(X)))
  } else {
    stop("Error: Unknown value for parameter sym.")
  }
  colnames(res) <- colnames(X)
  rownames(res) <- rownames(X)
  return(res)
}




#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
#' Extract the legend (object) of a given ggplot plot object
#' @param a.gplot A ggplot object
#' @return ggplot legend object
get.legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}




### Plot Functions -------------------------------------------------------------



#' Mean-Standard Deviation Plot
#'
#' Plots the marker-specific standard deviation over marker-specific means
#' for log2-normalized read counts. Bounds for the definition of positive
#' and negative controls can be highlighted by vertical lines (if given).
#'
#' @param raw Raw read count matrix (rows = genes, cols = samples).
#' @param tZero \emph{Optional.} Lower bound for negative controls.
#' @param tPoor \emph{Optional.} Upper bound for negative controls.
#' @param tWell \emph{Optional.} Lower bound for positive controls.
#' @param title Plot title
#' @param xlim.max Upper limit of the x-axis
#' @param ylim.max Upper limit of the y-axis
#'
#' @return ggplot object of the mean-sd plot
#' @export
#'
#' @examples
#' counts <- matrix(rnbinom(100, mu=3, size=0.01), nrow=20)
#' plot.mean.sd(counts)
plot.mean.sd <- function(raw, tZero, tPoor, tWell, title, xlim.max=NA, ylim.max=NA) {
  df <- data.frame(data.mean= rowMeans(log2(raw + 1)),
                   data.var = apply(log2(raw + 1), 1, stats::sd))
  
  if(missing(tZero)) {
    vline.t.zero <- ggplot2::geom_blank()
  } else {
    vline.t.zero <- ggplot2::geom_vline(xintercept = log2(tZero+1), color="#5851b8", linetype="twodash")
  }
  if(missing(tPoor)) {
    vline.t.poor <- ggplot2::geom_blank()
  } else {
    vline.t.poor <- ggplot2::geom_vline(xintercept = log2(tPoor+1), color="#5851b8", linetype="twodash")
  }
  if(missing(tWell)) {
    vline.t.well <- ggplot2::geom_blank()
  } else {
    vline.t.well <- ggplot2::geom_vline(xintercept = log2(tWell+1), color="#E7298A",  linetype="longdash")
  }
  
  if(missing(title)) {
    title <- ggplot2::geom_blank()
  } else {
    title <- ggplot2::ggtitle(title)
  }
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x=data.mean,y=data.var)) +
    ggplot2::geom_point(alpha=.25) +
    ggplot2::xlab("Mean") +
    ggplot2::ylab("Standard deviation") +
    vline.t.zero +
    vline.t.poor +
    vline.t.well +
    title +
    ggplot2::xlim(0, xlim.max) +
    ggplot2::ylim(0, ylim.max) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::theme_classic()
  return(p)
}



#' Read count histogram
#'
#' Plots a histogram of the marker-specific mean of log2-transformed
#' read counts. Bounds for the definition of positive
#' and negative controls can be highlighted by vertical lines (if given).
#'
#' @param raw Raw read count matrix (rows = genes, cols = samples).
#' @param tZero \emph{Optional.} Lower bound for negative controls.
#' @param tPoor \emph{Optional.} Upper bound for negative controls.
#' @param tWell \emph{Optional.} Lower bound for positive controls.
#' @param title Plot title
#' @param xlim.max Upper limit of the x-axis
#'
#' @return ggplot object of the read count histogram
#' @export
#'
#' @examples
#' counts <- matrix(rnbinom(100, mu=3, size=0.01), nrow=20)
#' plot.countHist(counts)
plot.countHist <- function(raw, binwidth=0.1, tZero, tPoor, tWell, title, xlim.max=NA, ylim.max=NA) {
  if(missing(tZero)) {
    vline.t.zero <- ggplot2::geom_blank()
  } else {
    vline.t.zero <- ggplot2::geom_vline(xintercept = log2(tZero+1), color="#5851b8", linetype="twodash")
  }
  if(missing(tPoor)) {
    vline.t.poor <- ggplot2::geom_blank()
  } else {
    vline.t.poor <- ggplot2::geom_vline(xintercept = log2(tPoor+1), color="#5851b8", linetype="twodash")
  }
  if(missing(tWell)) {
    vline.t.well <- ggplot2::geom_blank()
  } else {
    vline.t.well <- ggplot2::geom_vline(xintercept = log2(tWell+1), color="#E7298A",  linetype="longdash")
  }
  
  if(missing(title)) {
    title <- ggplot2::geom_blank()
  } else {
    title <- ggplot2::ggtitle(title)
  }
  
  df <- data.frame(log.expression = log2(rowMeans(raw) + 1))
  p <- ggplot2::ggplot(df, ggplot2::aes(x=log.expression)) +
    ggplot2::geom_histogram(binwidth = binwidth, color="black", fill="black") +
    ggplot2::xlab("Mean") +
    ggplot2::ylab("Frequency") +
    vline.t.zero +
    vline.t.poor +
    vline.t.well +
    title +
    ggplot2::xlim(0, xlim.max) +
    ggplot2::ylim(0, ylim.max) +
    ggplot2::theme_classic()
  return(p)
}



#' Scatter plot DANA assessment metrics
#'
plot.DANA.metrics <- function(metrics, label.size=3, label.repel=FALSE) {
  if (label.repel) {
    if (!requireNamespace("ggrepel", quietly = TRUE)) {
      stop("Package \"ggrepel\" needed for this function to work. Please install it.",
           call. = FALSE
      )
    }
  }
  
  scipen.prev <- getOption("scipen")
  options(scipen=4) # force non-scientific notation of x axis
  
  df <- data.frame(method=rownames(metrics),
                   cc = metrics$cc,
                   mcr = metrics$mcr)
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x=mcr, y=cc, label=method)) +
    ggplot2::geom_point(alpha=.75)
  
  if(label.repel) {
    p <- p + ggrepel::geom_text_repel(ggplot2::aes(label = method), size=label.size)
  } else {
    p <- p + ggplot2::geom_text(size=label.size, hjust=0, vjust=0)
  }
  
  p <- p +
    ggplot2::theme_classic() +
    ggplot2::xlab("Relative reduction of handling effects") +
    ggplot2::ylab("Biological signal preservation") +
    # ggplot2::geom_text_repel(aes(label = method), size=3) +
    ggplot2::scale_x_continuous(labels = scales::percent, limits = c(0,1)) +
    ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0,1))
  
  options(scipen=scipen.prev)  # reset scipen option
  return(p)
}



#' Plot a heatmap graph for a given correlation matrix using ggplot. 
#' miRNA clusters are highlighted.
#' @param corr A correlation matrix
#' @param clusters Clustering structure for the genes in corr
#' @param title Title of the plot. Default is title="", which gives no title
#' @param coords Coordinates data frame. If given, chromosomes are highlighted.
#' @param threshold Only values of the correlation matrix greater than the 
#' threshold in absolute value are plotted. Default is corr.threshold=0 which
#' corresponds to no thresholding.
#' @param subset If a value is given, plots only a given subset of the 
#' correlations. Subsets are given as an array of indexes.
#' @param graph Logical whether to only plot graph structure, 
#' set all correlations to 1. 
ggplot.corr <- function(corr, clusters, title, coords=NULL, threshold=FALSE, subset=FALSE, graph=FALSE){
  # reduce plot to a subset of the data
  if(any(subset)) {
    corr <- corr[subset, subset]
  }
  
  # Threshold matrix
  if(threshold) {
    corr[abs(corr) < threshold] <- 0
  }
  
  # Convert to graph 
  if(graph) {
    corr <- 1*(corr != 0)
  }
  
  if(missing(title)) {
    title <- ggplot2::geom_blank()
  } else {
    title <- ggplot2::ggtitle(title)
  }
  
  # Set corr below diagonal to NA
  corr[lower.tri(corr)] <- NA
  # Set diagonal to NA
  diag(corr) <- NA
  
  ## Plot using ggplot
  data <- expand.grid(X=rownames(corr), Y=colnames(corr))[upper.tri(corr, diag=FALSE),]
  data$Z <- as.vector(corr[upper.tri(corr, diag=FALSE)])
  
  clustered <- corr
  clustered[] <- 0 # Set all to false
  
  if(is.null(coords)) {
    diag(clustered) <- 1
    # color the miRNA clusters below the diagonal
    for (clust in clusters) {
      # find genes in clust in corr
      idx <- na.omit(match(unlist(clust), colnames(clustered)))
      # set edges in clusters to bright color in the corrplot
      clustered[idx, idx][lower.tri(clustered[idx,idx])] <- 1
    }
  } else {
    # coordinates data frame given -> highlight chromosomes
    genes <- rownames(clustered)
    chromosomes <- data.frame(chr=unique(coords[genes, "chr"]))
    chromosomes$val <- rep(c(1,2), ceiling(length(chromosomes$chr)/2))[1:length(chromosomes$chr)]
    rownames(chromosomes) <- chromosomes$chr
    for(i in 1:length(chromosomes$chr)) {
      idx <- which(coords[genes, "chr"] == chromosomes$chr[i])
      diag(clustered)[idx] <- chromosomes$val[i]
    }
    for (clust in clusters) {
      # find genes present in cluster clust in clustered
      idx <- na.omit(match(unlist(clust), colnames(clustered)))
      # set edges in clusters to specific color in the corrplot
      clustered[idx, idx][lower.tri(clustered[idx,idx])] <- chromosomes[coords[rownames(clustered[idx,idx])[1], "chr"], "val"]
    }
  }
  
  data.clust <- expand.grid(X=rownames(clustered), Y=colnames(clustered))[lower.tri(clustered, diag=TRUE),]
  data.clust$Z <- as.vector(clustered[lower.tri(clustered, diag=TRUE)])
  
  
  # data <- data %>% mutate(text = paste0("x: ", X, "\n", "y: ", Y, "\n", "Cor: ",round(Z,2)))
  p <- ggplot() + 
    geom_tile(data=data.clust, aes(X, Y, fill=factor(Z))) +
    scale_fill_manual(name="Clusters",
                      # values = c("0"="gray95", "1"="black", "2"="darkslategray4"), 
                      values = c("0"="white", "1"="black", "2"="black"), 
                      labels=c("un-clustered", "clustered", "clustered"),
                      guide=FALSE) +
    title +
    new_scale_fill() +
    geom_tile(data=data, aes(X, Y, fill= Z)) +
    scale_x_discrete(limits = rev(levels(as.factor(data$X)))) +
    # scale_fill_gradient2(name = "Cor", midpoint=0, low = "blue", mid = "white",
    #                      high = "red", limit = c(-1, 1), space = "Lab") +
    scale_fill_gradientn(colours = c("blue4", "blue", "white", "red", " red4"),
                         name = "Cor", limit = c(-1, 1), space = "Lab") +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank()) +
    theme(aspect.ratio=1) + 
    theme(panel.border = element_rect(colour = "gray70", fill=NA, size=.2, linetype="dashed"))
  return(p)
}



#' Plot an interactive heatmap for a given correlation matrix using ggplot. 
#' miRNA clusters are highlighted.
#' @param corr A correlation matrix
#' @param clusters Clustering structure for the genes in corr
#' @param title Title of the plot. Default is title="", which gives no title
#' @param coords Coordinates data frame. If given, chromosomes are highlighted
#' and information is displayed as mouse over text.
#' @param threshold Only values of the correlation matrix greater than the 
#' threshold in absolute value are plotted. Default is corr.threshold=0 which
#' corresponds to no thresholding.
#' @param subset If a value is given, plots only a given subset of the 
#' correlations. Subsets are given as an array of indexes.
#' @param graph Logical whether to only plot graph structure, 
#' set all correlations to 1. 
iggplot.corr <- function(corr, clusters, title="", coords=NULL, threshold=FALSE, subset=FALSE, graph=FALSE){
  # reduce plot to a subset of the data
  if(any(subset)) {
    corr <- corr[subset, subset]
  }
  
  # Threshold matrix
  if(threshold) {
    corr[abs(corr) < threshold] <- 0
  }
  
  # Convert to graph 
  if(graph) {
    corr <- 1*(corr != 0)
  }
  
  ## Value Guide
  # [-1,1] correlations
  # -2 Diagonal/cluster 1
  # -3 Clustered
  # -4 Un-Clustered
  # -5 Diagonal/cluster 2
  
  # Set corr below diagonal to -4 (un-clustered)
  corr[upper.tri(corr)] <- -4
  if(!is.null(coords)) {
    # coordinates data frame given -> highlight chromosomes
    genes <- rownames(corr)
    chromosomes <- data.frame(chr=unique(coords[genes, "chr"]))
    chromosomes$val <- rep(c(-2,-5), ceiling(length(chromosomes$chr)/2))[1:length(chromosomes$chr)]
    rownames(chromosomes) <- chromosomes$chr
    for(i in 1:length(chromosomes$chr)) {
      idx <- which(coords[genes, "chr"] == chromosomes$chr[i])
      diag(corr)[idx] <- chromosomes$val[i]
    }
    for (clust in clusters) {
      # find genes present in cluster clust in corr
      idx <- na.omit(match(unlist(clust), colnames(corr)))
      # set edges in clusters to specific color in the corrplot
      corr[idx, idx][upper.tri(corr[idx,idx])] <- chromosomes[coords[rownames(corr[idx,idx])[1], "chr"], "val"]
    }
  } else { 
    # Set diagonal to -2
    diag(corr) <- -2# color the miRNA clusters below the diagonal
    for (clust in clusters) {
      # find genes present in cluster clust in corr
      idx <- na.omit(match(unlist(clust), colnames(corr)))
      # set edges in clusters to bright color in the corrplot
      corr[idx, idx][upper.tri(corr[idx,idx])] <- -5
    }
  }
  
  
  
  ## Plot using ggplot
  data <- expand.grid(X=rownames(corr), Y=colnames(corr))
  data$Z <- as.vector(corr)
  
  ## Text for tooltips in interactive plot
  text <- corr
  if(is.null(coords)) {  # do not show chromosome
    text[lower.tri(text)] <- as.matrix(paste0("x: ", data$X, "\ny: ", data$Y, "\nCor: ",round(data$Z,3)), nrow=dim(text)[1])[lower.tri(text)]
    text[upper.tri(text)] <- as.matrix(paste0("x: ", data$X, "\ny: ", data$Y, "\n", ifelse(corr!=-4, "clustered", "un-clustered")), nrow=dim(text)[1])[upper.tri(text)]
    diag(text) <- paste0("x: ", rownames(text), "\ny: ", colnames(text))
  } else {  # show chromosome
    text[lower.tri(text)] <- as.matrix(paste0("x: ", coords[as.character(data$X), "chr"], " ", data$X, "\ny: ", coords[as.character(data$Y), "chr"], " ", data$Y, "\nCor: ",round(data$Z,3)), nrow=dim(text)[1])[lower.tri(text)]
    text[upper.tri(text)] <- as.matrix(paste0("x: ", coords[as.character(data$X), "chr"], " ", data$X, "\ny: ", coords[as.character(data$Y), "chr"], " ", data$Y, "\n", ifelse(corr!=-4, "clustered", "un-clustered")), nrow=dim(text)[1])[upper.tri(text)]
    diag(text) <- paste0("x: ", coords[rownames(text), "chr"], " ", rownames(text), "\ny: ", coords[rownames(text), "chr"], " ", colnames(text))
  }
  
  data$text <- as.vector(text)
  
  
  # data <- data %>% mutate(text = paste0("x: ", X, "\n", "y: ", Y, "\n", "Cor: ",round(Z,2)))
  p <- ggplot(data, aes(X, Y, fill= Z, text=text)) + 
    geom_tile() +
    # scale_x_discrete(limits = rev(levels(as.factor(data$X)))) +
    scale_y_discrete(limits = rev(levels(as.factor(data$Y)))) +
    scale_fill_gradientn(colours = c("black","grey95", "green", "darkslategray4", "blue", "white", "red"), limits = c(-5, 1)) +
    ggtitle(title) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          legend.position="none")
  p <- ggplotly(p, tooltip="text", width=600, height=600)
  return(p)
}




#' Plot a scatter plot comparison of two correlations
#' @param corr1 First correlation matrix
#' @param corr2 Second correlation matrix
#' @param clusters MiRNA clustering for the genes. Scatter plot will be reduced
#' to the members of miRNA clusters, if clusters is non-NULL.
#' @param title Title of the plot. Default is title="", which gives no title
#' @param xlab x-Axis label
#' @param ylab y-Axis label
#' @param subset If a value is given, plots only a given subset of the 
#' precisions Subsets are given as an array on indices.
#' @param tri Specifies if lower ("lower") or upper ("upper") triangle 
#' of the correlations is used. Default is tri="upper"
#' @param corr.threshold Only values of the correlation matrix greater than the 
#' threshold in absolute value are plotted. Default is corr.threshold=0 which
#' corresponds to no thresholding
#' @param ccc Concondance correlation coefficient. Is displayed on the lower 
#' right of the plot. This function computes the ccc between strictly positive 
#' correlations if no value is given.
#' @return Plot object for the scatter plot
plot.corr.scatter <- function (corr1, corr2, clusters=NULL, title="", xlab="", 
                               ylab="", subset=FALSE, tri="upper", threshold=0,
                               ccc=NULL) {
  if(any(dim(corr1) != dim(corr2))) {
    warning("Warning: Shapes of correlation matrices in plot.corr.scatter do not match! 
            Reducing correlation matrices to all mutual genes.\n")
    genes <- intersect(rownames(corr1), rownames(corr2))
    corr1 <- corr1[genes, genes]
    corr2 <- corr2[genes, genes]
  }
  if((tri!="lower")&&(tri!="upper")) {
    stop("Error: Unknown value of parameter tri in plot.corr.scatter")
  }
  
  # reduce plot to a subset of the data
  if(any(subset)) {
    precision <- precision[subset, subset]
  }
  
  # Threshold matrix
  if(threshold) {
    corr1[abs(corr1) < threshold] <- 0
    corr2[abs(corr2) < threshold] <- 0
  }
  
  # values for the plot
  valuesX <- c()
  valuesY <- c()
  
  ## Scatter plot of the full graph (within and in-between clusters)
  if(is.null(clusters)) {
    # indices of elements to consider
    idx <- corr1 | corr2
    if(tri=="lower") {
      idx <- idx & lower.tri(corr1)
    } else {  # case "upper"
      idx <- idx & upper.tri(corr1)
    }
    valuesX <- as.vector(corr1[idx])
    valuesY <- as.vector(corr2[idx])
    ## Scatter plot of within clusters co-expressions only
  } else {
    for(clust in clusters) {
      # find genes in clust in corr
      clust.genes <- colnames(corr1)[na.omit(match(unlist(clust), colnames(corr1)))]
      # elements that are present in any of the two part. correlation matrices
      idx <- corr1[clust.genes, clust.genes] | corr2[clust.genes, clust.genes]
      if(tri=="lower") {
        idx <- idx & lower.tri(idx)
      } else {  # case "upper"
        idx <- idx & upper.tri(idx)
      }
      valuesX <- c(valuesX, corr1[clust.genes, clust.genes][idx])
      valuesY <- c(valuesY, corr2[clust.genes, clust.genes][idx])
    }
  }
  
  # concordance correlation coefficient for the partial correlations
  if(is.null(ccc)) {
    correlation <- CCC(valuesX, valuesY)$rho.c$est
  } else {
    correlation <- ccc    
  }
  
  df <- data.frame(x=valuesX, y=valuesY)
  
  p <- ggplot(df,aes(x=x,y=y)) + 
    geom_point(alpha=.4) +
    xlim(ifelse(is.null(clusters), -1, 0),1) +
    ylim(ifelse(is.null(clusters), -1, 0),1) +
    labs(title=title, x=xlab, y=ylab) + 
    geom_abline(intercept = 0, slope = 1) +
    annotate("label", alpha=0.5, x=1, y=ifelse(is.null(clusters), -1, 0.02), hjust=1, vjust=0, 
             label=paste('cc =',round(correlation,3))) + 
    theme(aspect.ratio=1) +
    theme_classic()
  return(p)
}




#' Scatter plot for essential normPA result metrics
#' 
plot.result.statistics <- function(data.mposcor, norm.mposcor, norm.ccocorr, norm.names) {
  options(scipen=4) # force non-scientific notation of x axis
  library(ggplot2)
  library(ggrepel)
  
  # data frame for plot
  mposcorr.reduction <- (-norm.mposcor + data.mposcor) / data.mposcor
  df <- data.frame(x=mposcorr.reduction, 
                   y=norm.ccocorr,
                   label=norm.names)
  
  p.stats <- ggplot(df,aes(x=x,y=y, label=label)) + 
    geom_point(alpha=1) + 
    xlab("Rel. Noise Reduction (mcr-)") + 
    ylab("Biological Signal Preservation (cc+)") + 
    geom_text_repel(aes(label = label), size=3) +
    # geom_text(aes(label=label), size=3, nudge_y=0.02) +
    # xlim(c(min(c(0,mposcorr.reduction)),1.1*max(mposcorr.reduction))) +
    # scale_y_continuous(limits=c(0,1.05),breaks=c(0,0.25,0.5,0.75,1.0), labels=c(0,0.25,0.50,0.75,1)) +
    scale_x_continuous(labels = scales::percent, limits = c(min(0,mposcorr.reduction),1),
                       breaks = scales::pretty_breaks(n = 6)) +
    scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
    theme(aspect.ratio=1)
  return(p.stats)
}



#' Volcano plot generated for the results of a DEA
#'
#' @param top.table Resulting top.table of a differential expression analysis 
#' using voom. 
#' @param alpha p-value cutoff (significance level). 
#' @param q.values Logical if q.values (adjusted p-values using the 
#' Benjamini-Hochberg method) should be used instead of p-values. 
plot.DE.volcano <- function(top.table, alpha, q.values=FALSE, title="") {
  if(q.values){  # use q values (adjusted p-values, Benjamini-Hochberg)
    DEA.res <- data.frame(cbind(top.table$logFC, top.table$adj.P.Val))
    colnames(DEA.res) <- c("FC.log2", "q.value")
    DEA.res <- DEA.res[order(DEA.res[, "q.value"]), ]
    cat("number of DE genes:", sum(DEA.res[, "q.value"] < alpha), "\n")
    
    ## volcano plot for fold change and p-values
    par(mfrow=c(1,1))
    mask <- with(DEA.res, q.value < alpha)
    cols <- ifelse(mask,"red", "black")
    with(DEA.res, plot(FC.log2, -log10(q.value), cex = .5, pch = 16,col = cols, 
                       xlim = 1.05*c(-max(abs(FC.log2)),max(abs(FC.log2))),
                       ylim = 1.05*c(0, max(-log10(q.value))),
                       xlab = "Mean Difference", ylab = "-log10(q-value)",
                       main = title))
    abline(h = -log10(alpha), lty = 2)
  } else {  # use non-adjusted p values
    DEA.res <- data.frame(cbind(top.table$logFC, top.table$P.Value))
    colnames(DEA.res) <- c("FC.log2", "p.value")
    DEA.res <- DEA.res[order(DEA.res[, "p.value"]), ]
    cat("number of DE genes:", sum(DEA.res[, "p.value"] < alpha), "\n")
    
    ## volcano plot for fold change and p-values
    par(mfrow=c(1,1))
    mask <- with(DEA.res, p.value < alpha)
    cols <- ifelse(mask,"red", "black")
    with(DEA.res, plot(FC.log2, -log10(p.value), cex = .5, pch = 16,col = cols, 
                       xlim = 1.05*c(-max(abs(FC.log2)),max(abs(FC.log2))),
                       ylim = 1.05*c(0, max(-log10(p.value))),
                       xlab = "Mean Difference", ylab = "-log10(p-value)",
                       main = title))
    abline(h = -log10(alpha), lty = 2)
  }
}



#' Comcordance at the top plot
#'
#' Function for generating concordance at the top plot, which compares the 
#' concordance of the p-values for given differential expression analyzes 
#' to an assumed truth (a benchmark).
#'
#' @param DEA result for the data under study
#' @param truth.DEA DEA result for the assumed truth (the gold standard)
#' @param title Plot title
#' @param maxrank Optionally specify the maximum size of top-ranked items 
#' that you want to plot.
#' @param subset vector of a subset of genes/markers for this analysis
#'
#' @return a list of values about TPR, FPR, FDR, FNR
#'
plot.CAT <- function(DEA, truth.DEA, title="", maxrank=100, subset=NULL){
  # Subset DEAs to a given set of miRNAs
  if (!is.null(subset)) {
    DEA <- DEA[subset, ]
    truth.DEA <- truth.DEA[subset, ]
  }
  # Reduce Data to named p.values
  truth <- truth.DEA$P.Value
  names(truth) <- truth.DEA$genes
  compare <- list()
  for(i in 1:length(DEA)) {
    compare <- append(compare, list(DEA[[i]]$P.Value))
    names(compare[[i]]) <- DEA[[i]]$genes
  }
  names(compare) <- names(DEA)
  
  catplots <- list()
  for(i in 1:length(compare)) {
    catplots <- append(catplots, 
                       list(CATplot(compare[[i]], truth, 
                                    maxrank = maxrank, make.plot=FALSE)))
  }
  names(catplots) <- names(compare)
  
  data_cat <- data.frame(x = numeric(), y = numeric(), curve = character())
  for(i in 1:length(catplots)){
    y = catplots[[i]][,2]
    x = 1:length(y)
    data_cat = rbind(data_cat, data.frame(x, y, curve = names(catplots)[i]))
  }
  data_cat$curve = factor(data_cat$curve)
  
  return(ggplot(data_cat, aes(x, y, color = curve)) +
           theme(legend.title = element_blank()) +       
           geom_line(size=.75) +
           ylab("Rate of Agreement with Benchmark") +
           xlab("Significance Rank") +
           theme(legend.title=element_blank()) +
           ggtitle(title) +
           ylim(c(0,1)) + 
           theme_classic())
}



### Normalizaiton Functions ----------------------------------------------------



## Trimmed Mean of M-values (TMM)
#' https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
#' page 14
#' cpm() could compute normalized count per million based on DGE result
#' We compute cor() between cpm() and our normalized benchmark data and find
#' it to be 1
norm.TMM = function(raw, groups) {
  dat.DGE = DGEList(counts = matrix(raw, ncol = length(groups)),
                    group = factor(groups),
                    genes = rownames(raw))
  d = calcNormFactors(dat.DGE, method = "TMM")
  scaling.factor = d$samples$norm.factors * d$samples$lib.size / 1e6
  dat.normed = t(t(raw)/scaling.factor)
  # scaling.factor = d$samples$norm.factors * d$samples$lib.size /
  #   mean(dat.DGE$samples$lib.size)
  # dat.normed = round(t(t(raw)/scaling.factor))
  return(list(dat.normed = dat.normed,
              scaling.factor = scaling.factor,
              dge.normed = d))
}


## Total Count (TC)
#' https://academic.oup.com/bib/article/14/6/671/189645
norm.TC = function(raw, groups) {
  dat.DGE = DGEList(counts = matrix(raw, ncol = length(groups)),
                    group = factor(groups),
                    genes = rownames(raw))
  scaling.factor = dat.DGE$samples$lib.size/1e6
  dat.normed = t(t(raw)/scaling.factor)
  # scaling.factor <- dat.DGE$samples$lib.size / mean(dat.DGE$samples$lib.size)
  # dat.normed <- round(t(t(raw)/scaling.factor))
  return(list(dat.normed = dat.normed,
              scaling.factor = scaling.factor))
}


## Upper Quartile (UQ)

#' calcNormFactors could also obtain the UQ normalization factor
#' We compute cor() between its UQ and our UQ normalized benchmark data
#' and find it to be 1
norm.UQ = function(raw, groups) {
  dat.DGE = DGEList(counts = matrix(raw, ncol = length(groups)),
                    group = factor(groups),
                    genes = rownames(raw))
  q.factor = apply(dat.DGE$counts, 2, function(x) quantile(x[x != 0], probs = 0.75))
  scaling.factor = q.factor/1e6
  dat.normed = t(t(raw)/scaling.factor)
  # scaling.factor = q.factor/mean(dat.DGE$samples$lib.size)
  # dat.normed = round(t(t(raw)/scaling.factor))
  return(list(dat.normed = dat.normed,
              scaling.factor = scaling.factor))
}


## Median (Med)

#' https://academic.oup.com/bib/article/14/6/671/189645
norm.med = function(raw, groups) {
  dat.DGE = DGEList(counts = matrix(raw, ncol = length(groups)),
                    group = factor(groups),
                    genes = rownames(raw))
  m.factor = apply(dat.DGE$counts, 2, function(x) median(x[x != 0]))
  scaling.factor = m.factor/1e6
  dat.normed = t(t(raw)/scaling.factor)
  # scaling.factor = m.factor/mean(dat.DGE$samples$lib.size)
  # dat.normed = round(t(t(raw)/scaling.factor))
  return(list(dat.normed = dat.normed,
              scaling.factor = scaling.factor))
}


## DESeq

#' https://bioconductor.org/packages/release/bioc/vignettes/DESeq/inst/doc/DESeq.pdf
#' page 4
norm.DESeq = function(raw, groups) {
  condition = factor(groups)
  dat.DGE = estimateSizeFactors(newCountDataSet(round(raw), condition))
  scaling.factor = sizeFactors(dat.DGE)
  dat.normed = counts(dat.DGE, normalized = T)
  return(list(dat.normed = dat.normed,
              scaling.factor = scaling.factor))
}


## PoissonSeq

#' https://cran.r-project.org/web/packages/PoissonSeq/PoissonSeq.pdf
norm.PoissonSeq = function(raw, groups = NULL) {
  scaling.factor = PS.Est.Depth(raw)
  dat.normed = t(t(raw)/scaling.factor)
  return(list(dat.normed = dat.normed,
              scaling.factor = scaling.factor))
}


## Quantile Normalization (QN)

#' http://jtleek.com/genstats/inst/doc/02_05_normalization.html
norm.QN = function(raw, groups = NULL, filter = FALSE) {
  if (filter == TRUE) {
    raw = log2(raw + 1)
    raw = raw[rowMeans(raw) > 2, ]
  } else {
    raw = log2(raw + 1)
  }
  dat.log.normed = normalizeQuantiles(as.matrix(raw))
  dat.normed = 2^dat.log.normed - 1
  colnames(dat.normed) = colnames(raw)
  rownames(dat.normed) = rownames(raw)
  return(list(dat.normed = dat.normed))
}


## Remove Unwanted Variation (RUV)

#' http://jtleek.com/svaseq/zebrafish.html
#' http://www.bioconductor.org/packages/devel/bioc/vignettes/RUVSeq/inst/doc/RUVSeq.pdf
#' calibrator genes include tag "cali"
norm.RUV <- function(raw, groups, method = c("RUVg", "RUVs", "RUVr")) {
  filter <- apply(raw, 1, function(x) length(x[x > 5]) >= 2)
  dat.ruv <- raw[filter, ]
  genes <- rownames(dat.ruv)
  condition <- factor(groups)
  set <- newSeqExpressionSet(as.matrix(dat.ruv),
                             phenoData = data.frame(condition,
                                                    row.names = colnames(dat.ruv)))
  design <- model.matrix(~ condition, data = data.frame(condition,
                                                        row.names = colnames(dat.ruv)))
  y <- DGEList(counts = counts(set), group = condition)
  y <- calcNormFactors(y, method = "upperquartile")
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  fit <- glmFit(y, design)
  lrt <- glmLRT(fit, coef = 2)
  top <- topTags(lrt, n = nrow(set))$table
  spikes <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:0.15*nrow(raw)]))]
  
  if (method == "RUVg") {
    t <- RUVg(set, spikes, k = 1)
    dat.normed <- normCounts(t)
    return(list(dat.normed = dat.normed,
                adjust.factor = t$W))
  }else if (method == "RUVs") {
    differences <- makeGroups(condition)
    controls <- rownames(dat.ruv)
    t <- RUVs(set, controls, k = 1, differences)
    dat.normed <- normCounts(t)
    return(list(dat.normed = dat.normed,
                adjust.factor = t$W))
  }else if (method == "RUVr") {
    design <- model.matrix(~ condition, data = pData(set))
    y <- DGEList(counts = counts(set), group = condition)
    y <- calcNormFactors(y, method = "upperquartile")
    y <- estimateGLMCommonDisp(y, design)
    y <- estimateGLMTagwiseDisp(y, design)
    fit <- glmFit(y, design)
    res <- residuals(fit, type = "deviance")
    setUQ <- betweenLaneNormalization(set, which = "upper")
    controls <- rownames(dat.ruv)
    t <- RUVr(setUQ, controls, k = 1, res)
    dat.normed <- normCounts(t)
    
    return(list(dat.normed = dat.normed,
                adjust.factor = t$W))
  }
}





#' Applies normalization methods to given data set
#'
#' @details
#' Supported normalization methods are
#' \itemize{
#'   \item Total Count (TC)
#'   \item Upper Quartile (UQ)
#'   \item Median (med)
#'   \item Trimmed Median of Means (TMM)
#'   \item DESeq
#'   \item PoissonSeq
#'   \item Quantile Normalization (QN)
#'   \itme Remove Unwanted Variation (RUV)
#' }
#'
#' @param X Un-normalized data set of shape n x p, where n is the number of
#' samples and p is the number of markers
#' @param groups Vector of length n that maps the samples to sample groups
#' @param name Name of the data set
#' @param apply.TC Boolean if TC normalization is applied. Respective parameters
#' for all other methods.
normalize <- function(X, groups, name="X",
                      apply.TC=TRUE, apply.UQ=TRUE, apply.med=TRUE,
                      apply.TMM=TRUE, apply.DESeq=TRUE, apply.PoissonSeq=TRUE,
                      apply.QN=TRUE, apply.RUV=TRUE) {
  X.normalized <- list()
  
  if(apply.TMM) {
    X.normalized <- append(X.normalized, list(norm.TMM(X, groups)$dat.normed))
    names(X.normalized)[length(X.normalized)] <- paste(name, ".TMM", sep="")
  }
  if(apply.DESeq) {
    X.normalized <- append(X.normalized, list(norm.DESeq(X, groups)$dat.normed))
    names(X.normalized)[length(X.normalized)] <- paste(name, ".DESeq", sep="")
  }
  
  if(apply.PoissonSeq) {
    X.normalized <- append(X.normalized, list(norm.PoissonSeq(X, groups)$dat.normed))
    names(X.normalized)[length(X.normalized)] <- paste(name, ".PoissonSeq", sep="")
  }
  if(apply.TC) {
    X.normalized <- append(X.normalized, list(norm.TC(X, groups)$dat.normed))
    names(X.normalized)[length(X.normalized)] <- paste(name, ".TC", sep="")
  }
  if(apply.med) {
    X.normalized <- append(X.normalized, list(norm.med(X, groups)$dat.normed))
    names(X.normalized)[length(X.normalized)] <- paste(name, ".Med", sep="")
  }
  if(apply.RUV) {
    X.normalized <- append(X.normalized, list(norm.RUV(X, groups, method = "RUVg")$dat.normed))
    names(X.normalized)[length(X.normalized)] <- paste(name, ".RUVg", sep="")
    X.normalized <- append(X.normalized, list(norm.RUV(X, groups, method = "RUVs")$dat.normed))
    names(X.normalized)[length(X.normalized)] <- paste(name, ".RUVs", sep="")
    X.normalized <- append(X.normalized, list(norm.RUV(X, groups, method = "RUVr")$dat.normed))
    names(X.normalized)[length(X.normalized)] <- paste(name, ".RUVr", sep="")
  }
  if(apply.UQ) {
    X.normalized <- append(X.normalized, list(norm.UQ(X, groups)$dat.normed))
    names(X.normalized)[length(X.normalized)] <- paste(name, ".UQ", sep="")
  }
  if(apply.QN) {
    X.normalized <- append(X.normalized, list(norm.QN(X, groups)$dat.normed))
    names(X.normalized)[length(X.normalized)] <- paste(name, ".QN", sep="")
  }
  
  return(X.normalized)
}



### Definition of Control Genes ------------------------------------------------




#' Compute miRNA (polycistronic) clustering for given set of genes
#' @param genes List of genes for which the miRNA clustering is computed
#' @param coords Coordinate data set used for the location of the genes on the genome
#' @param threshold Maximum distance threshold in base-pairs for a pair of miRNAs
#' for them being considered as clustered.
#' @param single.element.clusters Boolean specifying if clusters are generated for
#' un-clustered genes, thus generating clusters with only a single element.
cluster.genes <- function(genes, coords, threshold=10000, single.element.clusters = TRUE) {
  chromosomes <- unique(coords[genes, "chr"])
  # Create empty list of clusters
  clusters <- list()
  
  # loop over the chromosomes and order/sort the genes on each chromosome by their location
  cluster.id <- 1
  for (chr in chromosomes) {
    # find genes on chromosome chrm
    genesOnChrm <- genes[coords[genes, "chr"] == chr]
    if (length(genesOnChrm) <= 2) {
      next  # next chromosome if there are less than two genes on the current chromosome
    }
    
    # local ordering of the genes on that chromosome
    localOrdering <- order(coords[genesOnChrm, "start"])
    genesOnChrm <- genesOnChrm[localOrdering]
    
    # initialize the first current cluster with the first gene
    curCluster <- genesOnChrm[1]
    # loop over all other genes on the chromosome
    for (i in 2:length(genesOnChrm)) {
      # if the distance of the current gene to the previous gene is less than threshold, add it to curCluster
      if(abs(coords[genesOnChrm[i], "start"] - coords[genesOnChrm[i-1], "end"]) <= threshold) {
        curCluster <- c(curCluster, genesOnChrm[i])
      } else {  # otherwise "close" the previous cluster, add it to the list of clusters and start a new cluster
        if(single.element.clusters || (length(curCluster) > 1)) {  # add only if the cluster contains more than 1 element
          clusters <- append(clusters, list(curCluster))
          names(clusters) <- c(names(clusters)[1:(length(clusters)-1)], paste("cluster-", curCluster[1], "(", length(curCluster), ")", sep=""))
          cluster.id <- cluster.id + 1
        }
        curCluster <- c(genesOnChrm[i])
      }
    }
    # append the last cluster to the list of clusters
    if(single.element.clusters || (length(curCluster) > 1)) {  # add only if the cluster contains more than 1 element
      clusters <- append(clusters, list(curCluster))
      names(clusters) <- c(names(clusters)[1:(length(clusters)-1)], paste("cluster-", curCluster[1], "(", length(curCluster), ")", sep=""))
      cluster.id <- cluster.id + 1
    }
  } # end loop over chr
  
  # return resulting list of clusters
  return(clusters)
}




#' Preprocessing for the normalization performance analysis. Defines positive  
#' and negative Control miRNAs 
#' 
#' You must either provide clusters (as a nested list of miRNA names) or a 
#' coordinates data frame with a cluster threshold. The coordinates data frame
#' must contain columns "chr", "start", and "end".
#' "chr" gives the chromosome for each gene, "start" and "end" are the starting 
#' and ending base-pair location of the gene, respectively. The rownames must
#' correspond to miRNA names.
#' 
#' @param data.RC Read count data matrix. Shape n x p, where n is the number 
#' of samples and p is the number of markers
#' @param t.zero Cutoff for zero-expressed genes
#' @param t.poor Cutoff for poorly-expressed genes. Poorly-expressed genes are
#' defined as genes with read count in [t.zero, t.poor]
#' @param t.well Cutoff for well-expressed genes. Well-expressed genes are 
#' defined as genes with read count in [t.well, inf)
#' @param t.cluster  Maximum distance threshold in base-pairs for a pair of 
#' miRNAs for them being considered as clustered. Default is 50000
#' @param coords A data.frame with one row for each marker with corresponding
#' rownames to the colnames in X. Must contain columns "chr", "start", and "end".
#' "chr" gives the chromosome for each gene, "start" and "end" are the starting 
#' and ending base-pair location of the gene, respectively.
#' @param clusters If clusters are already defined. They can be provided to 
#' this function by stecifying the clusters object (nested list).
define.controls <- function(data.RC,
                            t.zero,
                            t.poor,
                            t.well,
                            t.cluster=10000,
                            coords=NULL,
                            clusters=NULL) {
  genes <- rownames(data.RC)
  
  ## Perform miRNA clustering if no pre-defined clusters are given
  if(is.null(clusters)) {
    if(is.null(coords) || is.null(t.cluster)) {
      # If clustering is performed the coords data frame must be given
      stop("Error in define.controls: If no clusters are provided, 
           you must provide a coordinate data frame and a clustering threshold")
    }
    # Check if the coords data frame contains the columns "chr", "start" and "end"
    if(!all(c("chr", "start", "end") %in% colnames(coords))) {
      stop("Error in normPA: coords data frame does not contain required columns")
    }
    # only consider genes that are present in the coordinate data frame
    genes.in.coords <- rownames(coords)[na.omit(match(genes, rownames(coords)))]
    # sort genes on chromosomes
    genes.in.coords <- sort.on.chromosome(genes.in.coords, coords)
    clusters <- cluster.genes(genes.in.coords,
                              coords,
                              t.cluster,
                              single.element.clusters=TRUE)
  }
  
  ## positive control miRNAs
  genes.pos <- c()
  data.means <- rowMeans(data.RC)
  for(clust in clusters) {
    clust.genes <- unlist(clust)
    clust.genes.range <- names(which(data.means[clust.genes] >= t.well))
    if(length(clust.genes.range) >= 2) {
      genes.pos <- c(genes.pos, clust.genes.range)
    }
  }
  cat("Number of positive control markers with RC in [", t.well, 
      ", inf): ", length(genes.pos), "\n" , sep="")
  
  
  ## negative control miRNAs
  genes.neg <- c()
  data.means <- rowMeans(data.RC)
  for(clust in clusters) {
    clust.genes <- unlist(clust)
    clust.genes.range <- names(which((data.means[clust.genes] <= t.poor) &
                                       (data.means[clust.genes] >= t.zero)))
    if(length(clust.genes.range)!=0) {
      genes.neg <- c(genes.neg, clust.genes.range[length(clust.genes.range)])
    }
  }
  cat("Number of negative control markers with RC in [", t.zero, ", ", 
      t.poor, "]: ", length(genes.neg), "\n" , sep="")
  
  
  return(list(genes.pos=genes.pos,
              genes.neg=genes.neg,
              clusters=clusters))
}



### Precision and Partial Correlation Estimation -------------------------------



#' Pre-Processing for precision estimation 
#' 
#' Pre-processes a given data set for the precision estimation schemes. Includes:
#'   - transposes the data matrix
#'   - removes genes with read count 0
#'   - optional: removes genes with unique number of read counts less than or 
#'   equal to the given threshold unique.genes.threshold. These could cause 
#'   problems for the estimation in mb.estimate.
#'   - modified log2 transform
#'   - scaling to zero mean and variance 1.
#' @param X Data matrix (read counts), size (genes x samples): p x n
#' @param genes Subset of genes. If given, X is reduced to this subset
#' @param remove.zero.count Logical, if genes with read count 0 are removed
#' @param unique.genes.threshold Number of unique read counts. Genes with less 
#' unique read counts are removed. Default is unique.genes.threshold=FALSE (no filtering).
#' @param transform Type of data transformation: no ("n"), modified log2 ("log2"),
#' or cube root ("cube") transform.
#' @param scale.data Logical if the data is scaled column-wise to mean 0 and var 1.
#' @return Pre-processed data matrix with dimensions n x p
preprocess.precision <- function( X, 
                                  genes=NULL,
                                  remove.zero.count=TRUE, 
                                  unique.genes.threshold=FALSE,
                                  transform="log2",
                                  scale.data=TRUE) {
  # Apply genes subset
  if(!is.null(genes)) {
    genes.in.X <- rownames(X)[na.omit(match(genes, rownames(X)))]
    X <- X[genes.in.X, ]
  }
  
  # transpose data
  X <- t(X)
  
  # remove genes with read count 0
  if(remove.zero.count) X <- X[, !(colSums(abs(X))==0)]
  
  if(unique.genes.threshold) {
    # Remove genes with RC != 0 for only a single gene
    X <- X[, apply(X, 2, function(y) length(unique(y)) > unique.genes.threshold)]
  }
  
  # modified log2 transform
  # if(log2.transform) X <- log2(X+1)
  if(transform == "log2"){
    X <- log2(X+1)
  } else if(transform =="n") {
    # no tranform, do nothing
  } else if(transform == "cube") {
    X <- X^(1/3)
  } else {
    stop("Error: Unknown transformation in preprocessing.")
  }
  
  # scale the data by centering and dividing by standard deviation
  # Do not scale columns with standard deviation 0.
  if(scale.data) X <- apply(X, 2, function (y) (y-mean(y)) / sd(y) ^ as.logical(sd(y)))
  return(X)
}



#' Precision estimation using neighborhood selection (see Meinshausen, Buehlmann (2006))
#' @param X Design matrix. Size (samples x parameters): n x p
#' @param criterion Tuning parameter selection criterion. 
#' Available methods: Cross Validantion ("cv"), Akaike information criterion ("aic"),
#' Bayesian information criterion ("bic"), Adaptive Validation ("av"),
#' This paramtere is ignored iff a list of tuning parameters lambda is given. 
#' @param lambda List of p tuning parameters used for the p regressions in the 
#' neighborhood selection. If lambda=NULL, the tuning parameter selection 
#' criterion specified by criterion is used.
#' @return \describe{
#'   \item{prec}{Estimated precision matrix. Shape p x p}
#'   \item{corr}{Partial correlation matrix computed from prec. Shape p x p}
#'   \item{graph}{Co-expression graph computed from prec. }
#'   \item{lambda}{List of tuning parameters used for estimation of prec.}
#' }
mb.estimate <- function(X, criterion="bic", lambda=NULL) {
  library(glmnet)  # for lasso regression
  # Constants - Number of samples and parameters
  n <- nrow(X)  # num samples 
  p <- ncol(X)  # num parameters
  
  if(!is.null(lambda)) {
    criterion <- "lambda"
  }
  if(any(lambda<0)){
    stop("Error in mb.estimate: Tuning parameters lambda must be non-negative!")
  }
  
  betas <- matrix(-1, p, p)  # augmented betas
  tuning <- rep(NA, p)  # tuning parameters lambda
  # compute betas using the lasso estimator
  if(criterion == "cv") {
    # compute betas using the lasso estimator with cross-validation
    for (i in 1:p) {
      betas[-i, i] <- coef(cv.glmnet(x=X[, -i], y=X[, i], intercept=FALSE))[-1, ]
    }
  } else if(criterion == "aic") {
    # compute betas using the lasso estimator with Akaike information criterion
    for (i in 1:p) {
      lasso.fit <- glmnet(x=X[, -i], y=X[, i], intercept=FALSE)
      # AIC criterion
      tLL <- lasso.fit$nulldev- deviance(lasso.fit)  # log-likelihood
      dof <- lasso.fit$df  # degrees-of-freedom
      # lasso.AIC <- -tLL+2*k+2*k*(k+1)/(n-k-1)
      lasso.AIC <- -tLL+2*dof
      AIC.optimal <- which.min(lasso.AIC)
      betas[-i, i] <- lasso.fit$beta[, AIC.optimal]
      tuning[i] <- lasso.fit$lambda[AIC.optimal]
    }
  } else if(criterion == "bic") {
    # compute betas using the lasso estimator with Bayesian information criterion
    for (i in 1:p) {
      lasso.fit <- glmnet(x=X[, -i], y=X[, i], intercept=FALSE)
      # BIC criterion
      tLL <- lasso.fit$nulldev - deviance(lasso.fit)  # log-likelihood
      dof <- lasso.fit$df  # degrees-of-freedom
      # lasso.AIC <- -tLL+2*k+2*k*(k+1)/(n-k-1)
      lasso.BIC <- log(n)*dof - tLL
      BIC.optimal <- which.min(lasso.BIC)
      betas[-i, i] <- lasso.fit$beta[, BIC.optimal]
      tuning[i] <- lasso.fit$lambda[BIC.optimal]
    }
  } else if(criterion == "av") {
    # compute betas using the lasso estimator with adaptive validation    
    # adaptive validation "constant" parameter
    C.bar <- 0.75
    for (i in 1:p) {
      lasso.fit <- glmnet(x=X[, -i], y=X[, i], intercept=FALSE)
      beta <- as.matrix(lasso.fit$beta)  # lasso solution path
      lasso.lambda <- lasso.fit$lambda  # tuning parameter grid
      
      j <- 1
      found.lambda <- FALSE
      # Note that the ordering of the tuning parameters is different from the paper
      # Therefore, go trhough the lambdas from 1:N, not in reverse
      for(j in 1:length(lasso.lambda)) {
        for(k in 1:j) {
          lhs <- max(abs(beta[, j] - beta[, k])) / 
            (lasso.lambda[j] + lasso.lambda[k]) - C.bar
          if(lhs > 0) {
            found.lambda <- TRUE
            break
          }
        }
        if(found.lambda) break
      }
      AV.optimal <- j
      betas[-i, i] <- lasso.fit$beta[, AV.optimal]
      tuning[i] <- lasso.fit$lambda[AV.optimal]
    }
  } else if(criterion=="lambda") {
    # compute betas using the given tuning parameters lambda
    for (i in 1:p) {
      betas[-i, i] <- coef(glmnet(x=X[, -i], y=X[, i], intercept=FALSE, lambda=lambda[i]))[-1, ]
      tuning[i] <- lambda[i]
    }
  } else {
    stop(paste("Error in mb.estimate: Unknown tuning parameter criterion", criterion, "!"))
  }
  
  # inverse covriance matrix Theta
  Theta <- matrix(0, p, p)
  
  # compute diagonal
  for (i in 1:p) {
    Theta[i, i] <- n / ( sum((X[, i] - X[, -i] %*% betas[-i, i]) ^ 2) )
  }
  
  # compute full inverse covariance matrix
  for (i in 1:p) {
    for (j in 1:p) {
      Theta[i, j] <- -0.5 * (Theta[i, i] * betas[j, i] + 
                               Theta[j, j] * betas[i, j])
    }
  }
  
  colnames(Theta) <- rownames(Theta) <- colnames(X)
  partial.correlations <- prec2corr(Theta)
  
  return (list(prec=Theta, 
               corr=partial.correlations, 
               graph=to.graph(Theta),
               lambda=tuning))
}



#' Wrapper for the neighborhood selection using the huge package
#'
#' @param X Data matrix. Size (samples x genes): n x p
#' @param criterion Tuning parameter selection criterion. Available: "ric", "stars"
#' @param lambda A sequence of decreasing positive numbers to control the
#' regularization. Typical usage is to leave the input lambda = NULL 
#' and have the program compute its own lambda sequence based on 
#' nlambda and lambda.min.ratio.
#' @param nlambda The number of regularization/thresholding parameters. 
#' Default nlambda=10
#' @param lambda.min.ratio It is the smallest value for lambda, as a fraction 
#' of the upperbound (MAX) of the regularization/thresholding parameter which 
#' makes all estimates equal to 0. Default is lambda.min.ratio=0.1
#' @param data.name Name of the data set for which the estimation is performed
#' @return \describe{
#'   \item{prec}{Estimated precision matrix. Shape p x p}
#'   \item{corr}{Partial correlation matrix computed from prec. Shape p x p}
#'   \item{graph}{Co-expression graph computed from prec. }
#'   \item{lambda}{List of tuning parameters used for estimation of prec.}
#'   \item{huge.output}{Output of huge estimation.}
#' }
huge.mb.estimate <- function( X, criterion, lambda=NULL, nlambda=10, lambda.min.ratio=0.1, data.name="" ) {
  # output estimation text, if a data.name is given
  if(data.name != "") cat("Conducting huge.mb precision estimation for ",
                          data.name, "....", sep="")
  
  # Check for plausibility of selected tuning parameter calibration method
  if(!any(c(criterion=="ric", criterion=="stars"))) {
    stop(paste("Error in huge.mb.estimate: Unknown tuning parameter criterion", criterion, "!"))
  }
  
  if(is.null(lambda)) {  #no tuning parameter given
    # Estimate for range of tuning parameters
    out.huge <- huge(X, 
                     method="mb", 
                     nlambda=nlambda, 
                     lambda.min.ratio=lambda.min.ratio,
                     verbose=FALSE)
    # use huge.select to select the best tuning parameter
    out.huge <- huge.select(out.huge,criterion=criterion, verbose=FALSE)
    # refit the data useing the selected tuning parameter
    out.huge <- huge(X, 
                     method="mb", 
                     lambda=out.huge$opt.lambda,
                     verbose=FALSE)
  } else {
    # Estimate for range of tuning parameters
    out.huge <- huge(X, 
                     method="mb", 
                     lambda=lambda)
  }
  # Store results
  # TODO: This is a hack...
  # Store precisions as negative of beta
  # set diagonal to -1 to be able to transform to partial correlations
  # Not sure why this is not done by huge...
  precision <- -as.matrix(out.huge$beta[[1]])
  diag(precision) <- -1
  colnames(precision) <- rownames(precision) <- colnames(X)
  
  # used tuning parameter
  lambda <- as.numeric(out.huge$lambda)
  # partial correlations
  partial.correlations <- prec2corr(precision)  
  # graph
  graph <- to.graph(precision)
  
  # finish output
  if(data.name != "") cat("done\n", sep="")
  
  return( list(prec=precision, 
               corr=partial.correlations,
               graph=graph, 
               lambda=lambda, 
               huge.output=out.huge) )
}



#' Wrapper for the graphical lasso using the huge package
#'
#' @param X Data matrix. Size (samples x genes): n x p
#' @param criterion Tuning parameter selection criterion.  Available: "ric", 
#' "stars", and "ebic"
#' @param lambda A sequence of decreasing positive numbers to control the
#' regularization. Typical usage is to leave the input lambda = NULL 
#' and have the program compute its own lambda sequence based on 
#' nlambda and lambda.min.ratio.
#' @param nlambda The number of regularization/thresholding parameters. 
#' Default nlambda=10
#' @param lambda.min.ratio It is the smallest value for lambda, as a fraction 
#' of the upperbound (MAX) of the regularization/thresholding parameter which 
#' makes all estimates equal to 0. Default is lambda.min.ratio=0.1
#' @param data.name Name of the data set for which the estimation is performed
#' @return \describe{
#'   \item{prec}{Estimated precision matrix. Shape p x p}
#'   \item{corr}{Partial correlation matrix computed from prec. Shape p x p}
#'   \item{graph}{Co-expression graph computed from prec. }
#'   \item{lambda}{List of tuning parameters used for estimation of prec.}
#'   \item{huge.output}{Output of huge estimation.}
#' }
huge.glasso.estimate <- function( X, criterion, lambda=NULL, nlambda=10, lambda.min.ratio=0.1, data.name="" ) {
  # output estimation text, if a data.name is given
  if(data.name != "") cat("Conducting huge.glasso precision estimation for ",
                          data.name, "....", sep="")
  
  # Check for plausibility of selected tuning parameter calibration method
  if(!any(c(criterion=="ric", criterion=="stars", criterion=="ebic"))) {
    stop(paste("Error in huge.mb.estimate: Unknown tuning parameter criterion", criterion, "!"))
  }
  
  if(is.null(lambda)) {  #no tuning parameter given
    # Estimate for range of tuning parameters
    out.huge <- huge(X, 
                     method="glasso", 
                     nlambda=nlambda, 
                     lambda.min.ratio=lambda.min.ratio,
                     verbose=FALSE)
    # use huge.select to select the best tuning parameter
    out.huge <- huge.select(out.huge, criterion=criterion,
                            verbose=FALSE)
    # refit the data useing the selected tuning parameter
    out.huge <- huge(X, 
                     method="glasso", 
                     lambda=out.huge$opt.lambda,
                     verbose=FALSE)
  } else {
    # Estimate for range of tuning parameters
    out.huge <- huge(X, 
                     method="glasso", 
                     lambda=lambda)
  }
  
  # Precisions
  precision <- as.matrix(out.huge$icov[[1]])
  colnames(precision) <- rownames(precision) <- colnames(X)
  # Partial correlations
  corr <- prec2corr(precision)
  # graph
  graph <- to.graph(precision)
  
  # used tuning parameter
  lambda <- as.numeric(out.huge$lambda)
  
  # finish output
  if(data.name != "") cat("done\n", sep="")
  
  return( list(prec=precision, 
               corr=corr,
               graph=graph, 
               lambda=lambda, 
               huge.output=out.huge) )
}



#' Precision estimation using the FastGGM method. Paper: "FastGGM: An Efficient 
#' Algorithm for the Inference of Gaussian Graphical Model in Biological Networks"
#' by Wang et al. (2016).
#' This function provides a wrapper for the method `FastGGM` provided by the 
#' authors online at http://www.pitt.edu/~wec47/fastGGM.html
#' The default tuning parameter calibration method provided by the package is used.
#' @param X Design matrix. Size (samples x parameters): n x p
#' @param data.name Name of the data set for which the estimation is performed
#' @return \describe{
#'   \item{prec}{Estimated precision matrix. Shape p x p}
#'   \item{corr}{Partial correlation matrix computed from prec. Shape p x p}
#'   \item{graph}{Co-expression graph computed from prec. }
#'   \item{huge.output}{Output of FastGGM estimation.}
#' }
fastggm.estimate <- function(X, data.name="" ) {
  # output estimation text, if a data.name is given
  if(data.name != "") cat("Conducting fastggm precision estimation for ",
                          data.name, "....", sep="")
  
  fastggm.out <- FastGGM(X)
  
  # Results
  precision <- fastggm.out$precision
  rownames(precision) <- colnames(precision) <- colnames(X)
  corr <- fastggm.out$partialCor
  rownames(corr) <- colnames(corr) <- colnames(X)
  
  # finish output
  if(data.name != "") cat("done\n", sep="")
  
  return( list(prec=precision, 
               corr=corr,
               graph=to.graph(precision), 
               fastggm.output=fastggm.out) )
}




#' Estimates precision and partial correlation matrices for a single data matrix
#' @param data.RC Data set matrix of read counts. Size: p x n (genes x samples)
#' @param preprocess.transform Used transformation "log2" or "cube" for the preprocessing
#' @param name Name/Title of the data set
#' @param method Method used for the precision estimation.Available: "mb",
#' "huge.mb", "glasso", and "fastggm"
#' @param genes Subset of genes. If given, data.RC is reduced to this subset
#' @param tuning.criterion Used tuning parameter criterion in the precision estimation 
#' method. Options are "cv", "aic", "bic", and "av" for the mb method.
#' @return List of estimated precisions and partial correlations for the 
#' positive and negative control sets
estimate.precision.single <- function(data.RC, 
                                      preprocess.transform, 
                                      name, 
                                      method,
                                      genes=NULL, 
                                      tuning.criterion=NULL,
                                      lambda=NULL) {
  # Preprocessing
  X <- preprocess.precision(data.RC,
                            genes=genes,
                            remove.zero.count=TRUE, 
                            unique.genes.threshold=1,
                            transform=preprocess.transform,
                            scale.data=TRUE) 
  
  ## Perform precision estimation form positive and negative controls
  set.seed(42)  # for reprodicubility
  if(method=="mb") {
    # Perform the mb estimation
    result <- mb.estimate(X, criterion=tuning.criterion, lambda=lambda)
  } else if(method=="huge.mb") {
    # Perform the huge.mb estimation
    result <- huge.mb.estimate(X, criterion=tuning.criterion, lambda=lambda)
  } else if(method=="glasso") {
    # Perform the glasso estimation
    result <- huge.glasso.estimate(X, criterion=tuning.criterion, lambda=lambda)
  } else if(method=="fastggm") {
    # Perform the fastggm estimation
    result <- fastggm.estimate(X, data.name=name)
  }else {
    stop(paste("Error: Unknown value", method, "for parameter method in function estimate.precision.single"))
  }
  
  return(list(name=name,
              prec=result$prec, 
              corr=result$corr,
              X=X,
              lambda.pos=result$lambda))
}



### Marginal Correlation Estimation --------------------------------------------



#' Pre-Processing for marginal correlation estimation 
#' 
#' Pre-processes a given data set for the covariance estimation schemes. Includes:
#'   - transposes the data matrix
#'   - removes genes with read count 0
#'   - optional: removes genes with unique number of read counts less than or 
#'   equal to the given threshold unique.genes.threshold. These could cause 
#'   problems for the estimation in mb.estimate.
#'   - modified log2 transform
#' @param X Data matrix (read counts), size (genes x samples): p x n
#' @param genes Subset of genes. If given, X is reduced to this subset
#' @param remove.zero.count Logical, if genes with read count 0 are removed
#' @param unique.genes.threshold Number of unique read counts. Genes with less 
#' unique read counts are removed. Default is unique.genes.threshold=FALSE (no filtering).
#' @param transform Type of data transformation: no ("n"), modified log2 ("log2"),
#' or cube root ("cube") transform.
#' @return Pre-processed data matrix with dimensions n x p
preprocess.correlation <- function( X, 
                                    genes=NULL,
                                    remove.zero.count=TRUE, 
                                    unique.genes.threshold=FALSE,
                                    transform="log2") {
  # Apply genes subset
  if(!is.null(genes)) {
    genes.in.X <- rownames(X)[na.omit(match(genes, rownames(X)))]
    X <- X[genes.in.X, ]
  }
  
  # transpose data
  X <- t(X)
  
  # remove genes with read count 0
  if(remove.zero.count) X <- X[, !(colSums(abs(X))==0)]
  
  if(unique.genes.threshold) {
    # Remove genes with RC != 0 for only a single gene
    X <- X[, apply(X, 2, function(y) length(unique(y)) > unique.genes.threshold)]
  }
  
  # modified log2 transform
  # if(log2.transform) X <- log2(X+1)
  if(transform == "log2"){
    X <- log2(X+1)
  } else if(transform =="n") {
    # no tranform, do nothing
  } else if(transform == "cube") {
    X <- X^(1/3)
  } else {
    stop("Error: Unknown transformation in preprocessing.")
  }
  
  return(X)
}



#' Estimates marginal correlation matrices for a read count single data matrix
#' @param data.RC Data set matrix of read counts. Size: p x n (genes x samples)
#' @param preprocess.transform Used transformation "log2" or "cube" for the preprocessing
#' @param name Name/Title of the data set
#' @param method Method used for the correlation estimation. Available: "pearson"
#' and "spearman"
#' @param genes Subset of genes. If given, data.RC is reduced to this subset
#' @return List of estimated precisions and partial correlations for the 
#' positive and negative control sets
estimate.correlation.single <- function(data.RC, 
                                        preprocess.transform, 
                                        name, 
                                        method,
                                        genes=NULL) {
  # Preprocessing
  X <- preprocess.correlation(data.RC,
                              genes=genes,
                              remove.zero.count=TRUE, 
                              unique.genes.threshold=1,
                              transform=preprocess.transform) 
  
  ## Perform precision estimation form positive and negative controls
  set.seed(42)  # for reprodicubility
  if(method=="pearson") {
    # Perform the pearson correlation estimation
    result <- cor(X, method="pearson")
  } else if(method=="spearman") {
    # Perform the spearman correlation estimation
    result <- cor(X, method="spearman")
  } else {
    stop(paste("Error: Unknown value", method, "for parameter method in function estimate.correlation.single"))
  }
  
  return(list(name=name,
              corr=result,
              X=X))
}



### Correlation Estimation Wrapper ---------------------------------------------



#' Estimates correlations for positive and negative controls in a data set.
#' User can choose partial or marginal correlations for each set of controls
#' individually.
#' @param data.RC Data set matrix (read counts)
#' @param genes.pos Positive control genes
#' @param genes.neg Negative control genes
#' @param preprocess.transform Used transformation "log2" or "cube" for the preprocessing
#' @param name Name/Title of the data set
#' @param method.pos Method used for the correlation estimation for positive controls. 
#' Available: "mb", "huge.mb", "glasso", "fastggm", "spearman", and "pearson".
#' @param method.neg Method used for the correlation estimation for negative controls. 
#' Available: "mb", "huge.mb", "glasso", "fastggm", "spearman", and "pearson".
#' @param tuning.criterion.pos Used tuning parameter criterion for the correlation
#' estimation method method.pos. Options are "cv", "aic", "bic", and "av" for the "mb" method.
#' @param tuning.criterion.neg Used tuning parameter criterion for the correlation
#' estimation method method.neg. Options are "cv", "aic", "bic", and "av" for the "mb" method.
#' @return List of estimated correlations for positive and negative controls
estimate.correlations <-  function(data.RC, 
                                   genes.pos, 
                                   genes.neg, 
                                   preprocess.transform, 
                                   name, 
                                   method.pos,
                                   method.neg,
                                   tuning.criterion.pos,
                                   tuning.criterion.neg) {
  cat("Estimating correlations for pos. and neg. controls for data set ", name, "...", sep="")
  
  ## Estimation for positive controls
  if(method.pos %in% c("mb", "huge.mb", "glasso", "fastggm")) {
    # estimate partial correlations
    result.pos <- estimate.precision.single(data.RC=data.RC,
                                            preprocess.transform=preprocess.transform,
                                            name=name,
                                            method=method.pos,
                                            genes=genes.pos,
                                            tuning.criterion=tuning.criterion.pos)
  } else if(method.pos %in% c("pearson", "spearman")) {
    # estimate marginal correlations
    result.pos <- estimate.correlation.single(data.RC=data.RC,
                                              preprocess.transform=preprocess.transform,
                                              name=name,
                                              method=method.pos,
                                              genes=genes.pos)
  } else {
    stop(paste("Error: Unknown value", method.pos, "for parameter method.pos in function estimate.correlations"))
  }
  
  ## Estimation for negative controls
  if(method.neg %in% c("mb", "huge.mb", "glasso", "fastggm")) {
    # estimate partial correlations
    result.neg <- estimate.precision.single(data.RC=data.RC,
                                            preprocess.transform=preprocess.transform,
                                            name=name,
                                            method=method.neg,
                                            genes=genes.neg,
                                            tuning.criterion=tuning.criterion.neg)
  } else if(method.neg %in% c("pearson", "spearman")) {
    # estimate marginal correlations
    result.neg <- estimate.correlation.single(data.RC=data.RC,
                                              preprocess.transform=preprocess.transform,
                                              name=name,
                                              method=method.neg,
                                              genes=genes.neg)
  } else {
    stop(paste("Error: Unknown value", method.neg, "for parameter method.neg in function estimate.correlations"))
  }
  
  cat("done\n")
  return(list(name=name,
              pos=list(name=name,
                       corr=result.pos$corr,
                       X=result.pos$X),
              neg=list(name=name,
                       corr=result.neg$corr,
                       X=result.neg$X)))
}



### Assessment  Metrics --------------------------------------------------------



#' Compute the mean correlation reduction in negative controls
compute.mcr <- function(rawCor, normCor) {
  mean.rawCor <- sum(pmax(rawCor[upper.tri(rawCor)], 0)) / sum(upper.tri(rawCor))
  mean.normCor <- sum(pmax(normCor[upper.tri(normCor)], 0)) / sum(upper.tri(normCor))
  mcr <- (mean.rawCor - mean.normCor) / mean.rawCor
  return(list(mcr=mcr,
              rawCor=rawCor,
              normCor=normCor))
}



#' Compute the concordance correlation of clustered positive controls
compute.cc <- function(rawCor, normCor, clusters) {
  # remove genes that are not present in both models
  rawCor <- rawCor[!is.na(match(colnames(rawCor), colnames(normCor))),
                   !is.na(match(colnames(rawCor), colnames(normCor)))]
  normCor <- normCor[!is.na(match(colnames(normCor), colnames(rawCor))),
                     !is.na(match(colnames(normCor), colnames(rawCor)))]
  
  # only consider cluster with multiple genes
  clusters <- clusters[lengths(clusters) > 1]
  clusters.raw <- c()
  clusters.norm <- c()
  for (clust in clusters) {
    # only consider cluster genes that are positive controls
    clust.genes <- clust[stats::na.omit(match(colnames(rawCor), clust))]
    if (length(clust.genes) < 2) {
      next  # disregard clusters with less than 2 genes
    }
    # subset of correlations in the cluster
    rawCor.clust <- rawCor[clust.genes, clust.genes]
    clusters.raw <- c(clusters.raw, rawCor.clust[upper.tri(rawCor.clust)])
    normCor.clust <- normCor[clust.genes, clust.genes]
    clusters.norm <- c(clusters.norm, normCor.clust[upper.tri(normCor.clust)])
  }
  
  ## Compute the concordance correlation between nonzero partial correlations
  idx <- clusters.raw | clusters.norm
  if (sum(idx)>0) {
    cc <- DescTools::CCC(as.vector(clusters.raw[idx]),
                         as.vector(clusters.norm[idx]))$rho.c$est
  } else {
    cc <- 0
  }
  
  # For all cluster correlations (not just the non-zero ones)
  cc <- DescTools::CCC(as.vector(clusters.raw),as.vector(clusters.norm))$rho.c$est
  
  return(list(cc=cc,
              rawCor=rawCor,
              normCor=normCor))
}



#' Comparison of two models for both positive and negative controls. 
#' Includes scatter plot comparison and agreement computation.
#' @param rawModel First estimation model
#' @param normModel Second estimation model
#' @param clusters Clustering for all genes
#' @param coords Coordinates data frame. If given will be used to display
#' chromosomes in correlation plots
#' @param corr.threshold Only values of the correlation matrix greater than the 
#' threshold in absolute value are plotted. Default is corr.threshold=0 which
#' corresponds to no thresholding
#' @param display.results Boolean whether to display the resulting graphs and statistics
compare.models <- function(rawModel, 
                           normModel, 
                           clusters, 
                           coords=NULL,
                           corr.threshold=0, 
                           display.results=TRUE) {
  cat("Comparing models ", rawModel$name, " vs. ", normModel$name, "....", sep="")
  
  # Comparison of positive controls
  comparison.pos <- compute.cc(rawCor=rawModel$pos$corr, 
                               normCor=normModel$pos$corr, 
                               clusters=clusters)
  
  # Comparison of negative controls
  comparison.neg <- compute.mcr(rawCor=rawModel$neg$corr, normCor=normModel$neg$corr)
  
  ## Plot results
  
  if(display.results) {
    ## (clustered) graph plot
    # Positive controls
    p1 <- ggplot.corr(comparison.pos$rawCor, 
                      clusters=clusters, 
                      threshold=corr.threshold, 
                      title=paste("Corr - ", rawModel$name),
                      coords=coords)
    p2 <- ggplot.corr(comparison.pos$normCor, 
                      clusters=clusters, 
                      threshold=corr.threshold, 
                      title=paste("Corr - ", normModel$name),
                      coords=coords)
    grid.arrange(p1 + theme(legend.position="none"), 
                 p2 + theme(legend.position="none"),
                 get.legend(p1),
                 widths = c(2,2,1),
                 top="Positive Controls")
    
    ## scatter plot comparison
    par(mfrow=c(1,2))  # Scatter plots for pos controls
    p1 <- plot.corr.scatter(
      comparison.pos$rawCor, 
      comparison.pos$normCor,
      title="full graph",
      xlab=rawModel$name, 
      ylab=normModel$name, 
      threshold=corr.threshold)
    p2 <- plot.corr.scatter(
      comparison.pos$rawCor, 
      comparison.pos$normCor, 
      clusters=clusters,
      title="clusters",
      xlab=rawModel$name, 
      ylab=normModel$name, 
      threshold=corr.threshold)
    grid.arrange(p1, p2, ncol=2, top=paste(rawModel$name, " vs. ", normModel$name))
  }
  
  
  cat("done\n", sep="")
  
  
  ## Return results
  return(list(
    mcr=comparison.neg$mcr,
    cc=comparison.pos$cc,
    rawModel=rawModel,
    normModel=normModel,
    rawName=rawModel$name,
    normName=normModel$name,
    rawCorr.pos=comparison.pos$rawCor,
    rawCorr.neg=comparison.neg$rawCor,
    normCorr.pos=comparison.pos$normCor,
    normCorr.neg=comparison.neg$normCor
  ))
}



### Data-driven Normalization Assessment ---------------------------------------



#' Normalization Performance Analysis
#' 
#' This function performs an analysis of the performance of normalization 
#' methods applied to a data set. The normalization itself is not done by
#' `normPA`. Instead, the user must provide the data set as a matrix
#' and a list of normalized versions of this data matrix. 
#' Data matrices (normalized and un-normalized) must be
#' of shape p x n where p is the number of markers and n is the 
#' number of samples. The row names must correspond to the marker (gene) 
#' names. 
#' 
#' Two arrays of names of control markers must be provided in pos.controls
#' and neg.controls for positive controls and negative controls, respectively.
#' Furthermore, the used clusters for positive and negative controls
#' must be provided as a nested list of marker names. 
#' TODO Explain how the nested list looks like
#' A function for the definition of positive and negative controls and clusters
#' is provided by `normPA.preprocess`.
#' 
#' @param data.RC The original data matrix (read counts), 
#' size (genes x samples): p x n
#' @param data.norm List of normalized data matrices
#' @param pos.controls Array of names of positive control markers
#' @param neg.controls Array of names of negative control markers
#' @param clusters A list of lists of marker (gene) names. The top-level list 
#' holds the clusters. Each cluster is given by a list of marker names of 
#' all markers in the corresponding cluster.
#' @param coords Coordinates data frame. If given will be used to display
#' chromosomes in correlation plots.
#' @param case.name The name of the data set under study. Used for labels 
#' in plots and tables.
#' @param generate.plots Boolean if normPA generates result plots of the 
#' estimated models and displays an overview table of the result statistics.
#' @param preprocess.transform Pre-processing transformation that will be
#' applied to the read count data. Options are modified log2 
#' transformation ("log2"), cubic root transformation ("cube"), and 
#' no transformation ("n")
#' @param corr.method.pos Method used for the correlation estimation for positive controls. 
#' Available: "mb", "huge.mb", "glasso", "fastggm", "spearman", and "pearson".
#' @param corr.method.neg Method used for the correlation estimation for negative controls. 
#' Available: "mb", "huge.mb", "glasso", "fastggm", "spearman", and "pearson".
#' @param tuning.criterion.pos Used tuning parameter criterion for the correlation
#' estimation method method.pos. Options are "cv", "aic", "bic", and "av" for the "mb" method.
#' @param tuning.criterion.neg Used tuning parameter criterion for the correlation
#' estimation method method.neg. Options are "cv", "aic", "bic", and "av" for the "mb" method.
#' 
DANA <- function(data.RC, 
                 data.norm, 
                 pos.controls, 
                 neg.controls, 
                 clusters,
                 coords=NULL,
                 case.name="",
                 generate.plots=FALSE,
                 preprocess.transform="log2",
                 corr.method.pos="mb",
                 tuning.criterion.pos="bic",
                 corr.method.neg="pearson",
                 tuning.criterion.neg="") {
  # Model for un-normalized data
  data.model <- estimate.correlations(data.RC=data.RC,
                                      genes.pos=pos.controls,
                                      genes.neg=neg.controls,
                                      preprocess.transform=preprocess.transform,
                                      name=case.name,
                                      method.pos=corr.method.pos,
                                      method.neg=corr.method.neg,
                                      tuning.criterion.pos=tuning.criterion.pos,
                                      tuning.criterion.neg=tuning.criterion.neg)
  # Models for normalized data
  data.norm.models <- list()
  for(i in 1:length(data.norm)) {
    data.norm.models <- append(data.norm.models,
                               list(estimate.correlations(
                                 data.RC=data.norm[[i]],
                                 genes.pos=pos.controls,
                                 genes.neg=neg.controls,
                                 preprocess.transform=preprocess.transform,
                                 name=names(data.norm)[i],
                                 method.pos=corr.method.pos,
                                 method.neg=corr.method.neg,
                                 tuning.criterion.pos=tuning.criterion.pos,
                                 tuning.criterion.neg=tuning.criterion.neg)))
  }
  names(data.norm.models) <- names(data.norm)
  
  
  ## Compare Before/After Normalization
  
  
  comparison <- list()
  # # Compare un-normalized data with itself
  # comparison <- list(
  #   compare.models(data.model,
  #                  data.model,
  #                  clusters=clusters,
  #                  coords=coords,
  #                  corr.threshold=0,
  #                  display.results=generate.plots))
  for(y in data.norm.models) {
    comparison <-
      append(comparison,
             list(compare.models(rawModel=data.model,
                                 normModel=y,
                                 clusters=clusters,
                                 coords=coords,
                                 corr.threshold=0,
                                 display.results=generate.plots)))
  }
  names(comparison) <- names(data.norm)
  
  
  ## Summarize the results
  metrics <- data.frame(t(sapply(comparison, function(x) c(x$cc, x$mcr))))
  colnames(metrics) <- c("cc", "mcr")
  
  if(generate.plots) {
    # ASCII table plot of the resulting summary metrics
    stargazer(metrics, type="text", summary=FALSE, digits=4, 
              title="Comparison statistics", align=TRUE)
    
    ## Result statistics plot
    p.res <- plot.DANA.metrics(metrics, label.size=3, label.repel=TRUE)
    print(p.res)
  }
  
  return(list(data.model=data.model,
              norm.models=data.norm.models,
              metrics=metrics))
}



### Differential Expression Analysis -------------------------------------------



DE.voom <- function(data.RC, groups, plot=FALSE, plot.title="", adjust = NULL) {
  ## Differential expression Analysis using voom
  if(plot) par(mfrow=c(1,2), oma=c(0, 0, 2, 0))
  
  # Construct a DEList data frame
  x <- DGEList(counts=data.RC, group=groups, genes=rownames(data.RC))
  # keep <- filterByExpr(x, group=groups, min.count=0, min.total.count=10)
  # x <- x[keep,,]
  
  group <- factor(groups)
  
  if(is.null(adjust)) {
    design <- model.matrix(~ 0 + group)
    colnames(design) <- levels(group)
  } else {
    design <- model.matrix(~ 0 + group + adjust)
    colnames(design)[1:2] <- levels(group)
  }
  
  contr <- paste(levels(group)[2], "-", levels(group)[1])
  contr.matrix <- makeContrasts(contrasts=contr, 
                                levels = design)
  
  v <- voom(x, design, plot=plot)
  vfit <- lmFit(v, design)
  vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
  efit <- eBayes(vfit)
  
  # organize results into table
  top.table <- topTable(efit, sort.by = "P", n = Inf)
  
  if(plot) {
    plotSA(efit, main="Final model: Mean-variance trend")
    mtext(plot.title, outer = TRUE, cex = 1.5)
    par(mfrow=c(1,1))
  }
  return(top.table)
}



#' Comparison of a single DEA to a "truth" (a gold standard)
#'
#' Computing the true positive rate, false postive rate, false discovery rate, 
#' and false negative rate based on given golden standard.
#'
#' @param DEA result for the data under study
#' @param truth.DEA DEA result for the assumed truth (the gold standard)
#' @param alpha p-value cutoff (significance level). 
#' @param subset vector of a subset of genes/markers for this analysis
#'
#' @return a list of values about TPR, FPR, FDR, FNR
#'
compare.DEAs <- function(DEA, truth.DEA, alpha, subset=NULL) {
  # Subset DEAs to a given set of miRNAs
  if (!is.null(subset)) {
    DEA <- DEA[subset, ]
    truth.DEA <- truth.DEA[subset, ]
  }
  # compute DE genes
  truth.DE <- truth.DEA$genes[truth.DEA$P.Value < alpha]
  compare.DE <- DEA$genes[DEA$P.Value < alpha]
  
  stat.DE <- matrix("non-DE", nrow = nrow(truth.DEA), ncol = 2)
  rownames(stat.DE) <- sort(rownames(truth.DEA))
  colnames(stat.DE) <- c("Benchmark", "prediction")
  stat.DE[truth.DE, 1] <- "DE"
  stat.DE[compare.DE, 2] <- "DE"
  t <- table(prediction = stat.DE[,2], truth = stat.DE[,1])
  return(list(count = t,
              TPR = t[1,1]/sum(t[,1]),
              FPR = t[1,2]/sum(t[,2]),
              FDR = t[1,2]/sum(t[1,]),
              FNR = t[2,1]/sum(t[,1])))
}





























































