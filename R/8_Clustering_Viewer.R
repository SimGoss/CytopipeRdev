#' @title Plots the numbers of cells of each population (cluster or metacluster)
#'
#' @description This function aims to visualize the number of cells associated to each population. The populations can be clusters or metaclusters
#'
#' This representation displays the population in the X-axis and the total number of associated cells in the Y-axis.
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param population a character vector containing the identifiers of the population to use. By default, all the population are used
#' @param level a character value indicating the type of population plotted. Possible values are: "clusters", "metaclusters". By default, 'clusters' are used.
#' @param sort a boolean value indicating if population must be sorted by the number associated cluster
#' @param color.metadata a character value specifying the metadata used to color the barplot. By default, color.metadata is set to NULL (the barplot color is uniform)
#'
#' @return a ggplot2 object
#'
#' @export
#'

plotPopulationCounts <- function(CYTdata,
                                 population = NULL,
                                 level = c("clusters", "metaclusters"),
                                 sort = TRUE,
                                 color.metadata = NULL) {
  
  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.")}
  checkmate::qassert(population, c("0", "S*"))
  level = match.arg(level)
  checkmate::qassert(level, "S1")
  checkmate::qassert(sort, "B1")
  checkmate::qassert(color.metadata, c("0", "S1"))
  
  if (level == "clusters"){
    if (length(CYTdata@Clustering@clusters)==0) { stop("Error : Clustering step required before vizualization") }
    cellcount = CYTdata@Clustering@cellcount
  }
  if (level == "metaclusters"){
    if (length(CYTdata@Metaclustering@metaclusters)==0) { stop("Error : Metaclustering step required before vizualization") }
    cellcount = CYTdata@Metaclustering@cellcount
  }
  
  if (is.null(population)) { population = rownames(cellcount) }
  cellcount = cellcount[rownames(cellcount) %in% population, ]
  
  pop.effectif = apply(cellcount, 1, sum)
  
  if (!is.null(color.metadata)){
    
    if (!color.metadata %in% colnames(CYTdata@metadata)){ stop("Error : 'color.metadata' argument
    is not a biological condition present in metadata (", paste0(colnames(CYTdata@metadata), collapse=","), ")") }
    
    matrix.count = merge(CYTdata@metadata, t(cellcount), by = "row.names")
    matrix.count = plyr::ddply(matrix.count, color.metadata, function(x){
      return(colSums(x[,rownames(cellcount)]))
    })
    matrix.count = reshape2::melt(matrix.count, id = c(color.metadata))
    colnames(matrix.count) = c("condition", "population", "count")
    
    if (sort) {
      matrix.count$population = factor(matrix.count$population,
                                       levels = rownames(cellcount)[order(pop.effectif, decreasing=TRUE)])
    }
    
    plot <- ggplot2::ggplot(data = matrix.count,
                            ggplot2::aes(x = population,
                                         y = count,
                                         fill = condition)) +
      ggplot2::geom_bar(stat="identity", position = "stack", color="black") +
      viridis::scale_fill_viridis(option = "turbo", discrete = TRUE)
    
  } else {
    
    matrix.count = data.frame("count" = pop.effectif, "population" = rownames(cellcount))
    
    if (sort) {
      matrix.count$population = factor(matrix.count$population,
                                       levels = rownames(cellcount)[order(pop.effectif, decreasing=TRUE)])
    }
    
    plot <- ggplot2::ggplot(data = matrix.count,
                            aes(x = population,
                                y = count)) +
      ggplot2::geom_bar(stat="identity")
  }
  
  if (is.null(color.metadata)) {titre = paste("Number of cells per", level)
      
  } else {titre = paste("Number of cells per", level, "with", color.metadata, "proportion")}
  
  plot <- plot +
    ggplot2::ggtitle(titre) +
    ggplot2::xlab(level) + ggplot2::ylab("Number of cells") +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
                   legend.key = ggplot2::element_blank())
  
  return(plot)
}


#' @title Plots a diagram showing correlated markers
#'
#' @description This function aims to visualize the pairwise co-expression between all markers.
#' Each tile corresponds to the co-expression between two markers and is gradient-colored based on the Pearson or Spearman correlation
#'
#' @param CYTdata a CYTdata object
#' @param markers a character vector containing the name of the markers to analyse
#' @param population a character vector containing the identifiers of the population to use. By default, all the population are used
#' @param level a character value indicating the type of population plotted. Possible values are: "clusters", "metaclusters". By default, 'clusters' are used.
#' @param samples a character vector containing the name of the samples to consider
#' @param method a character value indicating the name of correlation method to use, "pearson" (default) or "spearman"
#' @param palette optionnal; a vector of color names to use for the gradient. 

#'
#' @return a list with:
#'  - corPlot: a ggplot2 object, the correlation diagram
#'  - corMatrix: the correlation matrix
#'
#' @export
#'
plotCorrDiagram <- function(CYTdata,
                            markers = NULL,
                            population = NULL,
                            level = c("clusters", "metaclusters"),
                            samples = NULL,
                            method = c("pearson","spearman"),
                            palette = NULL) { #c("yellow", "orange", "red", "brown")
  
  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { validObject(CYTdata) }
  
  markers = checkorderMarkers(CYTdata, markers, order=FALSE, checkDuplicates=TRUE)
  
  level = match.arg(level)
  checkmate::qassert(level, "S1")
  population = checkorderPopulation(CYTdata, population=population, level=level,
                                    order=TRUE, checkDuplicates=TRUE)
  samples = checkorderSamples(CYTdata, samples,
                              order=TRUE, checkDuplicates=TRUE)
  
  method = match.arg(method)
  checkmate::qassert(method, "S1")
  
  checkmate::qassert(palette, c("0", "S*"))
  
  if(level == "clusters") { popId = CYTdata@Clustering@clusters }
  else { popId = CYTdata@Metaclustering@metaclusters }
  
  data = subset(CYTdata@matrix.expression[,markers],
                popId %in% population & CYTdata@samples %in% samples)
  
  corMatrix = round(stats::cor(data, method = method), 3)
  # dist = stats::as.dist(1 - corMatrix) #CAUTION : this causes the data to not be in the same order as markers, which ultimately gives completely false graphs. 
  # hc = stats::hclust(dist)
  # corMatrix = corMatrix[hc$order, hc$order]
  corMatrix[upper.tri(corMatrix, diag = TRUE)] = NA
  
  # print(corMatrix)
  
  corMatrixmelted = reshape2::melt(corMatrix)
  
  # print(corMatrixmelted)
  
  plot <- ggplot2::ggplot(data = corMatrixmelted,
                          ggplot2::aes_string(x = "Var1",
                                              y = "Var2",
                                              fill = "value")) +
    ggplot2::ggtitle("Correlation diagram") +
    ggplot2::geom_tile(color = "white")
  
  
  if(!is.null(palette)){
    if(!all(areColors(palette))){
      stop("Error : 'palette' argument (", paste0(palette, collapse = ","),
           ") does not contain only hexadecimal color.)")
    }
    plot <- plot + ggplot2::scale_fill_gradientn(colours = palette, na.value = "white", name = paste(method, "correlation"))
  }
  else {
    plot <- plot + ggplot2::scale_fill_gradient2(low = "#20aaff", high = "orange", mid = "black",
                                                 midpoint = 0, limit = c(-1, 1), na.value = "white",
                                                 name = paste(method, "correlation"),
                                                 oob = scales::oob_squish)
  }
  
  plot <- plot +
    ggplot2::coord_fixed() +
    ggplot2::scale_y_discrete(position = "right") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 60, hjust = 1),
      axis.text.y = ggplot2::element_text(hjust = 0)) 
  
  return(list("corPlot" = plot,
              "corMatrix" = corMatrix))
}





