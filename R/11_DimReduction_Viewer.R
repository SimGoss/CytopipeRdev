
#' @title Visualize 2D reduced cell data according to density and metadata
#'
#' @description This function aims to visualize cells in the two dimensional embedded space after dimensionality reduction has been performed. 
#' This representation can be used on a CYTdata object for which a dimension reduction has been computed
#' Cells can be colored according to any kind of metadata stored in the CYTdata object, or according to density. 
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param samples a character vector containing the names of the biological samples to plot. Non included samples will still be ploted, but in grey. 
#' @param population a character vector containing the identifiers of the population to plot. By default, all the population are used. Non included populations will still be ploted, but in grey.
#' Please note that this filters the cells that are ploted but does not affect their colors, for that, please see plotDimRedPopulations
#' @param level a character value indicating the type of population filtered by the population argument. Possible values are: "clusters", "metaclusters". By default, 'clusters' are used.
#' @param plotType a character value indicating the type of coloring to apply. 
#' Possible values are "density" (default, the points are colors according to cell density as if in a 2D histogram) and "metadata" (the points are coloured according a specified metadata)
#' @param title a character value providing the title for the graph. Defaults to "Dimensionnality reduction view"
#' @param metadata a character vector containing the type of metadata to study, if plotType = "metadata".
#' @param axes a numeric of length 2 indicating the axes to display. Defaults to 1 and 2. 
#'
#' @return a ggplot2 object
#'
#' @export
#'

plotDimReduction <- function(CYTdata,
                             samples = NULL,
                             population = NULL,
                             level = c("clusters", "metaclusters"),
                             plotType = c("density", "metadata"),
                             title = "Dimensionnality reduction view",
                             metadata = "Timepoint",
                             axes = c(1, 2)) {
  
  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }
  
  checkmate::qassert(samples, c("0","S*"))
  checkmate::qassert(population, c("0","S*"))
  
  checkmate::qassert(plotType, "S1")
  plotType = match.arg(plotType)
  
  checkmate::qassert(title, "S1")
  
  checkmate::qassert(metadata, "S1")
  checkmate::qassert(axes, "N2")
  
  if (!metadata %in% colnames(CYTdata@metadata)){
    stop("Error : 'metadata' argument is not a metadata (e.g. 'Timepoint', 'Individual')")
  }
  
  data = CYTdata@DimReduction@coordinates[,axes]
  if (nrow(data)==0) { stop("Error : Dimensionnality reduction is null, please compute dimensionnality reduction algorithm before vizualization") }
  if (ncol(data)!=2) { stop("Error : Impossible to visualize, embedding not in 2D") }
  
  
  if (plotType == "density") {
    if (!is.null(samples) || !is.null(population)) {
      stop("Error : 'plotType' argument is set to 'density' and requires the entire dataset to be computed ('samples' and 'population' arguments must both be set to NULL)")
    }
    
    colnames(data) = c("x", "y") 
    plot <- ggplot2::ggplot(data, ggplot2::aes(x,y)) +
      ggpointdensity::geom_pointdensity(size=0.001) +
      viridis::scale_color_viridis() +
      ggplot2::labs(color = "Cell density")
    legendPosition = "right"
  }
  else {
    plot <- ggplot2::ggplot()
    samples = checkorderSamples(CYTdata, samples, order=TRUE, checkDuplicates=TRUE)
    data = cbind.data.frame(data, "spl" = CYTdata@samples)
    if (!is.null(population)) {
      level = match.arg(level)
      checkmate::qassert(level, "S1")
      population = checkorderPopulation(CYTdata, population, level=level, order=TRUE, checkDuplicates=TRUE)
      if (level=="clusters") { popId = CYTdata@Clustering@clusters }
      else { popId = CYTdata@Metaclustering@metaclusters }
      colorIndexes = (data$spl %in% samples) & (popId %in% population)
    }
    else { colorIndexes = (data$spl %in% samples) }
    
    dataGrey = subset(data, !colorIndexes)
    if (nrow(dataGrey)>0) {
      plot = plot +
        ggplot2::geom_point(data = dataGrey,
                            ggplot2::aes_string(x = colnames(dataGrey)[1],
                                                y = colnames(dataGrey)[2]),
                            color = "gray", size = 0.001) +
        ggnewscale::new_scale_color()
    }
    data = subset(data, colorIndexes)
    
    if (plotType == "metadata") {
      md = CYTdata@metadata
      md$spl = rownames(md)
      dims = colnames(data)
      data = merge(data, md, by = "spl")
      colnames(data)[colnames(data) == metadata] = "metadata"
      plot <- plot +
        ggplot2::geom_point(data = data,
                            ggplot2::aes_string(x = dims[1],
                                                y = dims[2],
                                                color = "metadata"),
                            size = 0.001) +
        ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(alpha = 1, size = 2), ncol = 1, title = metadata))
      legendPosition = "right"
    }
  }
  
  if (!is.null(CYTdata@DimReduction@optional_parameters$type)){
    labs = paste(CYTdata@DimReduction@optional_parameters$type, axes, sep = "")
  } else { labs = paste("dim", axes, sep = "")}
  
  plot <- plot +
    ggplot2::xlab(labs[1]) + ggplot2::ylab(labs[2]) +
    ggplot2::ggtitle(title) +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal()
  
  return(plot)
}

#' @title Visualize 2D reduced cell data according to marker expression
#'
#' @description This function aims to visualize cells in the embedded 2D space, by coloring them with the expression of a given marker.
#' This representation can be used on a CYTdata object for which a dimensionality reduction operation has been computed.
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param marker a character vector containing the name of the biological marker to plot
#' @param population a character vector containing the identifiers of the population to plot. By default, all the population are used. Non included populations will still be ploted, but in grey.
#' @param level a character value indicating the type of population filtered by the population argument. Possible values are: "clusters", "metaclusters". By default, 'clusters' are used.
#' @param samples a character vector containing the names of the biological samples to plot. Non included samples will still be ploted, but in grey. 
#' @param bounds a numeric providing the two quantiles that will be used to bind the expression of the marker. Defaults to c(0.05, 0.95).
#' @param relativeGradient a boolean. If TRUE, marker expression will be bound according to the marker's expression in the selected populations and/or samples. If FALSE (default), the entire cell population will be used. 
#' @param paletteGradient a character vector containing the colors to use for the gradient scale. Defaults to a blue/red gradient. 
#' @param axes a numeric of length 2 indicating the axes to display. Defaults to 1 and 2. 
#'
#' @return a ggplot2 object
#'
#' @export
#'

plotDimRedGradient <- function(CYTdata,
                               marker,
                               population = NULL,
                               level = c("clusters", "metaclusters"),
                               samples = NULL,
                               bounds = c(0.05, 0.95),
                               relativeGradient = FALSE,
                               paletteGradient = NULL,#c("yellow", "orange", "red", "brown")
                               axes = c(1, 2)) { 
  
  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }
  
  checkmate::qassert(marker, "S1")
  marker = checkorderMarkers(CYTdata, marker, order=TRUE, checkDuplicates=TRUE)
  checkmate::qassert(axes, "N2")
  
  checkmate::qassert(samples, c("0","S*"))
  samples = checkorderSamples(CYTdata, samples, order=TRUE, checkDuplicates=TRUE)
  if (!is.null(population)) {
    level = match.arg(level)
    checkmate::qassert(level, "S1")
    population = checkorderPopulation(CYTdata, population=population, level=level,
                                      order=TRUE, checkDuplicates=TRUE)
    if(level == "clusters") { popId = CYTdata@Clustering@clusters }
    else { popId = CYTdata@Metaclustering@metaclusters }
    colorIdx = (CYTdata@samples %in% samples) & (popId %in% population)
  }
  else { colorIdx = CYTdata@samples %in% samples }
  
  checkmate::qassert(bounds, "N2")
  if (!all(bounds>=0 & bounds<=1)) {
    stop("Error : 'bounds' argument is a vector of quantile bounds and must be two positive integer between 0 and 1")
  }
  checkmate::qassert(relativeGradient, "B1")
  
  checkmate::qassert(paletteGradient, c("0", "S*"))
  if(!is.null(paletteGradient)){
    areColors <- function(x) { sapply(x, function(X) {
      tryCatch(is.matrix(col2rgb(X)), error = function(e) FALSE) }) }
    if(!all(areColors(paletteGradient))){
      stop("Error : New gradient palette (", paste0(paletteGradient, collapse = ","),
           "), does not contain only hexadecimal color.)")
    }
  }
  else { paletteGradient = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9,"RdBu")))(85) }
  
  data = CYTdata@DimReduction@coordinates[,axes]
  if (nrow(data) == 0) { stop("Error : Dimensionnality reduction is null, please compute dimensionnality reduction algorithm before vizualization") }
  if (ncol(data)!=2) { stop("Error : Impossible to visualize, embedding not in 2D") }
  
  if (relativeGradient) {
    dataColored = subset(data, colorIdx)
    dataGray = subset(data, !colorIdx)
    dataColored$value = CYTdata@matrix.expression[colorIdx, marker]
    limits = stats::quantile(dataColored$value, probs = bounds)
    dataColored$value[dataColored$value < limits[1]] = limits[1]
    dataColored$value[dataColored$value > limits[2]] = limits[2]
  }
  else {
    data$value = CYTdata@matrix.expression[, marker]
    limits = stats::quantile(data$value, probs = bounds)
    data$value[data$value < limits[1]] = limits[1]
    data$value[data$value > limits[2]] = limits[2]
    dataColored = subset(data, colorIdx)
    dataGray = subset(data, !colorIdx)
  }
  
  plot <- ggplot2::ggplot()
  
  if (nrow(dataGray) != 0){
    plot <- plot + ggplot2::geom_point(data = dataGray,
                                       ggplot2::aes_string(x = colnames(dataGray)[1],
                                                           y = colnames(dataGray)[2]),
                                       col = "grey", size = 0.001)
  }
  
  
  if (!is.null(CYTdata@DimReduction@optional_parameters$type)){
    labs = paste(CYTdata@DimReduction@optional_parameters$type, axes, sep = "")
  } else { labs = paste("dim", axes, sep = "")}
  
  plot <- plot + ggplot2::geom_point(data = dataColored,
                                     ggplot2::aes_string(x = colnames(dataColored)[1],
                                                         y = colnames(dataColored)[2],
                                                         col = "value"),
                                     size = 0.001) +
    ggplot2::scale_color_gradientn(colours = paletteGradient) +
    ggplot2::scale_x_continuous(expand = c(0.01,0.01)) +
    ggplot2::scale_y_continuous(expand = c(0.01,0.01)) +
    ggplot2::labs(title = paste0(marker,"'s gradient representation"), col = marker) +
    ggplot2::xlab(labs[1]) + ggplot2::ylab(labs[2]) +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal()
  
  return(plot)
}

#' @title Visualize 2D reduced cell data according to clustering or metaclustering
#'
#' @description This function aims to visualize cells in the embedded 2D space, by coloring them with the expression of a given marker.
#' This representation can be used on a CYTdata object for which a dimensionality reduction operation has been computed.
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param population a character vector containing the identifiers of the population to plot. By default, all the population are used. Non included populations will still be ploted, but in grey.
#' @param level a character value indicating the type of population filtered by the population argument. Possible values are: "clusters", "metaclusters". By default, 'clusters' are used.
#' @param samples a character vector containing the names of the biological samples to plot. Non included samples will still be ploted, but in grey. 
#' @param printCentroid a boolean. If TRUE (default), displays the name of each population at that population's centroid position. 
#' @param axes a numeric of length 2 indicating the axes to display. Defaults to 1 and 2. 
#'
#' @return a ggplot2 object
#'
#' @export
#'

plotDimRedPopulations <- function(CYTdata,
                                  population = NULL,
                                  level = c("clusters", "metaclusters"),
                                  samples = NULL,
                                  printCentroid = TRUE,
                                  axes = c(1, 2)) {
  
  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }
  
  level = match.arg(level)
  checkmate::qassert(level, "S1")
  population = checkorderPopulation(CYTdata, population = population, level = level,
                                    order=TRUE, checkDuplicates=TRUE)
  samples = checkorderSamples(CYTdata, samples, order=TRUE)
  checkmate::qassert(printCentroid, "B1")
  checkmate::qassert(axes, "N2")
  
  if(level == "clusters") {
    palette = CYTdata@Clustering@palette
    popId = CYTdata@Clustering@clusters
  }
  else {
    palette = CYTdata@Metaclustering@palette
    popId = CYTdata@Metaclustering@metaclusters
  }
  
  data = CYTdata@DimReduction@coordinates[,axes]
  if (nrow(data) == 0) { stop("Error : Dimensionnality reduction is null, please compute dimensionnality reduction algorithm before vizualization") }
  if (ncol(data)!=2) { stop("Error : Impossible to visualize, embedding not in 2D") }
  data$popId = popId
  colorIdx = (CYTdata@samples %in% samples) & (popId %in% population)
  dataColored = subset(data, colorIdx)
  dataGray = subset(data, !colorIdx)
  
  plot <- ggplot2::ggplot()
  
  if (nrow(dataGray) > 0){
    plot <- plot +
      ggplot2::geom_point(data = dataGray,
                          ggplot2::aes_string(x = colnames(dataGray)[1],
                                              y = colnames(dataGray)[2]), size = 0.01,
                          #col = "#DCDCDC")
                          col = "grey")
  }
  
  cat("\nPlotting 2D-reduced data and using palette slot to color the following", level, ": \n - ", paste0(population, collapse = ", "))
  
  plot <- plot +
    ggplot2::geom_point(data = dataColored,
                        ggplot2::aes_string(x = colnames(dataColored)[1],
                                            y = colnames(dataColored)[2],
                                            col = "popId"),
                        size = 0.01) +
    ggplot2::scale_color_manual(values = palette) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(alpha = 1, size = 3), ncol = 1, title = "Population"))
  
  if (printCentroid) {
    centroids = plyr::ddply(dataColored, "popId",
                            function(x){ return(robustbase::colMedians(data.matrix(x[,colnames(dataColored)[1:2]]))) })
    plot <- plot +
      ggrepel::geom_text_repel(data = centroids,
                               ggplot2::aes_string(x = colnames(dataColored)[1],
                                                   y = colnames(dataColored)[2],
                                                   label ="popId"),
                               color = "black", size = 3, force = 5)
  }
  
  if (!is.null(CYTdata@DimReduction@optional_parameters$type)){
    labs = paste(CYTdata@DimReduction@optional_parameters$type, axes, sep = "")
  } else { labs = paste("dim", axes, sep = "")}
  
  plot <- plot +
    ggplot2::xlab(labs[1]) + ggplot2::ylab(labs[2]) +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal()
  
  return(plot)
}