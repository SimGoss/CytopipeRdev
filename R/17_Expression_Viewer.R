#' @title Plot a boxplot of marker expression
#'
#' @description This function plots a CYTdata object's expression data.
#' Warning : VERY slow. 
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param markers a character vector containing the names of biological markers to use. By default, all markers are plotted
#' @param rawData a boolean value specifying if data come from raw.matrix.expression slot of CYTdata object or if data come from matrix.expression slot.
#' By default, data from matrix.expression is used
#'
#' @return a ggplot2 object
#'
#' @export
#'

expression.boxplot <- function(CYTdata, markers = NULL, rawData = FALSE){
  
  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  checkmate::qassert(markers, c("0", "S*"))
  checkmate::qassert(rawData, "B1")
  
  if (rawData){
    if (ncol(CYTdata@raw.matrix.expression)>0){
      exprs = CYTdata@raw.matrix.expression
      title = "Raw ion counts"
    }
    else {
      message("Warning : raw.matrix.expression slot is empty, so data plotted in scatterplot are data from 'matrix.expression' dataframe")
      exprs = CYTdata@matrix.expression
      title = "Transformed ion counts"
    }
  }
  else {
    exprs = CYTdata@matrix.expression
    title = "Transformed ion counts"
  }
  
  if (!is.null(markers)){ exprs = exprs[,markers] }
  
  suppressMessages({ melt.exprs = reshape2::melt(exprs) })
  plot <- ggplot2::ggplot(melt.exprs, aes(y = value, x = variable)) +
    ggplot2::geom_boxplot(outlier.size=0.05) +
    ggplot2::ggtitle(title) + ggplot2::xlab("Markers") + ggplot2::ylab(title) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = element_text(hjust = 0.5)) +
    ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 0.5))
  
  return(plot)
}



#' @title Plot scatterplot of 2 marker's expression
#'
#' @description This function permits to plot a CYTdata object's expression data over scatterplot with 2 markers choosen
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param coupleMarkers a character vector providing the marker represented on each axis (axis x then axis y)
#' @param samples a character vector containing the names of biological samples to use. By default, all samples are used
#' @param level a character vector indicating whether to select "clusters" (default) or "metaclusters" for the population argument.
#' @param population a character vector indicating the identifiers of the clusters/metaclusters to compare. Defaults to all of them. 
#' @param rawData a boolean value specifying if data comes from raw.matrix.expression slot of CYTdata object or if data comes from matrix.expression slot. By default, data from matrix.expression is used
#' @param colorBy a character value providing the kind of coloration to apply to the plot. Possible values are : none (default), "clusters", "metaclusters", "density", "metadata". "density" will display a density plot. 
#' @param colorMetadata a character value providing the name of the metadata according to which the points should be colors if colorBy = "metadata".
#' @param neighboursDensity a numeric value providing the size of the grid for density calculations (effectively, the breaks of a 2D histogram). Defaults to 100. 
#' @param colorNone a character value providing the name or hexadecimal code for point color if colorBy = "none".
#' @param plotLinearRegression logical. If TRUE, a linear regression is performed and displayed on the plot. Please note no test is performed. Defaults to FALSE. 
#' @param downsamplingNumber a numeric value providing a downsampling target for displaying the points. Defaults to NULL, which is not recommended for large CYTdata objects. 
#' @param ylim a numeric of length 2 providing lower and upper bounds for the y axis
#' @param xlim a numeric of length 2 providing lower and upper bounds for the x axis
#'
#' @return a ggplot2 object
#'
#' @export
#'

markersScatterplot <- function(CYTdata,
                               coupleMarkers,
                               samples = NULL,
                               level = c("clusters", "metaclusters"),
                               population = NULL,
                               rawData = FALSE,
                               colorBy = c("none", "clusters", "metaclusters", "density", "metadata"),
                               colorMetadata = "Timepoint",
                               # colorDensity = NULL,
                               neighboursDensity = 100,
                               colorNone = NULL,
                               plotLinearRegression = FALSE,
                               downsamplingNumber = NULL,
                               ylim = NULL,
                               xlim = NULL){
  
  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }
  
  colorBy = match.arg(colorBy)
  checkmate::qassert(colorBy, "S1")
  
  if ((colorBy=="clusters" && length(CYTdata@Clustering@clusters)==0) ||
      (colorBy=="metaclusters" && length(CYTdata@Metaclustering@metaclusters)==0)) {
    stop("Error : 'colorBy' argument set to ", colorBy, ", but CYTdata's ", colorBy, "slot is empty.
         Please perform a previous clustering/metaclustering step before.")
  }
  
  checkmate::qassert(rawData, "B1")
  if (rawData){
    if (ncol(CYTdata@raw.matrix.expression)>0){ exprs = CYTdata@raw.matrix.expression }
    else { stop("Error : 'rawData' argument is set to TRUE but raw.matrix.expression slot is empty.") }
  }
  else { exprs = CYTdata@matrix.expression }
  
  checkmate::qassert(coupleMarkers, "S2")
  coupleMarkers = checkorderMarkers(CYTdata, coupleMarkers, order=FALSE, checkDuplicates=FALSE)
  data = cbind.data.frame(exprs[,coupleMarkers], "spl" = CYTdata@samples)
  colnames(data) = c("x", "y", "spl")
  
  checkmate::qassert(samples, c("0", "S*"))
  samples = checkorderSamples(CYTdata, samples, order=TRUE, checkDuplicates=TRUE)
  
  level = match.arg(level)
  checkmate::qassert(level, "S1")
  if (level=="clusters") { popId = CYTdata@Clustering@clusters }
  else { popId = CYTdata@Metaclustering@metaclusters }
  
  if (!is.null(population)) {
    if (length(popId)>0) {
      population = checkorderPopulation(CYTdata, population, level=level, order=TRUE, checkDuplicates=TRUE)
    }
    else {
      stop("Error : 'population' argument given and 'level' argument set to ", level, ". But CYTdata's ", level, " slot is empty.
           Please perform a previous clustering/metaclustering step before.")
    }
    subsetIdx = CYTdata@samples %in% samples & popId %in% population
    
  }
  else {
    subsetIdx = CYTdata@samples %in% samples
  }
  
  checkmate::qassert(downsamplingNumber, c("0", "N1"))
  if (!is.null(downsamplingNumber)) {
    
    subsetIdx = tapply(subsetIdx,
                       CYTdata@samples,
                       function(x) {
                         if (sum(x)<downsamplingNumber) { return(x) }
                         else {
                           y = rep(FALSE, length(x))
                           y[sample(which(x), downsamplingNumber, replace = FALSE)] = TRUE
                           return(y)
                         }
                       }) %>% unlist() %>% as.vector()
  }
  
  data = subset(data, subsetIdx)
  checkmate::qassert(plotLinearRegression, "B1")
  checkmate::qassert(ylim, c("0", "N2"))
  checkmate::qassert(xlim, c("0", "N2"))
  
  switch(colorBy,
         none = {
           checkmate::qassert(colorNone, c("0", "S*"))
           if (!is.null(colorNone)) { if(!areColors(colorNone)){ stop("Error : 'colorNone' argument is not a hexadecimal color.") } }
           else { colorNone = "#1E50BB" }
           plot <- ggplot2::ggplot(data) +
             ggplot2::geom_point(ggplot2::aes(x, y), size = 0.005, color = colorNone)
         },
         clusters = {
           data$Population = CYTdata@Clustering@clusters[subsetIdx]
           plot <- ggplot2::ggplot(data) +
             ggplot2::geom_point(ggplot2::aes(x, y, color = Population), size = 0.85) +
             ggplot2::scale_color_manual(values = CYTdata@Clustering@palette) +
             ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3)))
         },
         metaclusters = {
           data$Population = CYTdata@Metaclustering@metaclusters[subsetIdx]
           plot <- ggplot2::ggplot(data) +
             ggplot2::geom_point(ggplot2::aes(x, y, color = Population), size = 0.85) +
             ggplot2::scale_color_manual(values = CYTdata@Metaclustering@palette) +
             ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3)))
         },
         density = {
           checkmate::qassert(neighboursDensity, "N1")
           get_density <- function(x, y, n) {
             dens <- MASS::kde2d(x, y, n=n)
             ix <- findInterval(x, dens$x)
             iy <- findInterval(y, dens$y)
             ii <- cbind(ix, iy)
             return(dens$z[ii])
           }
           #checkmate::qassert(colorDensity, c("0", "S*"))
           #if(!is.null(colorDensity)){
           #   if(!all(areColors(colorDensity))){
           #       stop("Error : 'colorDensity' argument is a color palette (",
           #            paste0(colorDensity, collapse = ","), "), does not contain only hexadecimal color.") }
           #}
           #else { colorDensity = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11,'Spectral')))(11) }
           data$density <- get_density(data$x, data$y, n = neighboursDensity)
           plot <- ggplot2::ggplot(data) +
             ggplot2::geom_point(ggplot2::aes(x, y, color = density), size = 0.005) +
             #ggplot2::scale_color_manual(values = colorDensity) +
             viridis::scale_color_viridis() +
             ggplot2::labs(color = "Cell density")
         },
         metadata = {
           checkmate::qassert(colorMetadata, "S1")
           if (!colorMetadata %in% colnames(CYTdata@metadata)){
             stop("Error : 'colorMetadata' argument is not a metadata (ex : 'Timepoint', 'Individual')")
           }
           md = CYTdata@metadata
           md$spl = rownames(md)
           data = merge(md, data, by = "spl")
           data[[colorMetadata]] = droplevels(data[[colorMetadata]])
           colnames(data)[colnames(data)==colorMetadata] = "colorMetadata"
           plot <- ggplot2::ggplot(data = data) +
             ggplot2::geom_point(ggplot2::aes(x, y, color = colorMetadata), size = 0.005) +
             ggplot2::guides(color = ggplot2::guide_legend(title = colorMetadata, override.aes = list(size = 3)))
         })
  if (plotLinearRegression){
    plot <- plot +
      ggplot2::geom_smooth(ggplot2::aes(x, y, colour="linear", fill="linear"),
                           method="lm",
                           formula=y ~ x)
  }
  if (!is.null(ylim)){ plot <- plot + ggplot2::ylim(ylim[1], ylim[2]) }
  if (!is.null(xlim)){ plot <- plot + ggplot2::xlim(xlim[1], xlim[2]) }
  plot <- plot +
    #ggplot2::geom_hline(yintercept = asinh(20/5), lw = 1) +
    ggplot2::labs(title=paste(coupleMarkers[1], "vs", coupleMarkers[2], "comparison", sep=" ")) +
    ggplot2::theme_minimal() +
    ggplot2::ylab(coupleMarkers[2]) +
    ggplot2::xlab(coupleMarkers[1]) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   panel.background = ggplot2::element_blank(),
                   panel.border = ggplot2::element_rect(fill=NA),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   legend.position = "right")
  
  return(plot)
}

#' @title Plot distribution of marker expression
#'
#' @description This function aims to plot the distribution of marker expression for selected samples. Marker expression can be further divided according to any kind of metadata.
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param markers a character vector containing the names of a biological marker to use. Currently, only one of them can be used. 
#' @param samples a character vector containing the names of biological samples to use. By default, all samples are used
#' @param group.metadata a character value containing the names of the metadata to use to split the plot. Possible values are colnames of CYTdata's metadata slot. by default, it is set to NULL, which means that markers densities are not grouped by metadata condition.
#'
#' @return void
#'
#' @export
#'

markerDensity <- function(CYTdata,
                          markers,
                          samples = NULL,
                          group.metadata = NULL) {
  
  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.")}
  checkmate::qassert(group.metadata, c("0", "S1"))
  checkmate::qassert(markers, "S1")
  checkmate::qassert(samples, c("0", "S*"))
  
  if (is.null(markers)){ markers = CYTdata@markers }
  if(is.null(samples)){ samples = unique(CYTdata@samples) }
  data = cbind.data.frame(CYTdata@samples, CYTdata@matrix.expression[,markers])
  colnames(data) = c("sample", markers)
  data = subset(data, sample %in% samples)
  
  if(!is.null(group.metadata)){
    if (!group.metadata %in% colnames(CYTdata@metadata)){ stop("Error in metadata.densities function : 'group.metadata' argument is not a
    biological condition present in metadata (", paste0(colnames(CYTdata@metadata), collapse=","), ")") }
    md = CYTdata@metadata
    md$sample = rownames(md)
    data.merged = merge(data, md, by = "sample")
    
    for (mark in markers){
      data.sub = data.frame(mark.ex = data.merged[,mark],
                            group.ex = factor(data.merged[,group.metadata]))
      plot = ggplot2::ggplot(data.sub,
                             ggplot2::aes(x = mark.ex, y = group.ex, fill = group.ex)) +
        ggridges::geom_density_ridges() + ggridges::theme_ridges() +
        ggplot2::xlab(mark) + ggplot2::ggtitle(paste("Expression distribution for each", group.metadata, sep=" "))
    }
  }
  else{
    exprs = data[,markers]
    if (length(markers)==1){
      exprs.melted = data.frame("Markers" = rep(markers, length(exprs)), "Value" = exprs)
    }
    else {
      suppressMessages({ exprs.melted = reshape2::melt(exprs) })
      colnames(exprs.melted) = c("Markers", "Value")
    }
    plot <- ggplot2::ggplot(exprs.melted, aes(x=Value, color=Markers, fill=Markers)) +
      ggplot2::geom_density(alpha=0.6) +
      ggplot2::ylab("Distribution") + ggplot2::xlab("") +
      facet_wrap(~Markers, scales="free", ncol=8) +
      ggplot2::ggtitle(title) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "none",
        panel.spacing = unit(0.1, "lines"),
        strip.text.x = ggplot2::element_text(size = 9),
        axis.text.y = ggplot2::element_text(size=5.5),
        axis.text.x = ggplot2::element_text(size=5.5),
        strip.background = ggplot2::element_rect(colour="black", fill="white", size=0.5, linetype="solid"),
        plot.title = ggplot2::element_text(hjust = 0.5))
  }
  
  return(plot)
}

