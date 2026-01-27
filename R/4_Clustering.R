#' @title Import object of class 'Clustering' to put in CYTdata object
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param importedClusteringObject a S4 object of class 'Clustering' to put in CYTdata object
#' @param checkOverwrite a boolean value indicating whether to check if a Metaclustering has already been performed on the CYTdata object. 
#' If TRUE, it will check and ask for confirmation, if FALSE, any previous clustering will be overwritten by default. 
#' 
#' @return a S4 object of class 'CYTdata'
#'
#' @export
#'

importClustering <- function(CYTdata, importedClusteringObject, checkOverwrite=TRUE){
  
  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }
  
  if (class(importedClusteringObject)!="Clustering") {
    stop("Error : argument 'importedClusteringObject' must be a S4 object of class 'Clustering'.")
  }
  
  if (checkOverwrite && length(CYTdata@Clustering@clusters)!=0){
    reply <- readline(prompt="Clustering already performed, do you still want to continue and overwrite (yes or no): ")
    while (!tolower(reply) %in% c("yes", "y", "no", "n")) { reply <- readline(prompt="Reply not standard, please answer yes or no: ") }
    if (tolower(reply) %in% c("no", "n")){ stop("Function stopped, CYTdata object unchanged") }
  }
  
  message("\nUpdating CYTdata object...")
  CYTdata@Clustering = importedClusteringObject
  #CYTdata = MakeValid(CYTdata, verbose = TRUE)
  return(CYTdata)
}

#' @title Exports an  object of class 'Clustering' from a CYTdata object
#' 
#' @description Exports an  object of class 'Clustering' from a CYTdata object. Equivalent to "CYTdata@Clustering"
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#'
#' @return a S4 object of class 'Clustering'
#'
#' @export
#'

exportClustering <- function(CYTdata){
  validObject(CYTdata)
  return(CYTdata@Clustering)
}

#' @title Resets the 'Clustering' slot from a CYTdata object
#' 
#' @description The Clustering slot of a CYTdata object will be reset to an empty Clustering S4 object. 
#' The CYTdata object is validated (makeValid function) before and after removal. 
#' 
#' @param CYTdata a S4 object of class 'CYTdata'
#'
#' @return a S4 object of class 'CYTdata'
#'
#' @export
#'

removeClustering <- function(CYTdata) {
  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }
  
  message("\nUpdating CYTdata object...")
  CYTdata@Clustering = methods::new("Clustering")
  validObject(CYTdata)
  return(CYTdata)
}


#' @title Regroup cells with similar marker expressions using clustering techniques
#'
#' @description This function aims to regroup cells into clusters of cells with similar marker expression, based on a specific set of markers. 
#'
#' Several algorithms are available, but as of 01/2026, only SOM has been thoroughly tested :
#' - Self-Organizing Maps (SOM), a neural network-based clustering technique (som method of the kohonen package)
#' - FlowSOM, self-organizing maps with the FlowSOM package (FlowSOM method). Also performs Metaclustering.
#' - Phenograph, a graph-based clustering technique designed for high-dimensional single cells, with cytofkit2 (Rphenograph method, requires Python)
#' - SPADE (Spanning-tree Progression Analysis of Density-normalized Events), see Qiu et al. 2011
#' - DBSCAN, a commonly used density-based clustering algorithm used for various kind of data (dbscan package)
#' - Kmeans, which clusters the cells in n clusters based on the distance between the cells and the centroid of each cluster (stats package)
#' - Kmedoids, more robust version of Kmeans (pam function of the cluster package)
#' - clara, similar to Kmedoids, but deals with large datasets better (clara function of the cluster package)
#' - "PhenoGMM", "Mean_shift", "flowMeans" were planned, but not implemented yet. 
#' 
#' The whole set of cell markers or specific cell markers can be used during the clustering process.
#' 
#' The cell clustering can be performed on marker expression (default), or on the coordinates obtained with a dimension reduction method. 
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param markers a character vector providing the cell markers to use to generate the cluster data. By default, all markers are used
#' @param type a character value containing the type of clustering method to use. 
#' Possible values are: "SOM" (default), "FlowSOM", "Phenograph", "Spade", "DBSCAN", "Kmeans", "Kmedoids",  "CLARA" (, "PhenoGMM", "Mean_shift", "flowMeans").
#' @param scaledata boolean. If TRUE, centers (mean = 0) and scales (sd = 1) the expression of each marker
#' @param checkOverwrite a boolean value indicating whether to check if a Metaclustering has already been performed on the CYTdata object. 
#' If TRUE, it will check and ask for confirmation, if FALSE, any previous clustering will be overwritten by default. 
#' @param manifold a boolean value. If TRUE, dimension reduction coordinates are used instead of cell marker expression to perform the clustering. 
#' Requires a dimension reduction to have been performed and stored in the CYTdata object. Defaults to FALSE. 
#' @param seed a numeric value providing the random seed to use during stochastic operations
#' @param ... additional arguments passed on to method from R package.
#'
#' @return a S4 object of class 'CYTdata'
#'
#' @export

runClustering <- function(CYTdata,
                          markers = NULL,
                          type = c("SOM", "FlowSOM", "Phenograph", "Spade", "DBSCAN",
                                   "Kmeans", "Kmedoids", "CLARA",
                                   "PhenoGMM", "Mean_shift", "flowMeans"),
                          scaledata = FALSE,
                          checkOverwrite = TRUE,
                          manifold = FALSE,
                          seed = 42,
                          ...) {
  
  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }
  
  checkmate::qassert(markers, c(0,"S*"))
  type = match.arg(type)
  checkmate::qassert(type, "S1")
  checkmate::qassert(seed, "N1")
  
  if (checkOverwrite && length(CYTdata@Clustering@clusters)!=0){
    reply <- readline(prompt="Clustering already performed, do you still want to continue and overwrite (yes or no): ")
    while (!tolower(reply) %in% c("yes", "y", "no", "n")) { reply <- readline(prompt="Reply not standard, please answer yes or no: ") }
    if (tolower(reply) %in% c("no", "n")){ stop("Function stopped, CYTdata object unchanged") }
  }
  
  markers = checkorderMarkers(CYTdata, markers, order=TRUE, checkDuplicates=TRUE)
  
  
  if (manifold) {
    data = CYTdata@DimReduction@coordinates
  }
  else {
    data = CYTdata@matrix.expression[,markers, drop=FALSE]
    if (scaledata) { data = data %>% mutate_all(scale) }
  }
  
  cat("\nClustering data into cell clusters using ", type, " :")
  cat("\n\n - Algorithm is performed on the expression of the following markers : \n",
      paste0(markers, collapse=", "), "\n\n")
  parameters = append(list("type" = type,
                           "markers" = markers,
                           "seed" = seed),
                      list(...))
  
  switch (type,
          SOM = {
            t = system.time(somObject <- kohonen::som(as.matrix(data), ...))
            clusters <- somObject$unit.classif
            plots = list(somObject)
          },
          FlowSOM = {
            t = system.time(fSOMObject <- FlowSOM::FlowSOM(as.matrix(data), ...))
            clusters <- FlowSOM::GetClusters(fSOMObject)
            plots = list(fSOMObject)
          },
          Phenograph = {
            t = system.time(phenographObject <- cytofkit2::Rphenograph(data, ...))# phenograph_homemade(data, ...))
            clusters = phenographObject$membership
            plots = list(phenographObject)
          },
          Spade = {
            t = system.time(spadeObject <- spade_homemade(data, ...))
            clusters = spadeObject$clustering
            plots = list(spadeObject)
          },
          DBSCAN = {
            t = system.time(clusters <- dbscan::dbscan(data, ...)$cluster)
          },
          Kmeans = {
            t = system.time(clusters <- stats::kmeans(data, ...)$cluster)
          },
          Kmedoids = {
            t = system.time(clusters <- cluster::pam(data, ...)$clustering)
          },
          clara = {
            t = system.time(clusters <- cluster::clara(data, ...)$clustering)
          },
          PhenoGMM = {
            stop()
          },
          Mean_shift = {
            stop()
          },
          flowMeans = {
            stop()
          })
  
  cat("\n", type, "completed in", round(t[3], 3) ,"seconds.\n")
  
  if( length(clusters) != nrow(data) ){
    stop("Error : Clustering is not complete. Some events are not assigned.
            Please try other clustering methods.")
  }
  
  clusters = factor(as.vector(clusters))
  palette = rainbow(nlevels(clusters))
  names(palette) = levels(clusters)
  
  # Put into "Clustering" object
  cat("\n\nCreating Clustering object :")
  cat("\n - Computing cell cluster count matrix...")
  cellcount = compute.cellcount(clusters, CYTdata@samples)
  cat("\n - Computing cell cluster abundance matrix...")
  abundance = compute.abundance(cellcount)
  
  cat("Creating Clustering object.. \n")
  Clustering.object <- methods::new("Clustering",
                                    clusters = clusters,
                                    cellcount = cellcount,
                                    abundance = abundance,
                                    palette = palette,
                                    optional_parameters = parameters)
  
  if (type %in% c("FlowSOM", "SOM", "Phenograph", "Spade")) {
    Clustering.object@optional_plots = plots
  }
  
  cat("Updating CYTdata object.. \n")
  CYTdata@Clustering = Clustering.object
  
  if (type == "FlowSOM") {
    metaclusters = factor(as.vector(FlowSOM::GetMetaclusters(fSOMObject)))
    palette = rainbow(nlevels(metaclusters))
    names(palette) = levels(metaclusters)
    
    # Put into "Clustering" object
    cat("\n\nCreating Metaclustering object :")
    cat("\n - Computing cell metacluster count matrix...")
    cellcount = compute.cellcount(metaclusters, CYTdata@samples)
    cat("\n - Computing cell metacluster abundance matrix...")
    abundance = compute.abundance(cellcount)
    
    cat("Creating Metaclustering object.. \n")
    Metaclustering.object <- methods::new("Metaclustering",
                                          metaclusters = metaclusters,
                                          cellcount = cellcount,
                                          abundance = abundance,
                                          palette = palette,
                                          optional_parameters = parameters)
    cat("Updating CYTdata object.. \n")
    CYTdata@Metaclustering = Metaclustering.object
  }
  
  validObject(CYTdata)
  cat("Clustering done!\n")
  return(CYTdata)
}

#' @title Internal - Computes the number of cells for each population
#'
#' @description This function is used internally to computes the number of cells for each population.
#'
#' @param population a character vector providing the population to count the cells belonging to
#' @param samples a character vector providing for each cell the associated biological sample
#'
#' @return a dataframe containing the numbers of cells associated for each population for each sample (rownames = population / colnames = samples)
#'
#'
compute.cellcount <- function(population, samples) {
  
  checkmate::qassert(population, "F*")
  checkmate::qassert(samples, "F*")
  
  if(length(unique(samples))==1) {
    table = as.matrix(table(population))
    colnames(table) = unique(samples)
    return(as.data.frame(table))
  }
  else {
    table = sapply(levels(population), function(x){ return(table(samples[population == x])) })
    return(data.frame(t(table), check.names = FALSE))
  }
}

#' @title Internal - Computes the relative abundances for each cell population
#'
#' @description This function is used internally to computes the abundance of each population for each sample.
#'
#' @param cellcount a data.frame providing the numbers of cells associated to each population for each sample
#'
#' @return a data.frame containing the abundance of cells to each population for each sample
#'
compute.abundance <- function(cellcount) {
  
  checkmate::qassert(cellcount, "D*")
  
  abundance = apply(cellcount, 2, function(df){df / sum(df)}) * 100
  if(nrow(cellcount)==1){
    abundance = t(abundance)
    rownames(abundance) = rownames(cellcount)
  }
  abundance = data.frame(abundance, check.names = FALSE)
  return(abundance)
}





