#' @title Import object of class 'Metaclustering' to put in CYTdata object
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param importedMetaclusteringObject a S4 object of class 'Metaclustering' to put in CYTdata object
#' @param checkOverwrite a boolean value indicating whether to check if a Metaclustering has already been performed on the CYTdata object. 
#' If TRUE, it will check and ask for confirmation, if FALSE, any previous clustering will be overwritten by default. 
#' @param val a boolean value. If TRUE, performs validation of the CYTdata object (makeValid function) after importing the Metaclustering. 
#' @return a S4 object of class 'CYTdata'
#'
#' @export
#'

importMetaclustering <- function(CYTdata, importedMetaclusteringObject, checkOverwrite=TRUE, val  = TRUE){
  
  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }
  
  if (class(importedMetaclusteringObject)!="Metaclustering") {
    stop("Error : argument 'importedMetaclusteringObject' must be a S4 object of class 'Metaclustering'.")
  }
  
  if (checkOverwrite && length(CYTdata@Metaclustering@metaclusters)!=0){
    reply <- readline(prompt="Metaclustering already performed, do you still want to continue and overwrite (yes or no): ")
    while (!tolower(reply) %in% c("yes", "y", "no", "n")) { reply <- readline(prompt="Reply not standard, please answer yes or no: ") }
    if (tolower(reply) %in% c("no", "n")){ stop("Function stopped, CYTdata object unchanged") }
  }
  
  message("\nUpdating CYTdata object...")
  CYTdata@Metaclustering = importedMetaclusteringObject
  if (val) { CYTdata = MakeValid(CYTdata, verbose = TRUE) }
  return(CYTdata)
}


#' @title Resets the 'Metaclustering' slot from a CYTdata object
#' 
#' @description The metaclustering slot of a CYTdata object will be reset to an empty Metaclustering S4 object. 
#' The CYTdata object is validated (makeValid function) before and after removal. 
#' 
#' @param CYTdata a S4 object of class 'CYTdata'
#'
#' @return a S4 object of class 'CYTdata'
#'
#' @export
#'

removeMetaclustering <- function(CYTdata) {
  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }
  
  message("\nUpdating CYTdata object...")
  CYTdata@Metaclustering = methods::new("Metaclustering")
  validObject(CYTdata)
  return(CYTdata)
}


#' @title Internal - Rescales marker expression by quantile
#'
#' @description This function is used internally to rescale the marker expression by the quantile method.
#' Two quantiles are computed, and every expression value below the lower quantile or above the upper quantile is cropped to the respective quantile. Then, each expression value is rescaled from 0 to 1. 
#'
#' @param exprs a vector containing one marker expression
#' @param rescaleBounds a numeric value providing the value of the first quantile and last quantile
#'
#' @return a numeric containing quantile rescale marker expressions
#' 
#' @export
#'

rescale.exprs <- function(exprs,
                          rescaleBounds){
  
  checkmate::qassert(rescaleBounds, c("N2"))
  
  quantiles = stats::quantile(exprs, probs = rescaleBounds)
  low = quantiles[1]
  high = quantiles[2]
  
  exprs[exprs<low] = low
  exprs[exprs>high] = high
  
  exprs.scaled = scales::rescale(exprs, from = c(low, high), to = c(0,1))
  
  return(exprs.scaled)
}


#' @title Computes median or mean for a set of marker expression
#'
#' @description Computes median or mean expression for a set of given markers, per cluster or per metacluster
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param markers a vector object listing marker names
#' @param population a character value containing the identifiers of the population(s) (clusters or metaclusters) to study
#' @param level a character value indicating whether to compute MSI on the clusters or on the metaclusters. Possible values are: "clusters", "metaclusters". By default, 'clusters' are used.
#' @param typeMSI a character value to choose between "median" (default) and "mean".
#' @param rescaleBounds a numeric value providing the value of the first quantile and last quantile for the rescaled values, to be passed to rescale.exprs()
#' @param scale Logical. If TRUE, rescales the value of each marker's expression to the bounds indicated in rescaleBounds
#' 
#' @return a dataframe
#' 
#' @export
#'

computeMSI <- function(CYTdata,
                       markers = NULL,
                       population = NULL,
                       level = c("clusters", "metaclusters"),
                       typeMSI = c("median", "mean"),
                       rescaleBounds = c(0.05, 0.95),
                       scale = TRUE) {
  
  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }
  
  level = match.arg(level)
  checkmate::qassert(level, "S1")
  population = checkorderPopulation(CYTdata, population = population, level = level, order = TRUE, checkDuplicates = TRUE)
  markers = checkorderMarkers(CYTdata, markers, order = TRUE, checkDuplicates = TRUE)
  
  typeMSI = match.arg(typeMSI)
  checkmate::qassert(typeMSI, "S1")
  checkmate::qassert(scale, "B1")
  checkmate::qassert(rescaleBounds, c("N2"))

  if (level == "clusters"){ popId = CYTdata@Clustering@clusters }
  else { popId = CYTdata@Metaclustering@metaclusters }
  
  # Scale step's function
  rescale = function(x) { return(rescale.exprs(x, rescaleBounds)) }
  
  # Medoid step's function
  medoid = if (typeMSI == "median") robustbase::colMedians else colMeans
  
  Exprs = subset(CYTdata@matrix.expression[, markers, drop=FALSE], popId %in% population)
  
  if (scale) {
    
    Exprs = apply(Exprs, 2, rescale)
    
  }
  
  data = cbind.data.frame("popId" = popId[popId %in% population],
                          Exprs)
  
  MSI = plyr::ddply(data, "popId",
                    function(x) {
                      y = medoid(data.matrix(x[,-1]))
                      names(y) = markers
                      return(y)}) %>% remove_rownames() %>% column_to_rownames("popId")
  return(MSI)
}


#' @title Perform manual metaclustering of the clusters
#'
#' @description This function takes of list of metacluster with their associated clusters to integrate metaclustering to a CYTdata object. Useful if all clusters have been manually assigned to a specific cell population that can be stored in the Metaclustering slot of the CYTdata object. 
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param belongMetaclustering A list of manually set metaclusters with their associated cluster names.
#' For example : belongMetaclustering = list("Monocytes" = c("cluster1","cluster3","cluster4","cluster5"), Lymphocytes = c("cluster2","cluster6","cluster7"), etc)
#' @param completeWith a character value indicating how unassigned clusters should be treated. 
#' If "clusters", unassigned clusters will each be assigned to a metacluster named after them. 
#' If "metaclusters",unassigned clusters will keep their previous metacluster assignment (requires metaclustering to have been completed before).
#' @param checkOverwrite a boolean value indicating whether to check if a Metaclustering has already been performed on the CYTdata object. If TRUE, it will check and ask for confirmation, if FALSE, any previous clustering will be overwritten by default. 
#'
#' @return a S4 object of class 'CYTdata'
#'
#' @export
#'

manualMetaclustering <- function(CYTdata,
                                 belongMetaclustering,
                                 completeWith = c("clusters", "metaclusters"),
                                 checkOverwrite = TRUE) {
  
  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }
  
  checkmate::qassert(belongMetaclustering, "L*")
  completeWith = match.arg(completeWith)
  checkmate::qassert(completeWith, "S1")
  checkmate::qassert(checkOverwrite, "B1")
  
  names(belongMetaclustering) = as.character(names(belongMetaclustering))
  if(sum(sapply(names(belongMetaclustering), nchar)==0)>0) { stop("Error : 'belongMetaclustering' list argument contains empty names.")}
  belongMetaclustering = lapply(names(belongMetaclustering), function(x) {
    Y = unlist(belongMetaclustering[x])
    return(structure(rep(x, length(Y)), names = as.character(Y)))}) %>% unlist()
  
  if (length(CYTdata@Clustering@clusters)==0){
    stop("Error : Clustering slot is empty. Please compute clustering step before running metaclustering step.")
  }
  if(sum(is.na(belongMetaclustering))>0) {
    stop("Error : 'belongMetaclustering' argument contains NA values.
         Associated metaclusters names are missing for clusters : ",
         paste0(names(belongMetaclustering)[is.na(belongMetaclustering)], collapse=", "), ".")
  }
  if(sum(is.na(names(belongMetaclustering)))>0) {
    stop("Error : 'belongMetaclustering' argument's names contain NA values.
         It must be a fully named vector. Associated clusters names are missing
         for metaclusters : ",
         paste0(belongMetaclustering[is.na(names(belongMetaclustering))], collapse=", "), ".")
  }
  dupBelong = names(belongMetaclustering)[duplicated(names(belongMetaclustering))]
  if (length(dupBelong)>0) {
    stop("Error : 'belongMetaclustering' argument contain duplicated names (clusters).
         Several metaclusters are associated to same clusters ( ",
         paste0(dupBelong, collapse=", "), " ).")
  }
  naIdx = setdiff(names(belongMetaclustering), levels(CYTdata@Clustering@clusters))
  if (length(naIdx)>0) {
    stop("Error : 'belongMetaclustering' argument contain names non present in CYTdata clusters levels : ",
         paste0(naIdx, collapse = ", "))
  }
  naIdx = setdiff(levels(CYTdata@Clustering@clusters), names(belongMetaclustering))
  if (length(naIdx)==0) {
    message("Metaclustering identifiers do not need to be completed with other identifiers. belongMetaclustering permits to assign all the clusters.")
    completeWith = "clusters"
  }
  
  if (checkOverwrite && length(CYTdata@Metaclustering@metaclusters)!=0) {
    reply <- readline(prompt="Metaclustering already performed, do you still want to continue and overwrite (yes or no): ")
    while (!tolower(reply) %in% c("yes", "y", "no", "n")) { reply <- readline(prompt="Reply not standard, please answer yes or no: ") }
    if (tolower(reply) %in% c("no", "n")){ stop("Function stopped, CYTdata object unchanged") }
  }
  
  if (completeWith == "clusters") {
    metaclusters = do.call(dplyr::recode, c(list(CYTdata@Clustering@clusters), as.list(belongMetaclustering)))
    optional_parameters = list("manualMetaclustering" = "completeWithClusters")
  }
  else {
    if (length(CYTdata@Metaclustering@metaclusters)==0){
      stop("Error : manualMetaclustering function with 'completeWith' argument set to 'metaclusters' requires unempty Metaclustering slot, which is not the case.
           Please compute metaclustering step.")
    }
    
    
    metaclusters = as.character(CYTdata@Metaclustering@metaclusters)
    clusters = as.character(CYTdata@Clustering@clusters)
    metaclusters[clusters %in% names(belongMetaclustering)] = belongMetaclustering[clusters[clusters %in% names(belongMetaclustering)]]
    metaclusters = as.factor(metaclusters)
    optional_parameters = append(CYTdata@Metaclustering@optional_parameters, list("manualMetaclustering" = "completeWithMetaclusters"))
  }
  
  cellcount = compute.cellcount(metaclusters, CYTdata@samples)
  abundance = compute.abundance(cellcount)
  palette = structure(rainbow(nlevels(metaclusters)), names = levels(metaclusters))
  CYTdata@Metaclustering = methods::new("Metaclustering",
                                        metaclusters = metaclusters,
                                        cellcount = cellcount,
                                        abundance = abundance,
                                        palette = palette,
                                        optional_parameters = optional_parameters)
  
  validObject(CYTdata)
  return(CYTdata)
}


#' @title Performs hierarchical classification of the clusters (metaclustering) 
#'
#' @description This function aims to identify metaclusters, which are groups of clusters having similar expressions for selected markers. 
#' The mean or median marker expressions is computed for each cluster, and used to compute hierarchical clustering using stats::hclust
#' Several clustering criteria are available, please refer to the function 'hclust' of the 'stats' package for details
#' Data will automatically be scaled by quantiles, by default 0.05 and 0.95
#' 
#' If N_NewClusters > 0, a second level clustering will be performed, where pre-existing clusters will be clustered into N_NewClusters new clusters, stored in the clustering slot of the CYTdata object (previous clusters will thus be overwritten)
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param markers a character vector providing the marker names to use. By default, all markers are used
#' @param N_Metaclusters a numeric value providing the number of metaclusters to regroup the clusters into
#' @param N_NewClusters a numeric value providing the number of new clusters to regroup the previous clusters into
#' @param criterion a string value providing the agglomeration method to be use.
#' Possible values are: 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median' or 'centroid'
#' Please refer to the function 'hclust' of the 'stats' package for details
#' @param typeMSI a string value providing the type of operation to compute on the marker expression of each cluster ("median" or "mean"; "median" by default)
#' @param rescale.bounds a numeric vector providing the two quantiles (within [0,1]) that will be used as bounds for rescaling the data. Defaults to c(0.05, 0.95)
#' @param seed a numeric value providing the random seed to use during stochastic operations
#' @param checkOverwrite a boolean value indicating whether to check if a Metaclustering has already been performed on the CYTdata object. If TRUE, it will check and ask for confirmation, if FALSE, any previous clustering will be overwritten by default. 

#'
#' @return a S4 object of class 'CYTdata'
#'
#' @export
#'

run_Hierarchical_Metaclustering <- function(CYTdata,
                                            markers = NULL,
                                            N_Metaclusters,
                                            N_NewClusters = NULL,
                                            criterion = c("ward.D", "ward.D2", "single", "complete",
                                                          "average", "mcquitty", "median", "centroid"),
                                            typeMSI = c("median", "mean"),
                                            rescaleBounds = c(0.05, 0.95),
                                            #rescaleMetadata = NULL, #removed as it was removed from rescale.exprs
                                            seed = 42,
                                            checkOverwrite = TRUE) {
  
  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }
  
  checkmate::qassert(N_Metaclusters, "N1")
  checkmate::qassert(N_NewClusters, c("0", "N1"))
  criterion = match.arg(criterion)
  checkmate::qassert(criterion, "S1")
  checkmate::qassert(seed, "N1")
  checkmate::qassert(checkOverwrite, "B1")
  
  if (length(CYTdata@Clustering@clusters)==0){
    stop("Error : Clustering slot is empty. Please compute clustering step before running metaclustering step.") }
  
  if (N_Metaclusters<1) { stop("Error : argument 'N_Metaclusters' must be positive integer") }
  if (length(levels(CYTdata@Clustering@clusters)) < N_Metaclusters) {
    stop("Error : 'N_Metaclusters' argument must be smaller than the number of clusters to metacluster.")
  }
  if (checkOverwrite && length(CYTdata@Metaclustering@metaclusters)!=0){
    reply <- readline(prompt="Metaclustering already performed, do you still want to continue and overwrite (yes or no): ")
    while (!tolower(reply) %in% c("yes", "y", "no", "n")) { reply <- readline(prompt="Reply not standard, please answer yes or no: ") }
    if (tolower(reply) %in% c("no", "n")){ stop("Function stopped, CYTdata object unchanged") }
  }
  if (!is.null(N_NewClusters)) {
    if (N_Metaclusters > N_NewClusters) {
      stop("Error : argument 'N_NewClusters' must be smaller than 'N_Metaclusters'")
    }
    if (N_NewClusters < 1) {
      stop("Error : argument 'N_NewClusters' must be positive integer")
    }
    if (checkOverwrite && length(CYTdata@Clustering@clusters)!=0){
      reply <- readline(prompt="Two levels metaclustering is set and will change clusters identifiers, but clustering already performed.
                        Do you still want to continue and overwrite (yes or no): ")
      while (!tolower(reply) %in% c("yes", "y", "no", "n")) { reply <- readline(prompt="Reply not standard, please answer yes or no: ") }
      if (tolower(reply) %in% c("no", "n")){ stop("Function stopped, CYTdata object unchanged") }
    }
  }
  
  population = checkorderPopulation(CYTdata, population = NULL, level = "clusters")
  markers = checkorderMarkers(CYTdata, markers = markers, order = FALSE, checkDuplicates = TRUE)
  popId = CYTdata@Clustering@clusters
  
  data = computeMSI(CYTdata,
                    scale=TRUE,
                    markers = markers,
                    population = NULL,
                    level = "clusters",
                    typeMSI = typeMSI,
                    rescaleBounds = rescaleBounds)
  
  parameters = list("type" = "hierarchical_clustering",
                    "criterion" = criterion,
                    "seed" = seed,
                    "typeMSI" = typeMSI,
                    "rescaleBounds" = rescaleBounds)
  
  ### Hclust for clusters : clusters
  dataDist = stats::dist(data)
  dataDist[dataDist==0] = 0.000001
  hclustCols = stats::hclust(dataDist, method = criterion)
  isoMDSord = MASS::isoMDS(dataDist, k = 1)$points
  hclustCols = dendextend::rotate(hclustCols, rownames(isoMDSord)[order(isoMDSord)])
  labelCols = hclustCols$labels[hclustCols$order]
  
  # ### Hclust for clusters : markers
  # dataDist = stats::dist(t(data))
  # dataDist[dataDist==0] = 0.000001
  # hclustMarkers = stats::hclust(dataDist, method = criterion)
  # labelMarkers = hclustMarkers$labels[hclustMarkers$order]
  
  ### Metaclustering
  
  cat("\n\nGetting metaclusters.. \n")
  cutreeMetaclusters = stats::cutree(hclustCols, k = N_Metaclusters)
  newCutreeMetaclusters = match(cutreeMetaclusters, unique(cutreeMetaclusters[labelCols]))
  names(newCutreeMetaclusters) = names(cutreeMetaclusters)
  metaclustersBelong = newCutreeMetaclusters[rownames(data)]
  metaclusters = factor(as.vector(metaclustersBelong[popId]))
  
  # cat("Creating HierarchicalHeatmap optionnal list.. \n")
  # HierarchicalHeatmap <- list("data" = data,
  #                             "hierarchyCols" = hclustCols,
  #                             "labelCols" = labelCols,
  #                             "hierarchyMarkers" = hclustMarkers,
  #                             "labelMarkers" = labelMarkers)
  
  if(!is.null(N_NewClusters)){ ### 2nd level
    
    cat("\n\nGetting New Clusters.. \n")
    cutreeNewClusters = stats::cutree(hclustCols, k = N_NewClusters)
    newcutreeNewClusters = match(cutreeNewClusters, unique(cutreeNewClusters[labelCols]))
    names(newcutreeNewClusters) = names(cutreeNewClusters)
    NewClustersBelong = newcutreeNewClusters[rownames(data)]
    newclusters = factor(as.vector(NewClustersBelong[popId]))
    
    palette = rep("white",length(levels(newclusters)))
    names(palette) = levels(newclusters)
    
    # Put into "Clustering" object
    cat("\n - Computing cell cluster count matrix...")
    cellcount = compute.cellcount(newclusters, CYTdata@samples)
    cat("\n - Computing cell cluster abundance matrix...")
    abundance = compute.abundance(cellcount)
    
    cat("Creating clustering object.. \n")
    Clustering.object <- methods::new("Clustering",
                                      "clusters" = newclusters,
                                      "cellcount" = cellcount,
                                      "abundance" = abundance,
                                      "palette" = palette,
                                      "optional_parameters" = parameters)
    
    cat("Updating CYTdata object.. \n")
    CYTdata@Clustering = Clustering.object
    validObject(CYTdata)
    
    #HierarchicalHeatmap$subclusters = popId
  }
  
  palette = rainbow(N_Metaclusters)
  names(palette) = levels(metaclusters)
  
  # Put into "Metaclustering" object
  cat("\n - Computing cell metacluster count matrix...")
  cellcount = compute.cellcount(metaclusters, CYTdata@samples)
  cat("\n - Computing cell metacluster abundance matrix...")
  abundance = compute.abundance(cellcount)
  
  cat("Creating Metaclustering object.. \n")
  Metaclustering.object <- methods::new("Metaclustering",
                                        "metaclusters" = metaclusters,
                                        "cellcount" = cellcount,
                                        "abundance" = abundance,
                                        "palette" = palette,
                                        "optional_parameters" = parameters)#,
                                        #"optional_HierarchicalHeatmap" = HierarchicalHeatmap)
  
  cat("Updating CYTdata object.. \n")
  CYTdata@Metaclustering = Metaclustering.object
  validObject(CYTdata)
  return(CYTdata)
}



