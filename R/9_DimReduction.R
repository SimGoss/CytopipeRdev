#' @title Import an object of class 'DimReduction' to put into a CYTdata object
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param importedDimReductionObject a S4 object of class 'DimReduction' to put in CYTdata object
#' @param checkOverwrite a boolean, if TRUE (default), will pause and display a warning if a DimReduction is already present in the CYTdata and must be overwritten. 
#' 
#' @return a S4 object of class 'CYTdata'
#'
#' @export
#'

importDimReduction <- function(CYTdata, importedDimReductionObject, checkOverwrite=TRUE){
  
  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }
  
  if (class(importedDimReductionObject)!="DimReduction") {
    stop("Error : argument 'importedDimReductionObject' must be a S4 object of class 'DimReduction'.")
  }
  
  if (checkOverwrite && length(CYTdata@DimReduction@coordinates)!=0){
    reply <- readline(prompt="Dimensionnality reduction already performed, do you still want to continue and overwrite (yes or no): ")
    while (!tolower(reply) %in% c("yes", "y", "no", "n")) { reply <- readline(prompt="Reply not standard, please answer yes or no: ") }
    if (tolower(reply) %in% c("no", "n")){ stop("Function stopped, CYTdata object unchanged") }
  }
  
  message("\nUpdating CYTdata object...")
  CYTdata@DimReduction = importedDimReductionObject
  CYTdata = MakeValid(CYTdata, verbose = TRUE)
  return(CYTdata)
}

#' @title Export a S4 object of class 'DimReduction' from a CYTdata object
#' 
#' @description Strictly equivalent to "CYTdata@DimReduction"
#' 
#' @param CYTdata a S4 object of class 'CYTdata'
#'
#' @return a S4 object of class 'DimReduction'
#'
#' @export
#'

export.DimReduction <- function(CYTdata){
  validObject(CYTdata)
  return(CYTdata@DimReduction)
}

#' @title Empty the 'DimReduction' slot in a CYTdata object
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#'
#' @return a S4 object of class 'CYTdata'
#'
#' @export
#'

removeDimReduction <- function(CYTdata) {
  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }
  
  message("\nUpdating CYTdata object...")
  CYTdata@DimReduction = methods::new("DimReduction")
  validObject(CYTdata)
  return(CYTdata)
}

#' @title Dimensionality Reduction of data using dimension reduction techniques
#'
#' @description This function aims to generate coordinates of data in an 2D-reduced space, stored in a CYTdata object.
#'
#' The algorithm used available are :
#' - Principal component analysis (PCA), a linear dimension reduction technique.
#' - Uniform Manifold Approximation and Projection (UMAP), a non-linear dimension reduction technique
#' - t-distributed stochastic neighbor embedding (t-SNE), a non-linear dimension reduction technique
#' - LargeVis-like method (lvish), a non-linear dimension reduction technique
#' The whole set of cell markers or specific cell markers can be used during the dimensionality reduction process.
#' Note that for PCA, currently only the first two principal components are kept, with no way to check the percentage of variance explained by them. You may want to run the PCA yourself on CYTdata@matrix.expression[,markers] using prcomp while we integrate it better. 
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param markers a character vector providing the cell markers to use to generate the reduced data
#' @param type a character value containing the type of RD method to use. Possible values are: "UMAP", "tSNE", "lvish", "PCA" (default = UMAP)
#' @param seed a numeric value providing the random seed to use during stochastic operations
#' @param checkOverwrite a boolean, if TRUE (default), will pause and display a warning if a DimReduction is already present in the CYTdata and must be overwritten. 
#' @param ... additional arguments passed on to method from R package.
#' For PCA, please refer to prcomp method from stats package : https://rdrr.io/r/stats/stats-package.html
#' For UMAP, please refer to umap method from uwot package : https://cran.r-project.org/web/packages/uwot/uwot.pdf
#' For t-SNE, please refer to Rtsne method from Rtsne package : https://cran.r-project.org/web/packages/Rtsne/Rtsne.pdf
#' For lvish, please refer to lvish method from uwot package : https://cran.r-project.org/web/packages/uwot/uwot.pdf
#'
#' @return a S4 object of class 'CYTdata', to which has been added a S4 object of class "DimReduction". 
#' DimReduction contains "coordinates", a n*2 matrix of the n cells coordinates on the 2D reduced space, 
#' and "optional_parameters", a list of "type", the name of the dimension reduction method used, "markers", the list of markers used, and "seed", the seed for stochastic operations. 
#'
#' @export
#'

runDimReduction <- function(CYTdata,
                            markers = NULL,
                            type = c("UMAP", "tSNE", "lvish", "PCA"),
                            seed = 42,
                            checkOverwrite = TRUE,
                            ...){
  
  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }
  
  type = match.arg(type)
  checkmate::qassert(type, "S1")
  checkmate::qassert(markers, c(0,"S*"))
  checkmate::qassert(seed, "N1")
  list.parameters = list(...)
  
  if (checkOverwrite && length(CYTdata@DimReduction@coordinates)!=0){
    reply <- readline(prompt="Dimensionnality reduction already performed, do you still want to continue and overwrite (yes or no): ")
    while (!tolower(reply) %in% c("yes", "y", "no", "n")) { reply <- readline(prompt="Reply not standard, please answer yes or no: ") }
    if (tolower(reply) %in% c("no", "n")){ stop("Function stopped, CYTdata object unchanged") }
  }
  
  markers = checkorderMarkers(CYTdata, markers = markers, order = TRUE)
  data = CYTdata@matrix.expression[,markers]
  
  cat("\nGenerate manifold using", type, " :")
  cat("\n\n - Markers used to generate 2D-reduced data space : ", paste0(markers, collapse = ", "), "\n\n")

  set.seed(seed)
  
  switch(type,
         UMAP = {
           manifold <- uwot::umap(data, ...)
         },
         tSNE = {
           manifold <- Rtsne::Rtsne(data, ...)$Y
         },
         lvish = {
           manifold <- uwot::lvish(data, ...)
         },
         PCA = {
           PC <- stats::prcomp(data, ...) # Compute PCA
           manifold = PC$x[,1:2]
         })
  
  coordinates = data.frame(manifold)
  colnames(coordinates) = c("dim1", "dim2")
  
  # Put into "DimReduction" object
  cat("\n\nCreating DimReduction object and updating CYTdata object.. \n")
  DimReduction.object <- methods::new("DimReduction",
                                      coordinates = coordinates,
                                      optional_parameters = append(list.parameters,
                                                                   list("type" = type,
                                                                        "markers" = markers,
                                                                        "seed" = seed)))
  CYTdata@DimReduction = DimReduction.object
  validObject(CYTdata)
  return(CYTdata)
}










