#' @title CYTdata class definition
#'
#' @description The CYTdata object is a S4 object containing all cytometry expressions.
#'
#' @slot samples a character vector containing the names of the biological samples
#' @slot matrix.expression a data.frame containing the marker expressions of each cells (data to analyse)
#' @slot metadata a data.frame containing the metadata associated to each sample
#' @slot raw.matrix.expression a data.frame containing the marker expressions of cells from the whole dataset before transformation process (raw data)
#'
#' @slot Clustering a S4 object of class 'Clustering' containing the clustering step-related informations
#' @slot Metaclustering a S4 object of class 'Metaclustering' containing the metaclustering step-related informations
#' @slot DimReduction a S4 object of class 'DimReduction' containing the dimensionality reduction step-related informations
#'
#' @slot DiffAbundanceAnalysis a data.frame containing the statistics of comparison computed with 'run.DiffAbundanceAnalysis' function
#' @slot Kinetic a S4 object of class 'Kinetic'
#'
#'
#' @name CYTdata-class
#' @rdname CYTdata-class
#' @exportClass CYTdata
#'

CYTdata <- methods::setClass("CYTdata",
                             slots = c(samples = "factor",
                                       matrix.expression = "data.frame",
                                       metadata = "data.frame",
                                       raw.matrix.expression = "data.frame",
                                       
                                       Clustering = "Clustering",
                                       Metaclustering = "Metaclustering",
                                       DimReduction = "DimReduction",
                                       
                                       DiffAbundanceAnalysis = "data.frame",
                                       Kinetic = "Kinetic"),
                             
                             validity = function(object) {
                               
                               ### Validity of samples, matrix.expression, raw.matrix.expression
                               
                               ## Check samples vector
                               errLev = setdiff(levels(object@samples), unique(object@samples))
                               if (length(errLev)>0) {
                                 stop("Error in CYTdata object : : the samples factor vector contains
                                       levels not present in the vector (", paste0(errLev, collapse = ", "), ").
                                       Please drop absent levels using droplevels function.")
                               }
                               if (!is.factor(object@samples)) {
                                 stop("Error in CYTdata object : the samples vector is not a factor vector anymore.")
                               }
                               if(any(is.na(object@samples))){
                                 stop("Error in  CYTdata object: samples vector (",
                                      paste0(levels(object@samples), collapse = ","),
                                      ") contain NA values.")
                               }
                               if(length(object@samples)!=nrow(object@matrix.expression)){
                                 stop("Error in CYTdata object: The number of rows of the expression matrix: ",
                                      nrow(object@matrix.expression),
                                      ", is inconsistent with the length of samples vector: ",
                                      length(object@samples))
                               }
                               
                               
                               if(any(is.na(colnames(object@matrix.expression)))){
                                 stop("Error in CYTdata object: The column names of the expression matrix (",
                                      paste0(colnames(object@matrix.expression), collapse = ","),
                                      ") contain NA values.")
                               }
                               if(any(is.na(object@matrix.expression))){
                                 stop("Error in CYTdata object: The expression matrix contain NA values.")
                               }
                               doubleMarkers = duplicated(colnames(object@matrix.expression))
                               if(any(doubleMarkers)){
                                 stop("Error in CYTdata object: The column names of the raw expression matrix (",
                                      paste0(colnames(object@matrix.expression), collapse = ","),
                                      ") contain duplicated values (",
                                      paste0(colnames(object@matrix.expression)[doubleMarkers], collapse = ","),").")
                               }
                               
                               if (ncol(object@raw.matrix.expression)>0 &&
                                   !identical(colnames(object@matrix.expression), colnames(object@raw.matrix.expression))) {
                                 stop("Error in CYTdata object: The column names of the expression matrix (",
                                      paste0(colnames(object@matrix.expression), collapse = ","),
                                      ") are not exactly equal (same order) to the column names of the raw expression matrix (",
                                      paste0(colnames(object@raw.matrix.expression), collapse = ","),
                                      ")")
                               }
                               
                               
                               ### Validity of metadata
                               # dataframe with factor column, unique colnames, including 'Timepoint' and 'Individual'
                               # samples vector levels as rownames
                               
                               if (nrow(object@metadata) != 0){
                                 classes = sapply(object@metadata, class)
                                 names(classes) = colnames(object@metadata)
                                 if(any(classes!="factor")) {
                                   stop("Error in CYTdata object: Metadata columns must be factors (Current classes :",
                                        paste0(paste(names(classes), classes, sep = "_"),collapse = ", "), ")")
                                 }
                                 if (!identical(levels(object@samples), rownames(object@metadata))){ # levels and rownames are unique
                                   stop("Error in CYTdata object: Sample names (rownames) of the metadata (",
                                        paste0(rownames(object@metadata), collapse = ", "),
                                        ") are not identical (order included) with the sample levels (",
                                        paste0(levels(object@samples), collapse = ", "),
                                        ")")
                                 }
                                 if (length(unique(colnames(object@metadata)))!=ncol(object@metadata)) {
                                   stop("Error in CYTdata object: Colnames of the metadata are not unique.
                                         They must be unique because they represent different biological conditions")
                                 }
                                 if (!all(c("Timepoint", "Individual") %in% colnames(object@metadata))) {
                                   stop("Error in CYTdata object: Colnames of the metadata (",
                                        paste0(colnames(object@metadata), collapse = ",") ,
                                        ") must include 'Timepoint' and 'Individual' variable.")
                                 }
                               }
                               
                               ######## Validity of Clustering ##########
                               
                               errLev = setdiff(levels(object@Clustering@clusters), unique(object@Clustering@clusters))
                               if (length(errLev)>0) {
                                 stop("Error in CYTdata@Clustering object : : the clusters identifiers factor vector contains
                                       levels not present in the vector (", paste0(errLev, collapse = ", "), ").
                                       Please drop absent levels using droplevels function.")
                               }
                               
                               if (length(object@Clustering@clusters)!=0){
                                 
                                 if (!is.factor(object@Clustering@clusters)) {
                                   stop("Error in CYTdata@Clustering object : the clusters identifiers vector is not a factor vector anymore.")
                                 }
                                 
                                 if(any(is.na(object@Clustering@clusters))){
                                   stop("Error in  CYTdata@Clustering object: clusters identifiers vector (",
                                        paste0(levels(object@Clustering@clusters), collapse = ","),
                                        ") contain NA values.")
                                 }
                                 if (length(object@Clustering@clusters) != length(object@samples)) {
                                   stop("Error in CYTdata@Clustering object: The number of events in the samples vector,
                                         expression matrix (", length(object@samples), "), is inconsistent with the number
                                         of events in the clusters vector (", length(object@Clustering@clusters), ")")
                                 }
                                 
                                 ## Check clusters vector and palette vector
                                 if (!identical(levels(object@Clustering@clusters), names(object@Clustering@palette))){
                                   stop("Error in CYTdata@Clustering object: Clusters names (",
                                        paste0(levels(object@Clustering@clusters), collapse = ","),
                                        ") are not identical (order included) to the names stored of palette vector (",
                                        paste0(names(object@Clustering@palette), collapse = ","), ")")
                                 }
                                 if(!all(areColors(object@Clustering@palette))){
                                   stop("Error in CYTdata@Clustering object: Named palette vector (",
                                        paste0(object@Clustering@palette, collapse = ","),
                                        "), does not contain only hexadecimal color.)")
                                 }
                                 
                                 ##### Check clusters cellcount matrix
                                 if (ncol(object@Clustering@cellcount)==0 && nrow(object@Clustering@cellcount)==0){
                                   stop("Error in CYTdata@Clustering object : clusters identifiers given but no cellcount matrix computed.")
                                 }
                                 else {
                                   if (!identical(levels(object@Clustering@clusters), rownames(object@Clustering@cellcount))) {
                                     stop("Error in CYTdata@Clustering object : The rownames of cellcount matrix
                                        are not identical (order included) to the levels of clusters identifiers.")
                                   }
                                   if (!identical(levels(object@samples), colnames(object@Clustering@cellcount))) {
                                     stop("Error in CYTdata@Clustering object : The colnames of cellcount matrix
                                        are not identical (order included) to the levels of samples identifiers.")
                                   }
                                 }
                                 ##### Check clusters abundance matrix
                                 if (ncol(object@Clustering@abundance)==0 && nrow(object@Clustering@abundance)==0){
                                   stop("Error in CYTdata@Clustering object : clusters identifiers given but no abundance matrix computed.")
                                 }
                                 else {
                                   if (!identical(levels(object@Clustering@clusters), rownames(object@Clustering@abundance))) {
                                     stop("Error in CYTdata@Clustering object : The rownames of abundance matrix
                                        are not identical (order included) to the levels of clusters identifiers.")
                                   }
                                   if (!identical(levels(object@samples), colnames(object@Clustering@abundance))) {
                                     stop("Error in CYTdata@Clustering object : The rcolnames of abundance matrix
                                        are not identical (order included) to the levels of samples identifiers.")
                                   }
                                 }
                                 
                               }
                               else { ##### If clusters slot is empty, nothing should exist into Clustering object
                                 if (ncol(object@Clustering@cellcount)!=0 && nrow(object@Clustering@cellcount)!=0){
                                   stop("Error in CYTdata@Clustering object: clusters cellcount matrix is given but no
                                         clusters indentifiers given. Please remove it.")
                                 }
                                 if (ncol(object@Clustering@abundance)!=0 && nrow(object@Clustering@abundance)!=0){
                                   stop("Error in CYTdata@Clustering object: clusters abundance matrix is given but no
                                         clusters indentifiers given. Please remove it.")
                                 }
                                 if (length(object@Clustering@palette)!=0){
                                   stop("Error in CYTdata@Clustering object: palette vector is given but no
                                         clusters indentifiers given. Please remove it.")
                                 }
                                 if (length(object@Clustering@optional_parameters)!=0){
                                   stop("Error in CYTdata@Clustering object: optional_parameters list is given but no
                                          clusters indentifiers given. Please remove it.")
                                 }
                               }
                               
                               
                               ######## Validity of Metaclustering ##########
                               
                               errLev = setdiff(levels(object@Metaclustering@metaclusters), unique(object@Metaclustering@metaclusters))
                               if (length(errLev)>0) {
                                 stop("Error in CYTdata@Metaclustering object : the metaclusters identifiers factor vector contains
                                       levels not present in the vector (", paste0(errLev, collapse = ", "), ").
                                       Please drop absent levels using droplevels function.")
                               }
                               
                               if (length(object@Metaclustering@metaclusters)!=0){
                                 
                                 if (!is.factor(object@Metaclustering@metaclusters)) {
                                   stop("Error in CYTdata@Metaclustering object : the metaclusters identifiers vector is not a factor vector anymore.")
                                 }
                                 
                                 if(any(is.na(object@Metaclustering@metaclusters))){
                                   stop("Error in  CYTdata@Metaclustering object : metaclusters identifiers vector (",
                                        paste0(levels(object@Metaclustering@metaclusters), collapse = ","),
                                        ") contain NA values.")
                                 }
                                 if (length(object@Metaclustering@metaclusters) != length(object@samples)) {
                                   stop("Error in CYTdata@Metaclustering object : The number of events in the samples vector,
                                         expression matrix (", length(object@samples), "), is inconsistent with the number
                                         of events in the metaclusters vector (", length(object@Metaclustering@metaclusters), ")")
                                 }
                                 
                                 if (length(object@Clustering@clusters)==0) {
                                   stop("Error in CYTdata@Metaclustering object : metaclusters identifiers vector is given but no clusters identifiers
                                        vector given (Clustering@clusters slot is empty).")
                                 }
                                 
                                 
                                 ## Check metaclusters vector and palette vector
                                 if (!identical(levels(object@Metaclustering@metaclusters), names(object@Metaclustering@palette))){
                                   stop("Error in CYTdata@Metaclustering object: metaclusters names (",
                                        paste0(levels(object@Metaclustering@metaclusters), collapse = ","),
                                        ") are not identical (order included) to the names stored of palette vector (",
                                        paste0(names(object@Metaclustering@palette), collapse = ","), ")")
                                 }
                                 if(!all(areColors(object@Metaclustering@palette))){
                                   stop("Error in CYTdata@Metaclustering object: Named palette vector (",
                                        paste0(object@Metaclustering@palette, collapse = ","),
                                        "), does not contain only hexadecimal color.)")
                                 }
                                 
                                 ##### Check metaclusters cellcount matrix
                                 if (ncol(object@Metaclustering@cellcount)==0 && nrow(object@Metaclustering@cellcount)==0){
                                   stop("Error in CYTdata@Metaclustering object : metaclusters identifiers given but no cellcount matrix computed.")
                                 }
                                 else {
                                   if (!identical(levels(object@Metaclustering@metaclusters), rownames(object@Metaclustering@cellcount))) {
                                     stop("Error in CYTdata@Metaclustering object : The rownames of cellcount matrix
                                        are not identical (order included) to the levels of metaclusters identifiers.")
                                   }
                                   if (!identical(levels(object@samples), colnames(object@Metaclustering@cellcount))) {
                                     stop("Error in CYTdata@Metaclustering object : The colnames of cellcount matrix
                                        are not identical (order included) to the levels of samples identifiers.")
                                   }
                                 }
                                 ##### Check metaclusters abundance matrix
                                 if (ncol(object@Metaclustering@abundance)==0 && nrow(object@Metaclustering@abundance)==0){
                                   stop("Error in CYTdata@Metaclustering object : metaclusters identifiers given but no abundance matrix computed.")
                                 }
                                 else {
                                   if (!identical(levels(object@Metaclustering@metaclusters), rownames(object@Metaclustering@abundance))) {
                                     stop("Error in CYTdata@Metaclustering object : The rownames of abundance matrix
                                        are not identical (order included) to the levels of metaclusters identifiers.")
                                   }
                                   if (!identical(levels(object@samples), colnames(object@Metaclustering@abundance))) {
                                     stop("Error in CYTdata@Metaclustering object : The rcolnames of abundance matrix
                                        are not identical (order included) to the levels of samples identifiers.")
                                   }
                                 }
                                 
                               }
                               else { ##### If metaclusters slot is empty, nothing should exist into Metaclustering object
                                 if (ncol(object@Metaclustering@cellcount)!=0 && nrow(object@Metaclustering@cellcount)!=0){
                                   stop("Error in CYTdata@Metaclustering object: metaclusters cellcount matrix is given but no
                                         metaclusters indentifiers given. Please remove it.")
                                 }
                                 if (ncol(object@Metaclustering@abundance)!=0 && nrow(object@Metaclustering@abundance)!=0){
                                   stop("Error in CYTdata@Metaclustering object: metaclusters abundance matrix is given but no
                                         metaclusters indentifiers given. Please remove it.")
                                 }
                                 if (length(object@Metaclustering@palette)!=0){
                                   stop("Error in CYTdata@Metaclustering object: palette vector is given but no
                                         metaclusters indentifiers given. Please remove it.")
                                 }
                                 if (length(object@Metaclustering@optional_parameters)!=0){
                                   stop("Error in CYTdata@Metaclustering object: optional_parameters list is given but no
                                          metaclusters indentifiers given. Please remove it.")
                                 }
                               }
                               
                               
                               
                               ######## Validity of DimReduction ##########
                               
                               if (nrow(object@DimReduction@coordinates)!=0 && ncol(object@DimReduction@coordinates)!=0) {
                                 if(any(is.na(object@DimReduction@coordinates))){
                                   stop("Error in  CYTdata@DimReduction object: DimReduction coordinates dataframe contain NA values.")
                                 }
                                 if (nrow(object@DimReduction@coordinates) != length(object@samples)) {
                                   stop("Error in CYTdata@DimReduction object: The number of events in the samples vector/expression matrix (",
                                        length(object@samples), ", is inconsistent with the number of coordinates in the embedded space (",
                                        nrow(object@DimReduction@coordinates), ").")
                                 }
                                 if (any(is.na(colnames(object@DimReduction@coordinates)))) {
                                   stop("Error in  CYTdata@DimReduction object: DimReduction coordinates colnames contain NA values.")
                                 }
                                 if (ncol(object@DimReduction@coordinates)!=2) {
                                   stop("Error in  CYTdata@DimReduction object: DimReduction coordinates must be 2D dataframe ( ",
                                        ncol(object@DimReduction@coordinates), " dimensions here).")
                                 }
                               }
                               
                               # If no coordinates, nothing should exist
                               if (nrow(object@DimReduction@coordinates)==0 &&
                                   ncol(object@DimReduction@coordinates)==0 &&
                                   length(object@DimReduction@optional_parameters)!=0){
                                 stop("Error in CYTdata@DimReduction object: optional_parameters list is given but no
                                         DimReduction coordinates given. Please remove it.")
                                 
                               }
                               
                               return(TRUE)
                             }
)

#' @title Checks validity of all elements of a CYTdata object
#'
#' @description Checks validity of all elements of a CYTdata object (préciser les étapes ?)
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param verbose whether to display additional messages
#'
#' @return a S4 object of class 'CYTdata'
#'
#'@export
#'

MakeValid <- function(CYTdata, verbose = TRUE){
  
  checkmate::qassert(verbose, "B1")
  
  ### Validity of samples, matrix.expression, raw.matrix.expression
  
  ## Check samples vector
  errLev = setdiff(levels(CYTdata@samples), unique(CYTdata@samples))
  if (length(errLev)>0) {
    if (verbose) message("Warning in CYTdata object : the samples identifiers factor vector contains
      levels not present in the vector (", paste0(errLev, collapse = ", "), "). Dropping these absent levels..")
    CYTdata@samples = droplevels(CYTdata@samples)
  }
  if (!is.factor(CYTdata@samples)) {
    stop("Error in CYTdata object : the samples vector is not a factor vector anymore.")
  }
  if(any(is.na(CYTdata@samples))){
    stop("Error in  CYTdata object: samples vector (",
         paste0(levels(CYTdata@samples), collapse = ","),
         ") contain NA values.")
  }
  if(length(CYTdata@samples)!=nrow(CYTdata@matrix.expression)){
    stop("Error in CYTdata object: The number of rows of the expression matrix: ",
         nrow(CYTdata@matrix.expression),
         ", is inconsistent with the length of samples vector: ",
         length(CYTdata@samples))
  }
  
  
  if(any(is.na(colnames(CYTdata@matrix.expression)))){
    stop("Error in CYTdata object: The column names of the expression matrix (",
         paste0(colnames(CYTdata@matrix.expression), collapse = ","),
         ") contain NA values.")
  }
  if(any(is.na(CYTdata@matrix.expression))){
    stop("Error in CYTdata object: The expression matrix contain NA values.")
  }
  doubleMarkers = duplicated(colnames(CYTdata@matrix.expression))
  if(any(doubleMarkers)){
    stop("Error in CYTdata object: The column names of the raw expression matrix (",
         paste0(colnames(CYTdata@matrix.expression), collapse = ","),
         ") contain duplicated values (",
         paste0(colnames(CYTdata@matrix.expression)[doubleMarkers], collapse = ","),").")
  }
  
  if (ncol(CYTdata@raw.matrix.expression)>0 &&
      !identical(colnames(CYTdata@matrix.expression), colnames(CYTdata@raw.matrix.expression))) {
    stop("Error in CYTdata object: The column names of the expression matrix (",
         paste0(colnames(CYTdata@matrix.expression), collapse = ","),
         ") are not exactly equal (same order) to the column names of the raw expression matrix (",
         paste0(colnames(CYTdata@raw.matrix.expression), collapse = ","),
         ")")
  }
  
  
  ######## Validity of metadata ##########
  # dataframe with factor column, unique colnames, including 'Timepoint' and 'Individual'
  # samples vector levels as rownames
  
  if (nrow(CYTdata@metadata) != 0){
    classes = sapply(CYTdata@metadata, class)
    names(classes) = colnames(CYTdata@metadata)
    if(any(classes!="factor")) {
      if (verbose) message("Warning in CYTdata object: Metadata columns must be factors. But some columns are not factor yet ( ",
                           paste0(paste(names(classes)[classes!="factor"], classes, sep = "_"),collapse = ", "), ").
           Converting remaining metadata columns to factor..")
      notfactor <- function(x) return(!is.factor(x))
      factorByUnique <- function(x) return(factor(x, levels=unique(x)))
      CYTdata@metadata = dplyr::mutate_if(CYTdata@metadata, notfactor, factorByUnique)
    }
    
    if (setequal(levels(CYTdata@samples), rownames(CYTdata@metadata))) {
      if (!identical(levels(CYTdata@samples), rownames(CYTdata@metadata))) {
        if (verbose) message("Warning in CYTdata object: Sample names (rownames) of the metadata (",
                             paste0(rownames(CYTdata@metadata), collapse = ", "),
                             ") are identical with the sample levels (",
                             paste0(levels(CYTdata@samples), collapse = ", "),
                             ") but not in the same order. Ordering metadata rownames the same order than sample levels..")
        CYTdata@metadata = CYTdata@metadata[levels(CYTdata@samples),]
      }
    }
    else {
      stop("Error: Sample names (rownames) of the metadata (",
           paste0(rownames(CYTdata@metadata), collapse = ","),
           ") are not concordant with the sample levels (",
           paste0(levels(CYTdata@samples), collapse = ","),
           ").")
    }
    
    if (length(unique(colnames(CYTdata@metadata)))!=ncol(CYTdata@metadata)) {
      stop("Error in CYTdata object: Colnames of the metadata are not unique.
           They must be unique because they represent different biological conditions")
    }
    if (!all(c("Timepoint", "Individual") %in% colnames(CYTdata@metadata))) {
      stop("Error in CYTdata object: Colnames of the metadata (",
           paste0(colnames(CYTdata@metadata), collapse = ",") ,
           ") must include 'Timepoint' and 'Individual' variable.")
    }
  }
  
  ######## Validity of Clustering ##########
  
  errLev = setdiff(levels(CYTdata@Clustering@clusters), unique(CYTdata@Clustering@clusters))
  if (length(errLev)>0) {
    if (verbose) message("Warning in CYTdata@Clustering object : the clusters identifiers factor vector contains
      levels not present in the vector (", paste0(errLev, collapse = ", "), "). Dropping these absent levels..")
    CYTdata@Clustering@clusters = droplevels(CYTdata@Clustering@clusters)
  }
  
  if (length(CYTdata@Clustering@clusters)!=0){
    
    if (!is.factor(CYTdata@Clustering@clusters)) {
      stop("Error in CYTdata@Clustering object : the clusters identifiers vector is not a factor vector anymore.")
    }
    
    if(any(is.na(CYTdata@Clustering@clusters))){
      stop("Error in  CYTdata@Clustering object : clusters identifiers vector (",
           paste0(levels(CYTdata@Clustering@clusters), collapse = ","),
           ") contain NA values.")
    }
    if (length(CYTdata@Clustering@clusters) != length(CYTdata@samples)) {
      stop("Error in CYTdata@Clustering object : The number of events in the samples vector,
      expression matrix (", length(CYTdata@samples), "), is inconsistent with the number
           of events in the clusters vector (", length(CYTdata@Clustering@clusters), ")")
    }
    
    
    ##### Check clusters vector and palette vector
    if (length(CYTdata@Clustering@palette)==0) {
      if (verbose) message("Warning in CYTdata@Clustering object : clusters identifiers given but no palette given.
            Generating blank palette..")
      namesCs = levels(CYTdata@Clustering@clusters)
      blankPalette = rep("#FFFFFF", length(namesCs))
      names(blankPalette) = namesCs
      CYTdata@Clustering@palette = blankPalette
    }
    if (!identical(levels(CYTdata@Clustering@clusters), names(CYTdata@Clustering@palette))){
      stop("Error in CYTdata@Clustering object : Clusters names (",
           paste0(levels(CYTdata@Clustering@clusters), collapse = ","),
           ") are not identical (order included) to the names stored of palette vector (",
           paste0(names(CYTdata@Clustering@palette), collapse = ","), ")")
    }
    if(!all(areColors(CYTdata@Clustering@palette))){
      stop("Error in CYTdata@Clustering object : Named palette vector (",
           paste0(CYTdata@Clustering@palette, collapse = ","),
           "), does not contain only hexadecimal color.)")
    }
    
    ##### Check clusters cellcount matrix
    if (ncol(CYTdata@Clustering@cellcount)==0 && nrow(CYTdata@Clustering@cellcount)==0){
      if (verbose) message("Warning in CYTdata@Clustering object : clusters identifiers given but no cellcount matrix computed.
              Computing cellcount matrix (is empty)..")
      CYTdata@Clustering@cellcount = compute.cellcount(CYTdata@Clustering@clusters, CYTdata@samples)
    }
    else {
      if (!identical(levels(CYTdata@Clustering@clusters), rownames(CYTdata@Clustering@cellcount))) {
        stop("Error in CYTdata@Clustering object : The rownames of cellcount matrix are not identical (order included)
        to the levels of clusters identifiers. Computing cellcount matrix ..")
        CYTdata@Clustering@cellcount = compute.cellcount(CYTdata@Clustering@clusters, CYTdata@samples)
      }
      if (!identical(levels(CYTdata@samples), colnames(CYTdata@Clustering@cellcount))) {
        stop("Error in CYTdata@Clustering object : The colnames of cellcount matrix  are not identical (order included)
             to the levels of samples identifiers. Computing cellcount matrix ..")
        CYTdata@Clustering@cellcount = compute.cellcount(CYTdata@Clustering@clusters, CYTdata@samples)
      }
    }
    ##### Check clusters abundance matrix
    if (ncol(CYTdata@Clustering@abundance)==0 && nrow(CYTdata@Clustering@abundance)==0){
      if (verbose) message("Warning in CYTdata@Clustering object : clusters identifiers given but no abundance matrix computed.
              Computing abundance matrix (is empty)..")
      CYTdata@Clustering@abundance = compute.abundance(CYTdata@Clustering@cellcount)
    }
    else {
      if (!identical(levels(CYTdata@Clustering@clusters), rownames(CYTdata@Clustering@abundance))) {
        stop("Error in CYTdata@Clustering object : The rownames of abundance matrix are not identical (order included)
        to the levels of clusters identifiers. Computing abundance matrix ..")
        CYTdata@Clustering@abundance = compute.abundance(CYTdata@Clustering@cellcount)
      }
      if (!identical(levels(CYTdata@samples), colnames(CYTdata@Clustering@abundance))) {
        stop("Error in CYTdata@Clustering object : The colnames of abundance matrix  are not identical (order included)
             to the levels of samples identifiers. Computing abundance matrix ..")
        CYTdata@Clustering@abundance = compute.abundance(CYTdata@Clustering@cellcount)
      }
    }
  }
  else { ##### If clusters slot is empty, nothing should exist into Clustering object
    if (ncol(CYTdata@Clustering@cellcount)!=0 && nrow(CYTdata@Clustering@cellcount)!=0){
      if (verbose) message("Warning in CYTdata@Clustering object: clusters cellcount matrix is given but no clusters indentifiers given.
              Removing cellcount matrix (empty slot).")
      CYTdata@Clustering@cellcount = data.frame()
    }
    if (ncol(CYTdata@Clustering@abundance)!=0 && nrow(CYTdata@Clustering@abundance)!=0){
      if (verbose) message("Warning in CYTdata@Clustering object: clusters abundance matrix is given but no clusters indentifiers given.
              Removing abundance matrix (empty slot).")
      CYTdata@Clustering@abundance = data.frame()
    }
    if (length(CYTdata@Clustering@palette)!=0){
      if (verbose) message("Warning in CYTdata@Clustering object: palette vector is given but no clusters indentifiers given.
              Removing palette vector (empty slot).")
      CYTdata@Clustering@palette = c()
    }
    if (length(CYTdata@Clustering@optional_parameters)!=0){
      if (verbose) message("Warning in CYTdata@Clustering object: optional_parameters list is given but no clusters indentifiers given.
              Removing optional_parameters list (empty slot).")
      CYTdata@Clustering@optional_parameters = list()
    }
  }
  
  
  ######## Validity of Metaclustering ##########
  
  errLev = setdiff(levels(CYTdata@Metaclustering@metaclusters), unique(CYTdata@Metaclustering@metaclusters))
  if (length(errLev)>0) {
    if (verbose) message("Warning in CYTdata@Metaclustering object : the metaclusters identifiers factor vector contains
      levels not present in the vector (", paste0(errLev, collapse = ", "), "). Dropping these absent levels..")
    CYTdata@Metaclustering@metaclusters = droplevels(CYTdata@Metaclustering@metaclusters)
  }
  
  if (length(CYTdata@Metaclustering@metaclusters)!=0){
    
    if (!is.factor(CYTdata@Metaclustering@metaclusters)) {
      stop("Error in CYTdata@Metaclustering object : the metaclusters identifiers vector is not a factor vector anymore.")
    }
    
    if(any(is.na(CYTdata@Metaclustering@metaclusters))){
      stop("Error in  CYTdata@Metaclustering object : metaclusters identifiers vector (",
           paste0(levels(CYTdata@Metaclustering@metaclusters), collapse = ","),
           ") contain NA values.")
    }
    if (length(CYTdata@Metaclustering@metaclusters) != length(CYTdata@samples)) {
      stop("Error in CYTdata@Metaclustering object : The number of events in the samples vector,
      expression matrix (", length(CYTdata@samples), "), is inconsistent with the number
           of events in the metaclusters vector (", length(CYTdata@Metaclustering@metaclusters), ")")
    }
    
    ##### Check metaclusters vector and palette vector
    if (length(CYTdata@Metaclustering@palette)==0) {
      if (verbose) message("Warning in CYTdata@Metaclustering object : metaclusters identifiers given but no palette given.
            Generating blank palette..")
      namesCs = levels(CYTdata@Metaclustering@metaclusters)
      blankPalette = rep("#FFFFFF", length(namesCs))
      names(blankPalette) = namesCs
      CYTdata@Metaclustering@palette = blankPalette
    }
    if (!identical(levels(CYTdata@Metaclustering@metaclusters), names(CYTdata@Metaclustering@palette))){
      stop("Error in CYTdata@Metaclustering object : metaclusters names (",
           paste0(levels(CYTdata@Metaclustering@metaclusters), collapse = ","),
           ") are not identical (order included) to the names stored of palette vector (",
           paste0(names(CYTdata@Metaclustering@palette), collapse = ","), ")")
    }
    if(!all(areColors(CYTdata@Metaclustering@palette))){
      stop("Error in CYTdata@Metaclustering object : Named palette vector (",
           paste0(CYTdata@Metaclustering@palette, collapse = ","),
           "), does not contain only hexadecimal color.)")
    }
    
    ##### Check metaclusters cellcount matrix
    if (ncol(CYTdata@Metaclustering@cellcount)==0 && nrow(CYTdata@Metaclustering@cellcount)==0){
      if (verbose) message("Warning in CYTdata@Metaclustering object : metaclusters identifiers given but no cellcount matrix computed.
              Computing cellcount matrix (is empty)..")
      CYTdata@Metaclustering@cellcount = compute.cellcount(CYTdata@Metaclustering@metaclusters, CYTdata@samples)
    }
    else {
      if (!identical(levels(CYTdata@Metaclustering@metaclusters), rownames(CYTdata@Metaclustering@cellcount))) {
        if (verbose) message("Warning in CYTdata@Metaclustering object : The rownames of cellcount matrix are not identical (order included)
        to the levels of metaclusters identifiers. Computing cellcount matrix ..")
        CYTdata@Metaclustering@cellcount = compute.cellcount(CYTdata@Metaclustering@metaclusters, CYTdata@samples)
      }
      if (!identical(levels(CYTdata@samples), colnames(CYTdata@Metaclustering@cellcount))) {
        if (verbose) message("Warning in CYTdata@Metaclustering object : The colnames of cellcount matrix  are not identical (order included)
             to the levels of samples identifiers. Computing cellcount matrix ..")
        CYTdata@Metaclustering@cellcount = compute.cellcount(CYTdata@Metaclustering@metaclusters, CYTdata@samples)
      }
    }
    ##### Check metaclusters abundance matrix
    if (ncol(CYTdata@Metaclustering@abundance)==0 && nrow(CYTdata@Metaclustering@abundance)==0){
      if (verbose) message("Warning in CYTdata@Metaclustering object : metaclusters identifiers given but no abundance matrix computed.
              Computing abundance matrix (is empty)..")
      CYTdata@Metaclustering@abundance = compute.abundance(CYTdata@Metaclustering@cellcount)
    }
    else {
      if (!identical(levels(CYTdata@Metaclustering@metaclusters), rownames(CYTdata@Metaclustering@abundance))) {
        if (verbose) message("Warning in CYTdata@Metaclustering object : The rownames of abundance matrix are not identical (order included)
        to the levels of metaclusters identifiers. Computing abundance matrix ..")
        CYTdata@Metaclustering@abundance = compute.abundance(CYTdata@Metaclustering@cellcount)
      }
      if (!identical(levels(CYTdata@samples), colnames(CYTdata@Metaclustering@abundance))) {
        if (verbose) message("Warning in CYTdata@Metaclustering object : The colnames of abundance matrix  are not identical (order included)
             to the levels of samples identifiers. Computing abundance matrix ..")
        CYTdata@Metaclustering@abundance = compute.abundance(CYTdata@Metaclustering@cellcount)
      }
    }
    ##### Check clusters identifiers
    if (length(CYTdata@Clustering@clusters)==0) {
      if (verbose) message("Warning in CYTdata@Metaclustering object : metaclusters identifiers vector is given but no clusters identifiers
                           vector given (Clustering@clusters slot is empty). Clusters identifiers vectro will be copied from metaclusters identifiers vector")
      CYTdata@Clustering@clusters = CYTdata@Metaclustering@metaclusters
      CYTdata@Clustering@palette = CYTdata@Metaclustering@palette
      CYTdata@Clustering@cellcount = compute.cellcount(CYTdata@Clustering@clusters, CYTdata@samples)
      CYTdata@Clustering@abundance = compute.abundance(CYTdata@Clustering@cellcount)
    }
  }
  else { ##### If metaclusters slot is empty, nothing should exist into Metaclustering object
    if (ncol(CYTdata@Metaclustering@cellcount)!=0 && nrow(CYTdata@Metaclustering@cellcount)!=0){
      if (verbose) message("Warning in CYTdata@Metaclustering object: metaclusters cellcount matrix is given but no metaclusters indentifiers given.
              Removing cellcount matrix (empty slot).")
      CYTdata@Metaclustering@cellcount = data.frame()
    }
    if (ncol(CYTdata@Metaclustering@abundance)!=0 && nrow(CYTdata@Metaclustering@abundance)!=0){
      if (verbose) message("Warning in CYTdata@Metaclustering object: metaclusters abundance matrix is given but no metaclusters indentifiers given.
              Removing abundance matrix (empty slot).")
      CYTdata@Metaclustering@abundance = data.frame()
    }
    if (length(CYTdata@Metaclustering@palette)!=0){
      if (verbose) message("Warning in CYTdata@Metaclustering object: palette vector is given but no metaclusters indentifiers given.
              Removing palette vector (empty slot).")
      CYTdata@Metaclustering@palette = c()
    }
    if (length(CYTdata@Metaclustering@optional_parameters)!=0){
      if (verbose) message("Warning in CYTdata@Metaclustering object: optional_parameters list is given but no metaclusters indentifiers given.
              Removing optional_parameters list (empty slot).")
      CYTdata@Metaclustering@optional_parameters = list()
    }
  }
  
  ######## Validity of DimReduction ##########
  
  if (nrow(CYTdata@DimReduction@coordinates)!=0 && ncol(CYTdata@DimReduction@coordinates)!=0) {
    if(any(is.na(CYTdata@DimReduction@coordinates))){
      stop("Error in CYTdata@DimReduction object: DimReduction coordinates dataframe contain NA values.")
    }
    if (nrow(CYTdata@DimReduction@coordinates) != length(CYTdata@samples)) {
      stop("Error in CYTdata@DimReduction object: The number of events in the samples vector/expression matrix (",
           length(CYTdata@samples), ", is inconsistent with the number of coordinates in the embedded space (",
           nrow(CYTdata@DimReduction@coordinates), ").")
    }
    if (ncol(CYTdata@DimReduction@coordinates)!=2) {
      stop("Error in CYTdata@DimReduction object: DimReduction coordinates must be 2D dataframe ( ",
           ncol(CYTdata@DimReduction@coordinates), " dimensions here).")
    }
    if (any(is.na(colnames(CYTdata@DimReduction@coordinates)))) {
      if (verbose) message("Warning in CYTdata@DimReduction object: DimReduction coordinates colnames contain NA values.
              Setting it to default name (dim1, dim2)..")
      colnames(CYTdata@DimReduction@coordinates) = c("dim1", "dim2")
    }
  }
  
  # If no coordinates, nothing should exist
  if (nrow(CYTdata@DimReduction@coordinates)==0 &&
      ncol(CYTdata@DimReduction@coordinates)==0 &&
      length(CYTdata@DimReduction@optional_parameters)!=0){
    if (verbose) message("Warning in CYTdata@DimReduction object: optional_parameters list is given but no
            DimReduction coordinates given. Removing optional_parameters list (empty slot).")
    CYTdata@DimReduction@optional_parameters = list()
  }
  
  return(CYTdata)
}




#' @title Clustering class definition
#'
#' @description The Clustering object is a S4 object containing all the data concerning the Clustering step performed on a 'CYTdata' object
#'
#' @slot clusters a vector containing the identifiers of cell clusters as a factor
#' @slot cellcount a data.frame containing the number of cells associated to each cluster for each sample
#' @slot abundance a data.frame containing the percentage of cells associated to each cluster for each sample
#' @slot palette TO DO
#' #' @slot optional_parameters TO CHECK a vector containing the parameters used for the identification of the cell clusters (especially the argument of 'run.Clustering' function : seed, markers, architecture, type, etc.)
#' @slot optional_plots TO DO
#' 
#'
#' @name Clustering-class
#' @rdname Clustering-class
#' @exportClass Clustering

Clustering <- methods::setClass("Clustering",
                                slots = c(clusters = "factor",
                                          cellcount = "data.frame",
                                          abundance = "data.frame",
                                          palette = "vector",
                                          optional_parameters = "list",
                                          optional_plots = "list"))

#' @title Metaclustering class definition
#'
#' @description The Metaclustering object is a S4 object containing all the data concerning the Metaclustering step performed on a 'CYTdata' object
#'
#' @slot metaclusters a vector containing the identifiers of cell metaclusters (phenotypic families) as a factor
#' @slot cellcount a data.frame containing the number of cells associated to each metacluster for each sample
#' @slot abundance a data.frame containing the percentage of cells associated to each metacluster for each sample
#' @slot palette TO DO
#' @slot optional_HierarchicalHeatmap TO CHECK a list containing all the features about heatmap representing phenotypes of clusters and families
#' @slot optional_parameters TO CHECK a vector containing the parameters used for the identification of the cell clusters (especially the argument of 'run.Clustering' function : seed, markers, architecture, type, etc.)
#'
#' @name Metaclustering-class
#' @rdname Metaclustering-class
#' @exportClass Metaclustering
#'

Metaclustering <- methods::setClass("Metaclustering",
                                    slots = c(metaclusters = "factor",
                                              cellcount = "data.frame",
                                              abundance = "data.frame",
                                              palette = "vector",
                                              optional_parameters = "list",
                                              optional_HierarchicalHeatmap = "list"))

#' @title DimReduction class definition
#'
#' @description The DimReduction object is a S4 object containing all the data concerning the DimReduction step performed on a 'CYTdata' object
#'
#' @slot coordinates a data.frame containing the coordinates data points in the reduced data space
#' @slot optional_parameters a list containing the parameters used for DimReduction computation (especially the arguments of 'run.DimReduction' function : seed, markers, type, etc.)
#'
#' @name DimReduction-class
#' @rdname DimReduction-class
#' @exportClass DimReduction
#'

DimReduction <- methods::setClass("DimReduction",
                                  slots = c(coordinates = "data.frame",
                                            optional_parameters = "list"))

#' @title Kinetic class definition
#'
#' @description The Kinetic object is a S4 object containing all the data concerning a kinetic clustering step performed on a 'CYTdata' object
#'
#' @slot families a data.frame containing the identifiers of metaclusters kinetic families computed with 'run.Kinetic' function
#' @slot optional_parameters a list containing the parameters used for kinetic clustering (especially the arguments of 'run.Kinetic' function : N.kinetics, sample.subset, etc.)
#'
#' @name Kinetic-class
#' @rdname Kinetic-class
#' @exportClass Kinetic
#'

Kinetic <- methods::setClass("Kinetic",
                             slots = c(families = "factor",
                                       optional_parameters = "list"))

