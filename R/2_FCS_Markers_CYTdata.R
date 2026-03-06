#' @title Imports of cell expression profiles from FCS files
#'
#' @description This function aims to import acquired cell events into a CYTdata object.
#'
#' Input files must be FCS files.
#' Different transformations can be applied such as logicle, arcsinh. If yes, parameters are added "..."
#' Cell marker having technical or biological biaises can be excluded during the import.
#'
#' @param files a character vector specifying the path of the tab-separated or FCS files to load
#' @param channels a character vector providing the channels to import from FCS/text files. By default, only de common channels to all the FCS files are kept
#' @param transformList an S4 object of class transformList with features integrated (transform function, channels to transform, parameters). See Flowcore package : https://www.bioconductor.org/packages/release/bioc/manuals/flowCore/man/flowCore.pdf
#' @param rawData a boolean specifying if raW data should be also imported and stored into CYTdata's slot "rawExpression".
#' @param Ndownsampling an integer specifying the number of cells per sample for downsampling
#' @param ... argument to pass onto method read.FCS from Flowcore package
#' @return a S4 object of class 'CYTdata'
#'
#' @importFrom checkmate qassert
#' @importFrom flowCore read.FCS
#'
#' @export
#'

createCYTdata <- function(files,
                          channels,
                          transformList = NULL,
                          rawData = FALSE,
                          Ndownsampling = NULL,
                          ...){
  
  checkmate::qassert(files, "S*")
  checkmate::qassert(channels, "S*")
  if (!is.null(transformList) & class(transformList)!="transformList") {
    stop("Error : argument 'transformList' must be of class 'transformList' or NULL.") }
  checkmate::qassert(rawData, "B1")
  
  if (!is.null(transformList)){
    transformedChannels = names(transformList@transforms)
    cherr = setdiff(transformedChannels, channels)
    if (length(cherr)>0) {
      stop("Error : Some channels that should be transformed are not present in
           'channels' argument (non conformable channels : ",
           paste0(cherr, collapse=", ") ,
           "). Please reset 'channels' or 'transformList'  arguments.")
    }
  }
  
  files_extensions = files %>%
    basename() %>%
    strsplit(split="\\.") %>%
    lapply(FUN = function(x) return(x[-1])) %>%
    unlist()
  if(any(files_extensions != "fcs")){
    message("Warning : Files",
            paste0(basename(files)[files_extensions != "fcs"], collapse = ","),
            "are not FCS files. Removed from import.
            Please convert the file(s) into FCS format,
            import and concatenate with the object created here"
    )
    files = files[files_extensions == "fcs"]
  }
  
  cat("Starting import of FCS files..")
  exprsRaw = data.frame() # Final raw matrix.expression data.frame
  exprs = data.frame() # Final matrix.expression data.frame
  samples = c() # Final samples vector
  i = 1
  for (f in files) {
    
    cat(paste0("\n\nImporting ", basename(f)," file.. (", i, "/", length(files), ")"))
    i = i+1
    fcs = flowCore::read.FCS(f, ...)
    
    cherr = setdiff(channels, flowCore::colnames(fcs))
    if (length(cherr)>0) {
      stop("Error : Channels that should be kept are not present in ",
           basename(f), " file (", paste0(cherr, collapse=", "),
           "). Please reset 'channels' argument.")
    }
    
    ### Import raw data if specified
    if (rawData) {
      newExprsRaw = flowCore::exprs(fcs)[,channels]
      if (!is.null(Ndownsampling)) {
        newExprsRaw = newExprsRaw[sample(1:nrow(newExprsRaw), Ndownsampling, replace = F),]
      }
      
      if (nrow(exprsRaw)==0){
        exprsRaw = newExprsRaw
      } else {
        exprsRaw = rbind.data.frame(exprsRaw, newExprsRaw)
      }
    }
    
    ### Transform data if specified
    if (!is.null(transformList)){
      fcs = flowCore::transform(fcs, transformList)
    }
    
    ### Store data into "exprs" dataframe, get markers vector
    if (nrow(exprs)==0){
      exprs = flowCore::exprs(fcs)[,channels]
      if (!is.null(Ndownsampling)) {
        exprs = exprs[sample(1:nrow(exprs), Ndownsampling, replace = F),]
      }
      newSample = gsub(pattern = ".fcs", replacement = "", x = basename(f))
      samples = rep(newSample, nrow(exprs))
      samplesLevels = newSample
      cat("\n - Number of events for sample", basename(f), ":", nrow(exprs))
      FileMarkers = flowCore::markernames(fcs)[channels]
      FileMarkers = gsub(pattern = "-|\\.", replacement = "_", x = FileMarkers)
      if (any(is.na(FileMarkers))) {
        message("Removing channel(s) : ",
                paste0(channels[is.na(FileMarkers)], collapse=","),
                ". Because associated markers is NA value (first file).")
        channels = channels[!is.na(FileMarkers)]
      }
      
    }
    else {
      newExprs = flowCore::exprs(fcs)[,channels]
      if (!is.null(Ndownsampling)) {
        newExprs = newExprs[sample(1:nrow(newExprs), Ndownsampling, replace = F),]
      }
      newSample = gsub(pattern = ".fcs", replacement = "", x = basename(f))
      samples = c(samples, rep(newSample, nrow(newExprs)))
      samplesLevels = c(samplesLevels, newSample)
      cat("\n - Number of events for sample", basename(f), ":", nrow(newExprs))
      exprs = rbind.data.frame(exprs, newExprs)
    }
    
  }
  
  cat("\nTotal number of events :", nrow(exprs))
  
  ### Convert colnames (channels) to markers
  colnames(exprs) = FileMarkers
  rownames(exprs) = NULL
  samples = factor(samples, levels = samplesLevels)
  
  cat("\n\nCreating CYTdata object..")
  CYTdata <- methods::new("CYTdata",
                          matrix.expression = exprs,
                          samples = samples,
                          Clustering = Clustering(),
                          Metaclustering = Metaclustering(),
                          DimReduction = DimReduction(),
                          Kinetic = Kinetic())
  
  if (rawData){
    colnames(exprsRaw) = FileMarkers
    rownames(exprsRaw) = NULL
    print(exprsRaw)
    CYTdata@raw.matrix.expression = exprsRaw
  }
  
  validObject(CYTdata)
  return(CYTdata)
}

#' @title Internal - Creates of a FlowFrame object
#'
#' @description This function is used internally to create a flowframe object with the purpose of extracting marker intensities to FCS files.
#'
#' @param intensities a data.frame providing the cell profile intensities
#'
#' @return a flowframe object containing the marker expression
#'


createFlowframe <- function(intensities) {
  
  markers <- colnames(intensities)
  p <- c()
  description <- list()
  
  description[["$DATATYPE"]] <- "F"
  
  for (i in seq(1, ncol(intensities))) {
    name <- markers[i]
    min <- min(intensities[, i])
    max <- max(intensities[, i])
    range <- max - min + 1
    
    l <- matrix(c(name, name, range, min, max), nrow = 1)
    colnames(l) <- c("name", "desc", "range", "minRange", "maxRange")
    rownames(l) <- paste0("$P", i)
    p <- rbind(p, l)
    
    description[[paste("$P", i, "N", sep = "")]] <- name
    description[[paste("$P", i, "S", sep = "")]] <- name
    description[[paste("$P", i, "R", sep = "")]] <- toString(range)
    description[[paste("$P", i, "B", sep = "")]] <- "32"
    description[[paste("$P", i, "E", sep = "")]] <- "0,0"
  }
  
  intensities <- as.matrix(intensities)
  dataframe <- methods::as(data.frame(p), "AnnotatedDataFrame")
  
  flowframe <- suppressWarnings(flowCore::flowFrame(intensities, dataframe, description = description))
  
  return(flowframe)
}


#' @title Exports cell expression matrix to TSV or FCS files
#'
#' @description Exports cell expression matrix from a CYTdata object to a tab-separated or FCS files. One FCS file by sample (according to sample CYTdata sample vector)
#' Cell expression profiles can be exported for a set of given samples and for a set of given cell population
#'
#' @param CYTdata a CYTdata object
#' @param outputDir a character value providing the path of folder where the output files are stored
#' @param markers a character vector containing the markers of which the expression is to export.
#' @param samples a character vector containing the names of biological samples to export. By default, all samples are extracted.
#' @param population a character vector containing the identifiers of the population to export.
#' If set to NULL (default), the selection by population is not used (useful if no clustering/metaclustering step has been performed)
#' @param level a character value indicating the type of population extracted. Possible values are: "clusters", "metaclusters". By default, 'clusters' are used.
#' @param concatenate.samples If set to TRUE, sample expression are concatenated into the same FCS files. By default, one FCS is written file by sample.
#' @param extension a character being the type of file into which data is exported. Possible values are: "FCS", "TXT" (default =  "FCS")
#'
#' @return none
#'
#' @export
#'


export.CYTdata.expression <- function(CYTdata,
                                      outputDir = ".",
                                      markers = NULL,
                                      samples = NULL,
                                      population = NULL,
                                      level = c("clusters", "metaclusters"),
                                      concatenate.samples = FALSE,
                                      extension = c("FCS", "TXT")) {
  
  checkmate::qassert(outputDir, "S1")
  checkmate::qassert(markers, c("0", "S*"))
  checkmate::qassert(samples, c("0", "S*"))
  checkmate::qassert(population, c("0", "S*"))
  level = match.arg(level)
  checkmate::qassert(level, "S1")
  checkmate::qassert(concatenate.samples, "B1")
  extension = match.arg(extension)
  checkmate::qassert(extension, "S1")
  
  if (is.null(markers)) { markers = colnames(CYTdata@matrix.expression) }
  data = cbind.data.frame(CYTdata@matrix.expression[, markers], "sample" = CYTdata@samples)
  
  if (!is.null(population)){
    if (level == "clusters"){
      if (length(CYTdata@Clustering@clusters)==0) { stop("Error : Clustering step required before exportation. Please set 'population' argument to NULL") }
      popId = CYTdata@Clustering@clusters
    }
    if (level == "metaclusters"){
      if (length(CYTdata@Metaclustering@metaclusters)==0) { stop("Error : Metaclustering step required before exportation. Please set 'population' argument to NULL") }
      popId = CYTdata@Metaclustering@metaclusters
    }
    data = subset(data, popId %in% population)
    data$popId = NULL
  }
  
  if (is.null(samples)) { samples = unique(CYTdata@samples)}
  data = subset(data, sample %in% samples)
  
  if (concatenate.samples) {
    if (extension == "FCS"){
      data$sample <- as.numeric(factor(data$sample))
      flowFrame <- createFlowframe(data)
      filename = paste(rlang::as_string(CYTdata), "_expression_matrix.fcs")
      flowCore::write.FCS(flowFrame, filename = paste(outputDir, filename, sep="/"))
    }
    else {
      data$sample <- as.numeric(factor(data$sample))
      filename = paste(rlang::as_string(CYTdata), "_expression_matrix.txt")
      utils::write.table(data,
                         file = paste(outputDir, filename, sep="/"),
                         sep = "\t",
                         row.names = FALSE)
    }
  }
  else {
    for (spl in samples){
      data.spl = subset(data, sample==spl)
      data.spl$sample = NULL
      if (extension == "FCS"){
        flowFrame <- createFlowframe(data.spl)
        filename = paste0(spl, ".fcs")
        flowCore::write.FCS(flowFrame, filename = paste(outputDir, filename, sep="/"))
      }
      else {
        filename = paste0(spl, ".txt")
        utils::write.table(data.spl,
                           file = paste(outputDir, filename, sep="/"),
                           sep = "\t",
                           row.names = FALSE)
      }
    }
  }
}


#' @title Get the markers (and channels associated) common to a set of FCS or TXT files
#'
#' @description This function aims at importing marker and channel names from a set of FCS files
#'
#' Input files must be FCS or TXT files.
#' Only common channels and markers are returned
#'
#' @param files a character vector specifying the path of the tab-separated or FCS files to load
#' @param removeTime a boolean specifying whether to remove the Time column from the data
#' @param verbose a boolean specifying whether to return additionnal information
#'
#' @return a named character vector
#'
#' @export
#'

getChannelsMarkersFCS <- function(files, removeTime = TRUE, verbose = FALSE){
  
  checkmate::qassert(files, "S*")
  checkmate::qassert(removeTime, "B1")
  checkmate::qassert(verbose, "B1")
  
  files_extensions = files %>%
    basename() %>%
    strsplit(split="\\.") %>%
    lapply(FUN = function(x) return(x[-1])) %>%
    unlist()
  if(any(files_extensions != "fcs")){
    message("Error : Files",
            paste0(basename(files)[files_extensions != "fcs"], collapse = ","),
            "are not FCS files. Removed from import.
            Please convert the file(s) into FCS format,
            import and concatenate with the object created here"
    )
    files = files[files_extensions == "fcs"]
  }
  
  common.channels = c()
  ref.markers = c()
  for (i in 1:length(files)) {
    
    f = files[i]
    fcs = flowCore::read.FCS(f, which.lines = 1:2, truncate_max_range = FALSE)
    
    if(length(common.channels)==0 && length(ref.markers)==0) {
      cat("\n\n\n Checking, first, ", basename(f), " 's channels and markers :")
      common.channels = flowCore::colnames(fcs)
      ref.markers = flowCore::markernames(fcs)[common.channels]
      ref.markers = gsub(pattern = "-|\\.", replacement = "_", x = ref.markers)
      
      if (any(is.na(ref.markers))) {
        message("\n\nWarning : Removing channel(s) \n - ",
                paste0(common.channels[is.na(ref.markers)], collapse=","),
                ".\n Because associated markers is NA value (first file).")
        common.channels = common.channels[!is.na(ref.markers)]
        ref.markers = ref.markers[!is.na(ref.markers)]
      }
      
      if (verbose) {
        cat("\n\n - Channels :", paste0(common.channels, collapse = ", "))
        cat("\n\n - Markers :", paste0(ref.markers, collapse = ", "))
      }
    }
    else {
      if (verbose) { cat("\n\n Checking", basename(f), "'s channels and markers :") }
      newChannels = flowCore::colnames(fcs)
      newMarkers = flowCore::markernames(fcs)[common.channels]
      newMarkers = gsub(pattern = "-|\\.", replacement = "_", x = newMarkers)
      
      if (removeTime){ newChannels = newChannels[newChannels != "Time"]}
      
      if (setequal(common.channels, newChannels)) {
        if (verbose) { cat("\n\n - ", basename(f), "contains the same channels than previous files.") }
      }
      else {
        non.common.channels1 = setdiff(common.channels, newChannels)
        non.common.channels2 = setdiff(newChannels, common.channels)
        if (length(non.common.channels1)>0){
          message("\n\n - Warning : ", paste0(non.common.channels1, collapse = ", "),
                  " channel(s) are present in previous FCS files but not in ", basename(f), " file.")
        }
        if (length(non.common.channels2)>0){
          message("\n\n - Warning : ", paste0(non.common.channels2, collapse = ", "),
                  " channel(s) are present in ", basename(f), " file but not in previous FCS files.")
        }
        if (verbose) { cat("\n -> Keeping the common channels..") }
        common.channels = intersect(common.channels, newChannels)
      }
      
      if (setequal(ref.markers[common.channels], newMarkers[common.channels])) {
        if (verbose) { cat("\n - Concordant markers, for common channels, between ", basename(f),
                           " file and previous FCS file.\n\n") }
      }
      else {
        message("\n - Warning : Common channels but non concordant markers between ",
                basename(f), " file (",
                paste0(setdiff(newMarkers[common.channels], ref.markers[common.channels]), collapse = ", "),
                ") and previous FCS files (",
                paste0(setdiff(ref.markers[common.channels], newMarkers[common.channels]), collapse = ", "),
                ").")
      }
    }
  }
  if (verbose) {
    print(common.channels)
    cat("Returning common channels to all the files and associated markers according to very first file..\n\n")
  }
  
  res = list("Markers" = ref.markers[common.channels],
             "Channels" = common.channels)
  
  return(res)
}


#' @title Checks marker name conformity within a CYTdata object
#'
#' @description This function aims to check for marker name conformity
#' It will check whether a marker list is identical to the colnames of the CYTdata's marker expression matrix
#' If order = T, it will also order the marker according to colnames. 
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param markers a vector object listing marker names
#' @param order a boolean specifying whether to reorder the markers according to colnames of the CYTdata's matrix expression
#' @param checkDuplicates a boolean specifying whether to check for duplicates in the markers list
#'
#' @return a S4 object of class 'CYTdata'
#'
#'@export
#'

checkorderMarkers <- function(CYTdata, markers, order = TRUE, checkDuplicates = TRUE){
  
  checkmate::qassert(markers, c("0","S*"))
  checkmate::qassert(order, "B1")
  
  if (!is.null(markers) && length(markers)==0) {
    stop("Error : markers argument is an empty vector (length=0, but non NULL).")
  }
  
  checkmate::qassert(checkDuplicates, "B1")
  if (checkDuplicates) {
    markersdup = markers[duplicated(markers)]
    if (length(markersdup)>0) {
      stop("Error : markers argument contain duplicated values ( ",
           paste0(markersdup, collapse = ", "), " ). It must be vector of unique markers.")
    }
  }
  
  markersId = colnames(CYTdata@matrix.expression)
  if (is.null(markers)) { markers = markersId }
  else {
    markErr = setdiff(markers, markersId)
    if (length(markErr)>0) {
      stop("Error : 'markers' argument providing markers not present in
           matrix.expression's colnames (", paste0(markErr, collapse=", "), ")")
    }
    if (order) {
      markers = markersId[markersId %in% markers]
    }
  }
  return(markers)
}

#' @title Renames markers within a CYTdata object
#'
#' @description This function aims to rename cell markers stored within a CYTdata object.
#'
#' This function is interesting to remove the names of the fluorochromes or metals recorded during the acquisition process.
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param from a character vector providing the marker names to replace. By default, markers names from markers slot are replaced, in the order
#' @param to a character vector providing the new marker names to use
#' @param removeConjugate a boolean to indicate whether to remove the heavy metal name associated to the marker. Faster than specifying "to" manually, but specific. 
#' @param sep a character specifying the separator of the metal and the marker name. Ignored if removeConjugate = F
#'
#' @return a S4 object of class 'CYTdata'
#'
#'@export
#'

renameMarkers <- function(CYTdata, from=NULL, to=NULL, removeConjugate=FALSE, sep="_"){
  
  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }
  
  from = checkorderMarkers(CYTdata, markers=from, order=FALSE, checkDuplicates=TRUE)
  
  checkmate::qassert(removeConjugate, "B1")
  if (removeConjugate) {
    checkmate::qassert(sep, "S1")
    message("Warning : 'removeConjugate' argument is set to TRUE, 'to' argument is ignored.
            The first part (before undescore) of feature names given in 'from' argument are removed and features are renamed." )
    to = sapply(stringr::str_split(from, "_"), function(x){return(paste(x[-1], collapse=sep))})
  }
  else {
    checkmate::qassert(to, "S*")
    if (length(from)!=length(to)) {
      stop("Error : Length of argument 'from' (", length(from),
           ") and argument 'to' (", length(to), ") must be equal")
    }
  }
  
  cat("\n\nCurrent marker names are (in the order) :", paste0(colnames(CYTdata@matrix.expression), collapse = ", "))
  cat("\n\n\nThe following markers :")
  cat("\n\n - ", paste0(from, collapse=", "))
  cat("\n\n\nwill be renamed, in the order, by :")
  cat("\n - ", paste0(to, collapse=", "))
  
  # rename markers vector
  indexes = match(from, colnames(CYTdata@matrix.expression))
  colnames(CYTdata@matrix.expression)[indexes] = to
  if (ncol(CYTdata@raw.matrix.expression)>0) { colnames(CYTdata@raw.matrix.expression) = colnames(CYTdata@matrix.expression) }
  
  # Check if changed markers names are valid (no NA, no duplicated values)
  validObject(CYTdata)
  return(CYTdata)
}

#' @title TO TEST Removes duplicated columns in a marker expression matrix
#'
#' @description This function aims to remove duplicated marker columns by checking whther marker expression is identicial across any pair of columns. 
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param markers a character vector specifying a the name of the markers to test for duplication. If NULL, all markers will be tested. 
#'
#' @return a S4 object of class 'CYTdata'
#'
#'@export
#'

removeDuplicates <- function(CYTdata, markers = NULL) {
  
  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }
  
  markers = checkorderMarkers(CYTdata, markers = markers, order = TRUE, checkDuplicates = TRUE)
  if (length(markers)!=ncol(CYTdata@matrix.expression)) {
    message("\nCells will be removed if their expression duplicated for combination of following markers : ", paste0(markers, collapse=", "))
  }
  
  idxDuplicates = duplicated(CYTdata@matrix.expression[,markers])
  if (sum(idxDuplicates)==0) {
    message("\nNo duplicated data found, identical CYTdata object returned")
    return(CYTdata)
  }
  else {
    message("\n", sum(idxDuplicates), " cell(s) is/are duplicated. Removing it from CYTdata object..")
    gatedIdx = !idxDuplicates
    newmatrix.expression = subset(CYTdata@matrix.expression, gatedIdx)
    
    newsamples = CYTdata@samples[gatedIdx]
    newsamples = droplevels(newsamples)
    removedSpls = setdiff(levels(CYTdata@samples), levels(newsamples))
    if(length(removedSpls)>0){
      cat("\n\n - Gating operation remove following samples :", paste0(removedSpls, collapse = ", "))
    }
    
    cat("\n\n - Creating new CYTdata object..")
    newCYTdata = methods::new("CYTdata",
                              samples = newsamples,
                              matrix.expression = newmatrix.expression)
    if (ncol(CYTdata@raw.matrix.expression)>0) {
      newCYTdata@raw.matrix.expression = subset(CYTdata@raw.matrix.expression, gatedIdx)
    }
    
    if (nrow(CYTdata@metadata)>0) {
      newCYTdata@metadata = subset(CYTdata@metadata, rownames(CYTdata@metadata) %in% levels(newCYTdata@samples))
    }
    
    if(length(CYTdata@Clustering@clusters)>0){
      cat("\n\n - Updating Clustering slot")
      newClusters = CYTdata@Clustering@clusters[gatedIdx]
      newClusters = droplevels(newClusters)
      newpalette = CYTdata@Clustering@palette[levels(newClusters)]
      cat("\ncomputing new cell cluster count, abundance matrix...")
      newcellcount = compute.cellcount(newClusters, newsamples)
      newabundance = compute.abundance(newcellcount)
      cat("\nCreating new Clustering object")
      newCYTdata@Clustering = methods::new("Clustering",
                                           clusters = newClusters,
                                           cellcount = newcellcount,
                                           abundance = newabundance,
                                           palette = newpalette,
                                           optional_parameters = CYTdata@Clustering@optional_parameters)
      message("\n\nRemark : Clustering results are preserved during gating operation but it is recommended to the user to run a
        new clustering step with parameters adapted to the gated dataset")
    }
    if(length(CYTdata@Metaclustering@metaclusters)>0){
      cat("\n\n - Updating Metaclustering slot")
      newMetaclusters = CYTdata@Metaclustering@metaclusters[gatedIdx]
      newMetaclusters = droplevels(newMetaclusters)
      newpalette = CYTdata@Metaclustering@palette[levels(newMetaclusters)]
      cat("\ncomputing new cell metacluster count, abundance matrix...")
      newcellcount = compute.cellcount(newMetaclusters, newsamples)
      newabundance = compute.abundance(newcellcount)
      cat("\nCreating new Metaclustering object")
      newCYTdata@Metaclustering = methods::new("Metaclustering",
                                               metaclusters = newMetaclusters,
                                               cellcount = newcellcount,
                                               abundance = newabundance,
                                               palette = newpalette,
                                               optional_parameters = CYTdata@Metaclustering@optional_parameters)
    }
    if(nrow(CYTdata@DimReduction@coordinates)>0){
      cat("\n\n - Updating DimReduction slot")
      newcoordinates = subset(CYTdata@DimReduction@coordinates, gatedIdx)
      newCYTdata@DimReduction = methods::new("DimReduction",
                                             coordinates = newcoordinates,
                                             optional_parameters = CYTdata@DimReduction@optional_parameters)
    }
    validObject(newCYTdata)
    return(newCYTdata)
  }
}
