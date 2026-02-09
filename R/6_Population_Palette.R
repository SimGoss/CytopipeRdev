#' @title Check and order population names according to clusters/metaclusters
#'
#' @description This functions checks that names in a character vector are all present in a CYTdata's clustering or metaclustering object
#' If order = TRUE, it also orders the character vector according to the order of the levels of the CYTdata's clusters/metaclusters.
#' If checkDuplicates = TRUE, it also checks for duplicates in the character vector 
#' The function is useful, for example, to check that no mistakes were made when writing a list of cluster names.
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param population a character vector containing the identifiers of the population(s) (clusters or metaclusters) to check and potentially order
#' @param level a character value indicating whether to compare population to the clusters or the metaclusters. Possible values are: "clusters", "metaclusters". By default, 'clusters' are used.
#' @param order logical. If TRUE, population will be ordered according to the order of the levels of the clusters/metaclusters
#' @param checkDuplicates logical. If TRUE, will check for duplicates inside the given population character vector. 
#'
#' @return a character vector
#'
#'@export
#'

checkorderPopulation <- function(CYTdata, population, level, order = TRUE, checkDuplicates = TRUE) {
  
  checkmate::qassert(population, c("0","S*"))
  checkmate::qassert(order, "B1")
  checkmate::qassert(level, "S1")
  level = match.arg(level)
  
  if (!is.null(population) && length(population)==0) {
    stop("Error : 'population' argument is an empty vector (length=0, but non NULL).")
  }
  
  checkmate::qassert(checkDuplicates, "B1")
  if (checkDuplicates) {
    populationdup = population[duplicated(population)]
    if (length(populationdup)>0) {
      stop("Error : 'population' argument contains duplicated values (",
           paste0(populationdup, collapse = ", "),
           "). It must be a vector of unique population levels.")
    }
  }
  
  if (level == "clusters") {
    popId = CYTdata@Clustering@clusters
    if (length(popId)==0) {
      stop("Error : 'level' argument requires clustering step to be performed (Clustering@clusters slot is empty).")
    }
  }
  else {
    popId = CYTdata@Metaclustering@metaclusters
    if (length(popId)==0) {
      stop("Error : 'level' argument requires metaclustering step to be performed (Metaclustering@metaclusters slot is empty).")
    }
  }
  
  if (is.null(population)) { population = levels(popId) }
  else {
    popErr = setdiff(population, levels(popId))
    if (length(popErr)>0) {
      stop("Error : 'population' argument providing identifiers not present in ",
           level, " vector (", paste0(popErr, collapse=", "), ").")
    }
    if (order) {
      # order population vector according to levels
      population = levels(popId)[levels(popId) %in% population]
    }
  }
  
  return(population)
}

#' @title Renames clusters or metaclusters within a CYTdata object.
#'
#' @description This function aims to rename clusters or metaclusters stored within a CYTdata object. 
#' It can rename only a subset of the clusters/metaclusters, or all of them. 
#' If merge = TRUE, clusters/metaclusters which present the same name after renaming will be merged.
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param level a character value indicating whether to rename clusters or metaclusters. Possible values are: "clusters", "metaclusters". By default, 'clusters' are used.
#' @param from a character vector providing the cluster/metacluster names to replace. By default, all of them are used, in the level order.
#' @param to a character vector providing the new cluster/metacluster names to use
#' @param to logical. If TRUE, clusters/metaclusters with duplicated names will be merged. If FALSE (default), duplicated names will result in an error. 
#'
#' @return a S4 object of class 'CYTdata'
#'
#' @export
#'

renamePopulations <-  function(CYTdata, level = c("clusters", "metaclusters"), from = NULL, to, merge = FALSE) {
  
  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }
  
  level = match.arg(level)
  checkmate::qassert(level, "S1")
  from = checkorderPopulation(CYTdata=CYTdata, population=from, level=level, order=FALSE, checkDuplicates=TRUE)
  
  checkmate::qassert(to, "S*")
  checkmate::qassert(merge, "B1")
  
  to_cat = to
  if (length(from)!=length(to)) {
    if (length(to) == 1) { to = rep(to, length(from)) }
    else {
      stop("Error : Length of argument 'from' (", length(from),
           ") and argument 'to' (", length(to), ") must be equal")
    }
  }
  
  if (level=="clusters") {
    popId = CYTdata@Clustering@clusters
    exPalette = CYTdata@Clustering@palette
  }
  else {
    popId = CYTdata@Metaclustering@metaclusters
    exPalette = CYTdata@Metaclustering@palette
  }
  
  cat("\n\nCurrent ", level, " names are (in the order) :", paste0(levels(popId), collapse = ", "), "\n")
  cat("\n\nThe following ", level, " :")
  cat("\n - ", paste0(from, collapse=", "))
  cat("\n\nwill be renamed, in the same order, by :")
  cat("\n - ", paste0(to_cat, collapse=", "))
  
  exLevels = levels(popId)
  newLevels = levels(popId)
  newLevels[match(from, newLevels)] = to
  
  duplicatednewLevels = unique(newLevels[duplicated(newLevels)])
  if (length(duplicatednewLevels)>0){
    if (!merge) {
      stop("Error : After renaming, several ", level, " have the same name (",  paste0(duplicatednewLevels, collapse=", "),
           " ), either by renaming to an already existing and unchanged ", level, " name,
           or by duplicate in the 'to' argument, or both. CYTdata was unchanged (as 'merge' is set to FALSE).")
    }
    else {
      message("Warning : After renaming, several ", level, " have the same name (",  paste0(duplicatednewLevels, collapse=", "),
              " ), either by renaming to an already existing and unchanged ", level, " name,
           or by duplicate in the 'to' argument, or both. These ", level, " will be merged according to their names (as 'merge' is set to TRUE).")
    }
  }
  
  newpopId = plyr::mapvalues(popId, from = from, to = to)
  
  # message("After renaming,  ", level, " levels contained duplicated name(s) and were merged. As it concerns color palette : \n
  #             - If the name of merged ", level, " was already associated to an existing ", level, " name before renaming,
  #             the color of the existing ", level, " is the new color of newly renamed ", level, ".\n
  #             - If the name of merged ", level, " was a new one (argument 'to' contain this name several times).
  #             The color of the resulting  ", level, " is the one of the ", level, " contained in 'from' argument
  #             and which was ordered first (in 'from' argument) among all the ", level, " renamed to this name.")
  
  newpalette = c()
  for (pop in levels(newpopId)) {
    if (sum(newLevels == pop)>1){ # if duplicated
      if (pop %in% exLevels[newLevels == pop]) { v = exPalette[pop] }
      else { pal = exPalette[from[match(TRUE, to==pop)[1]]] }
    }
    else {
      if (pop %in% to) { pal = exPalette[from[to==pop]] }
      else { pal = exPalette[pop] }
    }
    newpalette = c(newpalette, pal)
  }
  names(newpalette) = levels(newpopId)
  
  if (level=="clusters") {
    CYTdata@Clustering@clusters = newpopId
    cat("\n - Computing new cell", level, "count and abundance matrix...")
    CYTdata@Clustering@cellcount = compute.cellcount(CYTdata@Clustering@clusters, CYTdata@samples)
    CYTdata@Clustering@abundance = compute.abundance(CYTdata@Clustering@cellcount)
    CYTdata@Clustering@palette = newpalette
  }
  else {
    CYTdata@Metaclustering@metaclusters = newpopId
    cat("\n - Computing new cell", level, "count and abundance matrix...")
    CYTdata@Metaclustering@cellcount = compute.cellcount(CYTdata@Metaclustering@metaclusters, CYTdata@samples)
    CYTdata@Metaclustering@abundance = compute.abundance(CYTdata@Metaclustering@cellcount)
    CYTdata@Metaclustering@palette = newpalette
  }
  
  validObject(CYTdata) # check if Clustering, Metaclustering slots are ok
  return(CYTdata)
}

#' @title Replaces a CYTdata's clustering slot with its metaclustering slot, or vice-versa
#' 
#' @description This functions copy-pastes the clustering slot into the metaclustering slot or the metaclustering slot into the clustering slot. It does not empty the copy-pasted slot. 
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param direction a string specifying the direction of the change, "clustersToMetaclusters" or "metaclustersToClusters"
#' @param checkOverwrite a boolean value indicating whether to check if a clustering/metaclustering has already been performed on the CYTdata object. 
#' 
#' @return a CYTdata object
#'
#' @export
#'

changeLevel <- function(CYTdata,
                        direction = c("clustersToMetaclusters", "metaclustersToClusters"),
                        checkOverwrite = TRUE) {
  
  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }
  
  direction = match.arg(direction)
  checkmate::qassert(direction, "S1")
  checkmate::qassert(checkOverwrite, "B1")
  
  if (direction == "clustersToMetaclusters") {
    if (length(CYTdata@Clustering@clusters)==0) {
      stop("Error : changeLevel with direction argument equal to 'clustersToMetaclusters' can not be performed on object with empty Clustering slot.
           Please perform Clusetring step")
    }
    if (checkOverwrite && length(CYTdata@Metaclustering@metaclusters)!=0){
      reply <- readline(prompt="Metaclustering slot is not empty, do you still want to continue and overwrite (yes or no): ")
      while (!tolower(reply) %in% c("yes", "y", "no", "n")) { reply <- readline(prompt="Reply not standard, please answer yes or no: ") }
      if (tolower(reply) %in% c("no", "n")){ stop("Function stopped, CYTdata object unchanged") }
    }
    CYTdata@Metaclustering = methods::new("Metaclustering",
                                          metaclusters = CYTdata@Clustering@clusters,
                                          cellcount = CYTdata@Clustering@cellcount,
                                          abundance = CYTdata@Clustering@abundance,
                                          palette = CYTdata@Clustering@palette,
                                          optional_parameters = append(CYTdata@Clustering@optional_parameters, list("changeLevel" = direction)))
  }
  else {
    if (length(CYTdata@Metaclustering@metaclusters)==0) {
      stop("Error : changeLevel with direction argument equal to 'metaclustersToClusters' can not be performed on object with empty Metaclustering slot.
           Please perform Metaclusetring step")
    }
    if (checkOverwrite && length(CYTdata@Clustering@clusters)!=0){
      reply <- readline(prompt="Clustering slot is not empty, do you still want to continue and overwrite (yes or no): ")
      while (!tolower(reply) %in% c("yes", "y", "no", "n")) { reply <- readline(prompt="Reply not standard, please answer yes or no: ") }
      if (tolower(reply) %in% c("no", "n")){ stop("Function stopped, CYTdata object unchanged") }
    }
    CYTdata@Clustering = methods::new("Clustering",
                                      clusters = CYTdata@Metaclustering@metaclusters,
                                      cellcount = CYTdata@Metaclustering@cellcount,
                                      abundance = CYTdata@Metaclustering@abundance,
                                      palette = CYTdata@Metaclustering@palette,
                                      optional_parameters = append(CYTdata@Metaclustering@optional_parameters, list("changeLevel" = direction)))
  }
  validObject(CYTdata)
  return(CYTdata)
}


#' @title Get the names of the clusters belonging to each metacluster. 
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param format a string specifying the returned data's format. Possible values are "list" (default) and "data.frame"
#'
#' @return a list or dataframe indicating the clusters associated with each metacluster
#'
#' @export
#'

getMetaclustersMembership <- function(CYTdata, format = c("list", "data.frame")){
  
  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }
  
  format = match.arg(format)
  checkmate::qassert(format, "S1")
  
  if (length(CYTdata@Metaclustering@metaclusters)==0 || length(CYTdata@Clustering@clusters)==0) {
    stop()
  }
  
  if (format=="list") {
    res = as.list(sapply(levels(CYTdata@Metaclustering@metaclusters),
                         FUN = function(mc){ return(unique(CYTdata@Clustering@clusters[CYTdata@Metaclustering@metaclusters==mc])) }))
    names(res) = levels(CYTdata@Metaclustering@metaclusters)
    return(res)
  }
  else {
    res = sapply(levels(CYTdata@Clustering@clusters),
                 FUN = function(cl){ return(unique(CYTdata@Metaclustering@metaclusters[CYTdata@Clustering@clusters==cl])) })
    res = data.frame("clusters" = levels(CYTdata@Clustering@clusters),
                     "metaclusters" = res)
    return(res)
  }
  
}


#' @title Change color palette associated with clusters/metaclusters
#'
#' @description This function aims to change the color palette associated with a specific set of clusters or metaclusters
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param level a character value indicating whether to change palette on clusters or metaclusters. Possible values are: "clusters", "metaclusters". By default, 'clusters' are used.
#' @param population a character vector containing the identifiers of the population(s) (clusters or metaclusters) for which to change the color. If NULL, (default) all clusters/metaclusters are used
#' @param autoColor logical. Whether to generate a color palette automatically (TRUE), or provide a custom one (FALSE)
#' @param autoColorType character string specifying the type of color palette to generate if autoColor is set to TRUE. Possible values are "rainbow" (default) and "blank", which will set all colors to white (#FFFFFF)
#' @param homemadePalette a named character vector containing the colors (hexadecimal code, or name) to be associated with each cluster/metacluster in population. Necessary if autoColor = FALSE. 
#'
#' @return a S4 object of class 'CYTdata'
#'
#'@export
#'

changePalettes <- function(CYTdata,
                           level = c("clusters", "metaclusters"),
                           autoColor = TRUE,
                           population = NULL,
                           autoColorType = c("rainbow","blank"),
                           homemadePalette = NULL) {
  
  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  #else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }
  
  level = match.arg(level)
  checkmate::qassert(level, "S1")
  checkmate::qassert(autoColor, "B1")
  
  if (!level %in%c("clusters", "metaclusters")) {
    stop("Error : 'level' argument must be either 'clusters' or 'metaclusters'")
  }
  if (level == "clusters") {
    popId = CYTdata@Clustering@clusters
    if (length(popId)==0) {
      stop("Error : 'level' argument require clustering step to be performed (Clustering@clusters slot is empty).")
    }
  }
  else {
    popId = CYTdata@Metaclustering@metaclusters
    if (length(popId)==0) {
      stop("Error : 'level' argument require metaclustering step to be performed (Metaclustering@metaclusters slot is empty).")
    }
  }
  
  if (autoColor) {
    autoColorType = match.arg(autoColorType)
    checkmate::qassert(autoColorType, "S1")
    population = checkorderPopulation(CYTdata, population = population, level = level, order = TRUE, checkDuplicates = TRUE)
    if (autoColorType=="blank") {
      if (level=="clusters") { CYTdata@Clustering@palette[population] = rep("#FFFFFF", length(population)) }
      else { CYTdata@Metaclustering@palette[population] = rep("#FFFFFF", length(population)) }
    }
    else {
      if (level=="clusters") { CYTdata@Clustering@palette[population] = rainbow(length(population)) }
      else { CYTdata@Metaclustering@palette[population] = rainbow(length(population)) }
    }
  }
  else {
    if (is.null(homemadePalette)) { stop("Error : 'autoColor' argument is set to FALSE but no 'homemadePalette' is given") }
    checkmate::qassert(homemadePalette, "S*")
    
    if(!all(areColors(homemadePalette))){
      stop("Error : 'homemadePalette' argument does not contain only hexadecimal color.)")
    }
    population = names(homemadePalette)
    if(sum(is.na(population))>0) {
      stop("Error : 'homemadePalette' argument's names contain NA values.
         It must be a fully named vector of hexadecimal color.
         Associated ", level, " names are missing for colors : ",
           paste0(homemadePalette[is.na(population)], collapse=","), ".")
    }
    populationdup = population[duplicated(population)]
    if (length(populationdup)>0) {
      stop("Error : 'homemadePalette' argument's names contain duplicated values (",
           paste0(populationdup, collapse = ", "),
           "). Names must be vector of unique ", level, " levels.")
    }
    
    if (level == "clusters") { CYTdata@Clustering@palette[population] = homemadePalette[population] }
    else { CYTdata@Metaclustering@palette[population] = homemadePalette[population] }
    
    popErr = setdiff(population, levels(popId))
    if (length(popErr)>0) {
      stop("Error : 'homemadePalette' argument's names providing identifiers not present in ",
           level, " vector (", paste0(popErr, collapse=", "), ").")
    }
  }
  
  nochangePop = setdiff(levels(popId),population)
  if(length(nochangePop)>0){
    cat("\n\n - The following ", level, " :", paste0(nochangePop, collapse=", "), " kept the same colors.")
    cat("\n\n - The following ", level, " :", paste0(population, collapse=", "), " changed colors.")
  }
  else { cat("\n\n - All the ", level, " changed colors.") }
  
  validObject(CYTdata)
  return(CYTdata)
}

#' @title Give clusters the same color as metaclusters
#'
#' @description This function  changes the color of the clusters to match that of their respective metaclusters. 
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param population a character vector containing the identifiers of the clusters for which to change the color. If NULL, (default) all clusters are used
#' @param addVariation logical. If TRUE, some brightness adjustments will be made to the metacluster colors to be able to distinguish the clusters while still retaining the hue of the metacluster. Defaults to FALSE.
#'
#' @return a S4 object of class 'CYTdata'
#'
#'@export
#'

heritatePalette <- function(CYTdata, population = NULL, addVariation = FALSE) {
  
  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }
  
  checkmate::qassert(addVariation, "B1")
  
  popId = CYTdata@Metaclustering@metaclusters
  if (length(CYTdata@Metaclustering@metaclusters)==0) {
    stop("Error : 'heritatePalette' function require metaclustering step to be performed (Metaclustering@metaclusters slot is empty).")
  }
  if (length(CYTdata@Clustering@clusters)==0) {
    stop("Error : 'heritatePalette' function require clustering step to be performed (Clustering@clusters slot is empty).")
  }
  palette = CYTdata@Clustering@palette
  
  metaclusters = checkorderPopulation(CYTdata, level = "metaclusters", order = TRUE, checkDuplicates = TRUE, population = population)
  for (mc in metaclusters) {
    cls = unique(CYTdata@Clustering@clusters[CYTdata@Metaclustering@metaclusters == mc])
    color = CYTdata@Metaclustering@palette[mc]
    n = length(cls)
    if (addVariation) {
      colors = grDevices::colorRampPalette(c(as.character(shades::brightness(color, max(0.4, 0.9-0.5*n/20))),
                                             as.character(shades::brightness(color, min(1.4, 0.9+0.5*n/20)))))(n)
    }
    else { colors = rep(color, n) }
    palette[cls] = colors
  }
  
  CYTdata@Clustering@palette = palette
  validObject(CYTdata)
  return(CYTdata)
}


#' @title Display the color palette associated with clusters or metaclusters as a plot
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param level a character value indicating whether to display clusters or metaclusters. Possible values are: "clusters" (default), "metaclusters".
#' @param population a character vector containing the identifiers of the clusters/metaclusters to display
#' @param ncol integer, number of columns in the graph. Defaults to one per 10 clusters/metaclusters
#' @param labelsize numeric, size of the population labels. Defaults to 5.  
#'
#' @return nothing, but displays a plot. 
#'
#' @export
#'
#'

showPalettes <- function(CYTdata,
                         level = c("clusters", "metaclusters"),
                         population = NULL,
                         ncol = NULL,
                         labelSize = 5) {
  
  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }
  
  level = match.arg(level)
  checkmate::qassert(level, "S1")
  checkmate::qassert(labelSize, "N1")
  checkmate::qassert(ncol, c(0,"N1"))
  if (labelSize<1) { stop("Error : argument 'labelSize' must be positive integer") }
  if (!is.null(ncol) && ncol<1) { stop("Error : argument 'ncol' must be positive integer") }
  population = checkorderPopulation(CYTdata, population = population, level = level, order = TRUE, checkDuplicates = TRUE) # Error messages
  
  print(population)
  if (is.null(ncol)) { ncol = length(population)%/%10 }
  if (level == "clusters"){ palette = CYTdata@Clustering@palette }
  else { palette = CYTdata@Metaclustering@palette }
  
  Q = length(population)%/%ncol
  R = length(population)%%ncol
  print(Q)
  X = rep(1:ncol, Q)
  Y = rep(1:Q, each=ncol)
  if (R!=0) {
    X = c(X, 1:R)
    Y = c(Y, rep(Q+1,R))
  }
  
  data = data.frame("pop" = names(palette),
                    "X" = X,
                    "Y" = Y)
  plot <- ggplot2::ggplot(data, mapping=ggplot2::aes(X, Y, fill= pop)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::geom_label(ggplot2::aes(X, Y, label= pop), fill = "white", size = labelSize) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   axis.title = ggplot2::element_blank(),
                   axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   axis.line = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   legend.position = "none")
  
  plot(plot)
}


#' @title Internal - check if characters in a vector are valid colors
#'
#' @param x a character vector of color names
#'
#' @return a logical vector
#'
#' @export
#'
#'

areColors <- function(x) {
  sapply(x, function(X) {
    tryCatch(is.matrix(col2rgb(X)),
             error = function(e) FALSE)
  })
}




