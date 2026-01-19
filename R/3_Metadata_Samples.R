#' @title Assigns meta-information about biological samples to a CYTdata object
#'
#' @description This function aims to attach meta-information to each biological sample present in a CYTdata object.
#'
#' Especially, the following meta-information of each sample can be specified for subsequent analyses.
#' - The biological individual
#' - The biological condition (groups, vaccines, studies, etc.)
#' - The timepoint
#' Timepoint and Individual data must be specified.
#'
#' @param CYTdata a CYTdata object
#' @param metadata a dataframe containing meta-information about the biological samples.
#' The columns must contain, at least, a column named "Timepoint" and an other named "Individual".
#' The rownames have to be the biological samples, thus the number rows has to be equal to the number of samples.
#'
#' @return a S4 object of class 'CYTdata'
#'
#'
#' @export
#'

importMetadata <- function(CYTdata, metadataDf, checkOverwrite = TRUE){
  
  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }
  
  checkmate::qassert(checkOverwrite, "B1")
  if (checkOverwrite && ncol(CYTdata@metadata)!=0){
    reply = readline(prompt="Metadata already present, do you still want to continue and overwrite (yes or no): ")
    while (!tolower(reply) %in% c("yes", "y", "no", "n")) { reply <- readline(prompt="Reply not standard, please answer yes or no: ") }
    if (tolower(reply) %in% c("no", "n")){ stop("Function stopped, CYTdata object unchanged") }
  }
  
  checkmate::qassert(metadataDf, "D*")
  message("Metadata are imported by data.frame. Checking format..")
  
  if (length(unique(colnames(metadataDf))) != length(colnames(metadataDf))) {
    stop("Error: Colnames are not unique (", paste0(colnames(metadataDf), collapse = ","), ")")
  }
  if (!setequal(levels(CYTdata@samples), rownames(metadataDf))){
    stop("Error: Sample names (rownames) of the metadata (", paste0(rownames(metadataDf), collapse = ","),
         ") are not concordant with the sample levels (", paste0(levels(CYTdata@samples), collapse = ","), ").")
  }
  if (!all(c("Individual", "Timepoint") %in% colnames(metadataDf))) {
    stop("Error: Colnames must contain 'Individual' and 'Timepoint' (But it contains : ", paste0(colnames(CYTdata@metadata), collapse = ","), ")")
  }
  
  cat("\n\nConverting metadata columns to factor..")
  factorByUnique <- function(x) { return(factor(x, levels=unique(x)))}
  metadataDf = dplyr::mutate_if(metadataDf, function(x) { return(!is.factor(x)) }, factorByUnique)
  
  cat("\n\nUpdating metadata...")
  CYTdata@metadata = metadataDf[levels(CYTdata@samples),]
  validObject(CYTdata)
  return(CYTdata)
}

#' @title Checks sample name conformity within a CYTdata object
#'
#' @description This function aims to check for sample name conformity
#' It will check whether a sample list is identical to levels of the CYTdata's samples factor
#' If order = T, it will also order the samples according to the CYTdata's samples factor
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param samples a vector object listing sample names
#' @param order a boolean specifying whether to reorder the samples according to the levels of the CYTdata's samples factor
#' @param checkDuplicates a boolean specifying whether to check for duplicates in the given sample list.
#'
#' @return a S4 object of class 'CYTdata'
#'
#'@export
#'

checkorderSamples <- function(CYTdata, samples, order = TRUE, checkDuplicates = TRUE){
  
  checkmate::qassert(samples, c("0","S*"))
  checkmate::qassert(order, "B1")
  
  if (!is.null(samples) && length(samples)==0) {
    stop("Error : 'samples' argument is an empty vector (length=0, but non NULL).")
  }
  
  checkmate::qassert(checkDuplicates, "B1")
  if (checkDuplicates) {
    samplesdup = samples[duplicated(samples)]
    if (length(samplesdup)>0) {
      stop("Error : 'samples' argument contain duplicated values ( ",
           paste0(samplesdup, collapse = ", "), " ). It must be vector of unique samples levels.")
    }
  }
  
  samplesId = levels(CYTdata@samples)
  if (is.null(samples)) { samples = samplesId }
  else {
    splErr = setdiff(samples, samplesId)
    if (length(splErr)>0) {
      stop("Error : 'samples' argument providing samples not present in
           CYTdata's samples factor vector")
    }
    if (order) {
      samples = samplesId[samplesId %in% samples]
    }
  }
  return(samples)
}


#' @title Renames samples within a CYTdata object
#'
#' @description This function aims to rename samples stored within a CYTdata object. It can also merge samples based on their new names. 
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param from a character vector providing the sample names to replace. By default, all sample names are replaced
#' @param to a character vector providing the new sample names to use
#' @param merge a boolean specifying whether to merge samples that have the same name after renaming. If set to FALSE, duplicated names will result in an error. 
#'
#' @return a S4 object of class 'CYTdata'
#'
#'@export
#'

renameSamples <- function(CYTdata, from = NULL, to, merge = FALSE){
  
  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }
  
  from = checkorderSamples(CYTdata, samples=from, order=FALSE, checkDuplicates=TRUE)
  
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
  
  cat("\n\nCurrent samples names are (in the order) :", paste0(levels(CYTdata@samples), collapse = ", "), "\n")
  cat("\n\nThe following samples :")
  cat("\n - ", paste0(from, collapse=", "))
  cat("\n\nwill be renamed, in the order, by :")
  cat("\n - ", paste0(to_cat, collapse=", "))
  
  exLevels = levels(CYTdata@samples)
  Exmetadata = CYTdata@metadata
  newLevels = levels(CYTdata@samples)
  newLevels[match(from, newLevels)] = to
  
  duplicatednewLevels = unique(newLevels[duplicated(newLevels)])
  if (length(duplicatednewLevels)>0){
    if (!merge) {
      stop("Error : After renaming, several samples have the same name (",  paste0(duplicatednewLevels, collapse=", "),
           " ), either by renaming to an already existing and unchanged sample name,
           or by duplicate in the 'to' argument, or both. CYTdata unchanged ('merge' set to FALSE).")
    }
    else {
      message("Warning : After renaming, several samples have the same name (",  paste0(duplicatednewLevels, collapse=", "),
              " ), either by renaming to an already existing and unchanged sample name,
           or by duplicate in the 'to' argument, or both. These samples will be merged according to their names ('merge' set to TRUE).")
    }
  }
  
  newsamples = plyr::mapvalues(CYTdata@samples, from = from, to = to)
  
  if (nrow(CYTdata@metadata)>0){
    
    message("After renaming, sample levels contained duplicated name(s) and were merged. As it concerns metadata : \n
              - If the name of merged samples was already associated to an existing sample name before renaming,
              the metadata of the existing sample are the new metadata of newly renamed sample(s).\n
              - If the name of merged samples was a new one (argument 'to' contain this name several times).
              The metadata of the resulting sample is the one of the sample contained in 'from' argument
              and which was ordered first (in 'from' argument) among all the samples renamed to this name.")
    
    newmetadata = data.frame()
    for (spl in levels(newsamples)) {
      if (sum(newLevels == spl)>1){ # if duplicated
        if (spl %in% exLevels[newLevels == spl]) { metadataRow = Exmetadata[spl,] }
        else { metadataRow = Exmetadata[from[match(TRUE, to==spl)[1]],] }
      }
      else {
        if (spl %in% to) { metadataRow = Exmetadata[from[to==spl],] }
        else { metadataRow = Exmetadata[spl,] }
      }
      newmetadata = rbind.data.frame(newmetadata, metadataRow)
    }
    rownames(newmetadata) = levels(newsamples)
    CYTdata@metadata = newmetadata
  }
  
  CYTdata@samples = newsamples
  
  if (length(CYTdata@Clustering@clusters)>0) {
    CYTdata@Clustering@cellcount = compute.cellcount(CYTdata@Clustering@clusters, CYTdata@samples)
    CYTdata@Clustering@abundance = compute.abundance(CYTdata@Clustering@cellcount)
  }
  if (length(CYTdata@Metaclustering@metaclusters)>0) {
    CYTdata@Metaclustering@cellcount = compute.cellcount(CYTdata@Metaclustering@metaclusters, CYTdata@samples)
    CYTdata@Metaclustering@abundance = compute.abundance(CYTdata@Metaclustering@cellcount)
  }
  
  # Check if changed sample names, new cellcount and abundance matrixes are valid
  validObject(CYTdata)
  return(CYTdata)
}


#' @title Get the sample names associated to a specific metadata
#' 
#'#' @description This function takes a given metadata (e.g. "Timepoint" or "Individual") and returns a list giving the sample names associated with each level of the condition, e.g, each sample associated to "Timepoint" = "BSL", "D1", etc.  

#' @param CYTdata a S4 object of class 'CYTdata'
#' @param metadataCondition a string specifying the metadata condition chosen. The string has to be the name of the column (in metadata data.frame) containing the condition chosen
#'
#' @return a list containing the samples for each value of condition chosen
#'
#' @export
#'

getSamplesMetadata <- function(CYTdata, metadataCondition){
  
  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }
  
  checkmate::qassert(metadataCondition, "S1")
  
  if (ncol(CYTdata@metadata)==0) {
    stop("Error : No metadata slot present in CYTdata object")
  }
  if (!metadataCondition %in% colnames(CYTdata@metadata)){
    stop("Error : 'metadataCondition' argument metadataCondition is not the
         name of a column present in metadata (", paste0(colnames(CYTdata@metadata), collapse=","), ")")
  }
  
  metacondVect = CYTdata@metadata[[metadataCondition]]
  res = lapply(levels(metacondVect),
               FUN = function(co){
                 spls = subset(CYTdata@metadata, metacondVect == co)
                 return (rownames(spls))
               })
  names(res) = levels(metacondVect)
  return(res)
}

#' @title Plot of the number of cells for each sample
#'
#' @description This function aims to visualize the number of cells associated to each sample. This representation displays the samples in the X-axis and the number of associated cells in the Y-axis.
#' Several statistics can be computed and shown.
#'
#' @details
#' The following statistic can be computed:
#' -'min' corresponds to the lowest number of cells within a data set
#' -'q25' corresponds to the number of cells separates the quantiles 25% within data set
#' -'median' corresponds to the number of cells separates the lower half from the upper half within data set
#' -'mean' corresponds to the number of cells quantity shared within data set
#' -'q75' corresponds to the number of cells separates the quantiles 75% within data set
#' -'max' corresponds to the largest number of cells within a data set
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param stats a character vector providing the statistics to display. Possible values are: 'min', 'median', 'mean', 'q75', 'max'
#' @param samples a character vector containing the names of biological samples to use. By default, all samples are used
#' @param sort a boolean value indicating if samples must be sorted by the number of cells
#'
#' @return a ggplot2 object
#'
#' @export
#'

samplesCellcounts <- function(CYTdata,
                              stats = c("min", "q25", "median", "mean", "q75", "max"),
                              samples = NULL,
                              sort = TRUE) {
  
  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.")}
  checkmate::qassert(stats, "S*")
  checkmate::qassert(samples, c("0", "S*"))
  checkmate::qassert(sort, "B1")
  
  data = cbind.data.frame("sample" = CYTdata@samples, CYTdata@matrix.expression)
  
  if (!is.null(samples)) { data = subset(data, sample %in% samples)}
  
  nb.cells = plyr::ddply(data, "sample", nrow)
  colnames(nb.cells) = c("samples", "counts")
  
  if (sort == TRUE) { nb.cells = nb.cells[order(nb.cells$counts, decreasing = TRUE),] }
  nb.cells$samples = factor(nb.cells$samples, levels = nb.cells$samples)
  
  mean = mean(nb.cells$counts)
  median = stats::median(nb.cells$counts)
  
  xscale.middle <- nrow(nb.cells)/2
  if(xscale.middle%%2==0){ xscale.middle <- xscale.middle + 0.5}
  
  plot <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = nb.cells, ggplot2::aes_string(x="samples", xend="samples", yend="counts"), y=0, color = "gray50") +
    ggplot2::geom_point(data = nb.cells, ggplot2::aes_string(x="samples", y="counts"), size=2)
  
  if("mean" %in% stats){
    mean = mean(nb.cells$counts)
    plot <- plot +
      ggplot2::geom_hline(yintercept=mean, color="green") +
      ggplot2::geom_text(data = nb.cells, ggplot2::aes_string(x="xscale.middle", y="mean"), label=paste0("mean: ", format(mean, scientific = FALSE, big.mark= "," )),color="green", vjust=0)
  }
  
  if("min" %in% stats){
    min = min(nb.cells$counts)
    plot <- plot +
      ggplot2::geom_hline(yintercept=min, color="blue") +
      ggplot2::geom_text(data = nb.cells, ggplot2::aes_string(x="xscale.middle", y="min"), label=paste0("min: ", format(min, scientific = FALSE, big.mark = ",")), color="blue", vjust=0)
  }
  
  if("q25" %in% stats){
    q25 = stats::quantile(nb.cells$counts, probs = 0.25)
    plot <- plot +
      ggplot2::geom_hline(yintercept=q25, color="purple") +
      ggplot2::geom_text(data = nb.cells, ggplot2::aes_string(x="xscale.middle", y="q25"), label=paste0("q25: ", format(q25, scientific = FALSE, big.mark = ",")), color="purple", vjust=0)
  }
  
  if("median" %in% stats){
    median = stats::median(nb.cells$counts)
    plot <- plot +
      ggplot2::geom_hline(yintercept=median, color="red") +
      ggplot2::geom_text(data = nb.cells, ggplot2::aes_string(x="xscale.middle", y="median"),label=paste0("median: ", format(median, scientific = FALSE, big.mark = ",")), color="red", vjust=0)
  }
  
  if("q75" %in% stats){
    q75 = stats::quantile(nb.cells$counts, probs = 0.75)
    plot <- plot +
      ggplot2::geom_hline(yintercept=q75, color="purple") +
      ggplot2::geom_text(data = nb.cells, ggplot2::aes_string(x="xscale.middle", y="q75"), label=paste0("q75: ", format(q75, scientific = FALSE, big.mark = ",")), color="purple", vjust=0)
  }
  
  if("max" %in% stats){
    max = max(nb.cells$counts)
    plot <- plot +
      ggplot2::geom_hline(yintercept=max, color="blue") +
      ggplot2::geom_text(data = nb.cells, ggplot2::aes_string(x="xscale.middle", y="max"), label=paste0("max: ", format(max, scientific = FALSE, big.mark = ",")), color="blue", vjust=0)
  }
  
  plot <- plot +
    ggplot2::scale_y_continuous(labels = function(x){format(x, scientific = FALSE, big.mark=",")}) +
    ggplot2::ylab("Number of cells") +
    ggplot2::xlab("Samples") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5, face = "italic"),
                   axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 15.5),
                   axis.title.x = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   panel.border = ggplot2::element_rect(fill=NA),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   legend.position = "none")
  
  return(plot)
}

