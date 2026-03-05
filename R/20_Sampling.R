
#' @title Density-based downsampling (Internal)
#' @description Downsample indexes representing sample's events with a density-based method (spade)
#'
#' @param indexes.to.downsample a numeric vector specifying the indexes of a sample's events in the dataset
#' @param target.size a numeric value specifying the number of events targeted as downsampled, is smaller than the number of indexes to downsample
#' @param whole.data.exprs a data.frame containing the entire dataset : all the markers' expression (used for density computation)
#' @param exclude.pctile a numeric value, between 0 and 1, specifying the quantile value used to determine the exclusion threshold for density values
#' @param compensate.uniform a boolean value : Density based downsampling doesn't permit us to get the exact size specified by "target.size". So, do we add a step to adjust randomly the donsampled indexes ?
#' @param ncores numeric, number of processor cores to use.
#'
#' @return downsampled.indexes a numeric vector specifying the indexes downsampled

downsampling.density <- function(indexes.to.downsample,
                                 target.size,
                                 whole.data.exprs,
                                 exclude.pctile,
                                 compensate.uniform,
                                 ncores){
  cat("- Exact size to downsample : ", target.size)
  
  data.exprs = whole.data.exprs[indexes.to.downsample,]
  
  ### SPADE density
  df.dists = parallelDist::parDist(x = as.matrix(data.exprs), method = "euclidean", threads = ncores) %>% as.matrix() %>% as.data.frame()
  idx.used.df.dists = sample(1:nrow(df.dists), 3*nrow(df.dists)/4)
  med.dists = apply(df.dists[idx.used.df.dists,], 1, min, na.rm=TRUE) %>% as.vector() %>% median()
  kernel.width = 20*med.dists
  boof.df.dists = (df.dists < kernel.width) %>% as.data.frame()
  density = apply(boof.df.dists, 1, sum, na.rm=TRUE) %>% as.vector()
  #################################
  
  data = cbind.data.frame("density" = density, "indexes" = indexes.to.downsample,  data.exprs)
  
  print(data)
  exclusion.boundary = stats::quantile(density, exclude.pctile, names=FALSE)
  
  ## First step of DS : Spade exclusion
  data = subset(data, density > exclusion.boundary) # Selection by exclude.pctile
  
  print(data)
  
  n.compensate = target.size - nrow(data)
  if (n.compensate > 0){
    cat("- Number of cells, after spade exclusion : ", nrow(data), ". So no more density downsampling. Using uniform sampling to compensate and sample the exact size..")
    remaining.idx = indexes.to.downsample[!indexes.to.downsample %in% data$indexes]
    compensate.idx = sample(remaining.idx, n.compensate)
    downsampled.indexes = c(data$indexes, compensate.idx)
  }
  else {
    ## Second step of DS : Target size and boundary estimation
    density.sorted = sort(data$density)
    cdf = rev(cumsum(1/rev(density.sorted)))
    target.boundary = target.size/cdf[1]
    if (target.boundary > density.sorted[1]) {
      targets = (target.size-seq(1,length(density.sorted))) / cdf
      target.boundary = targets[which.min(targets-density.sorted > 0)]
    }
    data = subset(data, target.boundary/density > stats::runif(nrow(data))) # Selection by ~target.size
    cat("- For a sample, the number of cells, after the finished downsampling (density) process, is : ", nrow(data))
    n.compensate = target.size - nrow(data)
    if (n.compensate > 0){
      cat("- Using uniform sampling to compensate and sample the exact size..")
      remaining.idx = indexes.to.downsample[!indexes.to.downsample %in% data$indexes]
      compensate.idx = sample(remaining.idx, n.compensate)
      downsampled.indexes = c(data$indexes, compensate.idx)
    }
    else if (n.compensate < 0){
      cat("- Using uniform sampling to compensate and sample the exact size..")
      n.compensate = abs(n.compensate)
      compensate.idx = sample(data$indexes, n.compensate)
      downsampled.indexes = data$indexes[!data$indexes %in% compensate.idx]
    }
  }
  return(downsampled.indexes)
}


#' @title Performs the downsampling of events
#' @description Perform the downsampling, using uniformly-based or density-based random selections, of events
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param type a string value specifying the type of downsampling used, "uniform" (default) or "density" (density-based). Density based downsampling may be slow. 
#' @param sizeType a string value specifying the way the number of downsampled events is determined ("bySplPercent" (default) to use a fraction of each sample's original size, "bySplSmall" to use the smallest sample's size, "bySplAbsolute" to use an absolute size per sample, "byOverallAbsolute" to use an absolute size for the whole dataset)
#' @param sizeAbsolute a numeric value specifying the target number of cells, either by sample (if sizeType = "bySplAbsolute") or overall (if sizeType = "byOverallAbsolute"). Defaults to 1000.
#' @param sizePercent a numeric value specifying the target fraction of sample size for downsampling (expected if sizeType = "bySplPercent"). Defaults to 0.2. 
#' @param seed a numeric value: the seed for RNG.
#' @param density.exclusion.pctile a numeric value, between 0 and 1, specifying the quantile value used to determine the exclusion threshold for density values. Defaults to 0.01
#' @param density.compensate logical. Density based downsampling doesn't allow getting the exact size specified by "target.size". If TRUE, a step is added to adjust randomly the downsampled indexes. Defaults to FALSE. 
#' @param density.ncores a numeric value specifying the number of logical processor used for parallel computation in the case of density based downsampling
#' 
#' @return a S4 object of class 'CYTdata'
#'
#' @export
#'

performDownsampling <- function(CYTdata,
                                type = c("uniform", "density"),
                                sizeType = c("bySplPercent", "bySplSmall",
                                             "bySplAbsolute", "byOverallAbsolute"),
                                sizeAbsolute = 1000,
                                sizePercent = 1/5,
                                seed = 42,
                                density.exclusion.pctile = 0.01,
                                density.compensate = FALSE,
                                density.ncores = 1){
  
  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }
  
  type = match.arg(type)
  checkmate::qassert(type, "S1")
  cat("\nDownsampling will be performed using", type, "sampling ('type' argument)")
  
  sizeType = match.arg(sizeType)
  checkmate::qassert(sizeType, "S1")
  
  checkmate::qassert(seed, "N1")
  set.seed(seed)
  switch(type,
         uniform = {
           downsampling.function <- function(indexes, size) {
             return(sample(indexes, size)) # replace = FALSE by default
           }
         },
         density = {
           checkmate::qassert(density.exclusion.pctile, "N1")
           checkmate::qassert(density.compensate, "B1")
           checkmate::qassert(density.ncores, "N1")
           list.density.parameters = list("density.exclusion.pctile" = density.exclusion.pctile,
                                          "density.compensate" = density.compensate,
                                          "density.ncores" = density.ncores)
           cat("\n\n - The downsampling density parameters are :",
               paste0(names(list.density.parameters), "=", list.density.parameters, collapse = ", "))
           downsampling.function <- function(indexes, size) {
             return(downsampling.density(indexes.to.downsample = indexes,
                                         target.size = size,
                                         whole.data.exprs = CYTdata@matrix.expression,
                                         exclude.pctile = density.exclusion.pctile,
                                         compensate.uniform = density.compensate,
                                         ncores = density.ncores) )}
         })
  
  checkmate::qassert(sizePercent, "N1")
  if (!(sizePercent>=0 && sizePercent<=1)) {
    stop("Error: 'sizePercent' argument is a proportion and must be a positive numeric between 0 and 1.")
  }
  checkmate::qassert(sizeAbsolute, "N1")
  if (sizeAbsolute<0) {
    stop("Error: 'sizeAbsolute' argument must be a positive integer.")
  }
  
  switch(sizeType,
         bySplPercent = {
           cat("\n'sizeType' argument set to 'bySplPercent': Downsampling will be performed independently on each samples.
               The resulting subsampled object will contain a different amount of cells by samples, equal to ",
               sizePercent*100, "% of sample size.")
           sizes = round(table(CYTdata@samples)*sizePercent)
           names(sizes) = levels(CYTdata@samples)
         },
         bySplSmall = {
           cat("\n'sizeType' argument set to 'bySplSmall': Downsampling will be performed equally across the samples.
               The resulting subsampled object will contain the same amount of cells in each samples equal to the smallest sample's size")
           sizes = rep(min(table(CYTdata@samples)), length(levels(CYTdata@samples)))
           names(sizes) = levels(CYTdata@samples)
         },
         bySplAbsolute = {
           cat("\n'sizeType' argument set to 'bySplAbsolute': Downsampling will be performed equally across the samples.
               The resulting subsampled object will contain the same amount of cells in each samples equal to ", sizeAbsolute, " cells.")
           sizes = rep(sizeAbsolute, length(levels(CYTdata@samples)))
           names(sizes) = levels(CYTdata@samples)
         },
         byOverallAbsolute = {
           cat("\n'sizeType' argument set to 'byOverallAbsolute': Downsampling will be performed on the whole dataset without distinguishing each sample.
               The resulting subsampled object will contain ", sizeAbsolute, " cells.")
           sizes = sizeAbsolute
         })
  
  if (length(sizes)==1) {
    DownsamplingIndexes = downsampling.function(indexes = 1:nrow(CYTdata@matrix.expression),
                                                size = sizes)
  }
  else {
    DownsamplingIndexes = c()
    for (spl in names(sizes)) {
      size = sizes[spl]
      splindexes = which(CYTdata@samples == spl)
      if (length(splindexes) <= size){
        #cat("\n - No downsampling for sample ", spl, "(enough events already or sample size smaller than downsampling target)")
        DSindexes = splindexes
      }
      else {
        #cat("\n - ", size, " events to downsample, over ", length(splindexes), ", for sample ", spl)
        DSindexes = downsampling.function(indexes = splindexes, size = size)
      }
      DownsamplingIndexes = c(DownsamplingIndexes, DSindexes)
      #cat("\n - Downsampling done for sample ", spl)
    }
  }
  
  cat("\n\nUpdating CYTdata object..\n")
  CYTdata@samples = CYTdata@samples[DownsamplingIndexes]
  CYTdata@samples = droplevels(CYTdata@samples)
  CYTdata@matrix.expression = CYTdata@matrix.expression[DownsamplingIndexes,]
  if (nrow(CYTdata@raw.matrix.expression)>0) {
    CYTdata@raw.matrix.expression = CYTdata@raw.matrix.expression[DownsamplingIndexes,]
  }
  if (length(CYTdata@Clustering@clusters)>0) {
    cat("\n - Downsampling clusters identifiers")
    CYTdata@Clustering@clusters = CYTdata@Clustering@clusters[DownsamplingIndexes]
    CYTdata@Clustering@clusters = droplevels(CYTdata@Clustering@clusters)
    CYTdata@Clustering@cellcount = compute.cellcount(CYTdata@Clustering@clusters,
                                                     CYTdata@samples)
    CYTdata@Clustering@abundance = compute.abundance(CYTdata@Clustering@cellcount)
    CYTdata@Clustering@palette = CYTdata@Clustering@palette[levels(CYTdata@Clustering@clusters)]
  }
  if (length(CYTdata@Metaclustering@metaclusters)>0) {
    cat("\n - Downsampling metaclusters identifiers")
    CYTdata@Metaclustering@metaclusters = CYTdata@Metaclustering@metaclusters[DownsamplingIndexes]
    CYTdata@Metaclustering@metaclusters = droplevels(CYTdata@Metaclustering@metaclusters)
    CYTdata@Metaclustering@cellcount = compute.cellcount(CYTdata@Metaclustering@metaclusters,
                                                         CYTdata@samples)
    CYTdata@Metaclustering@abundance = compute.abundance(CYTdata@Metaclustering@cellcount)
    CYTdata@Metaclustering@palette = CYTdata@Metaclustering@palette[levels(CYTdata@Metaclustering@metaclusters)]
  }
  cat("\n - Removing DimReduction slot if existing..")
  if (nrow(CYTdata@DimReduction@coordinates)>0) {
    CYTdata@DimReduction = methods::new("DimReduction",
                                        coordinates = data.frame(),
                                        optional_parameters = list())
  }
  validObject(CYTdata) # Check if Sampling.info object and changes on samples, matrix.expression are valid
  return(CYTdata)
}

#' @title Performs the upsampling of downsampled events
#'
#' @description This function aims to perform the upsampling of downsampled events events based on an existing CYTdata object and existing cell events stored in tab-separated or FCS files.
#'
#' Importantly, the identification of cell clusters must have been performed prior to this operation.
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param newCYTdata a S4 object of class 'CYTdata'
#' @param markers a character vector providing the markers to use to perform upsampling
#' @param type a character value containing the type of the method to apply. Possible values are: "KNN", "RF", "Logistic.regression" (default = KNN)
#'
#' @return a S4 object of class 'CYTdata'
#'
#' @export
#'

performUpsampling <- function(CYTdata,
                              newCYTdata,
                              markers = NULL,
                              type = c("KNN", "RF", "Logistic.regression")){
  
  # Check args
  checkmate::qassert(markers, c(0,"S*"))
  type = match.arg(type)
  checkmate::qassert(type, "S1")
  
  if (is.null(markers)){ markers = CYTdata@markers }
  markers = as.vector(markers)
  markers_interest = intersect(markers, newCYTdata@markers)
  if (length(markers_interest)<length(markers)) {
    stop("Error in performUpsampling : Some of the markers to use for Upsampling are not present in newCYTdata (",
         paste0(markers[!markers %in% markers_interest], collapse = ", "),")\n")
  }
  cat("Common markers used for Upsampling :", paste0(markers_interest, collapse=", "))
  
  exprs = CYTdata@matrix.expression
  new.exprs = newCYTdata@matrix.expression[, newCYTdata@markers %in% CYTdata@markers]
  
  raw.exprs = CYTdata@raw.matrix.expression
  new.raw.exprs = newCYTdata@raw.matrix.expression[, newCYTdata@markers %in% CYTdata@markers]
  
  overall.exprs = rbind.data.frame(exprs, new.exprs)
  overall.raw.exprs = rbind.data.frame(raw.exprs, new.raw.exprs)
  overall.clusters = c(CYTdata@Clustering@clusters, rep(NA, nrow(new.exprs)))
  overall.samples = c(CYTdata@samples, newCYTdata@samples)
  
  indexes.duplicated = duplicated(overall.exprs)
  
  overall.exprs = overall.exprs[indexes.duplicated,]
  overall.raw.exprs = overall.raw.exprs[indexes.duplicated,]
  overall.clusters = overall.clusters[indexes.duplicated]
  overall.samples = overall.samples[indexes.duplicated]
  
  data.clusters = cbind.data.frame(overall.exprs, "clusters" = overall.clusters)
  training.dataset = subset(data.clusters, !is.na(clusters))
  test.dataset = subset(data.clusters, is.na(clusters))
  test.dataset$clusters = NULL
  
  switch(type,
         KNN = {
           centroids = plyr::ddply(training.dataset, "clusters", function(x) {
             x$clusters <- NULL
             centers <- apply(x, 2, stats::median, rm.na = TRUE)
             return(centers)
           })
           upsampled.clusters <- FNN::knnx.index(centroids, test.dataset, k = 1, algorithm = "kd_tree")
         },
         RF = { stop("Not coded yet") },
         Logistic.regression = { stop("Not coded yet") })
  
  
  overall.clusters = c(CYTdata@Clustering@clusters, upsampled.clusters)
  
  # Put into "Clustering" object
  cat("\n\nCreating new Clustering object :")
  cat("\n - Computing cell cluster count matrix...")
  cellcount = compute.cell.count(clusters = overall.clusters, samples = overall.samples)
  cat("\n - Computing cell cluster abundance matrix...")
  abundance = compute.cell.abundance(count.matrix = cellcount)
  
  Clustering.object <- methods::new("Clustering",
                                    clusters = overall.clusters,
                                    cellcount = cellcount,
                                    abundance = abundance,
                                    parameters = append(CYTdata@Clustering@parameters, list("type.upsampling" = type,
                                                                                            "markers.upsampling" = markers)))
  cat("Creating new CYTdata object.. \n")
  CYTdata <- methods::new("CYTdata",
                          matrix.expression = overall.exprs,
                          markers = colnames(exprs),
                          samples = samples,
                          files = files,
                          raw.matrix.expression = overall.raw.exprs,
                          Clustering = Clustering.object)
  CYTdata@Clustering = Clustering.object
  validObject(CYTdata)
  return(CYTdata)
  
  CYTdata@matrix.expression.r <- matrix.expression.r[, colnames(matrix.expression.r) %in% colnames(downsampled.exp)]
  CYTdata@samples <- c(CYTdata@samples, upsampled.samples)
  CYTdata@identify.clusters <- c(CYTdata@identify.clusters, knn)
  
  message("computing cell cluster count matrix...")
  CYTdata@matrix.cell.count <- computeCellCounts(proj = CYTdata@matrix.expression,
                                                 clusters = CYTdata@identify.clusters,
                                                 samples = CYTdata@samples)
  
  message("computing cell cluster abundance matrix...")
  CYTdata@matrix.abundance <- computeClusterAbundances(count = CYTdata@matrix.cell.count)
  
  validObject(CYTdata)
  return(CYTdata)
}
