#' @title Check that no cluster is divided between two metaclusters
#'
#' @param CYTdata a CYTdata object with clustering and metaclustering performed
#' @param output a character vector indicating whether to display a message or outright stop the function if a cluster is divided between two metaclusters.
#'
#' @export
#'

check_HierarchyMetaclusters <- function(CYTdata, output = c("message", "stop")) {
  
  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }
  output = match.arg(output)
  checkmate::qassert(output, "S1")
  
  if (length(CYTdata@Metaclustering@metaclusters)!=0) {
    if (length(CYTdata@Clustering@clusters)==0) {
      stop("Error in CYTdata@Metaclustering object : metaclusters identifiers vector is given but no clusters identifiers
           vector given (Clustering@clusters slot is empty).")
    }
    else {
      clustersbelong = unlist(lapply(as.list(sapply(levels(CYTdata@Clustering@clusters),
                                                    function(cl) { return(unique(CYTdata@Metaclustering@metaclusters[CYTdata@Clustering@clusters == cl])) })), length))
      if (any(clustersbelong>1)) {
        if (output=="message") {
          message("Warning in CYTdata@Metaclustering object : The following clusters ( ",
                  paste0(names(clustersbelong)[clustersbelong>1], collapse = ", "), " ) are split in different metaclusters.
                  However, all the cell within a cluster should be in the same metacluster. As long as the hierarchy between clusters and metaclusters is not respected,
                  function using both identifiers (such as heatmaps clusters-metaclusters, etc.) will not work anymore.
                  Please check the potential duplicated clusters in metaclusters membership list generated with 'getMetaclustersMembership' function")
        }
        else {
          stop("Error in CYTdata@Metaclustering object : The following clusters ( ",
               paste0(names(clustersbelong)[clustersbelong>1], collapse = ", "), " ) are split in different metaclusters.
               However, for some functions (such as heatmaps clusters-metaclusters, etc.), hierarchy between clusters and metaclusters must be repected and all
               the cell within a cluster must be in the same metacluster.
               Please check the potential duplicated clusters in metaclusters membership list generated with 'getMetaclustersMembership' function")
        }
      }
    }
  }
  else {
    message("No hierarchy, Metaclustering slot is empty.")
  }
}



#Similar function to check_HierarchyMetaclusters. Paul's original function seemed to do more, but was likely left unfinished.

hierarchyMetaclusters <- function(CYTdata, checkOnly = FALSE, formatChecking = c("message", "stop")) {
  
        stop("Function hierarchyMetaclusters is deprecated. Please use check_HierarchyMetaclusters instead.")
  
}

#' @title  Computes the percentage of clusters with sufficient number of cells
#'
#' @description This function aims to compute and show cell clusters/metaclusters having a number of associated cells greater than a specific threshold.
#' Cluster cell count is computed for all samples at once. The optional heatmap, however, will show every cluster/sample pair.
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' level = c("clusters", "metaclusters")
#' @param level a character value indicating whether to display clusters or metaclusters. Possible values are: "clusters" (default), "metaclusters".
#' @param population a character vector containing the identifiers of the clusters/metaclusters to study. Defaults to all of them. 
#' @param thSize a numeric value providing the minimum number of cells needed for a cluster to be considered a small cluster (default value = 50)
#' @param generateHeatmap logical. If TRUE, also generates a heatmap with the clusters passing the threshold for each sample and in total. Defaults to FALSE
#'
#' @return a list containing QC size features:
#' - cellcount: a cellcount dataframe
#' - vectorScores: a boolean vector indicating for each cluster if it passes the threshold. 
#' - thSize: the threshold used
#' - score: the percentage of clusters that pass the threshold
#' - heatmap: a ggplot of the QC heatmap. 
#'
#' @export
#'

computePopulation_SizeQC <- function(CYTdata,
                                     population = NULL,
                                     level = c("clusters", "metaclusters"),
                                     thSize = 50,
                                     generateHeatmap = FALSE){
  
  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }
  
  level = match.arg(level)
  checkmate::qassert(level, "S1")
  population = checkorderPopulation(CYTdata=CYTdata, population=population, level=level, order=TRUE, checkDuplicates=TRUE)
  
  checkmate::qassert(thSize, "N1")
  if (thSize<=0) { stop("Error : 'thSize'  argument must be positive integer.") }
  checkmate::qassert(generateHeatmap, "B1")
  
  if (level == "clusters"){ cellcount = CYTdata@Clustering@cellcount[population,] }
  else { cellcount = CYTdata@Metaclustering@cellcount[population,] }
  
  effectif = apply(cellcount, 1, sum)
  vectorScores = as.logical(effectif >= thSize)
  score = mean(vectorScores)*100
  
  res = list("cellcount" = cellcount,
             "vectorScores" = vectorScores,
             "thSize" = thSize,
             "score" = score)
  
  if (generateHeatmap){
    cat("\n\n Generating heatmap of QCsize ..")
    data = reshape2::melt(data.matrix(cellcount))
    colnames(data) = c("groupId", "samples", "count")
    data$passed = TRUE
    data$passed[as.logical(data$count < thSize)] = FALSE
    data$count = NULL
    
    dataScore = cbind(population, "total.cells", vectorScores)
    colnames(dataScore) = c("groupId", "samples", "passed")
    
    data = rbind(data, dataScore)
    title = paste("Quality control on", level, "size", sep = " ")
    subtitle = paste("Percentage of", level, "having a sufficient number of cells =",
                     format(round(score), nsmall = 2), "%", sep=" ")
    
    heatmap <- ggplot2::ggplot() +
      ggplot2::geom_tile(data = data,
                         ggplot2::aes_string(x = "samples", y = "groupId", fill = "passed"),
                         colour = "black") +
      ggplot2::scale_fill_manual(values = c("TRUE" = "green", "FALSE" = "red")) +
      ggplot2::geom_vline(xintercept = (length(unique(data$samples)) - 0.5), colour = "black", size = 2) +
      ggplot2::scale_x_discrete(expand = c(0, 0)) + ggplot2::scale_y_discrete(expand = c(0, 0)) +
      ggplot2::xlab("samples") + ggplot2::ylab("clusters") +
      ggplot2::ggtitle(bquote(atop(.(title), atop(italic(.(subtitle)), "")))) +
      ggplot2::theme_minimal() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                     plot.background = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
                     axis.line = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank(),
                     panel.border = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     panel.grid.major = ggplot2::element_blank(),
                     legend.title = ggplot2::element_blank(),
                     legend.key = ggplot2::element_blank())
    res$heatmap = heatmap
  }
  
  cat("\n\n - Percentage of ", level, " having at least ", thSize ," cells = ",
      format(round(res$score, 2), nsmall = 2), "%\n\n")
  
  return(res)
}


#' @title  Computes the percentage of clusters with uniform phenotypes
#'
#' @description This function aims to identify and show cell clusters having a uniform phenotype, defined as an unimodal expression and a narrow range of expression for a selected set of markers
#'
#'  @details
#'
#' - Check unimodal expression: testing unimodal distribution of markers with Hartigan's dip test.
#' - Check the range of expression: check if the interquartile range of the specified markers is below a specified threshold. 
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param markers a character vector indicating the markers to use in the QC. 
#' By default, all of them are used, but it is recommended to use only the markers used for  clustering. 
#' @param population a character vector containing the identifiers of the clusters/metaclusters to compute the QC on. Defaults to all of them. 
#' @param level a character value indicating whether to compute the QC on clusters or metaclusters. Possible values are: "clusters" (default), "metaclusters".
#' @param thPvalue a numeric value providing the p-value threshold for the test (unimodal if pvalue > thPvalue). Defaults to 0.05
#' @param thIQR a numeric value providing the IQR (interquartile range) threshold to assume a distribution as narrow enough. Default to 2
#' @param thGood a numeric value providing the proportion of markers that need to pass the QC in order to consider that a cluster passes the QC
#' @param NbSim an integer. Number of replicates for the Monte Carlo simulation. Defaults to 2000.
#' @param generateHeatmap logical. If TRUE (default), a heatmap is generated with the QC. 
#' @param sortQCbyMarkers logical. If TRUE (default), markers will be sorted on the heatmap by the number of clusters that pass the QC for that marker.  
#' @param generateEcdf logical. If TRUE, an empirical cumulative distribution is generated with the QC, useful for comparing thGood values. CURRENTLY NON FUNCTIONNAL.  
#' @param parallelize logical. If TRUE, calculation of the test is done in parallel on several cores. If you get an "impossible to find function "%dopar%" error, try library(foreach)
#' @param nbCores an integer. Number of processor cores to use. If NULL, uses one less cores than available to leave computational power for using the computer. 
#' 
#' @return a list containing QC uniform features:
#' - Pvalue: a dataframe containing the pvalues for the dip test
#' - IQR: a dataframe containing interquartile range
#' - Uniform: a dataframe indicating which cluster/marker pair passes the QC for both uniformity and IQR
#' - goodScores: a named vector giving the percentage of markers passing the QC for each cluster
#' - Score: the overall percentage of clusters passing the QC
#' - thPvalue: the thPvalue parameter
#' - thIQR the thIQR parameter
#' - thGood: the thGood parameter
#' - heatmap: a ggplot of the QC heatmap. 
#' - ecdf: a ggplot of the empirical cumulative distribution
#'
#' @export
#'

computeQCuniform <- function(CYTdata,
                             markers = NULL,
                             population = NULL,
                             level = c("clusters", "metaclusters"),
                             thPvalue = 0.05,
                             thIQR = 2,
                             thGood = 0.75,
                             NbSim = 2000,
                             generateHeatmap = TRUE,
                             sortQCbyMarkers = TRUE,
                             generateEcdf = FALSE, # empirical cumulative distribution
                             parallelize = FALSE,
                             nbCores = NULL) {
  
  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }
  
  checkmate::qassert(thPvalue, "N1[0,1]")
  checkmate::qassert(thIQR, "N1(0,)")
  checkmate::qassert(thGood, "N1[0,1]")
  checkmate::qassert(NbSim, c("N1(0,)"))
  checkmate::qassert(generateHeatmap, "B1")
  checkmate::qassert(sortQCbyMarkers, "B1")
  checkmate::qassert(generateEcdf, "B1")
  checkmate::qassert(parallelize, "B1")
  checkmate::qassert(nbCores, c("0","N1(0,)"))
  

  level = match.arg(level)
  checkmate::qassert(level, "S1")
  population = checkorderPopulation(CYTdata, population=population, level=level,
                                    order=TRUE, checkDuplicates=TRUE)
  markers = checkorderMarkers(CYTdata, markers=markers, order=TRUE, checkDuplicates=TRUE)
  
  if (level == "clusters"){ popId = CYTdata@Clustering@clusters }
  else { popId = CYTdata@Metaclustering@metaclusters }
  
  getPvalue <- function(x) return(diptest::dip.test(x, B=NbSim)$p.value)
  
  if (length(population)==1) {
    cat("'population' argument has length equal to 1, representing a single population.
        \n Remark 1 : 'parallelize' argument is ignored here, no need to parallelize
        computation for a single population.
        \n Remark 2 : 'generateEcdf' argument is ignored as well.")
    data = CYTdata@matrix.expression[popId %in% population, markers]
    IQR = apply(data, 2,
                function(z){
                  quantiles = stats::quantile(z)
                  return(quantiles[4] - quantiles[2])})
    names(IQR) = markers
    pvalues = apply(data, 2, getPvalue)
    names(pvalues) = markers
    
    uniform = pvalues > thPvalue & IQR <= thIQR 
    names(uniform) = markers
    
    popScore = mean(uniform)*100
    score = ifelse(popScore>=thGood*100, 100, 0)
    
    res <- list("Pvalue" = IQR,
                "IQR" =  pvalues,
                "Uniform" = uniform,
                "goodScores" = popScore,
                "Score" = score,
                "thPvalue" = thPvalue,
                "thIQR" = thIQR,
                "thGood" = thGood)
    
    if (generateHeatmap){
      data = data.frame("passed" = as.vector(uniform),
                        "markers" = names(uniform))
      res$heatmap <- ggplot2::ggplot() +
        ggplot2::geom_tile(ggplot2::aes_string(x = "markers",
                                               y = 1,
                                               fill = "passed"),
                           colour = "black", size = 0.25) +
        ggplot2::ggtitle(label = paste("Uniform", level, "quality control", sep=" "),
                         subtitle = paste("Percentage of markers passing the QC (for ",
                                          level, " ", population, " ) : ", format(round(score, 2), nsmall = 2), "%")) +
        ggplot2::scale_fill_manual(values = c("TRUE" = "green", "FALSE" = "red")) +
        ggplot2::scale_x_discrete(expand = c(0, 0)) +
        ggplot2::scale_y_discrete(expand = c(0, 0)) +
        ggplot2::ylab(level) + ggplot2::xlab("Markers") +
        ggplot2::theme_minimal() +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                       plot.background  = ggplot2::element_blank(),
                       axis.text.x      = ggplot2::element_text(angle = 90, hjust = 1, vjust = 1),
                       axis.line        = ggplot2::element_blank(),
                       panel.background = ggplot2::element_blank(),
                       panel.border     = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank(),
                       panel.grid.major = ggplot2::element_blank())
    }
  }
  else {
    
    data = cbind.data.frame("popId" = popId, CYTdata@matrix.expression[,markers])
    
    IQR = plyr::ddply(subset(data, popId %in% population),
                      "popId",
                      function(y){
                        apply(y[,-1], 2, function(z){
                          quantiles = stats::quantile(z)
                          return(quantiles[4] - quantiles[2])
                        })
                      })
    rownames(IQR) = population
    IQR$popId = NULL
    print(IQR)
    
    if(parallelize){
      cat("'parallelize' argument set to TRUE, creating core clusters..")
      maxCores = parallel::detectCores() - 1
      if (is.null(nbCores)) nbCores = maxCores
      if (nbCores > maxCores) stop("Error : nbCores argument require more cores than available.")
      
      message("\n- Parallelization activated : Register the parallel backend using ",
              nbCores, " cores...")
      cl = parallel::makeCluster(nbCores)
      doParallel::registerDoParallel(cl)
      
      Pvalue = foreach(i = 1:length(population), .combine = 'rbind.data.frame') %dopar% {
        pop = population[i]
        subdata = data[data$popId == pop, -1]
        pvalues = apply(subdata, 2, FUN = getPvalue)
      }
    }
    else {
      Pvalue = data.frame()
      for (pop in population) {
        subdata = data[data$popId == pop, -1]
        pvalues = apply(subdata, 2, FUN = getPvalue)
        Pvalue = rbind.data.frame(Pvalue, pvalues)
      }
    }
    
    rownames(Pvalue) = population
    colnames(Pvalue) = colnames(data)[-1]
    print(Pvalue)
    Uniform = data.frame(Pvalue > thPvalue & IQR <= thIQR, 
                         check.names = FALSE)
    print(Uniform)
    rownames(Uniform) = population
    print(Uniform)
    
    popScores = apply(Uniform, 1, mean)*100
    overallPassed = as.logical(as.vector(popScores) >= thGood*100)
    score = mean(overallPassed)*100
    res <- list("Pvalue" = Pvalue,
                "IQR" =  IQR,
                "Uniform" = Uniform,
                "goodScores" = popScores,
                "Score" = score,
                "thPvalue" = thPvalue,
                "thIQR" = thIQR,
                "thGood" = thGood)
    
    thGood100 = thGood*100 #percentage
    
    if (generateHeatmap) {
      
      cat("\n\n Generating heatmap of QCuniform ..")
      
      data = cbind.data.frame(Uniform, "overallPassed" = overallPassed)
      
      if (sortQCbyMarkers) {
        markersQC = apply(data[,-ncol(data)], 2, mean)*100
        lev = c(names(markersQC)[order(markersQC, decreasing=TRUE)], "overallPassed")
      }
      else {
        lev = colnames(Uniform)
      }
      
      data = reshape2::melt(data.matrix(data))
      colnames(data) = c(level, "markers", "passed") # row must be cluster
      data$passed = as.logical(data$passed)
      data$markers = factor(data$markers, levels = lev)
      
      res$heatmap = ggplot2::ggplot(data = data) +
        ggplot2::geom_tile(ggplot2::aes_string(x = "markers",
                                               y = level,
                                               fill = "passed"),
                           colour = "black", size = 0.25) +
        ggplot2::ggtitle(label = paste("Uniform", level, "quality control", sep=" "),
                         subtitle = paste("Percentage of", level, "being considered as 'good'",
                                          level, ": ", format(round(score, 2), nsmall = 2), "%")) +
        ggplot2::scale_fill_manual(values = c("TRUE" = "green", "FALSE" = "red")) +
        ggplot2::scale_x_discrete(expand = c(0, 0)) +
        ggplot2::scale_y_discrete(expand = c(0, 0)) +
        ggplot2::ylab(level) + ggplot2::xlab("Markers") +
        ggplot2::theme_minimal() +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                       plot.background  = ggplot2::element_blank(),
                       axis.text.x      = ggplot2::element_text(angle = 90, hjust = 1, vjust = 1),
                       axis.line        = ggplot2::element_blank(),
                       panel.background = ggplot2::element_blank(),
                       panel.border     = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank(),
                       panel.grid.major = ggplot2::element_blank())
    }
    
    # if (generateEcdf) {
    #   
    #   cat("\n\n Generating empirical cumulative distribution..")
    #   
    #   res$ecdf = ggplot2::ggplot(data.frame("x" = popScores), #C'EST PAS CA QU'IL FAUT
    #                              ggplot2::aes(x)) +
    #     ggplot2::stat_ecdf(geom = "step", pad = FALSE) +
    #     ggplot2::geom_vline(xintercept = thGood100, colour = "red", size = 1, linetype = "longdash") +
    #     ggplot2::geom_text(ggplot2::aes(x = thGood100,
    #                                     y = 1,
    #                                     label = paste(thGood100, "%", sep = "")),
    #                        hjust = 1, colour = "red", size = 5) +
    #     ggplot2::geom_hline(yintercept = score, colour = "red", size = 1) +
    #     ggplot2::geom_text(ggplot2::aes(y = round(score, 2),
    #                                     x = 1,
    #                                     label = score),
    #                        vjust = -1, colour = "red", size = 5) +
    #     ggplot2::xlab(paste("Good", level, "threshold value", sep=" ")) +
    #     ggplot2::ylab(paste("Percentage of good", level, sep=" ")) +
    #     ggplot2::ggtitle(paste("Percentage of good", level, "over good",
    #                            level, "threshold value", sep=" ")) +
    #     ggplot2::theme_minimal() +
    #     ggplot2::theme(plot.title = ggplot2::element_text(size=15, hjust = 0.5),
    #                    axis.title = ggplot2::element_text(size = 18),
    #                    axis.text = ggplot2::element_text(size = 16))
    # }
  }
  
  cat("\n\n - Percentage of ", level, " having at least",
      thGood100 ,"% of marker expression distribution passing the QC =",
      format(round(res$Score, 2), nsmall = 2), "% \n\n")
  
  return(res)
}















