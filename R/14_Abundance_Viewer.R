#' @title Get the abundance of a specific cell population
#'
#' @description This function returns the relative or absolute abundance of specific cell populations (i.e., clusters or metaclusters)
#'
#' @param CYTdata a CYTdata object
#' @param population a character vector containing the identifiers of the population for which to return abundance or cell count. By default, all of them are used, which will return 100% if Yvalue is set to "percentage". 
#' @param level a character value indicating the type of population to study. Possible values are: "clusters" (default) and "metaclusters".
#' @param samples a character vector containing the names of biological samples to consider. By default, all samples are used
#' @param Yvalue a character value containing the parameters to use. Possible value are "percentage" (default) or "absolute" (i.e., cell count).
#'
#' @return a data frame with abundance per sample, along with the metadata of each sample
#'
#' @export

getRelativeAbundance <- function(CYTdata,
                                 population = NULL,
                                 level = c("clusters", "metaclusters"),
                                 samples = NULL,
                                 Yvalue = c("percentage", "absolute")){
  
  level = match.arg(level)
  checkmate::qassert(level, "S1")
  Yvalue = match.arg(Yvalue)
  checkmate::qassert(Yvalue, "S1")
  
  population = checkorderPopulation(CYTdata, population=population, level=level, order=TRUE, checkDuplicates=TRUE)
  samples = checkorderSamples(CYTdata, samples, order=TRUE, checkDuplicates=TRUE)
  
  if (level == "clusters"){ cellcount = CYTdata@Clustering@cellcount[, samples, drop=FALSE] }
  else { cellcount = CYTdata@Metaclustering@cellcount[, samples, drop=FALSE] }
  
  subcellcount =  colSums(cellcount[population,, drop=FALSE])
  
  if (Yvalue == "percentage") {
    cellcount = subcellcount / apply(cellcount, 2, sum) * 100
    Ylab = "Abundance of population"
  }
  else {
    cellcount = subcellcount
    Ylab = "Absolute count of population"
  }
  
  cellcount = reshape::melt(cellcount)
  data = merge(cellcount, CYTdata@metadata, by = "row.names")
  return(data)
}

#' @title Gets number of cells above an expression threshold for a given marker. 
#'
#' @description This function aims to compute cellcount or abundance across samples for a given population (one or several cluster or metacluster), counting only cells for which the expression of a given marker is above a given threshold (e.g., a positivity thressold for IdU+ cells. )
#' As of now, percentage results are computed over the total number of cells, not the total number of marker+ cells. 
#'
#' @param CYTdata a CYTdata object
#' @param marker a character vector containing the name of the marker to use. Must be included in CYTdata's expression matrix's colnames. 
#' @param th a numeric value giving the threshold to compare the marker expression to. If the data are transformed (e.g., an arcsinh), then the threshold must be transformed as well. 
#' @param population a character vector containing the identifiers of the clusters/metaclusters to study.
#' @param level a character value containing the levels to display in the plot, either "clusters" (default) or "metaclusters"
#' @param samples a character vector containing the names of the biological samples to use. By default, all samples are used.
#' @param Yvalue a character value containing the variable to be studied. Possible values are percentage (default) or absolute.
#'
#' @return a dataframe containing the number/percentage of cells per sample for the given pop, along with collated metadata
#'
#' @export
#'

getRelativeAbundanceThreshold = function(CYTdata, marker, th, population = NULL, level = c("clusters", "metaclusters"), 
                                         samples = NULL, Yvalue = c("percentage", "absolute")) 
{
  checkmate::qassert(marker, "S1")
  checkmate::qassert(th, "R1")
  level = match.arg(level)
  checkmate::qassert(level, "S1")
  Yvalue = match.arg(Yvalue)
  checkmate::qassert(Yvalue, "S1")
  
  population = checkorderPopulation(CYTdata, population = population, 
                                    level = level, order = TRUE, checkDuplicates = TRUE)
  samples = checkorderSamples(CYTdata, samples, order = TRUE, 
                              checkDuplicates = TRUE)
  
  if (level == "clusters") {
    cellcount = CYTdata@Clustering@cellcount[, samples, drop = FALSE] #Paul here counts everything (cellcount), then keeps only the population(s) of interest (subcellcount). If we gave the function population = levels(CYTdata@Metaclustering(metaclusters)), it will add up all metaclusters, and thus return 100 (if Yvalue = percentage)
    popCellID = CYT2304DS_MonoDC_lineplot@Clustering@clusters
  }
  else {
    cellcount = CYTdata@Metaclustering@cellcount[, samples, drop = FALSE]
    popCellID = CYT2304DS_MonoDC_lineplot@Metaclustering@metaclusters
  }
  
  cellcountMarkerPos = cellcount #Copying and emptying cellcount
  cellcountMarkerPos[,] = NA
  
  for (pop in rownames(cellcountMarkerPos)) { #Same logic as Paul. Computing the number of marker+ cells, then considering only the cells belonging to the population(s) of interest as defined in the function's arguments.  
    
    # cat(pop, "\n")
    for (samp in colnames(cellcountMarkerPos)) {
      
      cellcountMarkerPos[pop, samp] = (CYT2304DS_MonoDC_lineplot@matrix.expression[popCellID == pop & CYT2304DS_MonoDC_lineplot@samples == samp, marker] > th) %>%
        sum() #Select expression of the specified marker for the sample and the pop. Sum all above threshold. 
      
    }
  }
  
  subcellcount = colSums(cellcountMarkerPos[population, , drop = FALSE]) 
  if (Yvalue == "percentage") {
    cellcountMarkerPos = subcellcount/apply(cellcount, 2, sum) * 100 #Caution here. Using apply(cellcount, 2, sum) results in computing the number of cells positive to the specified marker divided by the entire population of cells for that sample. Using apply(cellcountMarkerPos, 2, sum) instead will give the proportion of marker+ cells of the interest population over the total number of marker+ cells. 
    Ylab = "Abundance of population"
  }
  else {
    cellcountMarkerPos = subcellcount
    Ylab = "Absolute count of population"
  }
  
  cellcountMarkerPos = reshape::melt(cellcountMarkerPos)
  data = merge(cellcountMarkerPos, CYTdata@metadata, by = "row.names")
  return(data)
  
}

#' @title Stats function for lineplots (internal)
#'
#' @description Computes pairwise abundance comparison between timepoints and a reference timepoint (usually the baseline) 
#' Currently (2026-03-02) bugged if less than three pairs to compute. 
#'
#' @param dataMelted A melted dataframe as used in plotLineplot3 
#' @param statPairs A list of character pairs indicating the timepoint pairs on which to compute tests
#' @param individuals A character vector containing the names of the individuals to consider.  
#' @param test.statistics A character string indicating the type of statistics to compute : "wilcox.test" (default) or "t.test"
#' @param paired Logical. If TRUE, performs a pairwise test.
#' @param corrMethod A character string, the correction method for the p.value. Possible values are "holm" (default), "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none", see p.adjust for details. 
#'
#' @return a dataframe containing the number/percentage of cells per sample for the given pop, along with collated metadata
#'
#' @export
#'


statsLineplot <- function(dataMelted,
                          statPairs,
                          individuals,
                          test.statistics = c("wilcox.test", "t.test"),
                          paired = FALSE, 
                          corrMethod = "holm"){ 
  
  checkmate::qassert(dataMelted, "D1")
  checkmate::qassert(statPairs, "L*")
  checkmate::qassert(individuals, "S*")
  checkmate::qassert(test.statistics, "S1")
  test.statistics = match.arg(test.statistics)
  checkmate::qassert(paired, "B1")
  checkmate::qassert(corrMethod, "S1")
  
  pDataList = list()
  
  for (pop in unique(dataMelted$L1)){
    
    if (test.statistics == "wilcox.test") {
      
      if (paired) {
        cat(paste("\nComputing paired wilcoxon test on the following pairs of timepoints:\n", paste(statPairs, collapse = ", "), "\n", sep = ""))
      } else {
        cat(paste("\nComputing unpaired wilcoxon test on the following pairs of timepoints:\n", paste(statPairs, collapse = ", "), "\n", sep = ""))
      }
      
      statsStatix = rstatix::wilcox_test(dataMelted[dataMelted$L1 == pop & dataMelted$Var2 %in% individuals,], 
                                         value~Var1, 
                                         comparisons = statPairs,
                                         paired = paired, 
                                         p.adjust.method = corrMethod)
    } else {
      
      if (paired) {
        cat(paste("\nComputing paired t-test on the following pairs of timepoints:\n", paste(statPairs, collapse = ", "), "\n", sep = ""))
      } else {
        cat(paste("\nComputing unpaired t-test on the following pairs of timepoints:\n", paste(statPairs, collapse = ", "), "\n", sep = ""))
      }
      
      statsStatix = rstatix::t_test(dataMelted[dataMelted$L1 == pop & dataMelted$Var2 %in% individuals,], 
                                    value~Var1, 
                                    comparisons = statPairs,
                                    paired = paired, 
                                    p.adjust.method = corrMethod)
    }
    
    statsStatix = tibble::column_to_rownames(statsStatix, "group2")
    pAdj = statsStatix$p.adj; names(pAdj) = row.names(statsStatix)
    pDataList[[pop]] = pAdj
    
  }
  
  #Checking for mistakes in the list names
  for (elem in pDataList) {if (!all(names(elem) == names(pDataList[[1]]))) {stop("ERROR: non indentical names in list of p values")}}
  
  #Formatting
  pDataF = do.call(rbind.data.frame, pDataList) %>% t()
  rownames(pDataF) = names(pDataList[[1]]); colnames(pDataF) = names(pDataList)
  
  return(pDataF)
  
}


#Author: Cyril ETIENNE
#Date : 2026-01-15

#' @title Plots cell population abundances across timepoints
#'
#' @description This function aims to visualize and compare the cell population abundances variation across timepoints for each biological condition using lineplot representations.
#'
#' Several sets of clusters or metaclusters can be displayed, but the plot usually gets very cluttered above ~6, unless displayType = "line".
#' The representation can be restricted to a specific set of samples.
#' Statistics can be computed for all comparisons.
#'
#' @param CYTdata a CYTdata object
#' @param population a character vector containing the identifiers of the clusters/metaclusters to display in the plot. By default, all clusters or metaclusters are used.
#' @param level a character value containing the levels to display in the plot, either "clusters" (default) or "metaclusters"
#' @param samples a character vector containing the names of the biological samples to use. By default, all samples are used.
#' @param title a character value containing the title of the plot
#' @param Yvalue a character value containing the variable to be studied. Possible values are percentage (default) or absolute.
#' @param stats a boolean value indicating whether statistical comparisons of each timepoint to the baseline should be performed.
#' @param test.statistics a character value providing the type of statistical test to use. Possible values are: 'wilcox.test' (default) or 't.test'
#' @param paired a boolean value indicating if a paired or unpaired comparison should be applied. Defaults to FALSE.
#' @param corrMethod a character value providing the type of correction for the p.value. Please refer to rstatix::wilcox_test or rstatix::t_test for details. Defaults to "holm". 
#' @param BSL a character value containing the name of the baseline sample against which the test will be performed for each timepoint. Defaults to "BSL2" (change ?)
#' @param displayType a character value containing the type of line plot. "errorbars" (default) plots error bars, "ribbon" plots a ribbon, "line" plots only the average. 
#' @param displayError a character value containing the boundaries to be used for displaying the error bars or the ribbon. 
#' "sd" (default) for +/- 1 standard deviation, "minmax" for the minimum and maximum. 
#' @param meanMed a character value indicating whether to display the mean ('mean', default) or the median ('median') in the lineplot. 
#' @param palette a character vector colors for the plots. if NULL (default), the CYTdata's colour palette will be used. 
#'
#' @return a ggplot2 object
#'
#' @export
#'

plotLineplot <- function(CYTdata,
                         population = NULL, #MC/subset name
                         level = c("clusters", "metaclusters"),
                         samples = NULL, 
                         title = NULL, #Title of the graph
                         Yvalue = c("percentage", "absolute"),
                         stats = TRUE,
                         test.statistics = c("wilcox.test", "t.test"),
                         paired = FALSE, #Whether to compute a paired test
                         corrMethod = "holm", #p value correction method
                         BSL = "BSL2", #Name of the baseline timepoint for wilcox tests
                         displayType = c("errorbars", "ribbon", "line"),
                         displayError = c("sd", "minmax"),
                         meanMed = c("mean", "median"), #display mean or median ? 
                         palette = NULL) { 
  
  #CHECKING VALIDITY OF ARGUMENTS
  
  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }
  
  if (ncol(CYTdata@metadata)==0) {
    stop("Error : Missing metadata slot for 'CYTdata' argument.")
  }
  
  level = match.arg(level) #checks if level matches what is defined in the function's definition
  checkmate::qassert(level, "S1") 
  
  Yvalue = match.arg(Yvalue)
  checkmate::qassert(Yvalue, "S1")
  
  test.statistics = match.arg(test.statistics)
  checkmate::qassert(test.statistics, "S1")
  
  displayType = match.arg(displayType)
  checkmate::qassert(displayType, "S1")
  
  displayError = match.arg(displayError)
  checkmate::qassert(displayError, "S1")
  
  meanMed = match.arg(meanMed)
  checkmate::qassert(meanMed, "S1")
  
  checkmate::qassert(BSL, "S1")
  checkmate::qassert(population, c("S*", "0"))
  checkmate::qassert(samples, c("S*", "0"))
  checkmate::qassert(stats, "B")
  checkmate::qassert(paired, "B")
  checkmate::qassert(corrMethod, "S1")
  checkmate::qassert(palette, c("S*", "0"))
  
  if (!is.null(title)) {checkmate::qassert(title, "S1")}
  
  #ASSIGNING DEFAULT VALUES
  
  if (is.null(title)) {title = paste("Population", Yvalue ,"evolution")}
  
  if (is.null(population) & level == "clusters") {population = CYTdata@Clustering@clusters %>% levels()}
  if (is.null(population) & level == "metaclusters") {population = CYTdata@Metaclustering@metaclusters %>% levels()}
  
  #EXTRACTING DATA
  
  popRelAbundance = NULL 
  
  for (pop in population){ 
    
    popRelAbundance[[pop]] = Bahamas::getRelativeAbundance(CYTdata,
                                                           population = pop,
                                                           samples = samples,
                                                           level = level, 
                                                           Yvalue = Yvalue) %>%
      reshape2::melt() %>% 
      reshape2::dcast(Individual~Timepoint) %>%
      tibble::column_to_rownames(var = "Individual")
    
    mean = colMeans(popRelAbundance[[pop]])
    max = lapply(popRelAbundance[[pop]], max, na.rm = TRUE)
    min = lapply(popRelAbundance[[pop]], min, na.rm = TRUE)
    sd   = popRelAbundance[[pop]] %>% as.matrix() %>% matrixStats::colSds()
    med   = popRelAbundance[[pop]] %>% as.matrix() %>% matrixStats::colMedians()
    
    popRelAbundance[[pop]]["mean",] = mean
    popRelAbundance[[pop]]["med",] = med
    popRelAbundance[[pop]]["sd",] = sd
    popRelAbundance[[pop]]["max",] = max
    popRelAbundance[[pop]]["min",] = min
    
    popRelAbundance[[pop]] = popRelAbundance[[pop]] %>% t()
    
  }
  
  
  #Creating dataframe for ggplot2
  
  dataMelted = reshape2::melt(popRelAbundance)
  
  data = purrr::reduce(.x = list(dataMelted[dataMelted[,"Var2"] == "mean",],
                                 dataMelted[dataMelted[,"Var2"] == "med",],
                                 dataMelted[dataMelted[,"Var2"] == "sd",],
                                 dataMelted[dataMelted[,"Var2"] == "min",],
                                 dataMelted[dataMelted[,"Var2"] == "max",]),
                       .f = merge,
                       by = c("Var1","L1"))
  
  data = data[,c(1,2,4,6,8,10,12)]
  colnames(data) = c("Timepoint","Population","mean","med","sd", "min", "max")
  
  if (meanMed == "median") {data$`meanMed` = data$med
  } else {data$`meanMed` = data$mean} #choosing between mean and median
  
  #COMPUTING STATISTICS
  
  if (stats) {  
    
    pairLevels = CYTdata@metadata$Timepoint %>% droplevels() %>% levels() %>% setdiff(BSL)
    statPairs = paste(BSL, pairLevels) %>% strsplit(" ")
    
    pDataF = statsLineplot(dataMelted = dataMelted,
                           statPairs = statPairs,
                           individuals = CYTdata@metadata$Individual %>% levels(),
                           test.statistics = test.statistics,
                           paired = paired,
                           corrMethod = corrMethod)
    
    #Defining symbols for dispalying the p.values on the graph
    
    pDataF = reshape2::melt(pDataF %>% as.matrix())
    colnames(pDataF) = c("Timepoint", "Population", "p")
    
    pDataF = pDataF[pDataF$p <= 0.1,]
    pDataF$`Symbol` = "***"
    pDataF$`Symbol`[pDataF$p > 0.001] = "**"
    pDataF$`Symbol`[pDataF$p > 0.01] = "*"
    pDataF$`Symbol`[pDataF$p > 0.05] = "°"
    
    #Defining the y position where the p.value's symbol will be displayed. 
    
    yPos = NULL
    
    for (i in 1:nrow(pDataF)){
      
      if (displayType == "errorbars" & displayError == "sd") {
        
        yPos[i] = data$meanMed[data$Population == pDataF$Population[i] & as.vector(data$Timepoint) == pDataF$Timepoint[i]] + 
          data$sd[data$Population == pDataF$Population[i] & as.vector(data$Timepoint) == pDataF$Timepoint[i]] #p displayed above meanMed + sd
        
      } else if (displayType == "errorbars" & displayError == "minmax") {
        
        yPos[i] = data$max[data$Population == pDataF$Population[i] & as.vector(data$Timepoint) == pDataF$Timepoint[i]] #p displayed above max
        
      } else if (displayType == "ribbon" | displayType == "line") {
        
        yPos[i] = data$meanMed[data$Population == pDataF$Population[i] & as.vector(data$Timepoint) == pDataF$Timepoint[i]]
        
      }
      
    }
    
    pDataF$`yPos` = yPos
    
  }  
  
  #CREATING THE GRAPH
  
  Lplot <- ggplot2::ggplot(data,
                           ggplot2::aes(y = meanMed,
                                        x = Timepoint,
                                        group = Population, color = Population)) +
    ggplot2::geom_line(linewidth=1) + 
    ggplot2::geom_point(size = 2)
  
  #Error bars or ribbon
  
  if (displayType == "errorbars" & displayError == "sd") {
    
    Lplot <- Lplot + 
      ggplot2::geom_errorbar(ggplot2::aes(ymin = meanMed-sd, ymax = meanMed+sd), 
                             width = .1)
    
  } else if (displayType == "errorbars" & displayError == "minmax") {
    
    Lplot <- Lplot + 
      ggplot2::geom_errorbar(ggplot2::aes(ymin = min, ymax = max), 
                             width = .1)
    
  } else if (displayType == "ribbon" & displayError == "sd") {
    
    Lplot <- Lplot + 
      ggplot2::geom_ribbon(ggplot2::aes(ymin=meanMed-sd, ymax=meanMed+sd), 
                           linetype=2, 
                           alpha=0.1) 
    
  } else if (displayType == "ribbon" & displayError == "minmax") {
    
    Lplot <- Lplot + 
      ggplot2::geom_ribbon(ggplot2::aes(ymin=min, ymax=max), 
                           linetype=2, 
                           alpha=0.1)  
    
  } 
  
  #Titles
  
  if (Yvalue == "percentage") {
    
    Lplot <- Lplot + ggplot2::labs(title = title,
                                   x = "Timepoint",
                                   y = "Relative abundance (%)")
    
  } else {
    
    Lplot <- Lplot + ggplot2::labs(title = title,
                                   x = "Timepoint",
                                   y = "Cell count")
    
  }
  
  #Statistical annotation
  
  if (stats) {  
    
    if (test.statistics == "wilcox.test" & paired) {
      
      annotateText = "Wilcoxon (paired)" #"Wilcoxon signed rank test (paired)"
      
    } else if (test.statistics == "wilcox.test" & !paired) {
      
      annotateText = "Wilcoxon (unpaired)" #"Wilcoxon rank sum test (unpaired)"
      
    } else if (test.statistics == "t.test" & paired) {
      
      annotateText = "t-Test (paired)" #"Student's t-Test (paired)"
      
    } else if (test.statistics == "t.test" & !paired) {
      
      annotateText = "t-Test (unpaired)" #"Student's t-Test (unpaired)"
      
    }
    
    
    Lplot <- Lplot + 
      ggplot2::geom_text(mapping = ggplot2::aes(x = Timepoint, y = yPos),
                         data = pDataF,
                         label = pDataF$Symbol,
                         size = 6,
                         position = ggplot2::position_nudge(x = 0.1, y = 0),#y = max(data$meanMed)/20),
                         hjust = "left",
                         show.legend = F) +
      ggplot2::annotate("text", x = Inf, y = Inf,
                        label = paste(annotateText, "\nTimepoint to ", BSL, ":\n*** p < 0.001\n** p < 0.01\n* p < 0.05\n° p < 0.1\nCorrection: ", corrMethod, sep = ""),
                        hjust = "left", vjust = "top")
    
  }
  
  #Theme
  
  if (is.null(palette)) {
    if (level == "clusters") {
      Lplot = Lplot + 
        ggplot2::scale_color_manual(values = CYTdata@Clustering@palette)
    } else {
      Lplot = Lplot + 
        ggplot2::scale_color_manual(values = CYTdata@Metaclustering@palette)
    }
  } else {
    Lplot = Lplot + 
      # ggplot2::scale_colour_brewer(palette = palette) #If a manual palette was set : override CYTdata palette
      ggplot2::scale_color_manual(values = palette) #If a manual palette was set : override CYTdata palette
  }
  
  Lplot <- Lplot + 
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme(legend.justification  = "bottom", legend.margin = ggplot2::margin(0.2,1,0.5,0.2,"cm")) #Top, right, bottom, left. Right at 1 cm to increase size and avoid annotate being cut off. 
  
  return(Lplot)
  
}


#Author: Cyril ETIENNE
#Date : 2026-02-05

#' @title Plots cell population abundances across timepoints using a lineplot representation, for cells above a marker expression threshold.
#'
#' @description Similar to plotLineplot3, this function aims to visualize and compare the cell population abundances variation across timepoints for each biological condition using lineplot representations.
#' However, it only considers cells that significantly express a given marker, i.e., that are above a given threshold. 
#' 
#' Several sets of clusters or metaclusters can be displayed, but the plot usually gets very cluttered above ~6, unless displayType = "line".
#' The representation can be restricted to a specific set of samples.
#' Statistics can be computed for all comparisons.
#'
#' @param CYTdata a CYTdata object
#' @param marker a character vector containing the name of the marker to use. Must be included in CYTdata's expression matrix's colnames. 
#' @param th a numeric value giving the threshold to compare the marker expression to. If the data are transformed (e.g., an arcsinh), then the threshold must be transformed as well. 
#' @param population a character vector containing the identifiers of the clusters/metaclusters to display in the plot. By default, all clusters or metaclusters are used.
#' @param level a character value containing the levels to display in the plot, either "clusters" (default) or "metaclusters"
#' @param samples a character vector containing the names of the biological samples to use. By default, all samples are used.
#' @param title a character value containing the title of the plot
#' @param Yvalue a character value containing the variable to be studied. Possible values are percentage (default) or absolute.
#' @param stats a boolean value indicating whether statistical comparisons of each timepoint to the baseline should be performed.
#' @param test.statistics a character value providing the type of statistical test to use. Possible values are: 'wilcox.test' (default) or 't.test'
#' @param paired a boolean value indicating if a paired or unpaired comparison should be applied. Defaults to FALSE.
#' @param corrMethod a character value providing the type of correction for the p.value. Please refer to rstatix::wilcox_test or rstatix::t_test for details. 
#' @param BSL a character value containing the name of the baseline sample against which the test will be performed for each timepoint. Defaults to "BSL2" (change ?)
#' @param displayType a character value containing the type of line plot. "errorbars" (default) plots error bars, "ribbon" plots a ribbon, "line" plots only the average. 
#' @param displayError a character value containing the boundaries to be used for displaying the error bars or the ribbon. 
#' "sd" (default) for +/- 1 standard deviation, "minmax" for the minimum and maximum. 
#' @param meanMed a character value indicating whether to display the mean ('mean', default) or the median ('median') in the lineplot. 
#' @param palette a character vector colors for the plots. if NULL (default), the CYTdata's colour palette will be used. 
#'
#' @return a ggplot2 object
#'
#' @export
#'

plotLineplotMarkerPos <- function(CYTdata,
                                  marker,
                                  th,
                                  population = NULL, #MC/subset name
                                  level = c("clusters", "metaclusters"),
                                  samples = NULL, 
                                  title = NULL, #Title of the graph
                                  Yvalue = c("percentage", "absolute"),
                                  stats = TRUE,
                                  test.statistics = c("wilcox.test", "t.test"),
                                  paired = FALSE, #Whether to compute a paired test
                                  corrMethod = "holm", #p value correction method
                                  BSL = "BSL2", #Name of the baseline timepoint for wilcox tests
                                  displayType = c("errorbars", "ribbon", "line"),
                                  displayError = c("sd", "minmax"),
                                  meanMed = c("mean", "median"), #display mean or median ? 
                                  palette = NULL) { 
  
  #CHECKING VALIDITY OF ARGUMENTS
  
  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }
  
  if (ncol(CYTdata@metadata)==0) {
    stop("Error : Missing metadata slot for 'CYTdata' argument.")
  }
  
  checkmate::qassert(marker, "S1")
  checkmate::qassert(th, "R1")
  
  level = match.arg(level) #checks if level matches what is defined in the function's definition
  checkmate::qassert(level, "S1") 
  
  Yvalue = match.arg(Yvalue)
  checkmate::qassert(Yvalue, "S1")
  
  test.statistics = match.arg(test.statistics)
  checkmate::qassert(test.statistics, "S1")
  
  displayType = match.arg(displayType)
  checkmate::qassert(displayType, "S1")
  
  displayError = match.arg(displayError)
  checkmate::qassert(displayError, "S1")
  
  meanMed = match.arg(meanMed)
  checkmate::qassert(meanMed, "S1")
  
  checkmate::qassert(BSL, "S1")
  checkmate::qassert(population, c("S*", "0"))
  checkmate::qassert(samples, c("S*", "0"))
  checkmate::qassert(stats, "B")
  checkmate::qassert(paired, "B")
  checkmate::qassert(corrMethod, "S1")
  checkmate::qassert(palette, c("S*", "0"))
  
  if (!is.null(title)) {checkmate::qassert(title, "S1")}
  
  #ASSIGNING DEFAULT VALUES
  
  if (is.null(title)) {title = paste("Population", Yvalue ,"evolution")}
  
  if (is.null(population) & level == "clusters") {population = CYTdata@Clustering@clusters %>% levels()}
  if (is.null(population) & level == "metaclusters") {population = CYTdata@Metaclustering@metaclusters %>% levels()}
  
  #EXTRACTING DATA
  
  popRelAbundance = NULL 
  
  for (pop in population){ 
    
    popRelAbundance[[pop]] = getRelativeAbundanceThreshold(CYTdata,
                                                           marker = marker,
                                                           th = th,
                                                           population = pop,
                                                           samples = samples,
                                                           level = level, 
                                                           Yvalue = Yvalue) %>%
      reshape2::melt() %>% 
      reshape2::dcast(Individual~Timepoint) %>%
      tibble::column_to_rownames(var = "Individual")
    
    mean = colMeans(popRelAbundance[[pop]])
    max = lapply(popRelAbundance[[pop]], max, na.rm = TRUE)
    min = lapply(popRelAbundance[[pop]], min, na.rm = TRUE)
    sd   = popRelAbundance[[pop]] %>% as.matrix() %>% matrixStats::colSds()
    med   = popRelAbundance[[pop]] %>% as.matrix() %>% matrixStats::colMedians()
    
    popRelAbundance[[pop]]["mean",] = mean
    popRelAbundance[[pop]]["med",] = med
    popRelAbundance[[pop]]["sd",] = sd
    popRelAbundance[[pop]]["max",] = max
    popRelAbundance[[pop]]["min",] = min
    
    popRelAbundance[[pop]] = popRelAbundance[[pop]] %>% t()
    
  }
  
  
  #Creating dataframe for ggplot2
  
  dataMelted = reshape2::melt(popRelAbundance)
  
  data = purrr::reduce(.x = list(dataMelted[dataMelted[,"Var2"] == "mean",],
                                 dataMelted[dataMelted[,"Var2"] == "med",],
                                 dataMelted[dataMelted[,"Var2"] == "sd",],
                                 dataMelted[dataMelted[,"Var2"] == "min",],
                                 dataMelted[dataMelted[,"Var2"] == "max",]),
                       .f = merge,
                       by = c("Var1","L1"))
  
  data = data[,c(1,2,4,6,8,10,12)]
  colnames(data) = c("Timepoint","Population","mean","med","sd", "min", "max")
  
  if (meanMed == "median") {data$`meanMed` = data$med
  } else {data$`meanMed` = data$mean} #choosing between mean and median
  
  #COMPUTING STATISTICS
  
  if (stats) {  
    
    #Creating data frame to store the p-value of each test
    
    pairLevels = CYTdata@metadata$Timepoint %>% droplevels() %>% levels() %>% setdiff(BSL)
    statPairs = paste(BSL, pairLevels) %>% strsplit(" ")
    
    pDataF = statsLineplot(dataMelted = dataMelted,
                           statPairs = statPairs,
                           individuals = CYTdata@metadata$Individual %>% levels(),
                           test.statistics = test.statistics,
                           paired = paired,
                           corrMethod = corrMethod)
    
    #Defining symbols for displaying the p.values on the graph
    
    pDataF = reshape2::melt(pDataF %>% as.matrix())
    colnames(pDataF) = c("Timepoint", "Population", "p")
    
    pDataF = pDataF[pDataF$p <= 0.1,]
    pDataF$`Symbol` = "***"
    pDataF$`Symbol`[pDataF$p > 0.001] = "**"
    pDataF$`Symbol`[pDataF$p > 0.01] = "*"
    pDataF$`Symbol`[pDataF$p > 0.05] = "°"
    
    #Defining the y position where the p.value's symbol will be displayed. 
    
    yPos = NULL
    
    for (i in 1:nrow(pDataF)){
      
      if (displayType == "errorbars" & displayError == "sd") {
        
        yPos[i] = data$meanMed[data$Population == pDataF$Population[i] & as.vector(data$Timepoint) == pDataF$Timepoint[i]] + 
          data$sd[data$Population == pDataF$Population[i] & as.vector(data$Timepoint) == pDataF$Timepoint[i]] #p displayed above mean + sd
        
      } else if (displayType == "errorbars" & displayError == "minmax") {
        
        yPos[i] = data$max[data$Population == pDataF$Population[i] & as.vector(data$Timepoint) == pDataF$Timepoint[i]] #p displayed above mean + max
        
      } else if (displayType == "ribbon" | displayType == "line") {
        
        yPos[i] = data$meanMed[data$Population == pDataF$Population[i] & as.vector(data$Timepoint) == pDataF$Timepoint[i]]
        
      }
      
    }
    
    pDataF$`yPos` = yPos
    
  }  
  
  #CREATING THE GRAPH
  
  #Base
  
  Lplot <- ggplot2::ggplot(data,
                           ggplot2::aes(y = meanMed,
                                        x = Timepoint,
                                        group = Population, color = Population)) +
    ggplot2::geom_line(linewidth=1) + 
    ggplot2::geom_point(size = 2)
  
  #Error bars or ribbon
  
  if (displayType == "errorbars" & displayError == "sd") {
    
    Lplot <- Lplot + 
      ggplot2::geom_errorbar(ggplot2::aes(ymin = meanMed-sd, ymax = meanMed+sd), 
                             width = .1)
    
  } else if (displayType == "errorbars" & displayError == "minmax") {
    
    Lplot <- Lplot + 
      ggplot2::geom_errorbar(ggplot2::aes(ymin = min, ymax = max), 
                             width = .1)
    
  } else if (displayType == "ribbon" & displayError == "sd") {
    
    Lplot <- Lplot + 
      ggplot2::geom_ribbon(ggplot2::aes(ymin=meanMed-sd, ymax=meanMed+sd), 
                           linetype=2, 
                           alpha=0.1) 
    
  } else if (displayType == "ribbon" & displayError == "minmax") {
    
    Lplot <- Lplot + 
      ggplot2::geom_ribbon(ggplot2::aes(ymin=min, ymax=max), 
                           linetype=2, 
                           alpha=0.1)  
    
  } 
  
  #Titles
  
  if (Yvalue == "percentage") {
    
    Lplot <- Lplot + ggplot2::labs(title = title,
                                   x = "Timepoint",
                                   y = "Relative abundance (%)")
    
  } else {
    
    Lplot <- Lplot + ggplot2::labs(title = title,
                                   x = "Timepoint",
                                   y = "Cell count")
    
  }
  
  #Statistical annotation
  
  if (stats) {  
    
    if (test.statistics == "wilcox.test" & paired) {
      
      annotateText = "Wilcoxon (paired)" #"Wilcoxon signed rank test (paired)"
      
    } else if (test.statistics == "wilcox.test" & !paired) {
      
      annotateText = "Wilcoxon (unpaired)" #"Wilcoxon rank sum test (unpaired)"
      
    } else if (test.statistics == "t.test" & paired) {
      
      annotateText = "t-Test (paired)" #"Student's t-Test (paired)"
      
    } else if (test.statistics == "t.test" & !paired) {
      
      annotateText = "t-Test (unpaired)" #"Student's t-Test (unpaired)"
      
    }
    
    
    Lplot <- Lplot + 
      ggplot2::geom_text(mapping = ggplot2::aes(x = Timepoint, y = yPos),
                         data = pDataF,
                         label = pDataF$Symbol,
                         size = 6,
                         position = ggplot2::position_nudge(x = 0.1, y = 0),#y = max(data$meanMed)/20),
                         hjust = "left",
                         show.legend = F) +
      ggplot2::annotate("text", x = Inf, y = Inf,
                        label = paste(annotateText, "\nTimepoint to ", BSL, ":\n*** p < 0.001\n** p < 0.01\n* p < 0.05\n° p < 0.1\nCorrection: ", corrMethod, sep = ""),
                        hjust = "left", vjust = "top")
    
  }
  
  #Theme
  
  if (is.null(palette)) {
    if (level == "clusters") {
      Lplot = Lplot + 
        ggplot2::scale_color_manual(values = CYTdata@Clustering@palette)
    } else {
      Lplot = Lplot + 
        ggplot2::scale_color_manual(values = CYTdata@Metaclustering@palette)
    }
  } else {
    Lplot = Lplot + 
      # ggplot2::scale_colour_brewer(palette = palette) #If a manual palette was set : override CYTdata palette
      ggplot2::scale_color_manual(values = palette) #If a manual palette was set : override CYTdata palette
  }
  
  Lplot <- Lplot + 
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme(legend.justification  = "bottom", legend.margin = ggplot2::margin(0.2,1,0.5,0.2,"cm")) #Top, right, bottom, left. Right at 1 cm to increase size and avoid annotate being cut off. 
  
  return(Lplot)
  
}

#' @title Lineplot of cell population abundances displaying each individual
#'
#' @description Similar to plotLineplot3, this function aims to visualize and compare the cell population abundances variation across timepoints for each biological condition using lineplot representations.
#' However, it does not display the mean and standard deviation, but the mean in bold and all other individuals with each an associated point shape. 
#' 
#' Several sets of clusters or metaclusters can theoretically be displayed, but it is recommended to use only one for plot readability.
#' The representation can be restricted to a specific set of samples.
#' Statistics can be computed for all comparisons.
#'
#' @param CYTdata a CYTdata object
#' @param population a character vector containing the identifiers of the clusters/metaclusters to display in the plot. By default, all clusters or metaclusters are used.
#' @param level a character value containing the levels to display in the plot, either "clusters" (default) or "metaclusters"
#' @param samples a character vector containing the names of the biological samples to use. By default, all samples are used.
#' @param title a character value containing the title of the plot.
#' @param Yvalue a character value containing the variable to be studied. Possible values are percentage (default) or absolute.
#' @param stats a boolean value indicating whether statistical comparisons of each timepoint to the baseline should be performed.
#' @param test.statistics a character value providing the type of statistical test to use. Possible values are: 'wilcox.test' (default) or 't.test'
#' @param paired a boolean value indicating if a paired or unpaired comparison should be applied. Defaults to FALSE.
#' @param corrMethod a character value providing the type of correction for the p.value. Please refer to rstatix::wilcox_test or rstatix::t_test for details. 
#' @param BSL a character value containing the name of the baseline sample against which the test will be performed for each timepoint. Defaults to "BSL2" (change ?)
#' @param meanMed a character value indicating whether to display the mean ('mean', default) or the median ('median') in the lineplot. 
#' @param palette a character vector colors for the plots. if NULL (default), the CYTdata's colour palette will be used. 
#'
#' @return a ggplot2 object
#'
#' @export
#'

plotLineplotInd <- function(CYTdata,
                            marker = NULL,
                            th = NULL,
                            population = NULL, #MC/subset name
                            level = c("clusters", "metaclusters"),
                            samples = NULL, 
                            title = NULL, #Title of the graph
                            Yvalue = c("percentage", "absolute"),
                            stats = TRUE,
                            test.statistics = c("wilcox.test", "t.test"),
                            paired = FALSE, #Whether to compute a paired test
                            corrMethod = "holm", #p value correction method
                            BSL = "BSL2", #Name of the baseline timepoint for wilcox tests
                            meanMed = c("mean", "median"), #display mean or median ? 
                            palette = NULL) { 
  
  #CHECKING VALIDITY OF ARGUMENTS
  
  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }
  
  if (ncol(CYTdata@metadata)==0) {
    stop("Error : Missing metadata slot for 'CYTdata' argument.")
  }
  
  level = match.arg(level) #checks if level matches what is defined in the function's definition
  checkmate::qassert(level, "S1") 
  
  Yvalue = match.arg(Yvalue)
  checkmate::qassert(Yvalue, "S1")
  
  test.statistics = match.arg(test.statistics)
  checkmate::qassert(test.statistics, "S1")
  
  meanMed = match.arg(meanMed)
  checkmate::qassert(meanMed, "S1")
  
  checkmate::qassert(BSL, "S1")
  checkmate::qassert(population, c("S*", "0"))
  checkmate::qassert(samples, c("S*", "0"))
  checkmate::qassert(stats, "B")
  checkmate::qassert(paired, "B")
  checkmate::qassert(corrMethod, "S1")
  checkmate::qassert(palette, c("S1", "0"))
  checkmate::qassert(marker, c("S1", "0"))
  checkmate::qassert(th, c("R1", "0"))
  
  if (!is.null(marker) & is.null(th)) {stop("ERROR: marker name provided with no threshold (th).")}
  if (is.null(marker) & !is.null(th)) {stop("ERROR: Threshold (th) provided with no marker name.")}
  
  if (!is.null(title)) {checkmate::qassert(title, "S1")}
  
  #ASSIGNING DEFAULT VALUES
  
  if (is.null(title)) {title = paste("Population", Yvalue ,"evolution")}
  
  if (is.null(population) & level == "clusters") {population = CYTdata@Clustering@clusters %>% levels()}
  if (is.null(population) & level == "metaclusters") {population = CYTdata@Metaclustering@metaclusters %>% levels()}
  
  #EXTRACTING DATA
  
  popRelAbundance = NULL 
  
  for (pop in population){ 
    
    if (!is.null(marker)) {
      
      popRelAbundance[[pop]] = getRelativeAbundanceThreshold(CYTdata,
                                                             marker = marker,
                                                             th = th,
                                                             population = pop,
                                                             samples = samples,
                                                             level = level, 
                                                             Yvalue = Yvalue) %>%
        reshape2::melt() %>% 
        reshape2::dcast(Individual~Timepoint) %>%
        tibble::column_to_rownames(var = "Individual")
      
    } else {
      
      popRelAbundance[[pop]] = Bahamas::getRelativeAbundance(CYTdata,
                                                             population = pop,
                                                             samples = samples,
                                                             level = level, 
                                                             Yvalue = Yvalue) %>%
        reshape2::melt() %>% 
        reshape2::dcast(Individual~Timepoint) %>%
        tibble::column_to_rownames(var = "Individual")
      
    }
    
    
    mean = colMeans(popRelAbundance[[pop]])
    med   = popRelAbundance[[pop]] %>% as.matrix() %>% matrixStats::colMedians()
    
    popRelAbundance[[pop]]["mean",] = mean
    popRelAbundance[[pop]]["med",] = med
    
    popRelAbundance[[pop]] = popRelAbundance[[pop]] %>% t()
    
  }
  
  #Creating dataframe for ggplot2
  
  data = reshape2::melt(popRelAbundance)
  
  #COMPUTING STATISTICS (split in two blocs for formatting reasons)
  
  if (stats) {  
    
    pairLevels = CYTdata@metadata$Timepoint %>% droplevels() %>% levels() %>% setdiff(BSL)
    statPairs = paste(BSL, pairLevels) %>% strsplit(" ")
    
    pDataF = statsLineplot(dataMelted = data,
                           statPairs = statPairs,
                           individuals = CYTdata@metadata$Individual %>% levels(),
                           test.statistics = test.statistics,
                           paired = paired,
                           corrMethod = corrMethod)
    
  }
  
  #Finishing formatting the dataframe
  
  colnames(data) = c("Timepoint","Individual","Value","Population")
  
  dataInd = data[data$Individual != "mean" & data$Individual != "med",]
  
  if (meanMed == "mean") {dataMeanMed = data[data$Individual == "mean",]
  } else {dataMeanMed = data[data$Individual == "med",]}
  
  dataInd$IndividualPop = paste(dataInd$Individual, dataInd$Population, sep = "_")
  dataMeanMed$Individual = dataMeanMed$Population
  
  if (stats) {
    
    #Defining symbols for displaying the p.values on the graph
    
    pDataF = reshape2::melt(pDataF %>% as.matrix())
    colnames(pDataF) = c("Timepoint", "Population", "p")
    
    pDataF = pDataF[pDataF$p <= 0.1,]
    pDataF$`Symbol` = "***"
    pDataF$`Symbol`[pDataF$p > 0.001] = "**"
    pDataF$`Symbol`[pDataF$p > 0.01] = "*"
    pDataF$`Symbol`[pDataF$p > 0.05] = "°"
    
    #Defining the y position where the p.value's symbol will be displayed. 
    
    yPos = NULL
    
    for (i in 1:nrow(pDataF)){
      
      yPos[i] = dataMeanMed$Value[which(dataMeanMed$Population == pDataF$Population[i] & as.vector(dataMeanMed$Timepoint) == pDataF$Timepoint[i])]
      
    }
    
    pDataF$`yPos` = yPos
    
  }  
  
  #CREATING THE GRAPH
  
  #Variables
  
  #Keeping only point shapes that are easy to differentiate, and in priority the easiest
  #It also possible to use unicode values
  pointShape = c(15,19,17,18,21,22,23,24,25,3,4,8,9,7,10,11,12,13,14)[1:(dataInd$Individual %>% droplevels() %>% levels() %>% length())]
  
  #Base
  
  Lplot = ggplot2::ggplot(dataInd, 
                          ggplot2::aes(x = Timepoint, y = Value, color = Population)) +
    ggplot2::geom_line(ggplot2::aes(group = IndividualPop),
                       alpha = .8) +
    ggplot2::geom_line(data = dataMeanMed, 
                       ggplot2::aes(group = Individual),
                       alpha = 1, 
                       size = 1.3) +
    ggplot2::geom_point(ggplot2::aes(shape = Individual), size = 2, fill = "white")
  
  #Titles
  
  if (Yvalue == "percentage") {
    
    Lplot <- Lplot + ggplot2::labs(title = title,
                                   x = "Timepoint",
                                   y = "Average relative abundance (%)")
    
  } else {
    
    Lplot <- Lplot + ggplot2::labs(title = title,
                                   x = "Timepoint",
                                   y = "Cell count")
    
  }
  
  #Statistical annotation
  
  if (stats) {  
    
    if (test.statistics == "wilcox.test" & paired) {
      
      annotateText = "Wilcoxon (paired)" #"Wilcoxon signed rank test (paired)"
      
    } else if (test.statistics == "wilcox.test" & !paired) {
      
      annotateText = "Wilcoxon (unpaired)" #"Wilcoxon rank sum test (unpaired)"
      
    } else if (test.statistics == "t.test" & paired) {
      
      annotateText = "t-Test (paired)" #"Student's t-Test (paired)"
      
    } else if (test.statistics == "t.test" & !paired) {
      
      annotateText = "t-Test (unpaired)" #"Student's t-Test (unpaired)"
      
    }
    
    
    Lplot <- Lplot + 
      ggplot2::geom_text(mapping = ggplot2::aes(x = Timepoint, y = yPos),
                         data = pDataF,
                         label = pDataF$Symbol,
                         size = 6,
                         position = ggplot2::position_nudge(x = 0.1, y = 0),
                         hjust = "left",
                         show.legend = F) +
      ggplot2::annotate("text", x = Inf, y = Inf,
                        label = paste(annotateText, "\nTimepoint to ", BSL, ":\n*** p < 0.001\n** p < 0.01\n* p < 0.05\n° p < 0.1\nCorrection: ", corrMethod, sep = ""),
                        hjust = "left", vjust = "top")
    
  }
  
  #Theme
  
  if (is.null(palette)) {
    if (level == "clusters") {
      Lplot = Lplot + 
        ggplot2::scale_color_manual(values = CYTdata@Clustering@palette)
    } else {
      Lplot = Lplot + 
        ggplot2::scale_color_manual(values = CYTdata@Metaclustering@palette)
    }
  } else {
    Lplot = Lplot + 
      # ggplot2::scale_colour_brewer(palette = palette) #If a manual palette was set : override CYTdata palette
      ggplot2::scale_color_manual(values = palette) #If a manual palette was set : override CYTdata palette
  }
  
  Lplot <- Lplot + 
    ggplot2::scale_shape_manual(values = pointShape) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme(legend.justification  = "bottom", legend.margin = ggplot2::margin(0.2,1,0.5,0.2,"cm")) #Top, right, bottom, left. Right at 1 cm to increase size and avoid annotate being cut off. 
  
  return(Lplot)
  
}


#' @title Plots cell cluster abundances using a boxplot representation
#'
#' @description This function aims to visualize and compare the cell cluster abundances for each biological condition using boxplot and jitter representations.
#'
#' The abundance of a specific cell cluster or a set of cell clusters can be displayed.
#' The representation can be restricted to a specific set of samples.
#' Moreover, boxplot can be constructed based on sample meta information.
#' Statistics can be computed for all comparisons.
#'
#' @param CYTdata a CYTdata object
#' @param population a character vector containing the names of the population to use.
#' @param level a character value containing the levels to display in the plot, either "clusters" (default) or "metaclusters"
#' @param samples a character vector containing the names of biological samples to use. By default, all samples are used
#' @param observationMetadata a character value containing the biological condition to use for the boxplot groups. 
#' The name must match a column of the CYTdata's metadata slot. By default, "Timepoint" is used.
#' @param groupMetadata a character value containing the biological condition (as for observationMetadata) to use for a secondary division of the boxplot groups. Defaults to NULL. 
#' @param pointMetadata a character value containing the biological condition (as for observationMetadata) to use for displaying point shape. Defaults to NULL. 
#' @param Yvalue a character value containing the quantitative variable to display. Possible values are "percentage" (the relative abundance, default) or "absolute." (the cell count)
#' @param computePval a boolean value indicating if statistics should be computed. 
#' @param test.statistics a character value providing the type of statistical test to use. Possible values are: 'wilcox.test' or 't.test'
#' @param paired a boolean value indicating if a paired or unpaired comparison should be applied
#' @param corrMethod a character value providing the type of correction for the p.value. Please refer to rstatix::wilcox_test or rstatix::t_test for details. 
#' 
#' @return a ggplot2 object
#'
#' @export
#'

plotBoxplot <- function(CYTdata,
                        population,
                        level = c("clusters", "metaclusters"),
                        samples = NULL,
                        observationMetadata = "Timepoint",
                        groupMetadata = NULL,
                        pointMetadata = NULL,
                        Yvalue = c("percentage", "absolute"),
                        computePval = FALSE,
                        test.statistics = c("wilcox.test", "t.test"),
                        paired = FALSE,
                        corrMethod = "holm") {
  
  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }
  
  if (ncol(CYTdata@metadata)==0) {
    stop("Error : Missing metadata slot for 'CYTdata' argument.")
  }
  checkmate::qassert(observationMetadata, "S1")
  checkmate::qassert(groupMetadata, c("0", "S1"))
  checkmate::qassert(pointMetadata, c("0", "S1"))
  if (!observationMetadata %in% colnames(CYTdata@metadata)){
    stop("Error : 'observationMetadata' argument is not a metadata (ex : 'Timepoint', 'Individual')")
  }
  if (!is.null(groupMetadata)){
    if (!groupMetadata %in% colnames(CYTdata@metadata)){
      stop("Error : 'groupMetadata' argument is not a metadata (ex : 'Timepoint', 'Individual')")
    }
    if (groupMetadata==observationMetadata){
      stop("Error : Arguments 'observationMetadata' and 'groupMetadata' are the same metadata")
    }
  }
  if (!is.null(pointMetadata)){
    if (!pointMetadata %in% colnames(CYTdata@metadata)){
      stop("Error : 'pointMetadata' argument is not a metadata (ex : 'Timepoint', 'Individual')")
    }
    if (pointMetadata==observationMetadata){
      stop("Error : Arguments 'observationMetadata' and 'pointMetadata' are the same metadata")
    }
  }
  
  level = match.arg(level)
  checkmate::qassert(level, "S1")
  Yvalue = match.arg(Yvalue)
  checkmate::qassert(Yvalue, "S1")
  data = getRelativeAbundance(CYTdata, population, level, samples, Yvalue)
  
  if (Yvalue == "percentage") { Ylab = "Abundance of population" }
  else { Ylab = "Absolute count of population" }
  
  checkmate::qassert(computePval, "B1")
  test.statistics = match.arg(test.statistics)
  checkmate::qassert(test.statistics, "S1")
  checkmate::qassert(paired, "B1")
  
  switch(test.statistics,
         wilcox.test = { test_fct <- rstatix::wilcox_test },
         t.test = { test_fct <- rstatix::t_test })
  
  data[[observationMetadata]] = droplevels(data[[observationMetadata]])
  colnames(data)[colnames(data)==observationMetadata] = "observationMetadata"
  if (!is.null(groupMetadata)){
    data[[groupMetadata]] = droplevels(data[[groupMetadata]])
    colnames(data)[colnames(data)==groupMetadata] = "groupMetadata"
  }
  if (!is.null(pointMetadata)){
    data[[pointMetadata]] = droplevels(data[[pointMetadata]])
    colnames(data)[colnames(data)==pointMetadata] = "pointMetadata"
    if (length(unique(data[[pointMetadata]]))>6) {
      message("Warning : 'pointMetadata' argument's corresponding metadata has more
            than 6 different values which is greater than the number of point shapes
            available in ggplot2/geom_point. pointMetadata not taken into account and
            argument set to NULL.")
      pointMetadata = NULL
    }
  }
  
  # return(data)
  
  if (!is.null(groupMetadata)){
    
    plot <- ggplot2::ggplot(data, ggplot2::aes(x = observationMetadata, y = value)) +
      ggplot2::geom_boxplot(ggplot2::aes(color = groupMetadata),
                            size=0.1, fatten=1, linewidth=1, outlier.shape=NA, width=0.3,
                            position=ggplot2::position_dodge(0.5)) +
      ggplot2::labs(fill = groupMetadata) +
      ggplot2::guides(color = ggplot2::guide_legend(title = groupMetadata, override.aes = list(size=3, shape=NA)))
    
    if (computePval) {
      statGroup = data %>%
        group_by(observationMetadata) %>%
        test_fct(value ~ groupMetadata, paired = paired, p.adjust.method = corrMethod) %>%
        rstatix::add_xy_position(x = "observationMetadata", dodge = 0.8) %>%
        rstatix::add_significance("p")
      print(statGroup)
      statGroup = statGroup %>%
        filter(p.adj.signif != "ns") %>%
        filter(!is.na(p.adj.signif))
      
      if (nrow(statGroup)>0) {
        plot <- plot + ggpubr::stat_pvalue_manual(data = statGroup,
                                                  label = "p = {p.adj.signif}",
                                                  hide.ns = F, color = "black", size = 5, tip.length = 0)
      }
      
    }
    
    if (!is.null(pointMetadata)){
      plot <- plot + ggplot2::geom_point(ggplot2::aes(group = groupMetadata, shape = pointMetadata),
                                         position = ggplot2::position_dodge(0.5),
                                         color = "black", size = 3) +
        ggplot2::scale_shape_manual(values=1:nlevels(data$pointMetadata)) +
        ggplot2::labs(shape = pointMetadata)
    } else {
      plot <- plot + ggplot2::geom_point(ggplot2::aes(group = groupMetadata),
                                         position = ggplot2::position_dodge(0.5),
                                         color = "black", shape = 16, size = 3)
    }
  } else {
    plot <- ggplot2::ggplot(data, ggplot2::aes(x = observationMetadata, y = value, color = observationMetadata)) +
      ggplot2::geom_boxplot(ggplot2::aes(color = observationMetadata), outlier.shape = NA, linewidth=1, size = 0.2, fatten = 1, width = 0.3, position=ggplot2::position_dodge(0.5)) +
      ggplot2::guides(color = ggplot2::guide_legend(title = observationMetadata, override.aes = list(size=3, shape=NA)))
    
    if (!is.null(pointMetadata)){
      plot <- plot + ggplot2::geom_point(ggplot2::aes(group = observationMetadata),
                                         position = ggplot2::position_dodge(0.5),
                                         color = "black", size = 5) +
        ggplot2::labs(shape = pointMetadata)
    } else {
      plot <- plot + ggplot2::geom_point(ggplot2::aes(group = observationMetadata),
                                         position = ggplot2::position_dodge(0.5),
                                         color = "black", shape = 16, size = 3)
    }
    
    if (computePval) {
      statObservation = data %>%
        test_fct(value ~ observationMetadata, paired = paired, p.adjust.method = corrMethod) %>%
        rstatix::add_xy_position(x = "observationMetadata") %>%
        rstatix::add_significance("p") %>%
        filter(p.adj.signif != "ns") %>%
        filter(!is.na(p.adj.signif))
      print(statObservation)
      
      if (nrow(statObservation)>0) {
        plot <- plot + ggpubr::stat_pvalue_manual(data = statObservation,
                                                  label = "p = {p.adj.signif}",
                                                  hide.ns = F, color = "black", size = 5, tip.length = 0)
      }
    }
  }
  
  plot <- plot +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.1))) +
    ggplot2::labs(title = "Boxplot representation",
                  subtitle = paste0(level, ": ", paste0(population, collapse = ", "))) +
    ggplot2::ylab(Ylab) +
    ggplot2::xlab(observationMetadata) +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5, face = "italic"),
                   aspect.ratio = 1,
                   panel.grid.minor =  ggplot2::element_blank(),
                   panel.grid.major =  ggplot2::element_blank(),
                   legend.position = "right")
  
  return(plot)
  
}


#' @title Plot pie charts of population abundance
#
#' @description This function aims to visualize the population composition of a CYTdata according to specified metadata (e.g., timepoint and individual), by plotting pie charts.
#' The pie charts are plotted on a 2D grid, with each axis showing one kind of metadata, allowing easy comparison of cell composition.  
#' Theoretically, the function should allow variations in the size of the pies to account for cell count. However, it is currently recommended to use value = "abundance" and to not set pieSize and byMeanSize to TRUE, as these parameters of the function are bugged.
#'
#' @param CYTdata a CYTdata object
#' @param samples a character vector containing the names of biological samples to use. By default, all samples are used
#' @param population a character vector containing the identifiers of the clusters/metaclusters to display in the plot. By default, all clusters or metaclusters are used.
#' @param value a character value containing the variable to be studied. Possible values are "abundance" (default, representing the relative abundance) or "count" (the cell count).
#' Please note that "count" is only relevant if pieSize is set to TRUE, and will make the pie vary in size depending on cell count. 
#' @param level a character value containing the levels to display in the plot, either "clusters" (default) or "metaclusters"
#' @param Xaxis a string being the metadata used to plot pie charts on the horizontal line. By default, "Timepoint" metadata is used.
#' @param Yaxis a string being the metadata used to plot pie charts on the vertical line. By default, "Individual" metadata is used.
#' @param NFSValues a named numeric with blood count (numération formule sanguine) per sample. Optional. If provided, the abundance data will be multiplied by the blood count for pie size adjustement.  
#' @param pieSize a boolean, controls pie size according to cell count. Bugged..
#' @param byMeanSize a boolean, controls pie size according to cell count. Bugged..
#' @param colorBar a character string indicating the color for the pie chart's borders. Defaults to black. 
#'
#' @return a ggplot2 object
#'
#' @export
#'

plotPiechart <- function(CYTdata,
                         samples = NULL,
                         population = NULL,
                         value = c("abundance", "count"),
                         level = c("clusters", "metaclusters"),
                         Xaxis = "Timepoint",
                         Yaxis = "Individual",
                         NFSvalues = NULL,
                         pieSize = FALSE, byMeanSize = FALSE,
                         colorBar = "#000000"){
  
  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }
  
  level = match.arg(level)
  checkmate::qassert(level, "S1")
  population = checkorderPopulation(CYTdata, population = population, level = level, order=TRUE, checkDuplicates=TRUE)
  samples = checkorderSamples(CYTdata, samples, order=TRUE)
  
  if (ncol(CYTdata@metadata)==0) {
    stop("Error : Missing metadata slot for 'CYTdata' argument.")
  }
  checkmate::qassert(Xaxis, "S1")
  checkmate::qassert(Yaxis, "S1")
  if (!Xaxis %in% colnames(CYTdata@metadata)){
    stop("Error : 'Xaxis' argument is not a metadata (ex : 'Timepoint', 'Individual')")
  }
  if (!Yaxis %in% colnames(CYTdata@metadata)){
    stop("Error : 'Yaxis' argument is not a metadata (ex : 'Timepoint', 'Individual')")
  }
  if (Xaxis == Yaxis){
    stop("Error : 'Xaxis' and 'Yaxis' arguments are the same metadata")
  }
  metadata = CYTdata@metadata[rownames(CYTdata@metadata) %in% samples,, drop=FALSE]
  Xlevels = levels(metadata[[Xaxis]])
  Ylevels = levels(metadata[[Yaxis]])
  comblevels = expand.grid(Xlevels, Ylevels)
  
  checkmate::qassert(NFSvalues, c("0", "N*"))
  checkmate::qassert(pieSize, "B1")
  checkmate::qassert(byMeanSize, "B1")
  checkmate::qassert(colorBar, "S1")
  if (!areColors(colorBar)) {
    stop("Error : 'colorBar' argument must be a hexadecimal color (ex : #FFFFFF).")
  }
  
  if (level == "metaclusters"){
    subOcject = CYTdata@Metaclustering
    lev = levels(subOcject@metaclusters)
  }
  else {
    subOcject = CYTdata@Clustering
    lev = levels(subOcject@clusters)
  }
  
  if (is.null(NFSvalues)) {
    value = match.arg(value)
    if (value=="abundance") {
      matrixCount = subOcject@abundance[population, samples, drop=FALSE]
    } else {
      matrixCount = subOcject@cellcount[population, samples, drop=FALSE]
    }
  }
  else {
    namesNFS = names(NFSvalues)
    if(sum(is.na(namesNFS))>0) {
      stop("Error : 'NFSvalues' argument's names contain NA values.
         It must be a fully named vector of numerical values.
         Associated ", level, " names are missing for NFS values : ",
           paste0(NFSvalues[is.na(namesNFS)], collapse=","), ".")
    }
    namesNFSdup = namesNFS[duplicated(namesNFS)]
    if (length(namesNFSdup)>0) {
      stop("Error : 'NFSvalues' argument's names contain duplicated values ( ",
           paste0(namesNFSdup, collapse = ", "),
           " ). Names must be vector of unique samples levels.")
    }
    
    splErr = setdiff(namesNFS, levels(CYTdata@samples))
    if (length(splErr)>0) {
      stop("Error : 'NFSvalues' argument's names providing identifiers not present in samples vector (",
           paste0(splErr, collapse=", "), ").")
    }
    splErr = setdiff(samples, namesNFS)
    if (length(splErr)>0) {
      stop("Error : 'NFSvalues' argument's names not providing all the samples needed by samples argument (",
           paste0(splErr, collapse=", "), ").")
    }
    matrixAbundance = subOcject@abundance[population, samples, drop=FALSE]
    matrixCount = data.matrix(matrixAbundance) %*% diag(NFSvalues[samples])
    matrixCount = data.frame(matrixCount)
    colnames(matrixCount) = samples
    rownames(matrixCount) = population
  }
  data = data.frame()
  for (i in 1:nrow(comblevels)) {
    spls = rownames(subset(metadata, (metadata[[Xaxis]] == comblevels[i,1]) & (metadata[[Yaxis]] == comblevels[i,2])))
    
    if (byMeanSize) {
      valueSpl = apply(matrixCount[,spls, drop=FALSE], 1, mean)
    } else {
      valueSpl = apply(matrixCount[,spls, drop=FALSE], 1, sum)
    }
    #if (length(spls)==1) { valueSpl = matrixCount[,spls] }
    #else { valueSpl = apply(matrixCount[,spls], 1, sum) }
    dataPie = data.frame("value" = valueSpl,
                         "group" = factor(population, levels = lev),
                         "Xaxis" = rep(comblevels[i,1], length(valueSpl)),
                         "Yaxis" = rep(comblevels[i,2], length(valueSpl)))
    if (pieSize) { dataPie$size = rep(sum(dataPie$value), nrow(dataPie)) }
    else { dataPie$size = rep(1, nrow(dataPie)) }
    data = rbind.data.frame(data, dataPie)
  }
  
  plot = ggplot2::ggplot(data,
                         ggplot2::aes(x= size/2, y = value, fill = group, width = size)) +
    ggplot2::geom_bar(position="fill", stat="identity", color = colorBar, width = 1) +
    ggplot2::scale_fill_manual(values = subOcject@palette) +
    ggplot2::facet_grid(factor(Yaxis, levels=Ylevels) ~ factor(Xaxis, levels=Xlevels)) +
    ggplot2::coord_polar("y", start = 0) +
    ggplot2::ylab("") + ggplot2::xlab("") +
    ggplot2::guides(fill=guide_legend(title=level)) +
    ggplot2::theme_classic() +
    ggplot2::theme(strip.text.x = ggplot2::element_text(size=12, face="bold.italic"),
                   strip.text.y = ggplot2::element_text(size=12, face="bold.italic"),
                   plot.title = ggplot2::element_text(hjust = 0.5),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5, face = "italic"),
                   axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   panel.grid  = ggplot2::element_blank(),
                   aspect.ratio = 1,
                   legend.position = "bottom")
  
  return(plot)
}


#' @title Plot a barplot of cell population
#'
#' @description This function aims to display a barplot showing the abundance of different cell populations across specific metadata. 
#'
#' @param CYTdata a CYTdata object.
#' @param population a character vector containing the names of the population to use.
#' @param level a character value containing the levels to display in the plot, either "clusters" (default) or "metaclusters"
#' @param samples a character vector containing the names of biological samples to use. By default, all samples are used
#' @param observationMetadata a character value containing the biological condition to use for the barplots, i.e., the column of the metadata dataframe. Defaults to "Timepoint".
#' @param groupMetadata a character value containing the biological condition to use for further subdivision of the barplots. Optional, defaults to NULL.
#' @param Yvalue a character value specifying whether the total cell count for each barplot should be normalized to 1 ("ABrelative", default; i.e., each barplot at the same height) or if actual cell count should be displayed with varying barplot height ("CCabsolute")
#'
#' @return a ggplot2 object
#'
#' @export
#'

plotStackedBarplot <- function(CYTdata,
                               population = NULL,
                               level = c("clusters", "metaclusters"),
                               samples = NULL,
                               observationMetadata = "Timepoint",
                               groupMetadata = NULL,
                               Yvalue = c("ABrelative", "CCabsolute")) {
  
  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }
  
  level = match.arg(level)
  checkmate::qassert(level, "S1")
  population = checkorderPopulation(CYTdata, population=population, level=level,
                                    order=TRUE, checkDuplicates=TRUE)
  samples = checkorderSamples(CYTdata, samples, order=TRUE, checkDuplicates=TRUE)
  
  if (ncol(CYTdata@metadata)==0) {
    stop("Error : Missing metadata slot for 'CYTdata' argument.")
  }
  checkmate::qassert(observationMetadata, "S1")
  checkmate::qassert(groupMetadata, c("0", "S1"))
  if (!observationMetadata %in% colnames(CYTdata@metadata)){
    stop("Error : 'observationMetadata' argument is not a metadata (ex : 'Timepoint', 'Individual')")
  }
  if (!is.null(groupMetadata)){
    if (!groupMetadata %in% colnames(CYTdata@metadata)){
      stop("Error : 'groupMetadata' argument is not a metadata (ex : 'Timepoint', 'Individual')")
    }
    if (groupMetadata==observationMetadata){
      stop("Error : Arguments 'observationMetadata' and 'groupMetadata' are the same metadata")
    }
  }
  
  Yvalue = match.arg(Yvalue)
  checkmate::qassert(Yvalue, "S1")
  
  if (level == "clusters"){
    data = CYTdata@Clustering@cellcount[population, samples, drop=FALSE]
    palette = CYTdata@Clustering@palette
    lev = levels(CYTdata@Clustering@clusters)
  }
  else {
    data = CYTdata@Metaclustering@cellcount[population, samples, drop=FALSE]
    palette = CYTdata@Metaclustering@palette
    lev = levels(CYTdata@Metaclustering@metaclusters)
  }
  
  data = suppressWarnings(reshape::melt(data.matrix(data)))
  colnames(data) = c("population", "samples", "value")
  md = CYTdata@metadata
  md$samples = rownames(md)
  data = merge(data, md, by = "samples")
  data[[observationMetadata]] = droplevels(data[[observationMetadata]])
  colnames(data)[colnames(data)==observationMetadata] = "observationMetadata"
  if (!is.null(groupMetadata)){
    data[[groupMetadata]] = droplevels(data[[groupMetadata]])
    colnames(data)[colnames(data)==groupMetadata] = "groupMetadata"
  }
  if (length(unique(data$observationMetadata))==1) {
    stop("Error : 'samples' argument procides only one level of 'observationMetadata' argument (",
         unique(data$observationMetadata), "). Impossible to build stacked area plot.")
  }
  data$population = factor(data$population, levels = lev)
  
  data %>% write.csv("StackedBarplot.csv")
  
  if (!is.null(groupMetadata)){
    data = data %>%
      dplyr::group_by(population, observationMetadata, groupMetadata) %>%
      dplyr::summarise(y = sum(value))
    
    if(Yvalue == "ABrelative") {
      data = data %>%
        dplyr::group_by(observationMetadata, groupMetadata) %>%
        dplyr::mutate(y = y / sum(y))
    }
    plot <- ggplot2::ggplot(data %>% arrange(observationMetadata),
                            ggplot2::aes(fill = population,
                                         y = y,
                                         x = groupMetadata))
    
    if(Yvalue == "ABrelative") {
      plot <- plot + ggplot2::geom_bar(position="fill", stat="identity", width = 0.75)
      Ylab = "Abundance (Relative)"
    }
    else {
      plot <- plot + ggplot2::geom_bar(position="stack", stat="identity", width = 0.75)
      Ylab =  "Absolute counts"
    }
    
    plot <- plot + facet_grid(~observationMetadata, switch = "x")
  }
  else {
    data = data %>%
      dplyr::group_by(population, observationMetadata) %>%
      dplyr::summarise(y = sum(value))
    if(Yvalue == "ABrelative") {
      data = data %>%
        dplyr::group_by(observationMetadata) %>%
        dplyr::mutate(y = y / sum(y))
    }
    plot <- ggplot2::ggplot(data,
                            ggplot2::aes(fill = population,
                                         y = y,
                                         x = observationMetadata))
    
    if(Yvalue == "ABrelative") {
      plot <- plot + ggplot2::geom_bar(position="fill", stat="identity")
      Ylab = "Abundance (Relative)"
    }
    else {
      plot <- plot + ggplot2::geom_bar(position="stack", stat="identity")
      Ylab =  "Absolute counts"
    }
    
  }
  
  plot <- plot +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::labs(title = "Stacked barplot representation",
                  fill = "Population") +
    ggplot2::ylab(Ylab) +
    #ggplot2::xlab(observationMetadata) +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size=55, face="bold"),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5, face = "italic"),
                   panel.grid.minor =  ggplot2::element_blank(),
                   panel.grid.major =  ggplot2::element_blank(),
                   text = ggplot2::element_text(size = 35),
                   axis.title.y = element_text(size = 40, face = "bold"),
                   axis.title.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(size = 35, face = "bold"),
                   axis.text.y = ggplot2::element_text(size = 35, face = "bold"),
                   legend.position = "right")
  if (!is.null(groupMetadata)){ plot <- plot + ggplot2::theme(strip.placement = "outside",
                                                              strip.text.x = element_text(size = 40, face = "bold"),
                                                              strip.background = element_rect(fill = NA, color = "white"),
                                                              panel.spacing = unit(1.4, "lines")) }
  
  return(plot)
}

#' @title Plot a stack area chart of cell population
#'
#' @description This function aims to display a stack area chart showing the abundance of different cell populations across specific metadata. 
#' It is particularly suited for visualization across timepoints. 
#'
#' @param CYTdata a CYTdata object.
#' @param population a character vector containing the names of the populations to use. By default, all clusters or metaclusters are used. 
#' @param level a character value containing the levels to display in the plot, either "clusters" (default) or "metaclusters"
#' @param samples a character vector containing the names of biological samples to use. By default, all samples are used
#' @param observationMetadata a character value containing the biological condition to use for the area chart, i.e., the column of the metadata dataframe. Defaults to "Timepoint".
#' @param Yvalue a character value specifying whether the total cell count for displayed metadata value should be normalized to 1 ("ABrelative", default) or if actual cell count should be displayed with varying plot height ("CCabsolute").
#'
#' @return a ggplot2 object
#'
#' @export
#'

plotStackedArea <- function(CYTdata,
                            population = NULL,
                            level = c("clusters", "metaclusters"),
                            samples = NULL,
                            observationMetadata = "Timepoint",
                            Yvalue = c("ABrelative", "CCabsolute")) {
  
  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }
  
  level = match.arg(level)
  checkmate::qassert(level, "S1")
  population = checkorderPopulation(CYTdata, population=population, level=level,
                                    order=TRUE, checkDuplicates=TRUE)
  samples = checkorderSamples(CYTdata, samples, order=TRUE, checkDuplicates=TRUE)
  
  if (ncol(CYTdata@metadata)==0) {
    stop("Error : Missing metadata slot for 'CYTdata' argument.")
  }
  checkmate::qassert(observationMetadata, "S1")
  if (!observationMetadata %in% colnames(CYTdata@metadata)){
    stop("Error : 'observationMetadata' argument is not a metadata (ex : 'Timepoint', 'Individual')")
  }
  
  Yvalue = match.arg(Yvalue)
  checkmate::qassert(Yvalue, "S1")
  
  if (level == "clusters"){
    data = CYTdata@Clustering@cellcount[population, samples, drop=FALSE]
    palette = CYTdata@Clustering@palette
    lev = levels(CYTdata@Clustering@clusters)
  }
  else {
    data = CYTdata@Metaclustering@cellcount[population, samples, drop=FALSE]
    palette = CYTdata@Metaclustering@palette
    lev = levels(CYTdata@Metaclustering@metaclusters)
  }
  
  data = suppressWarnings(reshape::melt(data.matrix(data)))
  colnames(data) = c("population", "samples", "value")
  md = CYTdata@metadata
  md$samples = rownames(md)
  data = merge(data, md, by = "samples")
  data[[observationMetadata]] = droplevels(data[[observationMetadata]])
  colnames(data)[colnames(data)==observationMetadata] = "observationMetadata"
  if (length(unique(data$observationMetadata))==1) {
    stop("Error : 'samples' argument procides only one level of 'observationMetadata' argument (",
         unique(data$observationMetadata), "). Impossible to build stacked area plot.")
  }
  data$population = factor(data$population, levels = lev)
  
  data = data %>%
    dplyr::group_by(population, observationMetadata) %>%
    dplyr::summarise(y = sum(value))
  if(Yvalue == "ABrelative") {
    data = data %>%
      dplyr::group_by(observationMetadata) %>%
      dplyr::mutate(y = y / sum(y))
    Ylab = "Abundance (Relative)"
  }
  else { Ylab =  "Absolute counts" }
  
  plot <- ggplot2::ggplot(data,
                          ggplot2::aes(x = observationMetadata,
                                       y = y)) +
    ggplot2::geom_area(ggplot2::aes(colour = population,
                                    group = population,
                                    fill = population),
                       alpha=0.6, size=1, colour="white") +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::labs(title = "Stacked Area representation",
                  fill = level) +
    ggplot2::ylab(Ylab) +
    ggplot2::xlab(observationMetadata) +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5, face = "italic"),
                   aspect.ratio = 1,
                   panel.grid.minor =  ggplot2::element_blank(),
                   panel.grid.major =  ggplot2::element_blank(),
                   legend.position = "right")
  return(plot)
}


