#' @title Compute Principal Variance Component Analysis (PVCA) to determine the proportion of variability explained by selected effects (metadata condition)
#'
#' @description This function is pvcaBatchAssess function from pvca package with some adaptations. This function assess the batch sources by fitting all "sources" as random effects
#' including two-way interaction terms in the Mixed Model(depends on lme4 package) to selected principal components, which were obtained from the original data correlation matrix.
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param markers a character vector providing the cell markers to use. By default, all markers are used
#' @param batch.factors a vector of factors that the mixed linear model will be fit on
#' @param threshold the percentile value of the minimum amount of the variabilities that the selected principal components need to explain
#'
#' @return a list containing the PVCA features and a barplot (ggplot object)
#'
#' @export

PVCA.batch <- function(CYTdata,
                       batch.factors,
                       threshold,
                       markers = NULL){
  
  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  checkmate::qassert(batch.factors, "S*")
  checkmate::qassert(markers, c("0", "S*"))
  checkmate::qassert(threshold, "N1")
  
  if (is.null(markers)) { markers = colnames(CYTdata@matrix.expression) }
  data.matrix = cbind(CYTdata@matrix.expression[,markers], "SPL" = CYTdata@samples)
  data.matrix = plyr::ddply(data.matrix, "SPL",
                            function(x) {
                              x$SPL = NULL
                              apply(x, 2, stats::median)
                            })
  rownames(data.matrix) = data.matrix$SPL
  data.matrix$SPL = NULL
  
  abatch <- Biobase::ExpressionSet(assayData = t(data.matrix),
                                   phenoData = new("AnnotatedDataFrame", data=CYTdata@metadata[rownames(data.matrix),]))
  
  theDataMatrix <- Biobase::exprs(abatch)
  dataRowN <- nrow(theDataMatrix)
  dataColN <- ncol(theDataMatrix)
  
  ########## Center the data (center rows) ##########
  theDataMatrixCentered = matrix(data = 0, nrow = dataRowN, ncol = dataColN)
  theDataMatrixCentered_transposed = apply(theDataMatrix, 1, scale, center = TRUE, scale = FALSE)
  theDataMatrixCentered = t(theDataMatrixCentered_transposed)
  
  ########## Compute correlation matrix &  Obtain eigenvalues ##########
  eigenData = theDataMatrixCentered %>% cor() %>% eigen()
  eigenValues = eigenData$values
  eigenVectorsMatrix = eigenData$vectors
  ev_n = length(eigenValues)
  eigenValuesSum = sum(eigenValues)
  percents_PCs = eigenValues /eigenValuesSum
  
  ##	Getting the experimental information ##
  expInfo <- Biobase::pData(abatch)[,batch.factors]
  exp_design <- as.data.frame(expInfo)
  expDesignRowN <- nrow(exp_design)
  expDesignColN <- ncol(exp_design)
  
  ########## Merge experimental file and eigenvectors for n components ##########
  my_counter_2 = 0
  my_sum_2 = 1
  for (i in ev_n:1){
    my_sum_2  = my_sum_2 - percents_PCs[i]
    if ((my_sum_2) <= threshold ){ my_counter_2 = my_counter_2 + 1 }
  }
  
  if (my_counter_2 < 3){ pc_n  = 3 }
  else { pc_n = my_counter_2 }
  
  ## pc_n is the number of principal components to model
  pc_data_matrix <- matrix(data = 0, nrow = (expDesignRowN*pc_n), ncol = 1)
  mycounter = 0
  for (i in 1:pc_n){
    for (j in 1:expDesignRowN){
      mycounter <- mycounter + 1
      pc_data_matrix[mycounter,1] = eigenVectorsMatrix[j,i]
    }
  }
  
  AAA <- exp_design[rep(1:expDesignRowN,pc_n),]
  Data <- cbind(AAA,pc_data_matrix)
  
  
  ####### Edit these variables according to your factors #######
  variables <-c (colnames(exp_design))
  for (i in 1:length(variables)) { Data$variables[i] <- as.factor(Data$variables[i]) }
  
  ########## Mixed linear model ##########
  op <- options(warn = (-1))
  effects_n = expDesignColN  + choose(expDesignColN, 2) + 1
  randomEffectsMatrix <- matrix(data = 0, nrow = pc_n, ncol = effects_n)
  
  ## 	Get model functions ##
  model.func <- c()
  index <- 1
  
  ##	level-1
  for (i in 1:length(variables)) {
    mod = paste("(1|", variables[i], ")",   sep="")
    model.func[index] = mod
    index = index + 1
  }
  
  ##	two-way interaction
  for (i in 1:(length(variables)-1)) {
    for (j in (i+1):length(variables)) {
      mod = paste("(1|", variables[i], ":", variables[j], ")",   sep="")
      model.func[index] = mod
      index = index + 1
    }
  }
  
  function.mods <- paste(model.func , collapse = " + ")
  
  ## 	Get random effects ##
  
  for (i in 1:pc_n){
    y = (((i-1)*expDesignRowN)+1)
    funct <- paste("pc_data_matrix", function.mods, sep =" ~ ")
    Rm1ML <- lme4::lmer(funct,
                        Data[y:(((i-1)*expDesignRowN)+expDesignRowN),],
                        REML = TRUE, verbose = FALSE, na.action = na.omit)
    randomEffects <- Rm1ML
    randomEffectsMatrix[i,] <- c(unlist(lme4::VarCorr(Rm1ML)),resid=sigma(Rm1ML)^2)
  }
  effectsNames <- c(names(lme4::getME(Rm1ML,"cnms")),"resid")
  
  ########## Standardize Variance ##########
  randomEffectsMatrixStdze <- matrix(data = 0, nrow = pc_n, ncol = effects_n)
  for (i in 1:pc_n){
    mySum = sum(randomEffectsMatrix[i,])
    for (j in 1:effects_n){ randomEffectsMatrixStdze[i,j] = randomEffectsMatrix[i,j]/mySum }
  }
  
  ########## Compute Weighted Proportions ##########
  
  randomEffectsMatrixWtProp <- matrix(data = 0, nrow = pc_n, ncol = effects_n)
  for (i in 1:pc_n){
    weight = eigenValues[i]/eigenValuesSum
    for (j in 1:effects_n){ randomEffectsMatrixWtProp[i,j] = randomEffectsMatrixStdze[i,j]*weight }
  }
  
  ########## Compute Weighted Ave Proportions ##########
  
  randomEffectsSums <- matrix(data = 0, nrow = 1, ncol = effects_n)
  randomEffectsSums <-colSums(randomEffectsMatrixWtProp)
  totalSum = sum(randomEffectsSums)
  randomEffectsMatrixWtAveProp <- matrix(data = 0, nrow = 1, ncol = effects_n)
  for (j in 1:effects_n){ randomEffectsMatrixWtAveProp[j] = randomEffectsSums[j]/totalSum }
  
  PVCA.obj = data.frame("dat"=as.vector(randomEffectsMatrixWtAveProp), "label"=effectsNames)
  
  plot = ggplot2::ggplot(PVCA.obj, aes(x=label, y=dat)) +
    ggplot2::geom_bar(stat="identity", color="blue", fill="blue") +
    ggplot2::ggtitle("PVCA estimation bar chart") +
    ggplot2::xlab("Effects") + ggplot2::ylab("Weighted average proportion variance") + ggplot2::ylim(0,1.1) +
    ggplot2::theme(legend.position="none") +
    ggplot2::theme_minimal()
  
  return(list("dat"=randomEffectsMatrixWtAveProp,
              "label"=effectsNames,
              "barplot"=plot))
}

#' @title function for the cyCombine batch correction workflow.
#'
#' @description Compute the batch correction on the data using the ComBat algorithm.
#'  Define a covariate, either as a character vector or name of tibble column.
#'  The covariate should preferable be the cell condition types but can be any column that infers heterogeneity in the data.
#'  The function assumes that the batch information is in the "batch" column and the data contains a "sample" column with sample information.
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param batch.metadata a character vector containing the names of metadata associated to batch acquisition (name of a CYTdata's metadata dataframe columns)
#' @param markers a character vector containing the names of biological markers to correct. By default, all Marker are corrected
#' @param norm_method a character being the normalization method. Should be either 'rank', 'scale' or 'qnorm'. Default: 'scale'. "rank" is recommended when combining data with heavy batch effects
#' @param ties.method a character being the method to handle ties, when using rank. Default: 'average'
#' @param seed a numeric being the seed to use when creating the SOM.
#' @param xdim a numeric being the x-dimension size of the SOM.
#' @param ydim a numeric being the y-dimension size of the SOM.
#' @param rlen a numeric being the number of times the data is presented to the SOM network. Consider a larger value, if results are not convincing (e.g. 100)
#' @param covar.metadata a character being the covariate ComBat should use. Must be the name of a specific metadata (name of a CYTdata's metadata dataframe columns)
#' @param parametric Default: TRUE. If TRUE, the parametric version of ComBat is used. If FALSE, the non-parametric version is used.
#'
#' @return a S4 object of class 'CYTdata'
#'
#' @export

run.cyCombine <- function(CYTdata,
                          batch.metadata,
                          markers = NULL,
                          norm_method = c("scale","rank","qnorm"),
                          ties.method = c("average", "first", "last", "random", "max", "min"),
                          seed = 92,
                          xdim, ydim, rlen,
                          covar.metadata = NULL,
                          parametric = TRUE){
  
  checkmate::qassert(batch.metadata, "S1")
  checkmate::qassert(markers, c("0", "S*"))
  if(is.null(markers)){ markers = colnames(CYTdata@matrix.expression) }
  markers = as.vector(markers)
  
  norm_method = match.arg(norm_method)
  checkmate::qassert(norm_method, "S1")
  ties.method = match.arg(ties.method)
  checkmate::qassert(ties.method, "S1")
  checkmate::qassert(xdim, "N1")
  checkmate::qassert(ydim, "N1")
  checkmate::qassert(rlen, "N1")
  checkmate::qassert(seed, "N1")
  
  checkmate::qassert(parametric, "B1")
  
  cat("\n\nPreparing Data for batch correction..")
  
  if (length(CYTdata@metadata)==0) { stop("Error preprocessing : Metadata slot of CYTdata object is empty but necessary.") }
  
  ### Prepare data
  
  data = cbind.data.frame(
    CYTdata@matrix.expression,
    "sample" = CYTdata@samples,
    "batch" = CYTdata@metadata[match(CYTdata@samples, rownames(CYTdata@metadata)), batch.metadata],
    "id" = 1:nrow(CYTdata@matrix.expression)
  )
  
  checkmate::qassert(covar.metadata, c("0", "S1"))
  if(!is.null(covar.metadata)){
    covar.metadata = CYTdata@metadata[match(CYTdata@samples, rownames(CYTdata@metadata)),covar.metadata]
    # vector of length = nrow(CYTdata@matrix.expression) with value for each cell of metadata condition given by covar.metadata
  }
  
  cat("\n\nPerforming batch correction..")
  
  ### Correct data
  data.corrected = cyCombine::batch_correct(data,
                                            markers = markers,
                                            covar = covar.metadata,
                                            parametric = parametric,
                                            norm_method = norm_method,
                                            ties.method = ties.method,
                                            xdim = xdim, ydim = ydim, rlen = rlen,
                                            seed = seed)
  
  cat("\n\nCreating new CYTdata object ..")
  newCYTdata = CYTdata
  newCYTdata@matrix.expression[,markers] = data.corrected[,markers]
  validObject(newCYTdata)
  
  return(newCYTdata)
}

#' @title Function for the CytoNorm batch correction workflow.
#'
#' @description Compute the batch correction on the data using the CytoNorm algorithm. This algorithm requires control samples for each batch acquisition.
#' No CYTdata object is required, only the file path of conbtrol data that will train the model and dnon control data that will be corrected
#' For more infomation, see CytoNorm package : https://github.com/saeyslab/CytoNorm, https://rdrr.io/github/saeyslab/CytoNorm/man/
#'
#' @param files.validation a character vector containing the paths of each data file which will be corrected
#' @param labels.validation a character vector containing the labels of batch acquisition variable for validation samples
#' @param files.training a character vector containing the paths of each data control file which will be used to train the model
#' @param labels.training labels.validation a character vector containing the labels of batch acquisition variable for validation samples
#' @param channels a character vector containing the names of channels to correct. By default, the common channels to the validation and training sets of files are selected
#' @param transformList.training Transformation list to pass to the flowCore to transform training data.
#' Default : arcsinh with cofactor 5 on all the common channels to the validation and training sets of files)
#' @param transformList.validation Transformation list to pass to the flowCore to transform validation data.
#' Default : same as transformList.training
#' @param transformList.reverse.validation Transformation list to pass to the flowCore to transform validation data after correction.
#' In an attempt to get raw corrected validation data. Default : reverse of transformList.validation transformation for the same channels
#' @param FlowSOM.params List with parameters to pass to the FlowSOM algorithm. Default = list(nCells = 1000000, xdim = 15, ydim = 15, nclus = 30, scale = FALSE).
#' @param normParams
#' @param tmpDir Directory to put the temporary files in. Default = "./tmp
#' @param outputDir Directory to put the temporary files in. Default = "."
#' @param prefix Prefix to put in front of the normalized file names. Default = "Norm_"
#' @param seed set.seed is called with this argument for reproducable results. Default = 1
#' @param clean Whether to remove the temporary files again at the end. Default = TRUE
#' @param verbose If TRUE, extra output is printed while running.
#' @param plot If TRUE, plots are saved to the current dir. Default = FALSE.
#' @param recompute If FALSE, will try to reuse previously saved FlowSOM model. If so, a warning message will be printed. Default = FALSE
#' @param write logical indicating whether the normalised samples should be written to new FCS files in a Normalized directory within outputDir, set to TRUE by default.
#' @param ... Additional arguments to pass to read.FCS
#'
#' @return the paths of corrected FCS files
#'
#' @export

run.CytoNorm <- function(files.validation,
                         labels.validation,
                         files.training,
                         labels.training,
                         channels,
                         transformList.training = NULL,
                         transformList.validation = NULL,
                         transformList.reverse.validation = NULL,
                         FlowSOM.params = list(nCells = 1000000,
                                               xdim = 15,
                                               ydim = 15,
                                               nClus = 30,
                                               scale = FALSE),
                         normParams = list(nQ = 101, goal = "mean"),
                         tmpDir = "./tmp",
                         outputDir = "./CytoNorm",
                         prefix = "CytoNorm_",
                         seed = 1,
                         clean = TRUE,
                         verbose = FALSE,
                         plot = FALSE,
                         recompute = FALSE,
                         write = TRUE,
                         ...){
  
  
  checkmate::qassert(files.validation, "S*")
  checkmate::qassert(labels.validation, "S*")
  checkmate::qassert(files.training, "S*")
  checkmate::qassert(labels.training, "S*")
  checkmate::qassert(channels, c("0","S*"))
  
  checkmate::qassert(FlowSOM.params, "L*")
  checkmate::qassert(normParams, "L*")
  checkmate::qassert(tmpDir, "S1")
  checkmate::qassert(outputDir, "S1")
  checkmate::qassert(prefix, "S1")
  checkmate::qassert(seed, "N1")
  checkmate::qassert(clean, "B1")
  checkmate::qassert(verbose, "B1")
  checkmate::qassert(plot, "B1")
  checkmate::qassert(recompute, "B1")
  checkmate::qassert(write, "B1")
  
  # cols_training = get.common.colnames(files.training)
  # cols_validation = get.common.colnames(files.validation)
  # common_channels = intersect(names(cols_validation), names(cols_training))
  # if (length(common_channels)==0){ stop("Error : No common channels between training and validation files set") }
  # if (!is.null(channels)){
  #   index.error = is.na(match(channels, common_channels))
  #   if (any(index.error)) {
  #     stop("Error : Some channels that should be normalized are not common to training and validation files set (non conformable channels : ",
  #          paste0(channels[index.error], collapse=", ") , ").") }
  # }
  # else {
  #   message("No channels given as input. All the channels common to training and validation files sets will be normalized")
  #   channels = common_channels
  # }
  # message("Markers to normalize in training set (ordered) : ", paste0(cols_training[channels], collapse=", "))
  # message("Markers to normalize in validation set (ordered) : ", paste0(cols_validation[channels], collapse=", "))
  
  
  if (!is.null(transformList.training)){
    if (class(transformList.training)!="transformList") {
      stop("Error : argument 'transformList.training' must be of class 'transformList' or NULL.")
    }
    channels.to.transform = names(transformList.training@transforms)
    index.error = is.na(match(channels.to.transform, channels))
    if (any(index.error)) {
      stop("Error : Some channels that should be transformed, for training step,
           are not present in the set of channels that will be normalized (non conformable channels : ",
           paste0(channels.to.transform[index.error], collapse=", ") , ").")
    }
  } else { transformList.training <- flowCore::transformList(channels, CytoNorm::cytofTransform) }
  
  if (!is.null(transformList.validation)){
    if (class(transformList.validation)!="transformList") {
      stop("Error : argument 'transformList.validation' must be of class 'transformList' or NULL.")
    }
    channels.to.transform = names(transformList.validation@transforms)
    index.error = is.na(match(channels.to.transform, channels))
    if (any(index.error)) {
      stop("Error : Some channels that should be transformed, for validation step,
           are not present in the set of channels that will be normalized (non conformable channels : ",
           paste0(channels.to.transform[index.error], collapse=", ") , ").")
    }
  } else { transformList.validation <- flowCore::transformList(channels, CytoNorm::cytofTransform) }
  
  if (!is.null(transformList.reverse.validation)){
    if (class(transformList.reverse.validation)!="transformList") {
      stop("Error : argument 'transformList.reverse.validation' must be of class 'transformList' or NULL.")
    }
    channels.to.transform = names(transformList.reverse.validation@transforms)
    index.error = is.na(match(channels.to.transform, channels))
    if (any(index.error)) {
      stop("Error : Some channels that should be reverse-transformed, for validation step,
           are not present in the set of channels that will be normalized (non conformable channels : ",
           paste0(channels.to.transform[index.error], collapse=", ") , ").")
    }
  } else { transformList.reverse.validation <- flowCore::transformList(channels, CytoNorm::cytofTransform.reverse) }
  
  model = CytoNorm::CytoNorm.train(files = files.training,
                                   labels = labels.training,
                                   channels = channels,
                                   transformList = transformList.training,
                                   outputDir = tmpDir,
                                   FlowSOM.params = FlowSOM.params,
                                   normMethod.train = CytoNorm::QuantileNorm.train,
                                   normParams = normParams,
                                   seed = seed,
                                   clean = clean,
                                   plot = plot,
                                   verbose = verbose,
                                   recompute = recompute,
                                   ...)
  
  res = CytoNorm::CytoNorm.normalize(model = model,
                                     files = files.validation,
                                     labels = labels.validation,
                                     transformList = transformList.validation,
                                     transformList.reverse = transformList.reverse.validation,
                                     normMethod.normalize = CytoNorm::QuantileNorm.normalize,
                                     outputDir = outputDir,
                                     prefix = prefix,
                                     clean = clean,
                                     verbose = verbose,
                                     write = write,
                                     ...)
  
  files_corr = outputDir %>% list.files(pattern = "fcs", full.names = TRUE)
  
  return(files_corr)
}

