#' @title Gate cells in a CYTdata object
#'
#' @description Gate cells in a CYTdata object according to clusters, metaclusters or samples and create new CYTdata object.
#' Dimension reduction coordinates are not kept
#' Gating can be performed according to a specific metadata by setting subset to the name of the samples corresponding to the desired metadata, in the form subset = rownames(filter(CYTdata@metadata, Cohort == "GC" & Timepoint %in% c("BSL", "D1", "D3", "D6")))
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param subset a character vector containing the identifier of the clusters, metaclusters or samples to keep
#' @param levels a character value indicating whether to gate according to "samples" (default), "clusters" or "metaclusters".
#' 
#' @return a S4 object of class 'CYTdata'
#'
#' @export
#'

CYTdataGating <- function(CYTdata,
                          subset,
                          level = c("samples", "clusters", "metaclusters")){
  
  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }
  
  level = match.arg(level)
  checkmate::qassert(level, "S1")
  
  checkmate::qassert(subset, "S*")
  
  switch(level,
         samples = {
           subset = checkorderSamples(CYTdata, subset, order = TRUE, checkDuplicates = TRUE)
           gatedIdx = CYTdata@samples %in% subset
         },
         clusters = {
           subset = checkorderPopulation(CYTdata, population = subset, level = level, order = TRUE, checkDuplicates = TRUE)
           gatedIdx = CYTdata@Clustering@clusters %in% subset
         },
         metaclusters = {
           subset = checkorderPopulation(CYTdata, population = subset, level = level, order = TRUE, checkDuplicates = TRUE)
           gatedIdx = CYTdata@Metaclustering@metaclusters %in% subset
         })
  
  cat("\nGating CYTdata object and creating one object containing cells belonging to the following", level, " :",
      paste0(subset, collapse=", "))
  newmatrix.expression = subset(CYTdata@matrix.expression, gatedIdx)
  cat("\n\n - Number of cells gated :", sum(gatedIdx))
  
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
    newCYTdata@metadata = CYTdata@metadata[rownames(CYTdata@metadata) %in% levels(newCYTdata@samples),] %>%
      mutate_all(.funs = function(x) { return(factor(x, levels = levels(x)[levels(x) %in% x])) })
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
                                         optional_parameters = append(CYTdata@Clustering@optional_parameters, list("gating" = TRUE)))
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
                                             optional_parameters = append(CYTdata@Metaclustering@optional_parameters, list("gating" = TRUE)))
    message("\n\nRemark : Metaclustering results are preserved during gating operation but it is recommended to the user to run a
        new Metaclustering step with parameters adapted to the gated dataset")
  }
  
  cat("\n\n - DimReduction and statistics related slots weren't preserved during gating operation. The user has to repeat these steps")
  validObject(newCYTdata)
  return(newCYTdata)
}


#' @title Merge two CYTdata objects
#'
#' @description Merge two CYTdata objects. 
#'
#' @param CYTdata1 a S4 object of class 'CYTdata'
#' @param CYTdata1Name a character vector containing the name for the first CYTdata object. Defaults to "CYTdata1"
#' @param CYTdata2 a S4 object of class 'CYTdata'
#' @param CYTdata2Name a character vector containing the name for the second CYTdata object. Defaults to "CYTdata2"
#' @param markers a character vector containing the name of the markers to keep in the cell expression matrix. Defaults to all common markers. 
#' @param checkDuplicates logical. If TRUE, checks for cell duplicates (i.e., cells with the exact same marker expression). May be slow. Defaults to FALSE. 
#' @param mergePopulation logical. If TRUE (default), clusters and metaclusters with identical names in both CYTdata objects will be merged. If FALSE, a suffix will be appended to all or identical (see renameforAll) clusters/metaclusters depending on their original CYTdata. 
#' @param force logical. If TRUE, CYTdata objects will be merged even if one of them contains no metadata slot. If FALSE (default) and one of their metadata slots is empty, will return an error. 
#' @param sameSplMetadata logical. If TRUE (default), both CYTdatas are expected to have the same sample names. If FALSE, they must have no samples in common. 
#' @param addMetadata logical. If TRUE (default), adds a new metadata column indicating the original CYTdata's name. Only works if sameSplMetadata = FALSE.
#' @param renameforAll logical. If TRUE (default) and mergePopulation = FALSE, all clusters/metaclusters will be appended the name of their original CYTdata, even if they have no cluster/metacluster with the same name in the other CYTdata. 
#' 
#' @return a S4 object of class 'CYTdata'
#'
#' @export
#'

CYTdataMerging <- function(CYTdata1, CYTdata1Name = "CYTdata1",
                           CYTdata2, CYTdata2Name = "CYTdata2",
                           markers = NULL,
                           checkDuplicates = FALSE,
                           mergePopulation = TRUE,
                           force = FALSE,
                           sameSplMetadata = TRUE,
                           addMetadata = TRUE, renameforAll = TRUE) {
  
  if (class(CYTdata1)!="CYTdata") { stop("Error : argument 'CYTdata1' a S4 object of class 'CYTdata'.") }
  else { CYTdata1 = MakeValid(CYTdata1, verbose = TRUE) }
  if (class(CYTdata2)!="CYTdata") { stop("Error : argument 'CYTdata2' a S4 object of class 'CYTdata'.") }
  else { CYTdata2 = MakeValid(CYTdata2, verbose = TRUE) }
  
  checkmate::qassert(CYTdata1Name, "S1")
  checkmate::qassert(CYTdata2Name, "S1")
  checkmate::qassert(markers, "S*")
  checkmate::qassert(checkDuplicates, "B1")
  checkmate::qassert(mergePopulation, "B1")
  checkmate::qassert(force, "B1")
  checkmate::qassert(sameSplMetadata, "B1")
  checkmate::qassert(renameforAll, "B1")
  checkmate::qassert(addMetadata, "B1")
  
  ######################## expression matrix ###############################
  
  if(is.null(markers)) { message("Common markers are kept") }
  markers1 = checkorderMarkers(CYTdata1, markers, order = FALSE)
  markers2 = checkorderMarkers(CYTdata2, markers, order = FALSE)
  markers = intersect(markers1, markers2)
  cat("\n\nMarkers kept during merging : ", paste0(markers, collapse=", "))
  newmatrix.expression = rbind.data.frame(CYTdata1@matrix.expression[,markers], CYTdata2@matrix.expression[,markers])
  
  if (sameSplMetadata) {
    if (!setequal(levels(CYTdata1@samples), levels(CYTdata2@samples))) {
      stop("Error : 'sameSplMetadata' argument set to TRUE but CYTdata1 and CYTdata2 have different samples ( CYTdata1 :",
           paste0(setdiff(levels(CYTdata1@samples), levels(CYTdata2@samples)), collapse=", "), ", CYTdata2 :",
           paste0(setdiff(levels(CYTdata2@samples), levels(CYTdata1@samples)), collapse=", "), ").
           Please rename these samples before merging the CYTdata objects.")
    }
    newsamples = c(CYTdata1@samples, CYTdata2@samples)
    newmetadata = cbind.data.frame(CYTdata1@metadata[levels(CYTdata1@samples),],
                                   CYTdata2@metadata[levels(CYTdata1@samples),])[,union(colnames(CYTdata1@metadata),
                                                                                        colnames(CYTdata2@metadata))]
  }
  else {
    ######################## samples ###############################
    
    commonSamples = intersect(levels(CYTdata1@samples), levels(CYTdata2@samples))
    if (length(commonSamples)>0) {
      stop("Error : CYTdata1 and CYTdata2 have common samples. Please rename these samples before merging the CYTdata objects. (", paste0(commonSamples, collapse=", "),").")
    }
    newsamples = c(CYTdata1@samples, CYTdata2@samples)
    
    ######################## metadata ###############################
    
    if (nrow(CYTdata1@metadata)>0 && nrow(CYTdata2@metadata)>0) {
      colmetadata = intersect(colnames(CYTdata1@metadata), colnames(CYTdata2@metadata))
      if (length(colmetadata)>0){
        cat("\n\nMerging common metadata.. (", paste0(colmetadata, collapse=", "), ")")
        newmetadata = rbind.data.frame(CYTdata1@metadata[,colmetadata],
                                       CYTdata2@metadata[,colmetadata])
        if (addMetadata) {
          newmetadata$CYTdata = factor(c(rep(CYTdata1Name,nlevels((CYTdata1@samples))),
                                         rep(CYTdata2Name,nlevels((CYTdata2@samples)))))
        }
      }
      else {
        if (addMetadata) {
          cat("\n\n - No common metadata. Only CYTdata identifiers are given as metadata ('addMetadata' argument set to TRUE).")
          newmetadata = data.frame("CYTdata" = factor(c(rep(CYTdata1Name,nlevels(CYTdata1@samples)),
                                                        rep(CYTdata2Name,nlevels(CYTdata2@samples)))))
          rownames(newmetadata) = c(levels(CYTdata1@samples),levels(CYTdata2@samples))
        }
        else {
          cat("\n\n - No common metadata. Resulting metadata slot is empty.")
          newmetadata = data.frame()
        }
      }
    }
    else {
      if (nrow(CYTdata1@metadata)>0){
        if (force) {
          message("Warning : CYTdata1's metadata slot is given but CYTdata2's metadata slot is empty.
                metadata slots are ignored ('force' argument to TRUE).")
        }
        else {
          stop("Error : CYTdata1's metadata slot is given but CYTdata2's metadata slot is empty.
         Please specify metadata in CYTdata2 object or set 'force' argument to TRUE in order to merge and ignore metadata.")
        }
      }
      if (nrow(CYTdata2@metadata)>0){
        if (force) {
          message("Warning : CYTdata1's metadata slot is given but CYTdata2's metadata slot is empty.
                metadata slots are ignored ('force' argument to TRUE).")
        }
        else {
          stop("Error : CYTdata2's metadata slot is given but CYTdata1's metadata slot is empty.
         Please specify metadata in CYTdata1 object or set 'force' argument to TRUE in order to merge and ignore metadata.")
        }
      }
      if (addMetadata) {
        cat("\n\nAdding CYTdata identifiers to resulting empty metadata slot ('addMetadata' argument set to TRUE).")
        newmetadata = data.frame("CYTdata" = factor(c(rep(CYTdata1Name,length(levels((CYTdata1@metadata)))),
                                                      rep(CYTdata2Name,length(levels((CYTdata2@metadata)))))))
      }
      else {
        cat("\n\n Missing metadata for at least one of the two CYTdata objects. Resulting metadata slot is empty.")
        newmetadata = data.frame()
      }
    }
  }
  
  ############################## Clustering data ###########################################
  
  if (length(CYTdata1@Clustering@clusters)>0 && length(CYTdata2@Clustering@clusters)>0) {
    commonClusters = intersect(levels(CYTdata1@Clustering@clusters), levels(CYTdata2@Clustering@clusters))
    newClusters = as.character(c(CYTdata1@Clustering@clusters, CYTdata2@Clustering@clusters))
    if (length(commonClusters)>0) {
      if (mergePopulation) {
        message("Warning : CYTdata1 and CYTdata2 contains clusters with the same name (", paste0(commonClusters, collapse=", "),
                "). All the cells belonging to such clusters are merged together into a common cluster ('mergePopulation' argument set to TRUE).")
      }
      else {
        message("Warning : CYTdata1 and CYTdata2 contains clusters with the same name.
                All the cells belonging to such clusters are not merged but a suffix is added to the original clusters' names. ('mergePopulation' argument set to FALSE).")
        Ids = c(rep(CYTdata1Name, length(CYTdata1@Clustering@clusters)),
                rep(CYTdata2Name, length(CYTdata2@Clustering@clusters)))
        if (renameforAll) {
          newClusters = paste(newClusters, Ids)
        }
        else {
          Ids = c(rep(CYTdata1Name, length(CYTdata1@Clustering@clusters)),
                  rep(CYTdata2Name, length(CYTdata2@Clustering@clusters)))
          newClusters[newClusters %in% commonClusters] = paste(newClusters[newClusters %in% commonClusters], Ids[newClusters %in% commonClusters])
        }
      }
    }
    newClusters = as.factor(newClusters)
    newClustersPalette = c()
    for (cl in levels(newClusters)){
      if (cl %in% levels(CYTdata1@Clustering@clusters)) { col = CYTdata1@Clustering@palette[cl] }
      else { col = CYTdata2@Clustering@palette[cl] }
      newClustersPalette = c(newClustersPalette, col)
    }
    names(newClustersPalette) = levels(newClusters)
  }
  else {
    if (length(CYTdata1@Clustering@clusters)>0) {
      if (force) {
        message("Warning : CYTdata1's Clustering@clusters slot is given but CYTdata2's Clustering@clusters slot is empty.
                Clustering slots are ignored and resulting slot is empty ('force' argument to TRUE).")
      }
      else {
        stop("Error : CYTdata1's Clustering@clusters slot is given but CYTdata2's Clustering@clusters slot is empty.
             Please specify Clustering slot in CYTdata2 object or set 'force' argument to TRUE in order to merge and ignore Clustering slot.")
      }
    }
    if (length(CYTdata2@Clustering@clusters)>0) {
      if (force) {
        message("Warning : CYTdata2's Clustering@clusters slot is given but CYTdata1's Clustering@clusters slot is empty.
                Clustering slots are ignored and resulting slot is empty ('force' argument to TRUE).")
      }
      else {
        stop("Error : CYTdata2's Clustering@clusters slot is given but CYTdata1's Clustering@clusters slot is empty.
             Please specify Clustering slot in CYTdata1 object or set 'force' argument to TRUE in order to merge and ignore Clustering slot.")
      }
    }
    newClusters = c()
  }
  
  ############################## Metaclustering data ###########################################
  
  if (length(CYTdata1@Metaclustering@metaclusters)>0 && length(CYTdata2@Metaclustering@metaclusters)>0) {
    commonMetaclusters = intersect(levels(CYTdata1@Metaclustering@metaclusters), levels(CYTdata2@Metaclustering@metaclusters))
    newMetaclusters = as.character(c(CYTdata1@Metaclustering@metaclusters, CYTdata2@Metaclustering@metaclusters))
    if (length(commonMetaclusters)>0) {
      if (mergePopulation) {
        message("Warning : CYTdata1 and CYTdata2 contains metaclusters with the same name (",
                paste0(commonMetaclusters, collapse=", "), "). All the cells belonging to such metaclusters are merged together into a common metacluster ('mergePopulation' argument set to TRUE).")
      }
      else {
        message("Warning : CYTdata1 and CYTdata2 contains metaclusters with the same name.
                All the cells belonging to such metaclusters are not merged but a suffix is added to the original metaclusters names ('mergePopulation' argument set to FALSE).")
        if (renameforAll) {
          newMetaclusters = paste(newMetaclusters, Ids)
        }
        else {
          Ids = c(rep(CYTdata1Name, length(CYTdata1@Metaclustering@metaclusters)),
                  rep(CYTdata2Name, length(CYTdata2@Metaclustering@metaclusters)))
          newMetaclusters[newMetaclusters %in% commonMetaclusters] = paste(newMetaclusters[newMetaclusters %in% commonMetaclusters], Ids[newMetaclusters %in% commonMetaclusters])
        }
        
      }
    }
    newMetaclusters = as.factor(newMetaclusters)
    newMetaclustersPalette = c()
    for (cl in levels(newMetaclusters)){
      if (cl %in% levels(CYTdata1@Metaclustering@metaclusters)) { col = CYTdata1@Metaclustering@palette[cl] }
      else { col = CYTdata2@Metaclustering@palette[cl] }
      newMetaclustersPalette = c(newMetaclustersPalette, col)
    }
    names(newMetaclustersPalette) = levels(newMetaclusters)
  }
  else {
    if (length(CYTdata1@Metaclustering@metaclusters)>0) {
      if (force) {
        message("Warning : CYTdata1's Metaclustering@metaclusters slot is given but CYTdata2's Metaclustering@metaclusters slot is empty.
                Metaclustering slots are ignored and resulting slot is empty ('force' argument to TRUE).")
      }
      else {
        stop("Error : CYTdata1's Metaclustering@metaclusters slot is given but CYTdata2's Metaclustering@metaclusters slot is empty.
             Please specify Metaclustering slot in CYTdata2 object or set 'force' argument to TRUE in order to merge and ignore Metaclustering slot.")
      }
    }
    if (length(CYTdata2@Metaclustering@metaclusters)>0) {
      if (force) {
        message("Warning : CYTdata2's Metaclustering@metaclusters slot is given but CYTdata1's Metaclustering@metaclusters slot is empty.
                Metaclustering slots are ignored and resulting slot is empty ('force' argument to TRUE).")
      }
      else {
        stop("Error : CYTdata2's Metaclustering@metaclusters slot is given but CYTdata1's Metaclustering@metaclusters slot is empty.
             Please specify Metaclustering slot in CYTdata1 object or set 'force' argument to TRUE in order to merge and ignore Metaclustering slot.")
      }
    }
    newMetaclusters = c()
  }
  
  ######################## raw expression matrix ###############################
  
  if (ncol(CYTdata1@raw.matrix.expression)>0 && length(CYTdata2@raw.matrix.expression)>0) {
    newraw.matrix.expression = rbind.data.frame(CYTdata1@raw.matrix.expression[,markers], CYTdata2@raw.matrix.expression[,markers])
  }
  else {
    if (ncol(CYTdata1@raw.matrix.expression)>0) {
      if (force) {
        message("Warning : CYTdata1's raw.matrix.expression slot is given but CYTdata2's raw.matrix.expression slot is empty.
                raw.matrix.expression slots are ignored and resulting slot is empty ('force' argument to TRUE).")
      }
      else {
        stop("Error : CYTdata1's raw.matrix.expression slot is given but CYTdata2's raw.matrix.expression slot is empty.
             Please specify raw.matrix.expression slot in CYTdata2 object or set 'force' argument to TRUE in order to merge
             and ignore raw.matrix.expression slot.")
      }
    }
    if (ncol(CYTdata2@raw.matrix.expression)>0) {
      if (force) {
        message("Warning : CYTdata2's raw.matrix.expression slot is given but CYTdata1's raw.matrix.expression slot is empty.
                raw.matrix.expression slots are ignored and resulting slot is empty ('force' argument to TRUE).")
      }
      else {
        stop("Error : CYTdata2's raw.matrix.expression slot is given but CYTdata1's raw.matrix.expression slot is empty.
             Please specify raw.matrix.expression slot in CYTdata1 object or set 'force' argument to TRUE in order to merge
             and ignore raw.matrix.expression slot.")
      }
    }
    newraw.matrix.expression = data.frame()
  }
  
  ######################## raw expression matrix ###############################
  
  if (checkDuplicates){
    idxsDuplicates = duplicated(newmatrix.expression)
    newmatrix.expression = newmatrix.expression[!idxsDuplicates,]
    
    newsamples = newsamples[!idxsDuplicates]
    newsamples = droplevels(newsamples)
    if (ncol(newmetadata)>0) {
      newmetadata = newmetadata[levels(newsamples),]
    }
    
    if (length(newClusters)>0) {
      newClusters = newClusters[!idxsDuplicates]
      newClusters = droplevels(newClusters)
      newClustersPalette  = newClustersPalette[levels(newClusters)]
    }
    
    if (length(newMetaclusters)>0) {
      newMetaclusters = newMetaclusters[!idxsDuplicates]
      newMetaclusters = droplevels(newMetaclusters)
      newMetaclustersPalette  = newMetaclustersPalette[levels(newMetaclusters)]
    }
    
    if (ncol(newraw.matrix.expression)>0) { newraw.matrix.expression = newraw.matrix.expression[!idxsDuplicates,] }
  }
  
  ########################  Creating new CYTdata object ###############################
  
  cat("\n\n - Creating new CYTdata object")
  newCYTdata = methods::new("CYTdata",
                            samples = newsamples,
                            matrix.expression = newmatrix.expression,
                            metadata = newmetadata)
  if (ncol(CYTdata1@raw.matrix.expression)>0 && length(CYTdata2@raw.matrix.expression)>0) { newCYTdata@raw.matrix.expression = newraw.matrix.expression }
  
  ########################  Creating Clustering, Metaclustering ###############################
  
  if (length(newClusters)>0) {
    cat("\ncomputing new cell cluster cellcount, abundance matrix...")
    cellcount = compute.cellcount(newClusters, newCYTdata@samples)
    abundance = compute.abundance(cellcount)
    optional_parameters = list("merging" = TRUE)
    if (length(CYTdata1@Clustering@optional_parameters)>0) {
      optional_parameters = append(optional_parameters,
                                   list("CYTdata1" = CYTdata1@Clustering@optional_parameters))
    }
    if (length(CYTdata2@Clustering@optional_parameters)>0) {
      optional_parameters = append(optional_parameters,
                                   list("CYTdata2" = CYTdata2@Clustering@optional_parameters))
    }
    cat("\n\n - Creating Clustering slot")
    newCYTdata@Clustering = methods::new("Clustering",
                                         clusters = newClusters,
                                         cellcount = cellcount,
                                         abundance = abundance,
                                         palette = newClustersPalette,
                                         optional_parameters = optional_parameters)
  }
  if (length(newMetaclusters)>0) {
    cat("\ncomputing new cell cluster cellcount, abundance matrix...")
    cellcount = compute.cellcount(newMetaclusters, newCYTdata@samples)
    abundance = compute.abundance(cellcount)
    optional_parameters = list("merging" = TRUE)
    if (length(CYTdata1@Metaclustering@optional_parameters)>0) {
      optional_parameters = append(optional_parameters,
                                   list("CYTdata1" = CYTdata1@Metaclustering@optional_parameters))
    }
    if (length(CYTdata2@Metaclustering@optional_parameters)>0) {
      optional_parameters = append(optional_parameters,
                                   list("CYTdata2" = CYTdata2@Metaclustering@optional_parameters))
    }
    cat("\n\n - Creating Metaclustering slot")
    newCYTdata@Metaclustering = methods::new("Metaclustering",
                                             metaclusters = newMetaclusters,
                                             cellcount = cellcount,
                                             abundance = abundance,
                                             palette = newMetaclustersPalette,
                                             optional_parameters = optional_parameters)
  }
  
  validObject(newCYTdata)
  return(newCYTdata)
}

