## script to merge gene modules from different WGCNA runs based on the overlap in genes

# Usage e.g.
# time Rscript /nfsdata/projects/jonatan/tools/wgcna-src/geneNetworkScripts/gene_module_merge2.R  --path_df_NWA /projects/jonatan/pub-perslab/18-liver-wgcna/tables/liver_perslab_int_wgcna1_cell_cluster_module_genes.csv.gz --colGeneWeights pkMs  --colGeneNames genes --minPropIntersect  0.75 --minWeightedCor  0.75 --minPropIntersect_iter_increase  1.01 --maxIter 25 --corMethod  pearson --dirOut /projects/jonatan/pub-perslab/18-liver-wgcna/ --prefixOut merge1 --RAMGbMax 250

# Algorithm:
# Do until k == maxIter or until there are no modules to merge:
# 1. for each module j, find the module i with which it has the highest gene overlap as 
#     a proportion of the total genes in j (number in [0:1])
# 2. if minPropIntersect and minWeightedCor conditions are met, 
#     before merging module j into i, check if module i itself might be 
#     merged into a third module m, and the intersecting genes between 
#     i and m are more correlated, weighted by the gene weights and 
#     the overall intersect size, than that of j and i

#TODO
# * renormalize gene weights?
# * compute correlations?

######################################################################
########################## FUNCTIONS #################################
######################################################################

library("optparse")


######################################################################
########################### OPTPARSE #################################
######################################################################

option_list <- list(
  make_option("--path_df_NWA", type="character",
              help = "path to dataframe from a gene network analysis run, in long format (i.e. one row per gene per celltype), containing 'cell_cluster', 'module', one or two gene name columns, and a column of numeric scores. I.e. one row per gene, e.g. rwgcna cell_cluster_module_genes.csv files"),  
  make_option("--colModule", type="character", default ="module",
              help = "nwa_df column with modules, [default %default]"),
  make_option("--colCellClust", type="character", default ="cell_cluster",
              help = "nwa_df column with modules, [default %default]"),
  make_option("--colGeneWeights", type="character", default ="pk.*M",
              help = "nwa_df column with gene weights, , [default %default]"),  
  make_option("--colGeneNames", type="character", default="genes",
              help ="string or regex for grepping nwa_df e.g. 'hgnc|symbol|gene_name' or 'ensembl', [default %default]"),
    make_option("--minPropIntersect", type="numeric", default = 0.7,
              help = "minimum proportion of module j intersecting with module i to discard j or merge j into i, [default %default]"),
  make_option("--minPropIntersect_iter_increase", type="numeric", default = 1.02,
              help = "At the end of each iteration, multiply the minPropIntersect threshold by some numeric constant >= 1 to help ensure convergence given that merging modules increases the probability that others will have a large overlap with them, [default %default]"),
  make_option("--minWeightedCor", type="numeric", default = 0.7,
              help = "if module j has at least minPropIntersect genes also in module i, set minimum correlation weighted by gene weight, to discard j or merge j into i, [default %default]"),
  make_option("--cellClusters_keep", type="character", default = NULL,
              help = "quoted vector of character, e.g. ''c('hepatocytes', 'stellate_cells')''. If provided, when comparing a module of this celltype with a module of another celltype, the script will always keep the former, overriding usual behaviour [default %default]"),
  make_option("--corMethod", type="character", default = "pearson",
              help = "pearson, spearman or kendall, [default %default]"),
  make_option("--mergeOrPrune", type="character", default = "prune",
              help = "'merge' to merge close modules or 'prune' to discard modules largely contained in others, [default %default]"),
  make_option("--maxIter", type="integer", default = 20,
              help = "maximum number of merge iterations , [default %default]"),
  make_option("--dirOut", type="character",
              help = "Outputs go to /tables and /RObjects subdirectories"),  
  make_option("--prefixOut", type="character", default = paste0(substr(gsub("-","",as.character(Sys.Date())),3,1000), "_", sample(x = 999, size = 1)),
              help = "Unique prefix for output files, [default %default]"),
  make_option("--doPlot", type="logical", default = F,
              help = "Whether or not to plot intersect and correlation matrices at each iteration, [default %default]"),
  make_option("--RAMGbMax", type="integer", default=250,
              help = "Upper limit on Gb RAM available. Taken into account when setting up parallel processes. [default %default]")
)


library("here")

######################################################################
######################### SOURCE FUNCTIONS ################################
######################################################################

source(here("perslab-sc-library", "utility_functions.R"))

######################################################################
############################ PARAMS ##################################
######################################################################

opt <- parse_args(OptionParser(option_list=option_list))

path_df_NWA <- opt$path_df_NWA 
colModule <- opt$colModule
colCellClust <- opt$colCellClust
colGeneWeights <- opt$colGeneWeights
colGeneNames <- opt$colGeneNames
maxIter <- opt$maxIter
dirOut <- opt$dirOut
prefixOut <- opt$prefixOut
corMethod <- opt$corMethod
mergeOrPrune = opt$mergeOrPrune
minPropIntersect <- opt$minPropIntersect
minPropIntersect_iter_increase <- opt$minPropIntersect_iter_increase
minWeightedCor <- opt$minWeightedCor
cellClusters_keep <- opt$cellClusters_keep
if (!is.null(cellClusters_keep)) cellClusters_keep <- eval(parse(text=cellClusters_keep))
doPlot <- opt$doPlot
RAMGbMax <- opt$RAMGbMax

######################################################################
########################### PACKAGES #################################
######################################################################

library("data.table")
library("magrittr")
library("Matrix")
library("parallel")
library("WGCNA")
library("readr")
if (doPlot) {
  library("corrplot")
  library("RColorBrewer")
}

######################################################################
############################# SET PARAMS #############################
######################################################################

options(stringsAsFactors = F, use="pairwise.complete.obs")

######################################################################
############################ CONSTANTS ###############################
######################################################################

# if specified output directory doesn't exist, create it 
if (!file.exists(dirOut)) {
  dir.create(dirOut) 
  message("dirOut not found, new one created")
}

dirPlots = paste0(dirOut,"plots/")
if (!file.exists(dirPlots)) dir.create(dirPlots) 

dirTables = paste0(dirOut,"tables/")
if (!file.exists(dirTables)) dir.create(dirTables)

dirRObjects = paste0(dirOut,"RObjects/")
if (!file.exists(dirRObjects)) dir.create(dirRObjects)

dirLog = paste0(dirOut,"log/")
if (!file.exists(dirLog)) dir.create(dirLog)

flagDate = substr(gsub("-","",as.character(Sys.Date())),3,1000)

randomSeed = 12345
set.seed(seed = randomSeed)

######################################################################
######################## LOAD MODULE DATA ############################
######################################################################

message("Loading gene module data")

dt_geneMod <- data.table::fread(path_df_NWA)

colGeneNames <- grep(colGeneNames, colnames(dt_geneMod), value=T)
colGeneWeights <- grep(colGeneWeights, colnames(dt_geneMod), value=T)#colnames(dt_geneMod)[which(sapply(X=colnames(dt_geneMod), FUN =function(colname) class(dt_geneMod[[colname]]))=="numeric")]

######################################################################
###################### INITIALISE MERGE LOOP #########################
######################################################################

dt_geneMod[["module_merged"]] <- dt_geneMod[[colModule]]
dt_geneMod[["cell_cluster_merged"]] <- dt_geneMod[[colCellClust]]

#dt_geneMod[[paste0(colGeneWeights, "_merged")]] <- dt_geneMod[[colGeneWeights]] 
#pathLog <- paste0(dirLog, prefixOut, "_mod_merge_log.txt")
iteration=1
modules_to_exclude <- c()

while (T) {
  
  message(paste0(prefixOut, " merge iteration ", iteration))  
  ######################################################################
  ########################### GET MODULES ##############################
  ######################################################################
  
  unique(dt_geneMod[["module_merged"]]) %>% sort %>% na.omit %>% as.character -> modules  
  
  modules <- modules[!modules %in% modules_to_exclude]
  
  modules <- modules[nchar(modules)>0]
  ######################################################################
  ######################## GET GENE WEIGHT VECTORS #####################
  ######################################################################
  # get a vector of gene weights, with gene as name, for each module
  
  message("Getting gene weight vectors")
  
  fun <- function(module) {
    dt_geneMod[[colGeneWeights]] %>% '['(dt_geneMod[["module_merged"]]==module) %>% na.omit %>% as.numeric -> geneWeights 
    dt_geneMod[[colGeneNames]] %>% '['(dt_geneMod[["module_merged"]]==module) %>% na.omit  %>% as.character -> names(geneWeights) 
    geneWeights <- geneWeights[nchar(names(geneWeights))>0]
    return(geneWeights)
  }
  
  list_iterable = list("X"=modules)
  list_geneWeights <- safeParallel(fun=fun,list_iterable=list_iterable)
  names(list_geneWeights) <- modules
  ######################################################################
  ############## Compute the mod j - mod i prop. overlap matrix ########
  ######################################################################
  
  # module * module matrix
  # for each column j, each row i gives the proportion of module j that overlaps with module i 
  
  message("Computing module-module gene intersect matrix")
  
  fun = function(geneWeights) {
    sapply(list_geneWeights, function(geneWeightsOther) {
      base::intersect(x=names(geneWeights), y=names(geneWeightsOther)) %>% length %>% '/'(length(geneWeights))
    }, simplify=T)
  }
  
  #list_iterable = list("X"=list_geneWeights)
  #mat_moduleGeneIntersect <- safeParallel(fun=fun,list_iterable=list_iterable, simplify = T)
  mat_moduleGeneIntersect <- sapply(FUN=fun,"X"=list_geneWeights, simplify = T)
  # set diagonal to zero
  mat_moduleGeneIntersect[col(mat_moduleGeneIntersect)==row(mat_moduleGeneIntersect)] <- 0
  
  ### Compute modules to exclude in the coming iterations ###
  # They have no big intersections with other modules as a prop of their genes
  logical_colNoBigIntersect <- apply(mat_moduleGeneIntersect, MARGIN=2, FUN=function(j) max(j) < 0.75*minPropIntersect)  
  # No other modules have a big intersection with them
  logical_rowNoBigIntersect <- apply(mat_moduleGeneIntersect, MARGIN=1, FUN=function(j) max(j) < 0.75*minPropIntersect)  
  logical_noBigIntersect <- logical_colNoBigIntersect | logical_rowNoBigIntersect 
  
  # exclude these from following iterations to save time
  modules_to_exclude <- c(modules_to_exclude, modules[!logical_noBigIntersect])
  
  ######################################################################
  ####################### Plot the overlap matrix ######################
  ######################################################################
  if (doPlot){
    message("Plotting gene intersect matrix")
    
    pdf(sprintf("%s%s_iter_%s_mod_intersect.pdf", dirPlots, prefixOut, iteration), 
        width=max(20,ncol(mat_moduleGeneIntersect) %/% 3),
        height=max(20,ncol(mat_moduleGeneIntersect) %/% 3))
    corrplot(corr = mat_moduleGeneIntersect, 
             method = "color",
             col = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")), bias = 1)(200),
             diag = F,
             is.corr=F,
             title=paste0(prefixOut, ": prop. of mod j intersecting mod i, iteration ", iteration),
             order = "original",#hclust",
             hclust.method = "average",
             addCoef.col = "black",
             tl.srt = 45,
             number.digits = 2L,
             number.cex = 0.5)
    
    invisible(dev.off())
  }
  ######################################################################
  ########### Compute the overlap weighted correlation matrix ##########
  ######################################################################
  
  message("Computing weighted intersect gene weight correlations for highly overlapping modules")
  
  mat_modIntersectCorrWeighted <- matrix(data=0, nrow=nrow(mat_moduleGeneIntersect), ncol=ncol(mat_moduleGeneIntersect))
  dimnames(mat_modIntersectCorrWeighted) <- dimnames(mat_moduleGeneIntersect)
  
  # traverse intersect matrix columns, find any where the max intersect is large
  idx_colBigIntersect <- apply(mat_moduleGeneIntersect, MARGIN=2, FUN=function(j) max(j) >= minPropIntersect)  %>% which
  # traverse columns, find the row that the intersect is large with
  idx_rowBigIntersect <- apply(mat_moduleGeneIntersect, MARGIN=2, which.max) %>% '['(idx_colBigIntersect)
   
  if (length(idx_colBigIntersect)==0) {
    message("no more modules to merge")
    break
  }
  
  # For columns and rows with big intersect, calculate the correlations
  for (k in 1:length(idx_colBigIntersect)) {
    j=idx_colBigIntersect[k]
    i=idx_rowBigIntersect[k]
    genesIntersect <- base::intersect(names(list_geneWeights[[j]]), names(list_geneWeights[[i]]))
    mat_modIntersectCorrWeighted[i,j] <- WGCNA::cor(x = list_geneWeights[[j]][match(genesIntersect, names(list_geneWeights[[j]]))], 
                                                     weights.x = list_geneWeights[[j]][match(genesIntersect, names(list_geneWeights[[j]]))],
                                                    # correlate the jth column
                                                     y= list_geneWeights[[i]][match(genesIntersect, names(list_geneWeights[[i]]))],
                                                     weights.y = list_geneWeights[[i]][match(genesIntersect, names(list_geneWeights[[i]]))],
                                                    # with the ith row
                                                     method = corMethod, 
                                                     quick=0.25, 
                                                     verbose=F, 
                                                     nThreads=30)
      
    
  }
  
  ######################################################################
  ############### Plot the overlap weighted correlation matrix #########
  ######################################################################
  if(doPlot) {
    message("Plotting gene intersect weighted correlation matrix")

    pdf(sprintf("%s%s_iter%s_mat_modIntersectCorrWeighted.pdf", dirPlots, prefixOut, iteration),
        width=max(20,ncol(mat_modIntersectCorrWeighted) %/% 3),
        height=max(20,ncol(mat_modIntersectCorrWeighted) %/% 3))
    corrplot(corr = mat_modIntersectCorrWeighted,
             method = "color",
             col = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")), bias = 1)(200),
             diag = F,
             is.corr=F,
             title=paste0(prefixOut, ": weighted correlation of mod j's intersection with mod i, iteration ", iteration),
             order = "original",#"hclust",
             #hclust.method = "average",
             addCoef.col = "black",
             tl.srt = 45,
             number.digits = 2L,
             number.cex = 0.5)

    invisible(dev.off())
  }
  ######################################################################
  ################### IDENTIFY MODULES TO MERGE ########################
  ######################################################################
  
  message("Identifying modules to merge")
  # 1. For each column in the intersection matrix, identify the maximum row value. If no values are over threshold, STOP
  
  # For each column j, find the index of the max row i of the intersect matrix, i!=j
  #idx_rowmax <- apply(mat_moduleGeneIntersect, MARGIN=2, which.max)
  
  # Find the corresponding values of the intersect
  #toprow_vals <- apply(mat_moduleGeneIntersect, MARGIN=2, max)
  
  # 2. For each maximum row value, check whether the column correspond to that row has any higher values.
  # Then check if the intersect meets the threshold
  # logical_colmerge <- toprow_vals >= sapply(idx_rowmax, function(j) {
  #   max(mat_moduleGeneIntersect[,j])
  # }) & toprow_vals >= minPropIntersect
  
  #   3. For each marked j:i pair, verify that the corresponding weighted correlation of the 
  # intersecting genes is > minWeightedCor. If not, do not merge
  # logical_highcorr <- sapply(which(logical_colmerge), function(j){
  #   mat_modIntersectCorrWeighted[idx_rowmax[j],j]>=minWeightedCor
  # })
  # 
  # logical_colmerge[which(logical_colmerge)][!logical_highcorr]
  
  logical2_colBigIntersectCorr <- sapply(1:length(idx_colBigIntersect), function(k){
    j = idx_colBigIntersect[k]
    i = idx_rowBigIntersect[k]
    # For each module pair j,i, before merging module j into i,
    # check if module i itself might be merged into a third module m, and 
    # the intersecting genes between i and m are more correlated, weighted
    # by the gene weights and the overall intersect size, than that of j and i
    
    idx_matchedMax <- which.max(mat_modIntersectCorrWeighted[,i])
    bool <- mat_modIntersectCorrWeighted[i,j]*mat_moduleGeneIntersect[i,j] >= max(mat_modIntersectCorrWeighted[idx_matchedMax,i])*mat_moduleGeneIntersect[idx_matchedMax,i]
    # is the intersect correlation sufficient?
    bool2 <- bool & mat_modIntersectCorrWeighted[i,j] >= minWeightedCor
  })

  if (!any(logical2_colBigIntersectCorr)) {
    message("no more modules to merge")
    break
  }
  
  ### Merge module j into i in dt_geneMod dataframe:
  
  # ii. for duplicated genes, add NA in dt_geneMod[["module_merged"]] for the j copy
  
  mods_to_merge <- modules[idx_colBigIntersect[logical2_colBigIntersectCorr]]
  mods_to_merge_into <- modules[idx_rowBigIntersect[logical2_colBigIntersectCorr]] 

  # Check whether to swap the two vectors to keep modules from favoured celltypes
  if (!is.null(cellClusters_keep)) {
    mods_to_merge_cellClusters = sapply(mods_to_merge, function(eachMod) {
      dt_geneMod[[colCellClust]][dt_geneMod[["module_merged"]]==eachMod][1]
      })
    
    if (any(mods_to_merge_cellClusters %in% cellClusters_keep)) {
      
      mods_to_merge_into_cellClusters = sapply(mods_to_merge_into, function(eachMod) {
        dt_geneMod[[colCellClust]][dt_geneMod[["module_merged"]]==eachMod][1]
      })
     
      # swap, where appropriate
      mods_to_merge <- ifelse(mods_to_merge_cellClusters %in% cellClusters_keep & 
                              !mods_to_merge_into_cellClusters %in% cellClusters_keep, mods_to_merge_into, mods_to_merge)
      mods_to_merge_into <- ifelse(mods_to_merge_cellClusters %in% cellClusters_keep & 
                              !mods_to_merge_into_cellClusters %in% cellClusters_keep, mods_to_merge, mods_to_merge_into)
      
    }
  }
  
  if (mergeOrPrune=="prune") mods_to_merge_into <- NULL 
  
  ######################################################################
  ###################### UPDATE dt_geneMod ##########################
  ######################################################################
  
  message("Updating module dataframe")
  
  for (i in 1:length(mods_to_merge)) {
    cell_cluster_to_merge_into <- if (mergeOrPrune=="merge") dt_geneMod[[colCellClust]][dt_geneMod[[colModule]]==mods_to_merge_into[i]][1] else NULL
    # module_merged column
    logical_to_merge <- dt_geneMod[["module_merged"]] %in% mods_to_merge[i] 
    logical_to_merge_into <- dt_geneMod[["module_merged"]] %in% mods_to_merge_into[i]

    dt_geneMod[["module_merged"]][logical_to_merge] <- if (mergeOrPrune=="merge") mods_to_merge_into[i] else NA
    dt_geneMod[["cell_cluster_merged"]][logical_to_merge] <- if (mergeOrPrune=="merge") cell_cluster_to_merge_into else NA
    
    if (mergeOrPrune=="merge") {
      # in module_merged column set duplicate genes in merged module to NA_character 
      logical_dup <- duplicated(dt_geneMod[["hgnc"]]) 
      logical_dup_merged <- logical_dup & logical_to_merge
      
      #logical2_isna <- is.na(dt_geneMod[["module_merged"]][logical_to_merge | logical_to_merge_into][logical2_dup])
      dt_geneMod[["module_merged"]][logical_dup_merged] <- NA_character_
      dt_geneMod[["cell_cluster_merged"]][logical_dup_merged] <- NA_character_ 
    }
  }
  
  # write out log 
  log_entry <- paste0(prefixOut, " merge iteration ", iteration, " with minPropIntersect = ", round(minPropIntersect,3), " and minWeightedCor = ", minWeightedCor)
  message(log_entry)
  #cat(log_entry, file = pathLog, append=T, sep = "\n")

  log_entry <- if (mergeOrPrune=="merge") { paste0("    merging ", mods_to_merge, " into ", mods_to_merge_into, collapse="\n") } else {
    paste0("    pruning ", mods_to_merge, " due to overlap with ", modules[idx_rowBigIntersect[logical2_colBigIntersectCorr]], collapse="\n") 
  }
  
  message(log_entry)
  #cat(log_entry, file = pathLog, append=T, sep = "\n")
  
  iteration = iteration+1
 
  if (iteration>maxIter) {
    message(paste0(maxIter, " iterations reached, stopping"))
  }
  
  #increment the minPropIntersect
  if (!is.null(minPropIntersect_iter_increase)) minPropIntersect <- min(1, minPropIntersect*minPropIntersect_iter_increase)
  
}

######################################################################
############################ WRAP UP #################################
######################################################################

data.table::fwrite(x = dt_geneMod, file= gzfile(paste0(gsub("\\.csv.gz", "", path_df_NWA), "_merged", ".csv.gz")), compress="gzip")

############################### FINISH ##################################

message("Script done!")
