# multi-WGCNA cross-dataset module preservation pipeline

# *inputs: 
#   * WGCNA outputs: cell_cluster_module_genes.csv dataframe
#   * one or more expression datasets in Seurat v2 or delim format
# *outputs:
#   * network connectivity preservation scores, value in [0:1]
#   * module preservation Z-scores and p-values


# usage: e.g.
## time Rscript /projects/jonatan/tools/wgcna-src/wgcna-toolbox/wgcna_preservation.R --path_dt_geneMod   "/projects/jonatan/tmp-epilepsy/tables/modMerge2_test_dt_geneMod_merged.csv.gz" --colGeneWeights  "pkIM" --vec_pathsDatExpr   'c("ep_ex" = "/projects/jonatan/tmp-epilepsy/RObjects/EP_ex_seurat_obj_filtered.RDS.gz")'  --vec_metadataIdentCols   'c("ep_ex" = "subtypesBoth")'  --list_list_vec_identLvlsMatch 'list("L6_Nr4a2"=list("ep_ex"=c("Exc_L5-6_THEMIS_DCSTAMP", "Exc_L5-6_THEMIS_CRABP1", "Exc_L5-6_THEMIS_FGF10")), "L3_Prss12"=list("ep_ex"=c("Exc_L3-4_RORB_CARM1P1")), "L5_Grin3a"=list("ep_ex"=c("Exc_L4-5_RORB_FOLH1B","Exc_L4-5_RORB_DAPK2", "Exc_L4-6_RORB_SEMA3E")), "L6_Syn3"=list("ep_ex"=c("Exc_L5-6_THEMIS_C1QL3")), "L2_3_Cux2"=list("ep_ex"=c("Exc_L2-3_LINC00507_FREM3")), "L2_Lamp5"=list("ep_ex"=c("Exc_L2-4_LINC00507_GLP2R","Exc_L2_LAMP5_LTK")), "L4_Rorb"=list("ep_ex"=c("Exc_L3-5_RORB_TWIST2", "Exc_L3-5_RORB_ESR1", "Exc_L3-5_RORB_COL22A1", "Exc_L3-5_RORB_FILIP1L")), "L6_tle4"=list("ep_ex"=c("Exc_L6_FEZF2_OR2T8","Exc_L6_FEZF2_SCUBE1", "Exc_L5-6_FEZF2_EFTUD1P1","Exc_L5-6_SLC17A7_IL15", "Exc_L5-6_FEZF2_ABO")), "L5_Htr2c"=list("ep_ex"=c("Exc_L4-6_FEZF2_IL26")))' --minGeneClusterSize 10 --minCellClusterSize 50 --colGeneNames "hgnc|symbol|gene_name_optimal" --dirOut    "/projects/jonatan/tmp-epilepsy/" --prefixOut    "ep_ex_4_multi_4" --dirTmp    "/scratch/tmp-wgcna/" --networkType  "signed hybrid" --corFnc    "cor" --dataOrganism    "hsapiens" --scaleCenterRegress T  --colMod module_merged --colCellClust cell_cluster_merged

######################################################################
########################## DEFINE FUNCTIONS ##########################
######################################################################

library("here")

source(here("perslab-sc-library", "utility_functions.R"))

######################################################################
########################### PACKAGES #################################
######################################################################

suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("WGCNA"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("Matrix"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("openxlsx"))
######################################################################
########################### OPTPARSE #################################
######################################################################

option_list <- list(
  make_option("--path_dt_geneMod", type="character",
              help = "path to dataframe from a gene network analysis run, in long format (i.e. one row per gene per celltype), containing 'cell_cluster', 'module', one or two gene name columns, and a column of numeric scores. I.e. one row per gene, e.g. rwgcna cell_cluster_module_genes.csv files"),  
  make_option("--colGeneNames", type="character", default="hgnc|symbol|gene_name_optimal",
              help ="WGCNA geneMod table column with gene names [default %default]"),
  make_option("--colGeneWeights", type="character",
              help = "WGCNA geneMod table column with  gene weights"), 
  make_option("--colMod", type="character", default="module_merged",
              help = "WGCNA geneMod table column with modules [default %default]"), 
  make_option("--colCellClust", type="character", default="cell_cluster_merged",
              help =" WGCNA geneMod table column with module cell cluster of origin, e.g. 'cell_cluster' or 'cell_cluster_merged', [default %default]"),
  make_option("--vec_pathsDatExpr", type="character",
              help = "Quoted vector of named strings, paths to normalized expression data in delimited text form (readable by data.table::fread) with gene names in the first column"),  
  make_option("--vec_pathsMetadata", type="character", 
              help = "Quoted vector of full paths to metadata in one of the standard (compressed) character separated formats. Should correspond to files in vec_pathsDatExpr. [default %default]"),  
  make_option("--vec_metadataIdentCols", type="character",
              help = "vector of characters to identify columns in metadata by which to split the expression data, named by the datExpr names, e.g. ''c(datExpr1='celltype', datExpr2='clust.res.1')''."),
  make_option("--list_list_vec_identLvlsMatch", type="character", default=NULL,
              help="a quoted list of named lists, named by 'reference' celltype, of vectors of components designating 'test' celltype(s), e.g. ''list('L6_Nr4a2'=list('ep_ex'=c('Exc_L5-6_THEMIS_DCSTAMP', 'Exc_L5-6_THEMIS_CRABP1', 'Exc_L5-6_THEMIS_FGF10')), 'L3_Prss12'=list('ep_ex'=c('Exc_L3-4_RORB_CARM1P1')), 'L5_Grin3a'=list('ep_ex'=c('Exc_L4-5_RORB_FOLH1B','Exc_L4-5_RORB_DAPK2', 'Exc_L4-6_RORB_SEMA3E')))''. If the argument is left as NULL, will match each level to all others. [default %default]"),       #help = "quoted list of named character vectors with named components. List names are labelling levels (e.g. 'tissue', 'celltype'), vector names are column in datExpr metadata, vector values are regex to match levels. Use NA in a vector to skip a dataset for that labelling. E.g. ''list('tissue' = c(tissue='.*', NA), 'sub_celltype'=c(tissue_cell_type = '.*', ClusterName = '.*'))''"),
  make_option("--minCellClusterSize", type="integer", default=50L,
              help="Minimum number of cells in a cell_cluster, integer, [default %default]."),
  make_option("--minpropModGenesInTestLvl", type="numeric", default=0.5,
              help="Minimum proportion of module genes in the test expression dataset [default %default]."),
  make_option("--minPropModGenesNon0inTestLvl", type="numeric", default=0.1,
              help="Minimum number of module genes not uniformly zero in test expression data subset (genes not in test level are counted as zeros). Numeric, [default %default]."),
  # make_option("--minGeneClusterSize", type="integer", default=10L,
  #             help="Minimum number of genes in a module, integer, [default %default]."),
   make_option("--dirOut", type="character",
              help = "Outputs go to /tables and /RObjects subdirectories"),  
  make_option("--prefixOut", type="character", default = paste0(substr(gsub("-","",as.character(Sys.Date())),3,1000), "_", sample(x = 999, size = 1)),
              help = "Unique prefix for output files, [default %default]"),
  make_option("--dirTmp", type="character", default = "/scratch/rkm916/tmp-wgcna/",
              help = "Outputs go to /tables and /RObjects subdirectories"),
  make_option("--networkType", type="character", default = "c('signed hybrid')",
              help = "for WGCNA modulePreservation function: one of 'signed', 'unsigned', ''c('signed hybrid')'',  [default %default]"),
  make_option("--corFnc", type="character", default = "cor",
              help = "for WGCNA modulePreservation function, 'cor' or 'bicor'"),
  make_option("--RAMGbMax", type="integer", default=250L,
              help = "Upper limit on Gb RAM available. Taken into account when setting up parallel processes. [default %default]"),
  make_option("--path_runLog", type="character", default=NULL,
              help = "Path to file to log the run and the git commit. If left as NULL, write to a file called runLog.text in the dirLog [default %default]")
)

######################################################################
########################### GET OPTIONS ##############################
######################################################################

opt <- parse_args(OptionParser(option_list=option_list))

path_dt_geneMod <- opt$path_dt_geneMod 
colGeneWeights <- opt$colGeneWeights
colMod <- opt$colMod
colCellClust <- opt$colCellClust
vec_pathsDatExpr <- eval(parse(text=opt$vec_pathsDatExpr))
vec_pathsMetadata <- eval(parse(text=opt$vec_pathsMetadata)) 
vec_metadataIdentCols <- eval(parse(text=opt$vec_metadataIdentCols)) # cannot be null
list_list_vec_identLvlsMatch <- opt$list_list_vec_identLvlsMatch
if (!is.null(list_list_vec_identLvlsMatch)) list_list_vec_identLvlsMatch <- eval(parse(text=list_list_vec_identLvlsMatch))
minCellClusterSize <- opt$minCellClusterSize
minpropModGenesInTestLvl <- opt$minpropModGenesInTestLvl
minPropModGenesNon0inTestLvl <- opt$minPropModGenesNon0inTestLvl
#minGeneClusterSize <- opt$minGeneClusterSize
colGeneNames <- opt$colGeneNames
dirOut <- opt$dirOut
prefixOut <- opt$prefixOut
dirTmp <- opt$dirTmp
networkType <- opt$networkType
if (grepl("hybrid", networkType)) networkType <- eval(parse(text=opt$networkType))
corFnc <- opt$corFnc
RAMGbMax <- opt$RAMGbMax
path_runLog <- opt$path_runLog

######################################################################
############################# SET PARAMS #############################
######################################################################

options(stringsAsFactors = F, 
        use="pairwise.complete.obs")

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
set.seed(randomSeed)

######################################################################
######################## LOAD MODULE DATA ############################
######################################################################

message("Loading data")

dt_geneMod <- fread(path_dt_geneMod)

list_metadata <- lapply(vec_pathsMetadata, fread)

list_datExpr <- lapply(vec_pathsDatExpr, function(pathDatExpr){
  dt_datExpr <- fread(pathDatExpr)
  if (any(duplicated(dt_datExpr[[1]]))) stop("datExpr has duplicate feature names, please ensure they are unique")
  datExpr <- as.matrix(dt_datExpr[,-1]) 
  rownames(datExpr) <- dt_datExpr[[1]]
  datExpr <- t(datExpr)
  return(datExpr)
})

names(list_datExpr) <- names(vec_pathsDatExpr)

######################################################################
####################### CHECK WHETHER REF LEVELS EXIST ###############
######################################################################

# check if the specified reference levels exist in the gene module df

if (!is.null(list_list_vec_identLvlsMatch)) {
  if (!all(names(list_list_vec_identLvlsMatch) %in% dt_geneMod[[colCellClust]]))  {
    missing = names(list_list_vec_identLvlsMatch)[!names(list_list_vec_identLvlsMatch) %in% dt_geneMod[[colCellClust]]] 
    if (length(missing)>0) stop(paste0(missing, " not found in gene module table ", collapse=" "))
  }
}

######################################################################
################### IF MATCHES NOT SPECIFIED, MATCH ALL TO ALL #######
######################################################################
#a quoted list of named lists, named by 'reference' celltype, of vectors of components 
# designating 'test' celltype(s). All , e.g. ''list('datExpr2'=list('n01'=c('microgl1',migrogl2','migrogl3'), 
# s02'=c('oligo1', 'oligo2','oligo3')), 'datExpr3'=list('microglia'=c('s01','s02'), 'oligo'=c('s03','s04')))''
if (is.null(list_list_vec_identLvlsMatch)){
  list_list_vec_identLvlsMatch <- list()
  for (lvlRef in sort(unique(dt_geneMod[[colCellClust]]))){ 
    if (nchar(lvlRef)>0 & !is.na(lvlRef)) {
      list_vec_identLvlsMatch <- list()
      for (datExprName in names(list_datExpr)) {    
        identCol <- vec_metadataIdentCols[[datExprName]]
        list_vec_identLvlsMatch[[datExprName]] <- sort(unique(list_metadata[[datExprName]][[identCol]]))
      }
    #list_list_vec_identLvlsMatch[[datExprName]]
      list_list_vec_identLvlsMatch[[lvlRef]] <- list_vec_identLvlsMatch
    }
  }
}

######################################################################
######### VERIFY THAT lvls_ref AND lvls_test EXIST IN METADATA #######
######################################################################

# check if the specified reference levels exist in the metadata annotation

if (!all(names(list_list_vec_identLvlsMatch) %in% list_metadata[[1]][[vec_metadataIdentCols[1]]]))  {
  missing = names(list_list_vec_identLvlsMatch)[!names(list_list_vec_identLvlsMatch) %in% list_metadata[[1]][[vec_metadataIdentCols[1]]]] 
  stop(paste0(missing, " not found in ", names(list_metadata)[1], " metadata", collapse=" "))
}

# check if the specified test levels exist in the metadata annotation

for (lvlRef in names(list_list_vec_identLvlsMatch)) {
  for (datExprTestName in names(list_list_vec_identLvlsMatch[[lvlRef]])) {
    if(!all(list_list_vec_identLvlsMatch[[lvlRef]][[datExprTestName]] %in% list_metadata[[datExprTestName]][[vec_metadataIdentCols[datExprTestName]]])) {
      missing = list_list_vec_identLvlsMatch[[lvlRef]][[datExprTestName]][!list_list_vec_identLvlsMatch[[lvlRef]][[datExprTestName]] %in% list_metadata[[datExprTestName]][vec_metadataIdentCols[datExprTestName]]]
      stop(paste0(missing, " not found in ", datExprTestName, " metadata", collapse=" "))
    }
  }
}


######################################################################
############################ INITIALISE ##############################
######################################################################

lvlRef_datExprTest_presNwOut <- lvlRef_datExprTest_presModsOut <- NULL

setwd(dir = dirTmp)

######################################################################
############### RUN CELLTYPE-CELLTYPE PRESERVATION ANALYSIS ##########
######################################################################

# loop over reference levels, i.e. reference cell clusters where WGCNA found modules

#for (lvlRef in names(list_list_vec_identLvlsMatch)) {

if (F) {
  fun1 = function(lvlRef)  {
  
    datExprTest_presNwOut <- list()
    
    ####################################################
    #### Get gene network coloring for cell_cluster ####
    ####################################################
    
    #idxDuplicateGenes <- duplicated(dt_geneMod[[colGeneNames]][dt_geneMod[[colCellClust]] == lvlRef])
    coloring <- dt_geneMod[[colMod]][dt_geneMod[[colCellClust]] == lvlRef]#[!idxDuplicateGenes]
    names(coloring) <- dt_geneMod[[colGeneNames]][dt_geneMod[[colCellClust]] == lvlRef]#[!idxDuplicateGenes]
    coloring <- coloring[!is.na(coloring)] 
    
    ####################################################
    ############ Prepare reference dataset #############
    ####################################################
    
    # Get metadata df
    metadataRef <- list_metadata[[1]]
    metadataColnameRef <- vec_metadataIdentCols[1]
    
    idxLvlRefCells <- as.integer(na.omit(match(metadataRef[[1]][metadataRef[[metadataColnameRef]] == lvlRef], 
                                                rownames(list_datExpr[[1]]))))
    
    if (length(idxLvlRefCells) < minCellClusterSize) {
      datExprTest_presNwOut <- NA_character_
      message(paste0(lvlRef, " skipped because it has < ", minCellClusterSize, " cells"))
      next
    } 
    
    idxGenes <- as.integer(na.omit(match(names(coloring), colnames(list_datExpr[[1]]))))
    datExprRefLvl <- list_datExpr[[1]][idxLvlRefCells, idxGenes]
    
    ######################################################################
    #### Filter out coloring genes which are missing in datExpr_ref ######
    ######################################################################
    
    coloring <- coloring[names(coloring) %in% colnames(datExprRefLvl)]
    
    ####################################################
    ############ Filter out small modules ##############
    ####################################################
    
    # modsTooSmall <- names(table(coloring))[table(coloring)<minGeneClusterSize]
    # coloring <- coloring[!coloring %in% modsTooSmall]
    
    ####################################################
    ############## Prepare test datasets ###############
    ####################################################
    
    for (datExprTestName in as.character(names(list_list_vec_identLvlsMatch[[lvlRef]]))) {
  
      message(paste0("Reference: ", lvlRef, ", testing preservation in subsets of ", datExprTestName))
      
      metadataTest <- list_metadata[[datExprTestName]]
      metadataColnameTest <- vec_metadataIdentCols[datExprTestName]
      
      idxGenes <- as.integer(na.omit(match(names(coloring), colnames(list_datExpr[[datExprTestName]]))))
        
      coloring_f = coloring[names(coloring) %in% colnames(list_datExpr[[datExprTestName]])]
      # # filter out modules missing completely in datExprTest
      # mods_absent <- names(table(coloring))[!names(table(coloring)) %in% names(table(coloring[names(coloring) %in% colnames(list_datExpr[[datExprTestName]])[idxGenes]]))]
      # coloring_f <- coloring
      # coloring_f[coloring %in% mods_absent] <- "grey"
      # 
      # # Filter out modules where fewer than half the genes present in datExprTest
      # vec_modPropGenesInTest <- table(coloring)/table(coloring[coloring %in% names(table(coloring_f))])
      # vec_modPropGenesInTestTooSmall<- names(vec_modPropGenesInTest[vec_modPropGenesInTest<0.5])
      # coloring_f <- coloring_f[!coloring_f %in% vec_modPropGenesInTestTooSmall]
  
      list_datExprTest <- list()
      
      # subset the test dataset by annotation level; the genes are the same
      for (lvlTest in as.character(list_list_vec_identLvlsMatch[[lvlRef]][[datExprTestName]])) {
        
        idxCells <- as.integer(na.omit(match(metadataTest[[1]][metadataTest[[metadataColnameTest]] == lvlTest], 
                           rownames(list_datExpr[[datExprTestName]]))))
        datExprLvlTest <- list_datExpr[[datExprTestName]][idxCells, idxGenes]
        
        # filter reference dataset so both dataset have same genes. (the genes are the same in all test lvl submatrices)
        datExprRefLvl = datExprRefLvl[,colnames(datExprRefLvl) %in% colnames(datExprLvlTest)]
    
        ######################################################################
        ################# CHECK SUBSETS FOR GOOD GENES & SAMPLES #############
        ######################################################################
        
        # We will record this in the module preservation loop
        
        vec_logicalGoodGenes <- goodSamplesGenes(datExprLvlTest)$goodGenes
        vec_logicalGoodSamples <- goodSamplesGenes(datExprLvlTest)$goodSamples
        
        # names(vec_logicalGoodGenes) <- colnames(datExprLvlTest)
        # names(vec_logicalGoodSamples) <- rownames(datExprLvlTest)
        # 
        # datExprTest_testLvlGoodSamples[[datExprTestName]][[lvlTest]] <- vec_logicalGoodSamples
        # datExprTest_testLvlGoodGenes[[datExprTestName]][[lvlTest]] <- datExprLvlTest
        # 
        datExprLvlTest <- datExprLvlTest[vec_logicalGoodSamples,vec_logicalGoodGenes]
        
        ######################################################################
        ################ Check that test subset has sufficient cells #########
        ######################################################################
        
        if (nrow(datExprLvlTest) < minCellClusterSize) {
          message(paste0(lvlRef, " preservation in ", lvlTest, " skipped because ", lvlTest, " has < ", minCellClusterSize, " cells"))
          next
        }
        list_datExprTest[[lvlTest]] <- datExprLvlTest
      }
      
      if (length(list_datExprTest)==0) {
        datExprTest_presNwOut[[datExprTestName]] <- NA_character_
        message(paste0(datExprTestName, " skipped because no test sets have < ", minCellClusterSize, " cells"))
        next
      } 
  
      ####################################################
      ################ Prepare multidata #################
      ####################################################
    
      # Prepare multidata
      multiData <- vector(mode="list",length = 1+length(list_datExprTest)) # 1 for reference network
      multiData[[1]] <- list("data"=datExprRefLvl)
      for (k in 1:(length(multiData)-1)){
        multiData[[k+1]] <- list("data"=list_datExprTest[[k]])
      }
      names(multiData)[1] <- lvlRef
      names(multiData)[2:length(multiData)] <- names(list_datExprTest)
      
      rm(list_datExprTest)
      
      ####################################################
      ####### Run preservationNetworkConnectivity ########
      ####################################################
  
      # This should test preservation in each test_lvl datExpr
      datExprTest_presNwOut[[datExprTestName]] <- 
        WGCNA::preservationNetworkConnectivity(multiExpr = multiData, 
                                               corFnc = corFnc, 
                                               corOptions = "use='pairwise.complete.obs'", 
                                               networkType = networkType, 
                                               power=8, 
                                               sampleLinks=T, 
                                               nLinks=50000, 
                                               blockSize=5000, 
                                               setSeed=randomSeed, 
                                               weightPower=2, 
                                               verbose=3, 
                                               indent=0)
      
      
    }
    return(datExprTest_presNwOut)
  }
  
  outfile = paste0(dirLog, prefixOut, "_networkPreservation_log.txt")
  list_iterable=list("X"=names(list_list_vec_identLvlsMatch))
  lvlRef_datExprTest_presNwOut <- NULL
  
  lvlRef_datExprTest_presNwOut <- tryCatch({
    safeParallel(fun=fun1, list_iterable = list_iterable, outfile=outfile)
    }, error=function(err) {
      warning(paste0("preservationNetworkConnectivity failed with the error ", err))
    })
  
  #lvlRef_datExprTest_presNwOut <- lapply(FUN=fun1,"X"=names(list_list_vec_identLvlsMatch))
  
  if (!is.null(lvlRef_datExprTest_presNwOut)) names(lvlRef_datExprTest_presNwOut) <- names(list_list_vec_identLvlsMatch)
}
######################################################################
################## RUN MODULE PRESERVATION ANALYSIS ##################
######################################################################

require("doParallel")

fun2 = function(lvlRef)  {
  
  datExprTest_presModsOut <- list()
  
  datExprTest_propModGenesInTestLvl <- list()
  datExprTest_propModGenesNon0inTestLvl <- list()
  datExprTest_testLvlGoodSamples <- list()
  datExprTest_testLvlGoodGenes <- list()
  
  ####################################################
  #### Get gene network coloring for cell_cluster ####
  ####################################################

  # get gene module allocations (coloring)
  coloring <- dt_geneMod[[colMod]][dt_geneMod[[colCellClust]] == lvlRef]#[!idxDuplicateGenes]
  names(coloring) <- dt_geneMod[[colGeneNames]][dt_geneMod[[colCellClust]] == lvlRef]#[!idxDuplicateGenes]
  coloring <- coloring[!is.na(coloring)]
  coloring <- coloring[!is.na(names(coloring))]
  coloring <- coloring[nchar(names(coloring))!=0]
  ####################################################
  ############ Prepare reference dataset #############
  ####################################################
  
  # Get metadata df for reference data
  metadataRef <- list_metadata[[1]]
  metadataColnameRef <- vec_metadataIdentCols[1]
  
  idxLvlRefCells <- as.integer(na.omit(match(metadataRef[[1]][(metadataRef[[metadataColnameRef]] == lvlRef)], 
                                             rownames(list_datExpr[[1]]))))
  

  # match 'coloring' genes to reference expression data column names (genes). All genes should be there and order should match!!
  idxGenes <- as.integer(na.omit(match(names(coloring), colnames(list_datExpr[[1]]))))
  datExprRefLvl <- list_datExpr[[1]][idxLvlRefCells, idxGenes]
  coloring <- coloring[names(coloring) %in% colnames(datExprRefLvl)]

  ####################################################
  ############## Prepare test datasets ###############
  ####################################################
  
  for (datExprTestName in as.character(names(list_list_vec_identLvlsMatch[[lvlRef]]))) {
    
    message(paste0("Reference: ", lvlRef, ", testing preservation in subsets of ", datExprTestName))
    
    metadataTest <- list_metadata[[datExprTestName]]
    metadataColnameTest <- vec_metadataIdentCols[datExprTestName]
    
    list_datExprTest <- list()
    
    datExprTest_propModGenesNon0inTestLvl[[datExprTestName]] <- list()
    datExprTest_propModGenesInTestLvl[[datExprTestName]] <- list()
    datExprTest_testLvlGoodSamples[[datExprTestName]] <- list()
    datExprTest_testLvlGoodGenes[[datExprTestName]] <- list()
    
    for (lvlTest in as.character(list_list_vec_identLvlsMatch[[lvlRef]][[datExprTestName]])) {

      idxCells <- as.integer(na.omit(match(metadataTest[[1]][(metadataTest[[metadataColnameTest]] == lvlTest)], 
                                           rownames(list_datExpr[[datExprTestName]]))))
      
      # this is going to be the same for each lvl within datExprTestName since the genes are the same
      # no need to subset by idxCells here
      datExprTest_propModGenesInTestLvl[[datExprTestName]][[lvlTest]] <- 
        sapply(unique(coloring), 
               function(module) {
                 sum(names(coloring)[coloring==module] %in% colnames(list_datExpr[[datExprTestName]])) / 
                   length(coloring[coloring==module])})
      
      datExprTest_propModGenesNon0inTestLvl[[datExprTestName]][[lvlTest]] <- 
        sapply(unique(coloring), 
                function(module) {
                  if (datExprTest_propModGenesInTestLvl[[datExprTestName]][[lvlTest]][module] < minpropModGenesInTestLvl) {
                    return(NA) # if not enough module genes are included in datExpr lvlRef sub 
                  } else {
                    # return the prop of mod genes with above-0 expression
                    logicalGenes = colnames(list_datExpr[[datExprTestName]]) %in% 
                      names(coloring[coloring==module])
                    return(
                      sum(colSums(list_datExpr[[datExprTestName]][idxCells, logicalGenes])>0) / 
                        length(coloring[coloring==module])) # this counts missing genes as zero
                  }}) 
      
      
      mods_absent = names(datExprTest_propModGenesNon0inTestLvl[[datExprTestName]][[lvlTest]]) %>% 
        '['(datExprTest_propModGenesInTestLvl[[datExprTestName]][[lvlTest]] < minpropModGenesInTestLvl |
              datExprTest_propModGenesNon0inTestLvl[[datExprTestName]][[lvlTest]] < minPropModGenesNon0inTestLvl) # will work even if one comparison is with NA
      
      coloring_f <- coloring
      coloring_f[coloring_f %in% mods_absent] <- "grey"
      
      if (all(coloring_f=="grey")) {
        message(paste0(lvlRef, " preservation in ", lvlTest, " skipped because no modules had sufficient number of genes and/or non-zero counts in ", lvlTest))
        next
      }
        
      idxGenes <- as.integer(na.omit(match(names(coloring_f), colnames(list_datExpr[[datExprTestName]]))))
      
      datExprLvlTest <- list_datExpr[[datExprTestName]][idxCells, idxGenes]
      
      ######################################################################
      ####################### Filter genes and samples #####################
      ######################################################################
      
      vec_logicalGoodGenes <- goodSamplesGenes(datExprLvlTest)$goodGenes
      vec_logicalGoodSamples <- goodSamplesGenes(datExprLvlTest)$goodSamples
      
      names(vec_logicalGoodGenes) <- colnames(datExprLvlTest)
      names(vec_logicalGoodSamples) <- rownames(datExprLvlTest)
        
      datExprTest_testLvlGoodSamples[[datExprTestName]][[lvlTest]] <- vec_logicalGoodSamples
      datExprTest_testLvlGoodGenes[[datExprTestName]][[lvlTest]] <- vec_logicalGoodGenes
      
      datExprLvlTest <- datExprLvlTest[vec_logicalGoodSamples,vec_logicalGoodGenes]
        
      ######################################################################
      ################ Check that test subset has sufficient cells #########
      ######################################################################
      
      if (nrow(datExprLvlTest) < minCellClusterSize) {
        message(paste0(lvlRef, " preservation in ", lvlTest, " skipped because ", lvlTest, " has < ", minCellClusterSize, " cells"))
        next
      }
      list_datExprTest[[lvlTest]] <- datExprLvlTest
    }
    
    if (length(list_datExprTest)==0) {
      datExprTest_presModsOut[[datExprTestName]] <- NA_character_
                             message(paste0(datExprTestName, " skipped because no test sets have < ", minCellClusterSize, " cells"))
                             next
    }
    
    ####################################################
    ################ Prepare multidata #################
    ####################################################
    
    # Prepare multidata
    multiData <- vector(mode="list",length = 1+length(list_datExprTest)) # 1 for reference network
    multiData[[1]] <- list("data"=datExprRefLvl)
    for (k in 1:(length(multiData)-1)){
      multiData[[k+1]] <- list("data"=list_datExprTest[[k]])
    }
    names(multiData)[1] <- lvlRef
    names(multiData)[2:length(multiData)] <- names(list_datExprTest)
    
    rm(list_datExprTest)

   ####################################################
   ############### Prepare multicolor #################
   ####################################################
   
   multiColor <- list(coloring_f)
   names(multiColor) <- lvlRef
   
   # Prepare for parallel computation (WGCNA multi threads)

   additionalGb = max(as.numeric(sapply(multiData, FUN = function(x) object.size(x), simplify = T)))/1024^3
   objSizeGb <- as.numeric(sum(sapply(ls(envir = .GlobalEnv), function(x) object.size(x=eval(parse(text=x)))))) / 1024^3
   nCores <- max(1, min(detectCores() %/% 3, RAMGbMax %/% (objSizeGb + additionalGb))-1)
   nCores <- min(nCores, 20)

   #disableWGCNAThreads()
   enableWGCNAThreads(nThreads = nCores)
   
   datExprTest_presModsOut[[datExprTestName]] <- tryCatch({
     modulePreservation(multiData=multiData,
                        multiColor=multiColor,
                        dataIsExpr = TRUE,
                        networkType = networkType, 
                        corFnc = corFnc,
                        corOptions = "use='pairwise.complete.obs'", NULL,#list(use="pairwise.complete.obs"),#NULL, #if (corFnc == "cor") list(use = 'p') else NULL,
                        referenceNetworks = 1, 
                        testNetworks = NULL,
                        nPermutations = 100, 
                        includekMEallInSummary = FALSE,
                        restrictSummaryForGeneralNetworks = TRUE,
                        calculateQvalue = FALSE,
                        randomSeed = randomSeed, 
                        maxGoldModuleSize = 1000, 
                        maxModuleSize = 1000, 
                        quickCor = 0, # set to one to speed up/risk errors and inaccuracies 
                        ccTupletSize = 2, 
                        calculateCor.kIMall = FALSE,
                        calculateClusterCoeff = T,
                        useInterpolation = FALSE, 
                        checkData = TRUE, 
                        greyName = NULL, 
                        savePermutedStatistics = TRUE, 
                        loadPermutedStatistics = FALSE, 
                        permutedStatisticsFile = "permutedStats-actualModules.RData",#if (useInterpolation) "permutedStats-intrModules.RData" else "permutedStats-actualModules.RData", 
                        plotInterpolation = TRUE, 
                        interpolationPlotFile = "modulePreservationInterpolationPlots.pdf", 
                        discardInvalidOutput = TRUE,
                        parallelCalculation = T,#TRUE,
                        verbose = 3, 
                        indent = 0)}, error = function(err) {
                          msg = paste0(lvlRef, ": modulePreservation failed with the following error: ", err, " - attempting with serial computation")
                          message(msg)
                          write(x = msg, file = outfile,
                                ncolumns = 1, append = T)
                          tryCatch({
                          modulePreservation(multiData=multiData,
                                            multiColor=multiColor,
                                            dataIsExpr = TRUE,
                                            networkType = networkType, 
                                            corFnc = corFnc,
                                            corOptions = "use='pairwise.complete.obs'", NULL,#list(use="pairwise.complete.obs"),#NULL, #if (corFnc == "cor") list(use = 'p') else NULL,
                                            referenceNetworks = 1, 
                                            testNetworks = NULL,
                                            nPermutations = 100, 
                                            includekMEallInSummary = FALSE,
                                            restrictSummaryForGeneralNetworks = TRUE,
                                            calculateQvalue = FALSE,
                                            randomSeed = randomSeed*2, 
                                            maxGoldModuleSize = 1000, 
                                            maxModuleSize = 1000, 
                                            quickCor = 0, # set to one to speed up/risk errors and inaccuracies 
                                            ccTupletSize = 2, 
                                            calculateCor.kIMall = FALSE,
                                            calculateClusterCoeff = T,
                                            useInterpolation = FALSE, 
                                            checkData = TRUE, 
                                            greyName = NULL, 
                                            savePermutedStatistics = TRUE, 
                                            loadPermutedStatistics = FALSE, 
                                            permutedStatisticsFile = "permutedStats-actualModules.RData",#if (useInterpolation) "permutedStats-intrModules.RData" else "permutedStats-actualModules.RData", 
                                            plotInterpolation = TRUE, 
                                            interpolationPlotFile = "modulePreservationInterpolationPlots.pdf", 
                                            discardInvalidOutput = TRUE,
                                            parallelCalculation = F,#TRUE,
                                            verbose = 3, 
                                            indent = 0)
                          }, error = function(err) {
                            msg = paste0(lvlRef, ": modulePreservation failed a second time the following error: ", err, " - returning NA_character_")
                            message(msg)
                            write(x = msg, file = outfile,
                                  ncolumns = 1, append = T)
                          return(NA_character_)
                          })
                        })
                           
  }
  return(list("propModGenesInTestLvl" = datExprTest_propModGenesInTestLvl,
              "propModGenesNon0inTestLvl" = datExprTest_propModGenesNon0inTestLvl,
              "goodSamples" = datExprTest_testLvlGoodSamples,
              "goodGenes" = datExprTest_testLvlGoodGenes,
              "presModsOut"=datExprTest_presModsOut))
}

outfile = paste0(dirLog, prefixOut, "_modulePreservation_log.txt")
#args=list("X"=names(list_list_vec_identLvlsMatch))
#lvlRef_datExprTest_presModsOut <- safeParallel(fun=fun2, args=args, outfile=outfile)
list_presModsOut <- lapply(FUN=fun2,"X"=names(list_list_vec_identLvlsMatch))

names(list_presModsOut) <- names(list_list_vec_identLvlsMatch)

datExprTest_testLvlGoodSamples = list()
datExprTest_testLvlGoodGenes <- list()
lvlRef_datExprTest_propModGenesInTestLvl <- list()
lvlRef_datExprTest_propModgenesNon0inTestLvl <- list()
lvlRef_datExprTest_presModsOut <- list()

for (refLvl in names(list_presModsOut)) {
  datExprTest_testLvlGoodSamples[[refLvl]] <- list_presModsOut[[refLvl]][["goodSamples"]]
  datExprTest_testLvlGoodGenes[[refLvl]] <- list_presModsOut[[refLvl]][["goodGenes"]]
  lvlRef_datExprTest_propModGenesInTestLvl[[refLvl]] <- list_presModsOut[[refLvl]][["propModGenesInTestLvl"]]
  lvlRef_datExprTest_propModgenesNon0inTestLvl[[refLvl]] <- list_presModsOut[[refLvl]][["propModGenesNon0inTestLvl"]] 
  lvlRef_datExprTest_presModsOut[[refLvl]] <- list_presModsOut[[refLvl]][["presModsOut"]]
}

######################################################################
################# REFORMAT MODULE PRESERVATION SCORES ################
######################################################################

lvlRef_df_datExprTestStats <- lapply(X=names(list_list_vec_identLvlsMatch), FUN=function(lvlRef) {
  
  if (all(is.na(lvlRef_datExprTest_presModsOut[[lvlRef]]))) return(NA_character_)
  
  list_datExprTestStats <- lapply(names(list_list_vec_identLvlsMatch[[lvlRef]]), function(datExprTestName) {
    
    if (all(is.na(lvlRef_datExprTest_presModsOut[[lvlRef]][[datExprTestName]])) |
        all(is.null(lvlRef_datExprTest_presModsOut[[lvlRef]][[datExprTestName]]))) return(NA_character_)

    modulePreservationOut = lvlRef_datExprTest_presModsOut[[lvlRef]][[datExprTestName]]
    # we should get a matrix with a column for each test set
    list_df_Zsummary.pres <- lapply("X"=names(modulePreservationOut[["preservation"]][["Z"]][[1]][-1]), 
                      FUN=function(lvlTest) {
                        df_presStats <- modulePreservationOut[["preservation"]][["Z"]][[1]][-1][[lvlTest]]                 
                        if(!all(is.na(df_presStats))) {
                          lvlTestCrop <- gsub("inColumnsAlsoPresentIn\\.", "", lvlTest)
                          df_Zsummary.pres <- data.frame(
                            "ref_dataset"=names(vec_pathsDatExpr)[1],
                            "ref_lvl"=lvlRef,
                            "test_dataset"=datExprTestName,
                            "test_lvl"=rep(lvlTestCrop, times= nrow(df_presStats)), 
                            "module" = rownames(df_presStats), 
                            "Z_stat"=df_presStats[["Zsummary.pres"]])
                 
                          return(df_Zsummary.pres)
                          } else {
                            return(NA_character_)
                          }
    })

    list_df_Zsummary.pres <- list_df_Zsummary.pres[!sapply(list_df_Zsummary.pres, function(df) all(is.na(df)))]
    df_Zsummary.pres <- if (length(list_df_Zsummary.pres)==1) list_df_Zsummary.pres[[1]] else Reduce(f=rbind, x=list_df_Zsummary.pres)
    
    #Z_stats_long <- reshape2::melt(data=Z_stats)
    #colnames(Z_stats_long) <- c("module","test_lvl","Zsummary.pres")
    # we should get a matrix with a column for each test set
    list_df_pBonf.pres <-lapply(names(modulePreservationOut[["preservation"]][["log.pBonf"]][[1]][-1]), FUN=function(lvlTest) {
      df_presStats <- modulePreservationOut[["preservation"]][["log.pBonf"]][[1]][-1][[lvlTest]] 
      if(!all(is.na(df_presStats))) {
        lvlTestCrop <- gsub("inColumnsAlsoPresentIn\\.", "", lvlTest)
        df_pBonf.pres <- data.frame(
          "ref_dataset"=names(vec_pathsDatExpr)[1],
          "ref_lvl"=lvlRef,
          "test_dataset"=datExprTestName,
          "test_lvl"=rep(lvlTestCrop, times= nrow(df_presStats)), 
          "module" = rownames(df_presStats), 
          "log.p.Bonfsummary.pres"=df_presStats[["log.p.Bonfsummary.pres"]])
        return(df_pBonf.pres)
      } else {
        return(NA_character_)
      }
    })
  
    list_df_pBonf.pres <- list_df_pBonf.pres[!sapply(list_df_pBonf.pres, function(df) all(is.na(df)))]
    df_pBonf.pres <- if (length(list_df_pBonf.pres)==1) list_df_pBonf.pres[[1]] else Reduce(f=rbind, x=list_df_pBonf.pres)
    
    # colnames(logpBonf_stats) <- gsub("inColumnsAlsoPresentIn\\.","",colnames(logpBonf_stats))
    # logpBonf_stats_long <- reshape2::melt(logpBonf_stats)[-c(1,2)]
    # colnames(logpBonf_stats_long) <- c("log.p.Bonfsummary.pres")

    df_presStats <- dplyr::full_join(df_Zsummary.pres, df_pBonf.pres, by=c("ref_dataset", "ref_lvl", "test_dataset", "test_lvl", "module"))
      
      
      # data.frame("ref_dataset" = names(vec_pathsDatExpr)[1],
      #                        "ref_lvl"= rep(lvlRef, times=nrow(Z_stats_long)),
      #                        "test_dataset"= rep(datExprTestName, times=nrow(Z_stats_long)),
      #                        Z_stats_long,
      #                        logpBonf_stats_long)
    return(df_presStats)
  })
  df_datExprTestStats <- if (length(list_datExprTestStats)==1) list_datExprTestStats[[1]] else Reduce(f=rbind, x=list_datExprTestStats)
  return(df_datExprTestStats)
})

names(lvlRef_df_datExprTestStats) <-names(list_list_vec_identLvlsMatch)
lvlRef_df_datExprTestStats <- lvlRef_df_datExprTestStats[!sapply(lvlRef_df_datExprTestStats, function(df) all(is.na(df)))]
df_modPresStats <- if (length(lvlRef_df_datExprTestStats)==1) lvlRef_df_datExprTestStats[[1]] else Reduce(f=rbind, x=lvlRef_df_datExprTestStats)
df_modPresStats <- df_modPresStats[!df_modPresStats[["module"]]=="gold",]

######################################################################
########### GET NETWORK CONNECTIVITY PRESERVATION STATS ##############
######################################################################
df_nwPresStats = NULL

if (!is.null(lvlRef_datExprTest_presNwOut)) {
  
  df_nwPresStats <- data.frame(ref_dataset=NA_character_, 
                           ref_lvl=NA_character_, 
                           test_dataset=NA_character_, 
                           test_lvl =NA_character_)
  
  
  k=0
  for (lvlRef in names(lvlRef_datExprTest_presNwOut)) { 
    datExprTest_presNwOut <- lvlRef_datExprTest_presNwOut[[lvlRef]]
    if (all(is.na(datExprTest_presNwOut))) next
      for (datExprTestName in names(datExprTest_presNwOut)) {
        presNwOut <- datExprTest_presNwOut[[datExprTestName]]
        if (all(is.na(presNwOut))) next
        for (i in 1:ncol(presNwOut[["pairwise"]])) {
          k=k+1
          df_nwPresStats[k,"ref_dataset"] <- names(vec_pathsDatExpr)[1]
          df_nwPresStats[k,"ref_lvl"] <- lvlRef
          df_nwPresStats[k,"test_dataset"] <- datExprTestName
          df_nwPresStats[k,"test_lvl"] <- list_list_vec_identLvlsMatch[[lvlRef]][[datExprTestName]][i]
          df_nwPresStats[k,"network_pairwise_min"] <- min(presNwOut[["pairwise"]][,i])
          df_nwPresStats[k,"network_pairwise_median"] <- median(presNwOut[["pairwise"]][,i])
          df_nwPresStats[k,"network_pairwise_max"] <- max(presNwOut[["pairwise"]][,i])
          df_nwPresStats[k,"network_pairwiseWeighted_min"] <- min(presNwOut[["pairwiseWeighted"]][,i])
          df_nwPresStats[k,"network_pairwiseWeighted_median"] <- median(presNwOut[["pairwiseWeighted"]][,i])
          df_nwPresStats[k,"network_pairwiseWeighted_max"] <- max(presNwOut[["pairwiseWeighted"]][,i])
      }
      # will come out as character vector
    }
  }
}

######################################################################
#################### JOIN PRESERVATION RESULTS #############################
######################################################################

if (!is.null(df_nwPresStats)){
  df_presStats <- dplyr::left_join(df_modPresStats, df_nwPresStats, by=c("ref_dataset","ref_lvl", "test_dataset", "test_lvl"))
} else {
  df_presStats = df_modPresStats
}

######################################################################
### ADD GOODGENES, GOODSAMPLES, PROP GENES NOT NA, PROP GENES NOT 0 ##
######################################################################

df_presStats[["propModGenesInTestLvl"]] <-
  df_presStats[["propModgenesNon0inTestLvl"]] <-
  df_presStats[["test_lvl_goodGenes"]] <-
  df_presStats[["test_lvl_goodGenes"]] <- rep(NA_real_,nrow(df_presStats))

for (lvlRef in names(lvlRef_datExprTest_propModGenesInTestLvl)) {
  for (datExprTest in names(lvlRef_datExprTest_propModGenesInTestLvl[[lvlRef]])) {
    for (lvlTest in names(lvlRef_datExprTest_propModGenesInTestLvl[[lvlRef]][[datExprTest]])) {
      for (module in names(lvlRef_datExprTest_propModGenesInTestLvl[[lvlRef]][[datExprTest]][[lvlTest]])) {
        
        logical_row <- df_presStats$test_dataset == datExprTest &
          df_presStats$ref_lvl == lvlRef &
          df_presStats$test_lvl == lvlTest &
          df_presStats$module == module
        
        if (sum(logical_row)==1) {
          
          df_presStats[["propModGenesInTestLvl"]][logical_row] <- lvlRef_datExprTest_propModGenesInTestLvl[[lvlRef]][[datExprTest]][[lvlTest]][module]
          df_presStats[["propModgenesNon0inTestLvl"]][logical_row] <- lvlRef_datExprTest_propModgenesNon0inTestLvl[[lvlRef]][[datExprTest]][[lvlTest]][module]
          df_presStats[["test_lvl_goodGenes"]][logical_row] <-
            sum(datExprTest_testLvlGoodGenes[[lvlRef]][[datExprTest]][[lvlTest]]) /
            length(datExprTest_testLvlGoodGenes[[lvlRef]][[datExprTest]][[lvlTest]])
          df_presStats[["test_lvl_goodSamples"]][logical_row] <-
            sum(datExprTest_testLvlGoodSamples[[lvlRef]][[datExprTest]][[lvlTest]]) /
            length(datExprTest_testLvlGoodSamples[[lvlRef]][[datExprTest]][[lvlTest]])
          
        } else if (sum(logical_row)>1) {
          
          stop(paste0("lvlRef:", lvlRef, ", datExprTest:",  datExprTest, ", lvlTest:", lvlTest, ", module: ", module, "produced
                      an error : more than one row matches"))
          
        } else if (sum(logical_row)==0) {
          # modulePreservation failed, but we can still return the other stats about number of common genes etc
          df_presStats <- rbind(df_presStats, df_presStats[nrow(df_presStats),]) # copy the last row
          
          df_presStats[["ref_dataset"]][nrow(df_presStats)] <- names(vec_pathsDatExpr)[1]
          df_presStats[["test_dataset"]][nrow(df_presStats)] <- datExprTest
          df_presStats[["ref_lvl"]][nrow(df_presStats)] <- lvlRef
          df_presStats[["test_lvl"]][nrow(df_presStats)] <- lvlTest
          df_presStats[["module"]][nrow(df_presStats)] <- module
          df_presStats[["propModGenesInTestLvl"]][nrow(df_presStats)] <- lvlRef_datExprTest_propModGenesInTestLvl[[lvlRef]][[datExprTest]][[lvlTest]][module]
          df_presStats[["propModgenesNon0inTestLvl"]][nrow(df_presStats)] <- lvlRef_datExprTest_propModgenesNon0inTestLvl[[lvlRef]][[datExprTest]][[lvlTest]][module]
          df_presStats[nrow(df_presStats), !colnames(df_presStats) %in%
                         c("ref_dataset",
                           "ref_lvl",
                           "test_dataset",
                           "test_lvl",
                           "module",
                           "propModGenesInTestLvl",
                           "propModgenesNon0inTestLvl",
                           "test_lvl_goodGenes",
                           "test_lvl_goodSamples")] <- NA_real_
          df_presStats[["test_lvl_goodGenes"]][nrow(df_presStats)] <-
            sum(datExprTest_testLvlGoodGenes[[lvlRef]][[datExprTest]][[lvlTest]]) /
            length(datExprTest_testLvlGoodGenes[[lvlRef]][[datExprTest]][[lvlTest]])
          df_presStats[["test_lvl_goodSamples"]][nrow(df_presStats)] <-
            sum(datExprTest_testLvlGoodSamples[[lvlRef]][[datExprTest]][[lvlTest]]) /
            length(datExprTest_testLvlGoodSamples[[lvlRef]][[datExprTest]][[lvlTest]])
        }
      }
    }
  }
}

######################################################################
############################ OUTPUTS #################################
######################################################################

fwrite(x=df_presStats, file = paste0(dirTables, prefixOut, "_df_preservationStats.csv"))
write.xlsx(x=df_presStats, file = paste0(dirTables, prefixOut, "_df_preservationStats.xlsx"))
saveRDS(object=datExprTest_testLvlGoodSamples, file=paste0(dirRObjects, prefixOut, "_list_list_vec_goodSamples.RDS.gz"), compress="gzip")
saveRDS(object=datExprTest_testLvlGoodGenes, file=paste0(dirRObjects, prefixOut, "_list_list_vec_goodGenes.RDS.gz"), compress="gzip")
############################ WRAP UP ################################

message("Script done!")
