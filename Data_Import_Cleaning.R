Data_Import_Cleaning <- function(){
  
  tryCatch({
    ### Sourcing tools from GEO ###
    
    source("http://www.bioconductor.org/biocLite.R")
    
    BiocManager::install("GenomicRanges")
    BiocManager::install("rtracklayer")
    BiocManager::install("usethis")
    library(GenomicRanges)
    library(rtracklayer)
    library(usethis)
    
    ### Installing package dependencies ###
    
    BiocManager::install("devtools")
    BiocManager::install("roxygen2")
    library(devtools)
    library(roxygen2)
  }, 
  
  # warning = function(warning_condition) {
  #   warning-handler-code
  # }, 
  
  error = function(e) {
    break()
  }, 
  
  finally={
    ### Importing the master table ###
    
    setwd(".")
    proj_set()
    ChIPSeqDataMaster <- read.csv(file = "GSAChIPSeqBenchmarkDatasetProfile.txt", sep = '\t', header = TRUE, quote = "")
    use_data(ChIPSeqDataMaster, internal = FALSE, overwrite = TRUE)
    
    ## The following code creates a list of samples for which we need to extract the BED files for.
    
    ChIPSeqSamples<-as.character(ChIPSeqDataMaster$GSM)
    
    ## Initializing list for storing BED files and the consecutive GRanges objects.
    Samples_in_BED = list()
    
    for(i in 1:length(ChIPSeqSamples))
    {
      Samples_in_BED[[i]] <- read.table(paste("./regen/",paste(eval(parse(text="ChIPSeqSamples[i]")),".bed", sep = ""), sep = ""), sep = "\t", header = FALSE)
      Samples_in_BED[[i]] <- Samples_in_BED[[i]][,1:3]
      colnames(Samples_in_BED[[i]]) <- c("chrom", "start", "end")
      Samples_in_BED[[i]] <- Samples_in_BED[[i]][order(Samples_in_BED[[i]]$chrom),]
      Samples_in_BED[[i]] <- GRanges(Samples_in_BED[[i]]$chrom, IRanges(Samples_in_BED[[i]]$`start`, Samples_in_BED[[i]]$`end`))
      genome(Samples_in_BED[[i]]) <- "hg19"
    }
    
    ## Saving BED files as GRanges objects ##
    Samples_in_BED <- GRangesList(Samples_in_BED)
    names(Samples_in_BED) <- ChIPSeqSamples
    use_data(Samples_in_BED, internal = FALSE, overwrite = TRUE)
  }
)}
  



