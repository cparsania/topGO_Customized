#library(psych)
#ibrary(topGO)
isTopGoPresent <- require(topGO)
if(!isTopGoPresent){
        cat("\n\n")
        cat("Installing top go package.....")
        Sys.sleep(time = 3)
        source("http://bioconductor.org/biocLite.R")
        biocLite("topGO")
        library("topGO")
}
#source("~/Documents/scripts/R_Script/topGo/topGOfunctions.R")
source("./R/topGOfunctions.R")
dbFilePath="./gene_association_files"

#################
# Run commmands #
#################
## getCommnadLineArguments
args <- commandArgs(TRUE)
validateCommnadLineArgs(args = args)

## get argument 1 which is species index given by user
userGivenSpeciesIndex <- as.integer(args[1])

#### prepare geneID2GO object required for topGO object 
geneID2GO <- getGeneToGoMapping(userGivenSpeciesIndex = userGivenSpeciesIndex)

## SUer arg 2 (geneset file path) to create list of genes given by user  
myInterestingGenes <- getUserGeneSet(args[2]) ## get genes from user 

## call processDataFunction for given genesets
for(i in seq_along(myInterestingGenes)){
        ## print message 
        cat("\n\n")
        cat("\t\t\t",rep("*",nchar(paste("\t\t\tperforming enrichment for ---> ", names(myInterestingGenes[i]) ," \t\t\t\n",sep=""))),"\t\t\t\n",sep="")
        cat("\t\t\tperforming enrichment for ---> ", names(myInterestingGenes[i]) ," \t\t\t\n")
        cat("\t\t\t",rep("*",nchar(paste("\t\t\tperforming enrichment for ---> ", names(myInterestingGenes[i]) ," \t\t\t\n",sep=""))),"\t\t\t\n",sep="")
        
        ## call function
        doProcess(myGeneSet = myInterestingGenes[[i]],backgroundGenes = geneID2GO,outFileName = paste(names(myInterestingGenes[i]),".txt",sep = ""))        
        
}
cat("\n\n")
cat("Go Enrichment Finished. Get your results from current dir.\n\n")
