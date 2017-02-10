#############
# Functions #
#############

getGeneToGoMapping <- function(userGivenSpeciesIndex){
        
        ## get file path
        file = getBgDBFile(speciesIndex = userGivenSpeciesIndex)
        
        ## read file 
        go_association_file <- readLines(file)
        
        ##seprate data and headers.
        isHeaderLine <- startsWith(go_association_file , "!")
        headerLine <- go_association_file[isHeaderLine]
        dataLines <- go_association_file[!isHeaderLine]
        
        ## extract taxon ID from organism specific line. 
        isOrganismInfo <- grepl("Organism ",headerLine,fixed = T)
        organismInfoLines <- headerLine[isOrganismInfo]
        
        ##reassign species index 
        userChoice <- userGivenSpeciesIndex
        
        ## seprate data as per taxon id 
        taxonToSelect <- makeUserChoiceForTaxonId(userChoice = userChoice)
        isTaxonSpecific <- grepl(paste("taxon:",taxonToSelect,sep = ""),dataLines,fixed = T)
        taxonSpecificLines <- dataLines[isTaxonSpecific]
        splitedLines <- strsplit(taxonSpecificLines,"\t",perl = TRUE) # split by one 
        asMat <- do.call(rbind,splitedLines)
        asDf <- as.data.frame(asMat)
        geneToGo <- asDf[,c(3,5)]   # column 3 belongs to gene_shortname and 5 belongs to go id      
        
        # convert to list object 
        geneID2GO <- convertManyToOne(geneToGo,1,2)        
        
        ##return
        geneID2GO
}

convertManyToOne <- function(data,indexOfUniqCol,indexOfRedundantCol,sep=", "){
        
        
        data[c(1:ncol(data))] = apply(data, 2, as.character) ## convert dataframe colum to char array. This is required to run aggregate fun properly 
        ag <- aggregate(data[indexOfRedundantCol], by=list(data[,indexOfUniqCol]), c) 
        
        # convert df to list form. 
        lst  <- ag[,2] 
        names(lst) <- ag[,1]
        return(lst)
}

## help message 
helpMessage <- function(CommandArgs){
        CommandArgs = "--help"
        if(CommandArgs == "--help"){
                cat("--help\n")
                cat("Run Commnad --> RScript toGO_customesed.R <species code> <geneset file>\n")  
                cat("******************\n")  
                cat("<species code>\n")  
                cat("******************\n")  
                allBgSpecies <- getBackGroundSpecies()
                cat(allBgSpecies,sep = "\n")
                cat("******************\n")  
                cat("<geneset file>\n")  
                cat("******************\n")  
                cat("Path of file containing one or more genesets. \nFile must be tab deliminated having one geneset in one column including header in first row\nGenes must be given in form of common name, if common name (like XBP1, areA etc...) not available then and then systamatic name acceptable (AN id , CAGL number) ")
                
        }
}

## validate user args 
validateCommnadLineArgs <- function(args){
        ## check total arguments match to required args or not 
        #args <- c("1","an_fang")
        requiredArgs <- 2
        totalArgs <- length(args)
        if(totalArgs == 1 & args[1] =="--help"){
                helpMessage()
                stop()
        }
        if(totalArgs != requiredArgs){
                cat("Two arguments required to run the script. Refer the help below \n\n")
                helpMessage()
                stop()
        }
        
        ## check individual args
        speciesCode <- as.integer(args[1])
        isGenesetFileExist <- file.exists(args[2])
        if(is.na(speciesCode)){
                cat("Species code not valid. It must be in between 1 to", length(getBackGroundSpecies()))
                helpMessage()
                stop()
        }else if(speciesCode < 1 || speciesCode >  length(getBackGroundSpecies())){
                cat("Species code not valid. It must be in between 1 to", length(getBackGroundSpecies()))
                helpMessage()
                stop()
        }else if(!isGenesetFileExist){
                cat("File doesn't exist. Give full path")
                helpMessage()
                stop()
        }else{
                cat("Given Arguments are OK. Performing Enrichment ......\n")
        }
        
        
}

## Set background genes file 
## NOTE : background geneset file can be downloaded from respective databases. 
getBgDBFile <-  function(speciesIndex){
        selectedSpeciesIndex = speciesIndex
        dbFiles=c(Candida=paste(dbFilePath,"gene_association.cgd",sep="/"),
                  Aspergillus=paste(dbFilePath,"gene_association.aspgd",sep="/"),
                  Saccharomyces=paste(dbFilePath,"gene_association.sgd",sep = "/"))
        
        if(selectedSpeciesIndex >= 1 & selectedSpeciesIndex <= 11){
                dbFile =dbFiles["Candida"]
        }else if(selectedSpeciesIndex > 11 & selectedSpeciesIndex <= 31){
                dbFile =dbFiles["Aspergillus"]
        }else if(selectedSpeciesIndex == 32){
                dbFile =dbFiles["Saccharomyces"]
        }
        return(dbFile)
}

## vector of all background species 
getBackGroundSpecies <-function(){
        totalSpecies <- c("1: Candida albicans SC5314,genome version: A22-s07-m01-r18, taxon: 5476",
                          "2: Candida glabrata CBS138, genome version: s02-m07-r13, taxon: 5478",
                          "3: Candida parapsilosis CDC317, genome version: s01-m03-r18, taxon: 5480",
                          "4: Candida dubliniensis CD36, genome version: s01-m02-r11, taxon: 42374",
                          "5: Candida orthopsilosis Co 90-125, taxon: 1136231",
                          "6: Debaryomyces hansenii CBS767, taxon: 284592",
                          "7: Candida albicans WO-1, taxon: 294748",
                          "8: Lodderomyces elongisporus NRLL YB-4239, taxon: 379508",
                          "9: Candida tropicalis MYA-3404, taxon: 294747",
                          "10: Candida lusitaniae ATCC 42720, taxon: 306902",
                          "11: Candida guilliermondii ATCC 6260, taxon: 294746",
                          "12: Aspergillus nidulans FGSC A4, genome version: s10-m04-r06, taxon: 162425",
                          "13: Aspergillus fumigatus Af293, genome version: s03-m05-r07, taxon: 746128",
                          "14: Aspergillus niger CBS 513.88, genome version: s01-m07-r03, taxon: 5061",
                          "15: Aspergillus oryzae RIB40, genome version: s01-m09-r04, taxon: 5062",
                          "16: Aspergillus terreus NIH2624, taxon: 341663",
                          "17: Aspergillus brasiliensis CBS 101740, taxon: 767769",
                          "18: Aspergillus zonatus, taxon: 41063",
                          "19: Aspergillus acidus CBS 106.47, taxon: 1137211",
                          "20: Aspergillus wentii DTO 134E9, taxon: 1073089",
                          "21: Aspergillus versicolor, taxon: 46472",
                          "22: Aspergillus carbonarius ITEM 5010, taxon: 602072",
                          "23: Aspergillus aculeatus ATCC16872, taxon: 690307",
                          "24: Aspergillus glaucus CBS 516.65, taxon: 1160497",
                          "25: Neosartorya fischeri NRRL 181, taxon: 331117",
                          "26: Aspergillus clavatus NRRL 1, taxon: 344612",
                          "27: Aspergillus flavus NRRL 3357, taxon: 332952",
                          "28: Aspergillus tubingensis CBS 134.48, taxon: 767770",
                          "29: Aspergillus sydowii, taxon: 75750",
                          "30: Aspergillus fumigatus A1163, taxon: 451804",
                          "31: Aspergillus kawachii, taxon: 1033177",
                          "32: Saccharomyces Cerevisiae, taxon: 559292")
        totalSpecies
}


## process data 
doProcess <- function(myGeneSet,backgroundGenes,outFileName){
        # i <- 1
        # myGeneSet = myInterestingGenes[[i]]
        # backgroundGenes = geneID2GO
        # outFileName = names(myInterestingGenes[i])
        
        geneNames <- names(backgroundGenes) ## all genes will be used as BG genes        
        geneList <- factor(as.integer(toupper(geneNames) %in% toupper(myGeneSet)))
        
        if(length(levels(geneList)) < 2){
                if(levels(geneList) == 0){
                        cat("\t --> Enrichment can not be performed \t\t\t\n")
                        cat("\t --> None of the gene present in the back ground data\t\t\t\n")
                        cat("\t --> Either you made wrong choice of organism or give differnt geneset \t\t\t\n")        
                }else if(levels(geneList) == 1){
                        cat("\t --> Enrichment can not be performed \t\t\t\n")
                        cat("\t --> All of the genes present in the back ground data\t\t\t\n")
                }
                stop()
        }
        names(geneList) <- geneNames
        str(geneList)
        
        ## prepare master object for GO enrichment 
        GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
        
        ## perform enrichment test 
        algo <- "elim"  ##  choice of the algorithms are classic,elim and weight
        stat <- "fisher"
        resultFisher.elim <- runTest(GOdata, algorithm = algo, statistic = stat)
        #resultFisher.weight <- runTest(GOdata, algorithm = "weight", statistic = "fisher")
        
        ## write Data 
        allRes <- GenTable(GOdata, elimFisher = resultFisher.elim, topNodes = 100,numChar=1000)
        #allRes <- GenTable(GOdata, classicFisher = resultFisher.classic, topNodes = 20)
        #allRes <- GenTable(GOdata, weightFisher = resultFisher.weight, topNodes )
        
        ## Get genes for enriched GO
        all_go_with_all_genes = genesInTerm(GOdata)
        all_go_with_my_genes = lapply(all_go_with_all_genes,function(x) x[x %in% myGeneSet] )
        enriched_go_with_my_genes <- lapply(all_go_with_my_genes[allRes[,1]], paste0, collapse = ",")
        allRes$genes <- unlist(enriched_go_with_my_genes)
        
        ## write output
        write.table(x = allRes,file = outFileName,row.names = F,quote = F,sep = "\t")
        
}

getUserGeneSet <- function(geneSetFile){
        
        userGeneSets <- read.table(geneSetFile,header = T,sep="\t",colClasses = "character",fill = TRUE, na.strings = c("","NA"))
        cat("# given genesets are ",dim(userGeneSets)[2] , "\n")
        userGeneSets <- as.list(userGeneSets) ## convert to list
        naRemovedGeneSets <- lapply(userGeneSets, function(elem) return(elem[!is.na(elem)])) ## remove NA
        geneSetLen <- lengths(naRemovedGeneSets) ## get length of each geneset
        minimumGenesRequiredForEnrichment <- 5
        for(i in seq_along(geneSetLen)){
                len <- geneSetLen[i]
                if(len <= minimumGenesRequiredForEnrichment){
                        cat("enrichment won't be performed for \"", names(geneSetLen[i]) ,"\" Minimum genes must be more than 5\n",sep="")
                }
        }
        geneSetsForEnrichment <- naRemovedGeneSets[geneSetLen > minimumGenesRequiredForEnrichment]
        geneSetsForEnrichment
}

# selectSpecies <- function(){
#         
#         speciesVector <- c("Candida","Aspergillus","Saccharomyces")
#         index <- seq_along(speciesVector)
#         names(index) <- speciesVector
#         cat("\n\n")
#         cat("\t\t\t*********************************\t\t\t\n")
#         cat("\t\t\t\tselect your species\t\t\t\n")
#         cat("\t\t\t*********************************\t\t\t\n")
#         
#         for(i in seq_along(speciesVector)){
#                 cat("\n\n")
#                 cat("\t\t\t",i,". For ", names(index[i]), " related species\t\t\t\n" , sep = "")
#         }
#         selectedSpecies <- readline()
#         speciesIndex <- as.integer(selectedSpecies)
#         cat("\t\t\t",rep(".",nchar(paste("\t\t\t you have selected ", names(index[speciesIndex]), "\t\t\t\n",sep=""))),"\t\t\t\n",sep="")
#         cat("\t\t\tyou have selected ", names(index[speciesIndex]), "\t\t\t\n")
#         cat("\t\t\t",rep(".",nchar(paste("\t\t\t you have selected ", names(index[speciesIndex]), "\t\t\t\n",sep=""))),"\t\t\t\n",sep="")
#         toBeReturned <- index[speciesIndex]
#         return(toBeReturned)
# }

makeUserChoiceForTaxonId <-  function(userChoice){
        ##input is char vactor having lines of organism information
        lineWithTaxonId <- getBackGroundSpecies()[userChoice]
        splitted <- strsplit(lineWithTaxonId ," ")[[1]] ## split by space
        selectedTaxon <- as.integer(splitted[length(splitted)])  ## get last elem of vector

        ## return
        selectedTaxon
}

# displayMessage <- function(displayText){
#         #displayText <- input        
#         cat("\n\n")
#         cat("\t\t\t*********************************\t\t\t\n")
#         cat("\t\t\tselect organism for GO enrichment\t\t\t\n")
#         cat("\t\t\t*********************************\t\t\t\n")
#         cat("\n\n")
#         cat(displayText,sep = "\n")
# }


