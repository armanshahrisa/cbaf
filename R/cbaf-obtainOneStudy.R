#' @title Obtaining the requested data for various subgroups of a cancer study.
#'
#' @description This function Obtains the requested data for the given genes across multiple subgroups of a cancer. It can
#' check whether or not all genes are included in subgroups of a cancer study and, if not, looks for the alternative gene names.
#'
#' @details
#' \tabular{lllll}{
#' Package: \tab cbaf \cr
#' Type: \tab Package \cr
#' Version: \tab 0.99.0 \cr
#' Date: \tab 2017-06-28 \cr
#' License: \tab Artistic-2.0 \cr
#' }
#'
#'
#'
#' @importFrom cgdsr CGDS getCancerStudies getCaseLists getGeneticProfiles getProfileData
#'
#' @importFrom BiocFileCache BiocFileCache bfcnew bfcquery bfcpath
#'
#'
#'
#' @usage obtainOneStudy(genesList, submissionName, studyName, desiredTechnique,
#' desiredCaseList = FALSE, validateGenes = TRUE)
#'
#' @param genesList a list that contains at least one gene group
#'
#' @param submissionName a character string containing name of interest. It is used for naming the process.
#'
#' @param studyName a character string showing the desired cancer name. It is an standard cancer study name that can be found on
#' cbioportal.org, such as 'Acute Myeloid Leukemia (TCGA, NEJM 2013)'.
#'
#' @param desiredTechnique a character string that is one of the following techniques: 'RNA-seq', 'microRNA-Seq', 'microarray.mRNA',
#' 'microarray.microRNA' or 'methylation'.
#'
#' @param desiredCaseList a numeric vector that contains the index of desired cancer subgroups, assuming the user knows index of
#' desired subgroups. If not, desiredCaseList must be set as 'none', function will show the available subgroups and ask the
#' user to enter the desired ones during the process. The default value is 'none'.
#'
#' @param validateGenes a logical value that, if set to be 'TRUE', function will check each cancer study to find whether
#' or not each gene has a record. If the given cancer doesn't have a record for specific gene, it checks for alternative gene
#' names that cbioportal might use instead of the given gene name.
#'
#'
#'
#' @return a list that contains the obtained data without further processing. Name of the list starts with 'obS' and contains
#' submissionName. Inside the list, there is one subgroup for every gene group, which itself contains one matrix for every study
#' subgroup. In addition, if validateGenes = TRUE, a secondary list containing gene validation results will be stored.
#' Name of the second list starts with 'vaS' and containes submissionName.
#'
#'
#'
#' @examples
#' genes <- list(K.demethylases = c("KDM1A", "KDM1B", "KDM2A"))
#'
#' # obtainOneStudy(genes, "test", "Breast Invasive Carcinoma (TCGA, Cell 2015)",
#' # "RNA-seq", desiredCaseList = c(3,4,5))
#'
#'
#'
#' @author Arman Shahrisa, \email{shahrisa.arman@hotmail.com} [maintainer, copyright holder]
#' @author Maryam Tahmasebi Birgani, \email{tahmasebi-ma@ajums.ac.ir}
#'
#' @export



#########################################################################
#########################################################################
########## Obtain the requested data for Subtypes of a Cancer ###########
#########################################################################
#########################################################################

obtainOneStudy <- function(genesList, submissionName, studyName, desiredTechnique, desiredCaseList = FALSE, validateGenes = TRUE){

  ##########################################################################
  ########## Prerequisites

  # Check genes

  if(!is.list(genesList)){

    stop("'genes' must be entered as a list containing at list one group of genes with descriptive group name for logistical purposes")

  }



  # Check submissionName

  if(!is.character(submissionName)){

    stop("'submissionName' must be entered as a character string for naming the process")

  }



  # cancer name

  if(!is.character(studyName)){

    stop("'studyName' must be entered as a character string")

  }



  # high-throughput data type

  if(is.character(desiredTechnique)){

    if(!(desiredTechnique %in% c("RNA-seq", "microRNA-Seq", "Microarray.mRNA", "Microarray.microRNA", "methylation")) | length(desiredTechnique)!= 1){

      stop("'desiredTechnique' must contain one of the following techniques: 'RNA-seq', 'microRNA-Seq', 'microarray.mRNA', 'microarray.microRNA' or

           'methylation'")

    }

  } else {

    stop("'desiredTechnique' must be entered as a character string describing a technique name")

  }



  # Choice of high-throughput data type

  if(desiredTechnique == "RNA-seq"){

    L2.characteristics <- c("mRNA Expression z-Scores (RNA Seq V2 RSEM)", "mRNA Expression z-Scores (RNA Seq RPKM)")

  } else if(desiredTechnique == "microRNA-Seq"){

    L2.characteristics <- c("microRNA expression Z-scores")

  } else if(desiredTechnique == "microarray.mRNA"){

    L2.characteristics <- c("mRNA Expression z-Scores (microarray)", "mRNA Expression z-Scores (U133 microarray only)", "mRNA expression z-scores (Illumina)", "mRNA expression Z-scores (all genes)", "mRNA Expression Z-Scores vs Normals", "mRNA Expression z-Scores (combined microarray)")

  } else if(desiredTechnique == "microarray.microRNA"){

    L2.characteristics <- c("mRNA Expression Z-Scores vs Normals", "mRNA Expression z-Scores (combined microarray)")

  } else if(desiredTechnique == "methylation"){

    L2.characteristics <- c("Methylation (HM450)", "Methylation (HM27)", "Methylation")

  } else{

    stop("desiredTechnique field can not be left empety. It should be chosen as 'RNA-seq', 'microRNA-Seq', 'microarray.mRNA', 'microarray.microRNA'or 'methylation'")

  }



  ##########################################################################
  ########## Decide whether functions should stops now!

  # Store the new parameteres

  newParameters <-list()

  newParameters$genesList <- genesList

  newParameters$submissionName <- submissionName

  newParameters$studyName <- studyName

  newParameters$desiredTechnique <- desiredTechnique

  if(is.logical(desiredCaseList)){

    newParameters$desiredCaseList <- 0

  } else{

    newParameters$desiredCaseList <- desiredCaseList

  }

  newParameters$validateGenes <- validateGenes

  newParameters$HaultOrder <- FALSE





  # Check wheather the requested data exists

  if(exists(paste("bfc_", submissionName, sep = ""))){

    bfc <- get(paste("bfc_", submissionName, sep = ""))

    if(nrow(bfcquery(bfc, "Parameters for obtainOneStudy()")) == 1){

      load(bfcpath(bfc, bfcquery(bfc, c("Parameters for obtainOneStudy()"))$rid))

      if(identical(oldParameters[-7], newParameters[-7])){

        continue <- FALSE

        newParameters$HaultOrder <- TRUE

        oldParameters <- newParameters

        save(oldParameters, file=bfc[[bfcquery(bfc, "Parameters for obtainOneStudy()")$rid]])

        assign(paste("bfc_", submissionName, sep = ""), bfc, envir = globalenv())

        print("--- Function 'obtainOneStudy()' was skipped: the requested data already exist ---")

      }else{

        continue <- TRUE

      }

    }else{

      continue <- TRUE

    }

  } else{

    continue <- TRUE

  }





  if(continue == TRUE){

    ##########################################################################
    ########## Set the function ready to work

    # Set cgdsr

    mycgds = CGDS("http://www.cbioportal.org/")

    # First step of procedure

    mycancerstudy = getCancerStudies(mycgds)[which(getCancerStudies(mycgds)[,2]==studyName),1]



    ##########################################################################
    ########## Core segment

    # Chosing the desired case lists

    if(desiredCaseList == FALSE && is.logical(desiredCaseList)){

      Choices <- getCaseLists(mycgds,mycancerstudy)[,2]

      print(paste(1:length(Choices), Choices, sep=". "))

      writeLines("")

      print(paste("Please enter the numeric index of desired case list(s) for ", studyName, ", seperated by comma. For instance, 1 and 2 must be enterd as: 1, 2.", sep=""))

      inputCases <- readline(prompt = "Enter the numeric index(es): ")

      inputCases <- as.numeric(unlist(strsplit(inputCases, ",")))

      if(is.character(inputCases)){

        stop("Desired case list(s) must contain numbers only")
      }

    } else {

      if(is.numeric(desiredCaseList)){

        inputCases <- desiredCaseList

      }

    }


    # Creating a vector which contains names of inputCases

    inputCases.names <- getCaseLists(mycgds,mycancerstudy)[inputCases ,2]





    # Finding the second characteristic of data in the cancer

    s.condition <- (getGeneticProfiles(mycgds,mycancerstudy)[,2])[getGeneticProfiles(mycgds,mycancerstudy)[,2] %in% L2.characteristics]

    mygeneticprofile = getGeneticProfiles(mycgds,mycancerstudy)[which(getGeneticProfiles(mycgds,mycancerstudy)[,2] == if(length(s.condition) >= 1){

      s.condition[1]

    } else if (length(s.condition) == 0){

      stop(paste(studyName, "lacks", desiredTechnique, "data!", sep=" "))

    }) ,1]


    # Shorten studyName - Temporarily inactive

    #  if(shortenStudyName == TRUE){

    #  studyName <- gsub(" ", ".", sapply(strsplit(as.character(studyName), split=" (", fixed=TRUE), function(x) (x[1])))

    #  }


    # Create parent list for storing final results in the global environment

    rawList <- list()

    # Creating child lists

    for(nname in 1:length(genesList)){

      rawList[[nname]] <- list()

      names(rawList)[nname] <- names(genesList)[nname]

    }

    # Creating a list for gene validation results

    if(validateGenes == TRUE){

      validationResult <- list()

      for(nname in 1:length(genesList)){

        validationResult[[nname]] <- "x"

        names(validationResult)[nname] <- names(genesList)[nname]

      }

    }





    # Report

    print(paste("***", "Obtaining the requested data for", submissionName, "***", sep = " "))

    # Creating progress bar

    obtainOneStudyProgressBar <- txtProgressBar(min = 0, max = length(inputCases), style = 3)



    ## Getting the required gene expresssion profile ...

    # 'for' control structure for obtaining data and calculating the requested parameters

    for(i in 1:length(inputCases)){

      # Determining name for list subset of study name

      groupName <- inputCases.names[i]

      # Correcting possible errors of list names

      groupName <- gsub(groupName, pattern = "\\+ ", replacement = " possitive ", ignore.case = TRUE)

      groupName <- gsub(groupName, pattern = "\\- ", replacement = " negative ", ignore.case = TRUE)

      groupName <- gsub(groupName, pattern = " ", replacement = "_", ignore.case = TRUE)



      # Setting the first characteristics of data according to the desired case list

      ind <- getCaseLists(mycgds,mycancerstudy)[inputCases[i] ,2]

      mycaselist = getCaseLists(mycgds,mycancerstudy)[inputCases[i] ,1]





      # obtaining data for every genegroup

      for(group in 1:length(genesList)){

        # Chose one group of genes

        genesNames <- genesList[[group]]



        # Obtaining Expression x-scores fore the requested genes

        # Assaign data to specific list member

        rawList[[group]][[i]] <- data.matrix(getProfileData(mycgds,genesNames[order(genesNames)],mygeneticprofile,mycaselist))

        names(rawList[[group]])[i] <- groupName

        # Find whether alternative gene names are used

        # Alter c.genes to be compatible with gene names in cBioPortal output

        alteredGeneNames <- sort(gsub("-", ".", genesNames))

        # Obtain name of genes that are absent in requested cancer

        absentGenes <- alteredGeneNames[!alteredGeneNames %in% colnames(rawList[[group]][[i]])]

        # For loop for determining changed genes

        alternativeGeneNames <- vector("character", length = length(absentGenes))

        # For loop

        for(ab in 1:length(absentGenes)){

          absentGeneProfileData <- colnames(data.matrix(getProfileData(mycgds, absentGenes[ab], mygeneticprofile,mycaselist)))

          # Check wheter gene has an alternative name or missed in the database

          if(length(absentGeneProfileData) == 1 & is.matrix(absentGeneProfileData)){

            alternativeGeneNames[ab] <- absentGeneProfileData

          } else{

            alternativeGeneNames[ab] <- "-"

          }

        }

        # Naming Alternative.gene.names

        names(alternativeGeneNames) <- absentGenes

        # Seperating genes with alternative names from those that are absent

        genesLackData <- alternativeGeneNames[alternativeGeneNames == "-"]

        genesWithData <- alternativeGeneNames[alternativeGeneNames != "-"]



        # modifying gene names containing an alternative name

        for(re in 1:length(genesWithData)){

          colnames(rawList[[group]][[i]])[colnames(rawList[[group]][[i]]) %in% genesWithData[re]] <- paste(genesWithData[re], " (", names(genesWithData[re]), ")", sep = "")

        }





        # validateGenes

        if(validateGenes == TRUE){

          # Empty validation matrix

          validationMatrix <- matrix(, ncol = length(genesNames), nrow = 1)

          # Naming empty matrix

          if(length(genesLackData) != 0){

            dimnames(validationMatrix) <- list(inputCases.names[i], c(colnames(rawList[[group]][[i]]), names(genesLackData)))

          } else{

            dimnames(validationMatrix) <- list(inputCases.names[i], colnames(rawList[[group]][[i]]))

          }



          # modifying gene names containing an alternative name

          for(re in 1:length(genesWithData)){

            colnames(validationMatrix)[colnames(validationMatrix) %in% genesWithData[re]] <- paste(genesWithData[re], " (", names(genesWithData[re]), ")", sep = "")

          }





          # Puting value for genes lacking data

          validationMatrix[,colnames(validationMatrix) %in% names(genesLackData)] <- "-"



          for(eval in 1:ncol(rawList[[group]][[i]])){

            ## Validating Genes

            # Correct those that are not found

            if(length(((rawList[[group]][[i]])[,eval])[!is.nan((rawList[[group]][[i]])[,eval])]) > 0 & all(!is.finite((rawList[[group]][[i]])[,eval])) &

               is.nan(mean(as.vector((rawList[[group]][[i]])[,eval])[abs((rawList[[group]][[i]])[,eval])], na.rm=TRUE))){

              validationMatrix[1, eval] <- "-"

            } else {

              validationMatrix[1, eval] <- "Found"

            }

          }

          # Storing the results in validation result

          validationMatrix <- validationMatrix[,sort(colnames(validationMatrix)), drop=FALSE]

          if(i == 1){

            validationResult[[group]]    <- validationMatrix

          } else if(i > 1){

            validationResult[[group]]    <- rbind(validationResult[[group]], validationMatrix)

          }

        }

      }

      # Update progressbar

      setTxtProgressBar(obtainOneStudyProgressBar, i)

    }

    # Closing progress bar

    close(obtainOneStudyProgressBar)




    ## bfc object

    # create bfc object in global environment

    if(!exists(paste("bfc_", submissionName, sep = ""))){

      bfc <- BiocFileCache(file.path(tempdir(), submissionName))

    } else{

      bfc <- get(paste("bfc_", submissionName, sep = ""))

    }



    # Store the obtained Data

    if(nrow(bfcquery(bfc, "Obtained Data for single study")) == 0){

      save(rawList, file=bfcnew(bfc, "Obtained Data for single study", ext="RData"))

    } else if(nrow(bfcquery(bfc, "Obtained Data for single study")) == 1){

      save(rawList, file=bfc[[bfcquery(bfc, "Obtained Data for single study")$rid]])

    }



    # Store the validation data

    if(validateGenes == TRUE){

      if(nrow(bfcquery(bfc, "Validation Data for single study")) == 0){

        save(validationResult, file=bfcnew(bfc, "Validation Data for single study", ext="RData"))

      } else if(nrow(bfcquery(bfc, "Validation Data for single study")) == 1){

        save(validationResult, file=bfc[[bfcquery(bfc, "Validation Data for single study")$rid]])

      }

    }

    # Store the parameters for this run

    oldParameters <- newParameters

    if(nrow(bfcquery(bfc, "Parameters for obtainOneStudy()")) == 0){

      save(oldParameters, file=bfcnew(bfc, "Parameters for obtainOneStudy()", ext="RData"))

    } else if(nrow(bfcquery(bfc, "Parameters for obtainOneStudy()")) == 1){

      save(oldParameters, file=bfc[[bfcquery(bfc, "Parameters for obtainOneStudy()")$rid]])

    }



    # Store bfc in global environmet

    assign(paste("bfc_", submissionName, sep = ""), bfc, envir = globalenv())

  }

}
