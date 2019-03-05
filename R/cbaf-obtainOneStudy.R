#' @title Obtain the requested data for various subgroups of a cancer study.
#'
#' @description This function Obtains the requested data for the given genes
#' across multiple subgroups of a cancer. It can check whether or not all genes
#' are included in subgroups of a cancer study and, if not, looks for the
#' alternative gene names.
#'
#' @details
#' \tabular{lllll}{
#' Package: \tab cbaf \cr
#' Type: \tab Package \cr
#' Version: \tab 1.5.2 \cr
#' Date: \tab 2019-03-05 \cr
#' License: \tab Artistic-2.0 \cr
#' }
#'
#'
#'
#' @importFrom cgdsr CGDS getCancerStudies getCaseLists getGeneticProfiles getProfileData
#'
#' @importFrom BiocFileCache BiocFileCache bfcnew bfcquery bfcpath
#'
#' @importFrom utils head setTxtProgressBar txtProgressBar
#'
#'
#'
#' @usage obtainOneStudy(genesList, submissionName, studyName, desiredTechnique,
#'   desiredCaseList = FALSE, validateGenes = TRUE)
#'
#'
#'
#' @param genesList a list that contains at least one gene group
#'
#' @param submissionName a character string containing name of interest. It is
#' used for naming the process.
#'
#' @param studyName a character string showing the desired cancer name. It is an
#'  standard cancer study name that can be found on cbioportal.org, such as
#'  \code{"Acute Myeloid Leukemia (TCGA, NEJM 2013)"}.
#'
#' @param desiredTechnique a character string that is one of the following
#' techniques: \code{"RNA-Seq"}, \code{"microRNA-Seq"}, \code{"microarray.mRNA"}
#' , \code{"microarray.microRNA"} or \code{"methylation"}.
#'
#' @param desiredCaseList a numeric vector that contains the index of desired
#' cancer subgroups, assuming the user knows index of desired subgroups. If not,
#'  desiredCaseList is set to \code{"none"}, function will show the available
#'  subgroups and ask the user to enter the desired ones during the
#'  process. The default value is \code{"none"}.
#'
#' @param validateGenes a logical value that, if set to be 'TRUE', causes the
#' function to check each cancer study to find whether or not each gene has a
#' record. If a cancer doesn't have a record for specific gene, function looks
#' for alternative gene names that cbioportal might use instead of the given
#' gene name.
#'
#'
#'
#' @return a BiocFileCach object that contains the obtained data without further
#'  processing. Name of the object is combination of `bfc_` and submissionName.
#'  Inside it, there is a section for the obtained data, which is stored as a
#'  list. At first level, this list is subdivided into diferent groups based on
#'  the list of genes that user has given the function, then each gene group
#'  itself contains one matrix for every study subgroup. Additonally, if
#'  validateGenes = TRUE, another section that contains gene validation results
#'  will be created in the BiocFileCach object.
#'
#'
#'
#' @examples
#' genes <- list(K.demethylases = c("KDM1A", "KDM1B", "KDM2A", "KDM2B", "KDM3A",
#'  "KDM3B", "JMJD1C", "KDM4A"), K.methyltransferases = c("SUV39H1", "SUV39H2",
#'  "EHMT1", "EHMT2", "SETDB1", "SETDB2", "KMT2A", "KMT2A"))
#'
#' obtainOneStudy(genes, "test", "Breast Invasive Carcinoma (TCGA, Cell 2015)",
#' "RNA-Seq", desiredCaseList = c(2,3,4,5))
#'
#'
#'
#' @author Arman Shahrisa, \email{shahrisa.arman@hotmail.com} [maintainer,
#' copyright holder]
#' @author Maryam Tahmasebi Birgani, \email{tahmasebi-ma@ajums.ac.ir}
#'
#' @export



################################################################################
################################################################################
############## Obtain the requested data for Subtypes of a Cancer ##############
################################################################################
################################################################################

obtainOneStudy <- function(

  genesList,

  submissionName,

  studyName,

  desiredTechnique,

  desiredCaseList = FALSE,

  validateGenes = TRUE

  ){

  ##############################################################################
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

    supported.techniques <- c("RNA-Seq",

                              "microRNA-Seq",

                              "Microarray.mRNA",

                              "Microarray.microRNA",

                              "methylation")

    if(!(desiredTechnique %in% supported.techniques) |

       length(desiredTechnique)!= 1){

      stop("'desiredTechnique' must contain one of the following techniques: 'RNA-Seq', 'microRNA-Seq', 'microarray.mRNA', 'microarray.microRNA' or 'methylation'")

    }else if(desiredTechnique %in% supported.techniques |

             length(desiredTechnique)== 1){

      if(desiredTechnique == "RNA-Seq"){

        L2.characteristics <-

          c("mRNA Expression z-Scores (RNA Seq V2 RSEM)",

            "mRNA Expression z-Scores (RNA Seq RPKM)")

      } else if(desiredTechnique == "microRNA-Seq"){

        L2.characteristics <-

          c("microRNA expression Z-scores")

      } else if(desiredTechnique == "microarray.mRNA"){

        L2.characteristics <-

          c("mRNA Expression z-Scores (microarray)",

            "mRNA Expression z-Scores (U133 microarray only)",

            "mRNA expression z-scores (Illumina)",

            "mRNA expression Z-scores (all genes)",

            "mRNA Expression Z-Scores vs Normals",

            "mRNA Expression z-Scores (combined microarray)")

      } else if(desiredTechnique == "microarray.microRNA"){

        L2.characteristics <-

          c("mRNA Expression Z-Scores vs Normals",

            "mRNA Expression z-Scores (combined microarray)")

      } else if(desiredTechnique == "methylation"){

        L2.characteristics <-

          c("Methylation (HM450)",

            "Methylation (HM27)",

            "Methylation")

      }

    }

  } else {

    stop("'desiredTechnique' must be entered as a character string describing a technique name")

  }



  # checking only the first member of desiredCaseList in case it is vector of

  # length > 1

  if(desiredCaseList[1] == TRUE |

     !is.logical(desiredCaseList) & !is.numeric(desiredCaseList)){

    stop("'desiredCaseList' must be set as FALSE or a be a numeric vector.")

  }



  # Check validateGenes

  if(!is.logical(validateGenes)){

    stop("'validateGenes' can only accept logical values: TRUE or FALSE .")

  }





  ##############################################################################
  ########## Decide whether function should stops now!

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





  # Check the database

  database <- system.file("extdata", submissionName, package="cbaf")



  # Remove old database

  if(dir.exists(database) & !(submissionName %in% c("test", "test2"))){

    creation.time <- file.info(database , extra_cols = FALSE)$ctime

    past.time <-

      as.numeric(difftime(Sys.time(), creation.time, units = c("days")))

    if(past.time >= 5){

      unlink(database, recursive = TRUE)

    }

  }



  # Check wheather the requested data exists

  if(dir.exists(database)){

    bfc <- BiocFileCache(

      file.path(system.file("extdata", package = "cbaf"), submissionName)

    )

    if(nrow(bfcquery(bfc, "Parameters for obtainOneStudy()")) == 1){

      oldParameters <- readRDS(

        bfcpath(bfc, bfcquery(bfc, c("Parameters for obtainOneStudy()"))$rid)

      )

      if(identical(oldParameters[-7], newParameters) |

         submissionName %in% c("test", "test2")){

        continue <- FALSE

        # Store the last parameter

        newParameters$lastRunStatus <- "skipped"

        oldParamObtainOneStudy <- newParameters

        saveRDS(

          oldParamObtainOneStudy,

          file=bfc[[bfcquery(bfc, "Parameters for obtainOneStudy()")$rid]]

        )

        if(submissionName %in% c("test", "test2")){

          message("--- 'test' and 'test2' databases contain sample data and therefore, are not changable. Please use a different submission name. ---")

        }

        message("--- Function 'obtainOneStudy()' was skipped: the requested data already exist ---")

      }else{

        continue <- TRUE

      }

    }else{

      continue <- TRUE

    }

  } else{

    continue <- TRUE

  }





  if(continue){

    ############################################################################
    ########## Set the function ready to work

    # Set cgdsr

    mycgds = CGDS("http://www.cbioportal.org/")

    # First step of procedure

    mycancerstudy =

      getCancerStudies(mycgds)[which(getCancerStudies(mycgds)[,2]==studyName),1]

    caseList <- getCaseLists(mycgds,mycancerstudy)



    ############################################################################
    ########## Core segment

    # Chosing the desired case lists

    if(is.logical(desiredCaseList)){

      Choices <- caseList[,2]

      print(paste(seq_along(Choices), Choices, sep = ". "))

      writeLines("")

      message("Please enter the numeric index of desired case list(s) for ", studyName, ", seperated by comma. For instance, 1 and 2 must be enterd as: 1, 2")

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

    inputCases.names <- caseList[inputCases ,2]





    # Finding the second characteristic of data in the cancer

    existing.L2.charac <-

      getGeneticProfiles(mycgds,mycancerstudy)[,2] %in% L2.characteristics


    s.condition <-

      (getGeneticProfiles(mycgds,mycancerstudy)[,2])[existing.L2.charac]


    s.condition <-

      if(length(s.condition) >= 1){

        s.condition[1]

      } else if(length(s.condition) == 0){

        stop(studyName, "lacks", desiredTechnique, "data!", sep=" ")

      }


    s.condition.idx <-

      which(getGeneticProfiles(mycgds,mycancerstudy)[,2] == s.condition)


    mygeneticprofile =

      getGeneticProfiles(mycgds,mycancerstudy)[s.condition.idx ,1]


    # Shorten studyName - Temporarily inactive

    #  if(shortenStudyName){

    #  studyName <- gsub(" ", ".", sapply(strsplit(as.character(studyName),

    #                    split=" (", fixed=TRUE), function(x) (x[1])))

    #  }



    # Create parent list for storing final results in the global environment

    rawList <- list()

    # Creating child lists

    for(nname in seq_along(genesList)){

      rawList[[nname]] <- list(); names(rawList)[nname] <-

        names(genesList)[nname]

    }

    # Creating a list for gene validation results

    if(validateGenes){

      validationResult <- list()

      for(nname in seq_along(genesList)){

        validationResult[[nname]] <- "x"; names(validationResult)[nname] <-

          names(genesList)[nname]

      }

    }

    # Create Empty List to fill with validation matrices

    validationMList <- vector("list", length(genesList)*length(inputCases))





    # Report

    message("***", " Obtaining the requested data for ", submissionName, " ***")

    # Creating progress bar

    obtainOneStudyProgressBar <-

      txtProgressBar(min = 0, max = length(inputCases), style = 3)



    ## Getting the required gene expresssion profile ...

    # 'for' control structure for obtaining data and calculating the parameters

    for(i in seq_along(inputCases)){

      # Determining name for list subset of study name

      groupName <- inputCases.names[i]

      # Correcting possible errors of list names

      groupName <- gsub(

        groupName, pattern = "\\+ ", replacement = " possitive ",

        ignore.case = TRUE

      )

      groupName <- gsub(

        groupName, pattern = "\\- ", replacement = " negative ",

        ignore.case = TRUE

      )



      # Setting the first characteristics of data according to the case list

      ind <- getCaseLists(mycgds,mycancerstudy)[inputCases[i] ,2]

      mycaselist = getCaseLists(mycgds,mycancerstudy)[inputCases[i] ,1]





      # obtaining data for every genegroup

      for(group in seq_along(genesList)){

        # Chose one group of genes

        genesNames <- genesList[[group]]

        numberOfGenes <- length(genesNames)



        # Obtaining Expression x-scores for the requested genes

        # Check number of genes first

        if(numberOfGenes <= 250){

          ProfileData <- getProfileData(

            mycgds, genesNames[order(genesNames)], mygeneticprofile, mycaselist

          )

        }else{

          # split genes in groups of 250 names

          operational_gene_number <- split(

            genesNames[order(genesNames)], ceiling(seq_len(numberOfGenes)/250)

            )


          # Create empty list for gene_matrices

          separated_results <- vector(

            "list", length = operational_gene_number

            )

          for(operational in seq_along(operational_gene_number)){

            separated_results[[operational]] <- getProfileData(

              mycgds,

              operational_gene_number[[operational]],

              mygeneticprofile,

              mycaselist

            )

          }


          # Merging data

          ProfileData <- do.call("cbind", separated_results)

          ProfileData <- ProfileData[,order(colnames(ProfileData))]

        }




        # Assaign data to specific list member

        rawList[[group]][[i]] <- data.matrix(ProfileData)

        names(rawList[[group]])[i] <- groupName

        # For convenience

        this.segment <- rawList[[group]][[i]]



        # Find whether alternative gene names are used

        # Alter c.genes to be compatible with gene names in cBioPortal output

        alteredGeneNames <- sort(gsub("-", ".", genesNames))

        # Obtain name of genes that are absent in requested cancer

        absentGenes <-

          alteredGeneNames[!alteredGeneNames %in% colnames(this.segment)]

        # For loop for determining changed genes


        if(length(absentGenes) != 0){

          alternativeGeneNames <-

            vector("character", length = length(absentGenes))

          # For loop

          for(ab in seq_along(absentGenes)){

            absent.gene.profile.data <- getProfileData(

              mycgds, absentGenes[ab], mygeneticprofile, mycaselist

            )


            alternative.gene.name <- colnames(

              data.matrix(absent.gene.profile.data)

            )


            # Check wheter gene has an alternative name or missed from the

            # database

            if(length(alternative.gene.name) == 1){

              alternativeGeneNames[ab] <- alternative.gene.name

            } else if(length(alternative.gene.name) == 0){

              alternativeGeneNames[ab] <- "-"

            }

          }

          # Naming Alternative.gene.names

          names(alternativeGeneNames) <- absentGenes

          # Seperating genes with alternative names from those that are absent

          genesLackData <- alternativeGeneNames[alternativeGeneNames == "-"]

          genesWithData <- alternativeGeneNames[alternativeGeneNames != "-"]



          # modifying gene names containing an alternative name

          for(re in seq_along(genesWithData)){

            colnames.idx <-

              colnames(rawList[[group]][[i]]) %in% genesWithData[re]


            colnames(rawList[[group]][[i]])[colnames.idx] <-

              paste0(genesWithData[re], " (", names(genesWithData[re]), ")")

          }


        }else{

          genesLackData <- NULL

          genesWithData <- NULL

        }





        # validateGenes

        if(validateGenes){

          # Empty validation matrix

          validationMatrix <- matrix(, ncol = ncol(this.segment), nrow = 1)

          # Naming empty matrix

          if(length(genesLackData) != 0){

            dimnames(validationMatrix) <- list(

              inputCases.names[i],

              c(colnames(this.segment), names(genesLackData))

            )

          } else{

            dimnames(validationMatrix) <-

              list(inputCases.names[i], colnames(this.segment))

          }



          # modifying gene names containing an alternative name

          if(length(genesWithData) != 0){

            for(re in seq_along(genesWithData)){

              colnames.idx <-

                colnames(validationMatrix) %in% genesWithData[re]


              colnames(validationMatrix)[colnames.idx] <-

                paste0(genesWithData[re], " (", names(genesWithData[re]), ")")

            }

          }





          # Puting value for genes lacking data

          colnames.idx <- colnames(validationMatrix) %in% names(genesLackData)

          validationMatrix[,colnames.idx] <- "-"



          for(eval in seq_len(ncol(this.segment))){

            loop.section <- (this.segment)[,eval]

            ## Validating Genes

            # Correct those that are not found

            if(length((loop.section)[!is.nan(loop.section)]) > 0 &

               all(!is.finite(loop.section)) &

               is.nan(mean(as.vector(loop.section)[abs(loop.section)],

                           na.rm=TRUE))){

              validationMatrix[1, eval] <- "-"

            } else {

              validationMatrix[1, eval] <- "Found"

            }

          }

          # Storing the results in validationMList

          validationMatrix <-

            validationMatrix[,sort(colnames(validationMatrix)), drop=FALSE]

          idx <- ((group-1)*length(inputCases))+i

          validationMList[[idx]] <- validationMatrix

        }

      }

      # Update progressbar

      setTxtProgressBar(obtainOneStudyProgressBar, i)

    }

    # Closing progress bar

    close(obtainOneStudyProgressBar)




    ## bfc object

    # create bfc object

    if(!dir.exists(database)){

      bfc <- BiocFileCache(

        file.path(system.file("extdata", package = "cbaf"), submissionName)

      )

    }



    # Store the obtained Data

    number.of.rows.obtained.data <-

      nrow(bfcquery(bfc, "Obtained data for single study"))

    if(number.of.rows.obtained.data == 0){

      saveRDS(

        rawList,

        file=bfcnew(bfc, "Obtained data for single study", ext="RDS")

      )

    } else if(number.of.rows.obtained.data == 1){

      saveRDS(

        rawList,

        file=bfc[[bfcquery(bfc, "Obtained data for single study")$rid]]

      )

    }



    # Fill the Validation Result

    for(mix in seq_along(genesList)){

      validationResult[[mix]] <- do.call(

        "rbind",

        validationMList[((mix-1)*length(inputCases))+seq_along(inputCases)]

      )

    }

    # Store the validation data

    if(validateGenes){


      number.of.rows.validaion.data <-

        nrow(bfcquery(bfc, "Validation data for single study"))


      if(number.of.rows.validaion.data == 0){

        saveRDS(

          validationResult,

          file=bfcnew(bfc, "Validation data for single study", ext="RDS")

        )

      } else if(number.of.rows.validaion.data == 1){

        saveRDS(

          validationResult,

          file=bfc[[bfcquery(bfc, "Validation data for single study")$rid]]

        )

      }

    }

    # Store the last parameter

    newParameters$lastRunStatus <- "succeeded"

    oldParamObtainOneStudy <- newParameters


    # Store the parameters for this run

    number.of.rows.parameters <-

      nrow(bfcquery(bfc, "Parameters for obtainOneStudy()"))

    if(number.of.rows.parameters == 0){

      saveRDS(

        oldParamObtainOneStudy,

        file=bfcnew(bfc, "Parameters for obtainOneStudy()", ext="RDS")

      )

    } else if(number.of.rows.parameters == 1){

      saveRDS(

        oldParamObtainOneStudy,

        file=bfc[[bfcquery(bfc, "Parameters for obtainOneStudy()")$rid]]

      )

    }

  }

}
