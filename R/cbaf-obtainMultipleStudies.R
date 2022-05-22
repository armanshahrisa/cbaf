#' @title Obtain the requested data for various cancer studies.
#'
#' @description This function Obtains the requested data for the given genes
#' across multiple cancer studies. It can check whether or not all genes are
#' included in cancer studies and, if not, looks for the alternative gene names.
#'
#' @details
#' \tabular{lllll}{
#' Package: \tab cbaf \cr
#' Type: \tab Package \cr
#' Version: \tab 1.18.0 \cr
#' Date: \tab 2022-04-24 \cr
#' License: \tab Artistic-2.0 \cr
#' }
#'
#'
#'
#' @importFrom cBioPortalData cBioPortal getStudies sampleLists molecularProfiles getDataByGenes
#'
#' @importFrom BiocFileCache BiocFileCache bfcnew bfcquery bfcpath
#'
#' @importFrom utils head setTxtProgressBar txtProgressBar
#'
#'
#'
#' @usage obtainMultipleStudies(genesList, submissionName, studiesNames,
#'   desiredTechnique, cancerCode = FALSE, validateGenes = TRUE)
#'
#'
#'
#' @param genesList a list that contains at least one gene group
#'
#' @param submissionName a character string containing name of interest. It is
#' used for naming the process.
#'
#' @param studiesNames a character vector or a matrix that containes desired
#' cancer names. The character vector containes standard
#' names of cancer studies that can be found on cbioportal.org, such as
#' \code{"Acute Myeloid Leukemia (TCGA, NEJM 2013)"}. Alternatively, a matrix can
#' used if users prefer user-defined cancer names. In this case, the first
#' column of matrix comprises the standard cancer names while the second column
#' must contain the desired cancer names.
#'
#' @param desiredTechnique a character string that is one of the following
#' techniques: \code{"RNA-Seq"}, \code{"RNA-SeqRTN"}, \code{"microRNA-Seq"},
#' \code{"microarray.mRNA"}, \code{"microarray.microRNA"} or
#' \code{"methylation"}.
#'
#' @param cancerCode a logical value that tells the function to use cbioportal
#' abbreviated cancer names instead of complete cancer names, if set to be
#' \code{TRUE}. For example, \code{"laml_tcga_pub"} is the abbreviated name for
#' \code{"Acute Myeloid Leukemia (TCGA, NEJM 2013)"}.
#'
#' @param validateGenes a logical value that, if set to be \code{TRUE}, function
#'  will check each cancer study to find whether or not each gene has a record.
#'  If a cancer study doesn't have a record for specific gene, it checks for
#'  alternative gene names that cbioportal might use instead of the given gene
#'  name.
#'
#'
#'
#' @return a BiocFileCach object that contains the obtained data without further
#'  processing. Name of the object is combination of \code{bfc_} and
#'  \code{submissionName}. Inside it, there is a section for the obtained data,
#'  which is stored as a list. At first level, this list is subdivided into
#'  diferent groups based on the list of genes that user has given the function,
#'   then each gene group itself contains one matrix for every cancer study.
#'   Additonally, if \code{validateGenes = TRUE}, another section that contains
#'   gene validation results will be created in the BiocFileCach object.
#'
#'
#'
#' @examples
#' genes <- list(K.demethylases = c("KDM1A", "KDM1B", "KDM2A", "KDM2B", "KDM3A",
#'  "KDM3B", "JMJD1C", "KDM4A"), K.methyltransferases = c("SUV39H1", "SUV39H2",
#'  "EHMT1", "EHMT2", "SETDB1", "SETDB2", "KMT2A", "KMT2A"))
#'
#' studies <- c("Acute Myeloid Leukemia (TCGA, Provisional)",
#' "Adrenocortical Carcinoma (TCGA, Provisional)",
#' "Bladder Urothelial Carcinoma (TCGA, Provisional)",
#' "Brain Lower Grade Glioma (TCGA, Provisional)",
#' "Breast Invasive Carcinoma (TCGA, Provisional)")
#'
#' obtainMultipleStudies(genes, "test2", studies, "RNA-Seq")
#'
#' @author Arman Shahrisa, \email{shahrisa.arman@hotmail.com} [maintainer,
#' copyright holder]
#' @author Maryam Tahmasebi Birgani, \email{tahmasebi-ma@ajums.ac.ir}
#'
#' @export



################################################################################
################################################################################
################ Obtain the requested data for multiple Cancers ################
################################################################################
################################################################################

obtainMultipleStudies <- function(

  genesList,

  submissionName,

  studiesNames,

  desiredTechnique,

  cancerCode = FALSE,

  validateGenes = TRUE

){

  ##############################################################################
  ########## Prerequisites

  # Check genes

  if(!is.list(genesList)){

    # stop("[obtainMultipleStudies] 'genes' must be a list that contains at list one group of genes")

    if(is.vector(genesList)){

      genesList <- list(a = genesList)

    }

  }



  # Check submissionName

  if(!is.character(submissionName)){

    stop("[obtainMultipleStudies] 'submissionName' must be a character string!")

  }



  # Check studiesNames

  if(!is.character(studiesNames)){

    stop("[obtainMultipleStudies] 'studiesNames' must be character vector!")

  }



  # Check desiredTechnique

  if(is.character(desiredTechnique)){

    supported.techniques <- c("RNA-Seq",

                              "RNA-SeqRTN",

                              "microRNA-Seq",

                              "Microarray.mRNA",

                              "Microarray.microRNA",

                              "methylation")

    if(!(desiredTechnique %in% supported.techniques) |

       length(desiredTechnique)!= 1){

      stop("[obtainMultipleStudies] 'desiredTechnique' must be either 'RNA-Seq', 'microRNA-Seq', 'microarray.mRNA', 'microarray.microRNA' or 'methylation' !")

    } else if(desiredTechnique %in% supported.techniques |

              length(desiredTechnique)== 1){

      if(desiredTechnique == "RNA-Seq"){

        L1.characteristics <- RNA.Seq_L1.terms

        L2.characteristics <- RNA.Seq_L2.terms

      } else if(desiredTechnique == "RNA-SeqRTN"){

        L1.characteristics <- RNA.Seq_L1.terms

        L2.characteristics <- RNA.Seq_rtn_L2.terms

      } else if(desiredTechnique == "microRNA-Seq"){

        L1.characteristics <- microRNA.Seq_L1.terms

        L2.characteristics <- microRNA.Seq_L2.terms

      } else if(desiredTechnique == "microarray.mRNA"){

        L1.characteristics <- microarray.with.mRNA_L1.terms

        L2.characteristics <- microarray.with.mRNA_L2.terms

      } else if(desiredTechnique == "microarray.microRNA"){

        L1.characteristics <- microarray.with.microRNA_L1.terms

        L2.characteristics <- microarray.with.microRNA_L2.terms

      } else if(desiredTechnique == "methylation"){

        L1.characteristics <- methylation_L1.terms

        L2.characteristics <- methylation_L2.terms

      }

    }

  } else {

    stop("[obtainMultipleStudies] 'desiredTechnique' must be a character string!")

  }



  # Check cancerCode

  if(!is.logical(cancerCode)){

    stop("[obtainMultipleStudies] 'cancerCode' must be either TRUE or FALSE!")

  }



  # Check validateGenes

  if(!is.logical(validateGenes)){

    stop("[obtainMultipleStudies] 'validateGenes' must be either TRUE or FALSE!")

  }





  ##############################################################################
  ########## Decide whether function should stops now!

  # Set cgdsr, stop if submissionName is either "test" or "test2" to
  # improving package test speed

  if(!(submissionName %in% c("test", "test2"))){

    cbio <- cBioPortal()

    #! mycgds = CGDS("http://www.cbioportal.org/")


    # Getting cancer names

    studies <- getStudies(cbio)

    supportedCancers <- studies

    #! supportedCancers_old <- getCancerStudies(mycgds)

  }


  # Check if all studiesNames are included in supportedCancers

  if(!(submissionName %in% c("test", "test2"))){

    if(any(studiesNames %in% supportedCancers$name)){

      studiesNames <- studiesNames[studiesNames %in% supportedCancers$name]

    }else{

      stop("[obtainMultipleStudies] None of the requested cancer studies is supported!")

    }

  }


  # Store the new parameteres

  newParameters <-list()

  newParameters$genesList <- genesList

  newParameters$submissionName <- submissionName

  newParameters$studyName <- studiesNames

  newParameters$desiredTechnique <- desiredTechnique

  newParameters$cancerCode <- cancerCode

  newParameters$validateGenes <- validateGenes





  # Check the database

  database <-

    system.file("extdata", submissionName, package="cbaf")



  # Remove old database

  if(dir.exists(database) & !(submissionName %in% c("test", "test2"))){

    creation.time <- file.info(database, extra_cols = FALSE)$ctime

    past.time <-

      as.numeric(difftime(Sys.time(), creation.time, units = c("days")))

    if(past.time >= 5){

      message("[obtainMultipleStudies] The downloaded data are outdated!")

      message("[obtainMultipleStudies] Removing the previous data.")

      unlink(database, recursive = TRUE)

    }

  }



  # Check wheather the requested data exists

  if(dir.exists(database)){

    bfc <- BiocFileCache(

      file.path(system.file("extdata", package = "cbaf"), submissionName),

      ask = FALSE

    )

    if(nrow(bfcquery(bfc, "Parameters for obtainMultipleStudies()")) == 1){

      oldParameters <- readRDS(

        bfcpath(

          bfc,

          bfcquery(bfc, c("Parameters for obtainMultipleStudies()"))$rid

        )

      )

      if(identical(oldParameters[-7], newParameters) |

         submissionName %in% c("test", "test2")){

        continue <- FALSE

        # Store the last parameter

        newParameters$lastRunStatus <- "skipped"

        oldParamObtainMultipleStudies <- newParameters

        saveRDS(

          oldParamObtainMultipleStudies,

          file=bfc[[bfcquery(bfc, "obtainMultipleStudies()")$rid]]

        )

        if(submissionName %in% c("test", "test2")){

          message("[obtainMultipleStudies] Please choose a name other than 'test' and 'test2'.")

        }

        message("[obtainMultipleStudies] The requested data already exist locally.")

        message("[obtainMultipleStudies] The function was haulted!")

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

    # Creating a vector for cancer names and subsequent list subset name

    if(!is.matrix(studiesNames)){

      studiesNames <- studiesNames[order(studiesNames)]

      studiesNamesMatrix <- studiesNames


      if(cancerCode){

        cancer.studies.idx <-

          which(supportedCancers$name %in% as.character(studiesNames))

        groupNames <- supportedCancers$studyId[cancer.studies.idx]

      } else if(!cancerCode){

        groupNames <- as.character(studiesNames)

      }


    } else if(is.matrix(studiesNames)){

      studiesNamesMatrix <- studiesNames[order(studiesNames[,1]),]

      studiesNames <- studiesNamesMatrix[,2]

      groupNames <- studiesNames

    }



    ############################################################################
    ########## Core segment

    # Report

    message("[obtainMultipleStudies] Downloading the required data.")



    # create progress bar

    obtainMultipleStudiesProgressBar <-

      txtProgressBar(min = 0, max = length(studiesNames), style = 3)



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

    validationMList <- vector("list", length(genesList)*length(studiesNames))

    # List of cancers with corrupted data

    cancersWithCorruptedData <- NULL

    # List of cancers lacking approprate data

    cancersLackingData <- NULL





    ## Getting the required gene expresssion profile ...

    # 'for' control structure for obtaining data and calculating the parameters

    for(c in seq_along(studiesNames)){

      # Determining name for list subset of study name

      groupName <- groupNames[c]

      # Correcting possible errors of list names

      groupName <- gsub(

        groupName, pattern = "\\+ ", replacement = " possitive ",

        ignore.case = TRUE

      )

      groupName <- gsub(

        groupName, pattern = "\\- ", replacement = " negative ",

        ignore.case = TRUE

      )





      # Find cancer abbreviated name

      CancerStudies.idx <-

        which(supportedCancers$name == as.character(studiesNames[c]))

      mycancerstudy = as.character(

        supportedCancers[CancerStudies.idx, "studyId"]

      )



      # Finding the first characteristics of data in the cancer

      AvailableDataTypes <- sampleLists(cbio, mycancerstudy)

      #! AvailableDataTypes_old <- getCaseLists(mycgds, mycancerstudy)

      if(length(AvailableDataTypes) > 1){

        # Available terms and priotizing based on cbaf-constants.R

        match_index_1 <- match(L1.characteristics, AvailableDataTypes$name)

        f.condition <- AvailableDataTypes$name[match_index_1]

        f.condition <- f.condition[!is.na(f.condition)]

        if(length(f.condition) >= 1){

          f.condition <- f.condition[1]

          f.condition.idx <- which(AvailableDataTypes$name == f.condition)

          mycaselist <- as.character(

            AvailableDataTypes[f.condition.idx ,"sampleListId"]

          )


          CancerPossessData <- TRUE

        } else if(length(f.condition) == 0){

          CancerPossessData <- FALSE

        }

        CancerPossessCorruptedData <- FALSE

      } else{

        CancerPossessCorruptedData <- TRUE

        # Updating cancers list

        if(is.null(cancersWithCorruptedData)){

          cancersWithCorruptedData <- studiesNames[c]

        }else{

          cancersWithCorruptedData <-
            c(cancersWithCorruptedData, studiesNames[c])

        }

      }





      # Finding the second characteristics of data in the cancer

      if(!CancerPossessCorruptedData){

        if(CancerPossessData){

          AvailableDataFormats <- molecularProfiles(cbio, mycancerstudy)

          #! AvailableDataFormats_old <- getGeneticProfiles(mycgds, mycancerstudy)

          match_index_2 <- match(L2.characteristics, AvailableDataFormats$name)

          s.condition <- AvailableDataFormats$name[match_index_2]

          s.condition <- s.condition[!is.na(s.condition)]



          if(length(s.condition) >= 1){

            s.condition <- s.condition[1]

            s.condition.idx <- which(AvailableDataFormats$name == s.condition)

            mygeneticprofile <- as.character(

              AvailableDataFormats[s.condition.idx ,"molecularProfileId"]

            )

            CancerPossessData <- TRUE

          } else if (length(s.condition) == 0){

            CancerPossessData <- FALSE

          }

        }

      }


      # Updating cancers list

      if(!CancerPossessData){

        if(is.null(cancersLackingData)){

          cancersLackingData <- studiesNames[c]

        }else{

          cancersLackingData <-
            c(cancersLackingData, studiesNames[c])

        }

      }



      # obtaining data for every genegroup

      for(group in seq_along(genesList)){

        # Chose one group of genes

        genesNames <- unique(genesList[[group]])

        numberOfGenes <- length(genesNames)

        # Merging four constitutive genes with data in case all genes are NA

        geneNames_plus_constitutive_genes <- c(genesNames, constitutive_genes)

        order_index <- order(geneNames_plus_constitutive_genes)

        ordered_genesNames <- geneNames_plus_constitutive_genes[order_index]

        number_Of_OrderedGenes <- length(ordered_genesNames)


        # Obtaining Expression z-scores for the requested genes

        # Check number of genes first

        if(number_Of_OrderedGenes <= 250){

          if(!CancerPossessCorruptedData){

            if(CancerPossessData){

              #! ProfileData_old <- getProfileData(

              #!  mycgds,

              #!  genesNames[order(genesNames)],

              #!  mygeneticprofile,

              #!  mycaselist

              #! )

              # Obtaining Gene Data

              Unprocessed_ProfileData_list <-

                getDataByGenes(

                  cbio,

                  studyId = mycancerstudy,

                  genes = ordered_genesNames,

                  by = "hugoGeneSymbol",

                  sampleListId = mycaselist,

                  molecularProfileIds = mygeneticprofile
                )

              # Extracting data.frame from List

              Unprocessed_ProfileData <- Unprocessed_ProfileData_list[[1]]

              # Subseting data.frame to contain the needed Columns

              Filtered_Unprocessed_ProfileData <-

                Unprocessed_ProfileData[,c("sampleId", "value", "hugoGeneSymbol")]

              # Spliting data.frame by sampleId

              patient_genes_list <- split(Filtered_Unprocessed_ProfileData,

                                          Filtered_Unprocessed_ProfileData$sampleId)

              # Making each table in the list a one column table

              for(hugo in seq_len(length(patient_genes_list))){


                present_gene_table <- patient_genes_list[[hugo]]

                present_gene_table_2 <-

                  present_gene_table[,"value", drop = FALSE]

                # Converting to matrix
                present_gene_matrix <- as.matrix(present_gene_table_2)

                # Giving gene names and patient id to the values
                colnames(present_gene_matrix) <-

                  names(patient_genes_list)[hugo]

                rownames(present_gene_matrix) <-

                  present_gene_table$hugoGeneSymbol

                patient_genes_list[[hugo]] <- t(present_gene_matrix)


              }


              # Generating old ProfileData format by collapsing the list
              ProfileData <- do.call(rbind, patient_genes_list)

              # Sorting the ProfileData by column and rown names
              ProfileData <- ProfileData[,order(colnames(ProfileData))]

              ProfileData <- ProfileData[order(rownames(ProfileData)),]

            }

          }

        }else{

          # split genes in groups of 250 names

          operational_gene_number <- split(

            ordered_genesNames, ceiling(seq_len(numberOfOrderedGenes)/250)

          )


          # Create empty list for gene_matrices

          separated_results <- vector(

            "list", length = length(operational_gene_number)

          )

          for(operational in seq_along(operational_gene_number)){

            #! separated_results[[operational]] <- getProfileData(

            #!   mycgds,

            #!   operational_gene_number[[operational]],

            #!   mygeneticprofile,

            #!   mycaselist

            #! )

            Unprocessed_ProfileData_list <-

              getDataByGenes(

                cbio,

                studyId = mycancerstudy,

                genes = operational_gene_number[[operational]],

                by = "hugoGeneSymbol",

                sampleListId = mycaselist,

                molecularProfileIds = mygeneticprofile
              )

            # Extracting data.frame from List

            Unprocessed_ProfileData <- Unprocessed_ProfileData_list[[1]]

            # Subseting data.frame to contain the needed Columns

            Filtered_Unprocessed_ProfileData <-

              Unprocessed_ProfileData[,c("sampleId", "value", "hugoGeneSymbol")]

            # Spliting data.frame by sampleId

            patient_genes_list <- split(Filtered_Unprocessed_ProfileData,

                                        Filtered_Unprocessed_ProfileData$sampleId)

            # Making each table in the list a one column table

            for(hugo in seq_len(length(patient_genes_list))){


              present_gene_table <- patient_genes_list[[hugo]]

              present_gene_table_2 <-

                present_gene_table[,"value", drop = FALSE]

              # Converting to matrix
              present_gene_matrix <- as.matrix(present_gene_table_2)

              # Giving gene names and patient id to the values
              colnames(present_gene_matrix) <-

                names(patient_genes_list)[hugo]

              rownames(present_gene_matrix) <-

                present_gene_table$hugoGeneSymbol

              patient_genes_list[[hugo]] <- t(present_gene_matrix)


            }


            # Generating old ProfileData format by collapsing the list
            ProfileData <- do.call(rbind, patient_genes_list)

            # Sorting the ProfileData by column and rown names
            ProfileData <- ProfileData[,order(colnames(ProfileData))]

            ProfileData <- ProfileData[order(rownames(ProfileData)),]

            separated_results[[operational]] <- ProfileData


          }


          # Merging data

          ProfileData <- do.call("cbind", separated_results)

          ProfileData <- ProfileData[,order(colnames(ProfileData))]

        }


        # Check if all requested genes are present and remove four constitutive genes

        presence_index <- colnames(ProfileData) %in% constitutive_genes

        number_of_present_constitutive_genes <- sum(presence_index)

        # Determine the first constitutive gene for gene validation

        first_constitutive_gene <- (colnames(ProfileData)[presence_index])[1]


        if(ncol(ProfileData) <= number_of_present_constitutive_genes){

          stop("None of the requested genes is in database!")

        }else{

          ProfileData <-

            ProfileData[,!colnames(ProfileData) %in% constitutive_genes]

        }






        if(!CancerPossessCorruptedData){

          if(CancerPossessData){

            # Assaign data to specific list member

            rawList[[group]][[c]] <- data.matrix(ProfileData)

            names(rawList[[group]])[c] <- groupName

            # For convenience

            this.segment <- rawList[[group]][[c]]



            # Find whether alternative gene names are used

            # Alter c.genes to be compatible with gene names in cBioPortal
            # output

            alteredGeneNames <- sort(gsub("-", ".", genesNames))

            # Obtain unique name of genes that are absent in requested cancer

            presence_index_2 <-

              unique(alteredGeneNames) %in% colnames(this.segment)

            absentGenes <- alteredGeneNames[!presence_index_2]

            # For loop for determining changed genes

            if(length(absentGenes) != 0){

              alternativeGeneNames <-

                vector("character", length = length(absentGenes))

              # For loop

              for(ab in seq_along(absentGenes)){

                #! absent.gene.profile.data <- getProfileData(

                #!   mycgds, absentGenes[ab], mygeneticprofile, mycaselist

                #!)

                OneAbsentGene_OneContutiveGenes <-

                  c(absentGenes[ab], first_constitutive_gene)


                Unprocessed_ProfileData_list <-

                  getDataByGenes(

                    cbio,

                    studyId = mycancerstudy,

                    genes = OneAbsentGene_OneContutiveGenes,

                    by = "hugoGeneSymbol",

                    sampleListId = mycaselist,

                    molecularProfileIds = mygeneticprofile
                  )

                # Extracting data.frame from List

                Unprocessed_ProfileData <- Unprocessed_ProfileData_list[[1]]

                # Subseting data.frame to contain the needed Columns

                Filtered_Unprocessed_ProfileData <-

                  Unprocessed_ProfileData[,c("sampleId", "value", "hugoGeneSymbol")]

                # Spliting data.frame by sampleId

                patient_genes_list <- split(Filtered_Unprocessed_ProfileData,

                                            Filtered_Unprocessed_ProfileData$sampleId)

                # Making each table in the list a one column table

                for(hugo in seq_len(length(patient_genes_list))){


                  present_gene_table <- patient_genes_list[[hugo]]

                  present_gene_table_2 <-

                    present_gene_table[,"value", drop = FALSE]

                  # Converting to matrix
                  present_gene_matrix <- as.matrix(present_gene_table_2)

                  # Giving gene names and patient id to the values
                  colnames(present_gene_matrix) <-

                    names(patient_genes_list)[hugo]

                  rownames(present_gene_matrix) <-

                    present_gene_table$hugoGeneSymbol

                  patient_genes_list[[hugo]] <- t(present_gene_matrix)


                }


                # Generating old ProfileData format by collapsing the list
                absent.gene.profile.data <- do.call(rbind, patient_genes_list)

                # Sorting the ProfileData by column and rown names
                absent.gene.profile.data <-

                  absent.gene.profile.data[,order(colnames(absent.gene.profile.data)), drop = FALSE]

                absent.gene.profile.data <-

                  absent.gene.profile.data[order(rownames(absent.gene.profile.data)),, drop = FALSE]



                absentGeneProfileData <- colnames(

                  data.matrix(absent.gene.profile.data)

                )

                absentGeneProfileData_pure <-

                  absentGeneProfileData[!absentGeneProfileData %in% constitutive_genes]


                # Check whether gene has an alternative name or missed from the

                # database

                if(length(absentGeneProfileData_pure) == 1){

                  alternativeGeneNames[ab] <- absentGeneProfileData_pure

                } else if(length(absentGeneProfileData_pure) == 0){

                  alternativeGeneNames[ab] <- "-"

                }

              }

              # Naming Alternative.gene.names

              names(alternativeGeneNames) <- absentGenes

              # Seperating genes with alternative names from those that lack

              genesLackData <- alternativeGeneNames[alternativeGeneNames == "-"]

              genesWithData <- alternativeGeneNames[alternativeGeneNames != "-"]



              # modifying gene names containing an alternative name

              for(re in seq_along(genesWithData)){

                colnames.idx <-

                  colnames(rawList[[group]][[c]]) %in% genesWithData[re]


                colnames(rawList[[group]][[c]])[colnames.idx] <-

                  paste0(genesWithData[re], " (", names(genesWithData[re]), ")")

              }


            }else{

              genesLackData <- NULL

              genesWithData <- NULL

            }





            # validateGenes

            if(validateGenes){

              # Empty validation matrix

              validationMatrix <- matrix(, ncol = numberOfGenes, nrow = 1)

              # Naming empty matrix

              if(length(genesLackData) != 0){

                dimnames(validationMatrix) <- list(

                  groupNames[c],

                  c(colnames(this.segment), names(genesLackData))

                )

              } else{

                dimnames(validationMatrix) <-

                  list(groupNames[c], colnames(this.segment))

              }



              # modifying gene names containing an alternative name

              if(length(genesWithData) != 0){

                for(re in seq_along(genesWithData)){

                  colnames.idx <-

                    colnames(validationMatrix) %in% genesWithData[re]


                  colnames(validationMatrix)[colnames.idx] <-

                    paste0(genesWithData[re],
                           " (",
                           names(genesWithData[re]), ")")

                }

              }





              # Puting value for genes lacking data

              colnames.idx <-
                colnames(validationMatrix) %in% names(genesLackData)

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

              idx <- ((group-1)*length(studiesNames))+c

              validationMList[[idx]] <- validationMatrix

            }

          }

        }

      }

      # Update progressbar

      setTxtProgressBar(obtainMultipleStudiesProgressBar, c)

    }

    # Closing progress bar

    close(obtainMultipleStudiesProgressBar)





    # Print studies with Corrupted or absent data

    if(!is.null(cancersWithCorruptedData)){

      if(length(cancersWithCorruptedData) == 1){

        message("[obtainMultipleStudies] Corrupted / under modification data:")

      } else{

        message("[obtainMultipleStudies] Corrupted / being modification data:")

      }

      print(cancersWithCorruptedData)

    }



    if(!is.null(cancersLackingData)){

      if(length(cancersLackingData) == 1){

        message("[obtainMultipleStudies] Following study lacks ", desiredTechnique, " Data:")

      } else{

        message("[obtainMultipleStudies] Following study contains corrupted Data:")

      }

      print(cancersLackingData)

    }





    # Check if any cancer has data

    StudiesWithResults_idx <-
      !(studiesNames %in% c(cancersWithCorruptedData, cancersLackingData))

    StudiesWithResults <- studiesNames[StudiesWithResults_idx]

    if( length(StudiesWithResults) < 1 ){

      stop("[obtainMultipleStudies] No cancer study exists with ", desiredTechnique, " data or data for all studies are completely corrupted!")

    }





    ## bfc object

    # create bfc object

    if(!dir.exists(database)){

      bfc <- BiocFileCache(

        file.path(system.file("extdata", package = "cbaf"), submissionName),

        ask = FALSE

      )

    }



    # Remove NA objects from List objects

    for(MainList in seq_along(rawList)){

      rawList[[MainList]] <-

        Filter(function(a) any(!is.na(a)), rawList[[MainList]])

    }


    validationMList <- Filter(function(a) any(!is.na(a)), validationMList)





    # Store the obtained Data

    number.of.rows.obtained.data <-

      nrow(bfcquery(bfc, "Obtained data for multiple studies"))

    if(number.of.rows.obtained.data == 0){

      saveRDS(

        rawList,

        file=bfcnew(bfc, "Obtained data for multiple studies", ext="RDS")

      )

    } else if(number.of.rows.obtained.data == 1){

      saveRDS(

        rawList,

        file=bfc[[bfcquery(bfc, "Obtained data for multiple studies")$rid]]

      )

    }





    # Fill the Validation Result

    for(mix in seq_along(genesList)){

      validationResult[[mix]] <- do.call(

        "rbind",

        validationMList[((mix-1)*length(StudiesWithResults))+seq_along(StudiesWithResults)]

      )

    }

    # Store the validation data

    if(validateGenes){


      number.of.rows.validaion.data <-

        nrow(bfcquery(bfc, "Validation data for multiple studies"))


      if(number.of.rows.validaion.data == 0){

        saveRDS(

          validationResult,

          file=bfcnew(bfc, "Validation data for multiple studies", ext="RDS")

        )

      } else if(number.of.rows.validaion.data == 1){

        saveRDS(

          validationResult,

          file=bfc[[bfcquery(bfc, "Validation data for multiple studies")$rid]]

        )

      }

    }

    # Store the last parameter

    newParameters$lastRunStatus <- "succeeded"

    oldParamObtainMultipleStudies <- newParameters

    CorrectCnacer_idx <-

      !(studiesNames %in% c(cancersWithCorruptedData, cancersLackingData))

    newParameters$studyName <- studiesNames[CorrectCnacer_idx]


    # Store the parameters for this run

    number.of.rows.parameters <-

      nrow(bfcquery(bfc, "Parameters for obtainMultipleStudies()"))

    if(number.of.rows.parameters == 0){

      saveRDS(

        oldParamObtainMultipleStudies,

        file=bfcnew(bfc, "Parameters for obtainMultipleStudies()", ext="RDS")

      )

    } else if(number.of.rows.parameters == 1){

      saveRDS(

        oldParamObtainMultipleStudies,

        file=bfc[[bfcquery(bfc, "Parameters for obtainMultipleStudies()")$rid]]

      )

    }

    # message("[obtainMultipleStudies] Finished.")

  }

}
