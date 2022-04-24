#' @title Check which Data types are available for each cancer study.
#'
#' @description This function checks all the cancer studies that are registered
#' in 'cbioportal.org' to examine whether or not they contain RNA-Seq,
#' microRNA-Seq, microarray(mRNA), microarray(miRNA) and methylation data.
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
#' @importFrom cBioPortalData cBioPortal getStudies sampleLists molecularProfiles
#'
#' @importFrom openxlsx createWorkbook addWorksheet writeData saveWorkbook
#'
#' @importFrom utils head setTxtProgressBar txtProgressBar
#'
#'
#'
#' @usage availableData(excelFileName)
#'
#' @param excelFileName a character string that is required to name the output
#' and, if requested, excel file.
#'
#'
#' @return An excel file that contains all the cancer studies versus available
#' data types
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
################ Dataset availability in all cBioportal Cancers ################
################################################################################
################################################################################

availableData <- function(excelFileName){


  ##############################################################################
  # Check input parameters

  if(!is.character(excelFileName)){

    stop("[availableData] 'excelFileName' must be a character string!")

  }



  if(file.exists(paste(excelFileName, ".xlsx", sep = ""))){

    message("[availableData] Warning! '", excelFileName, ".xlsx", "' already exists!")

    choiceYesNo <- readline(prompt = "[availableData] Overwrite the file? (yes/no): ")

    if(choiceYesNo == "yes"){

      # Remove the previous file

      file.remove(paste(excelFileName, ".xlsx", sep = ""))

      continue <- TRUE

    }else if(choiceYesNo == "no"){

      continue <- FALSE

    }else{

      stop("[availableData] please type 'yes' or 'no'!")

    }

  } else{

    continue <- TRUE

  }



  if(continue){

    ############################################################################
    # Prerequisites

    # Prerequisites for cBioportal

    cbio <- cBioPortal()

    studies <- getStudies(cbio)

    # Converting new format to the old format

    colnames(studies)[colnames(studies) == "studyId"] <- "cancer_study_id"

    list_of_studies <- studies

    #! mycgds = CGDS("http://www.cbioportal.org/")

    #! list_of_studies <- getCancerStudies(mycgds)



    message("[availableData] Checking all cancer studies")


    # create progress bar

    pb <- txtProgressBar(min = 0, max = nrow(list_of_studies)*2, style = 3)

    i <- 0



    ############################################################################
    ## core segment

    # looking for supported techniques at level 1

    list_of_available_data_L1 <- sapply(

      list_of_studies$cancer_study_id, function(cs, cbio_api) {


      # Obtain available techniques

      samp <- sampleLists(cbio_api, cs)


      # Converting new format to the old format

      colnames(samp)[colnames(samp) == "sampleListId"] <- "case_list_id"

      colnames(samp)[colnames(samp) == "name"] <- "case_list_name"

      available_options_1 <- samp

      #! available_options_1 <- getCaseLists(cgds, cs)



      # Update progressbar

      i <<- i + 1

      setTxtProgressBar(pb, i)


      if(length(available_options_1) > 1){

        if(any(colnames(available_options_1) == "case_list_name")){


          description <- available_options_1$case_list_name


          c(RNA.Seq = as.character(

            any(RNA.Seq_L1.terms %in% description)

          ),


          RNA.Seq.RTN = as.character(

            any(RNA.Seq_L1.terms %in% description)

          ),


          microRNA.Seq = as.character(

            any(microRNA.Seq_L1.terms %in% description)

          ),


          microarray_with_mRNA_data = as.character(

            any(microarray.with.mRNA_L1.terms %in% description)

          ),


          microarray_with_microRNA_data = as.character(

            any(microarray.with.microRNA_L1.terms %in% description)

          ),


          methylation = as.character(

            any(methylation_L1.terms %in% description))

          )


        } else{


          c(RNA.Seq = "FALSE",

            RNA.Seq.RTN = "FALSE",

            microRNA.Seq = "FALSE",

            microarray_with_mRNA_data = "FALSE",

            microarray_with_microRNA_data = "FALSE",

            methylation = "FALSE")


        }

      }else{

        c(RNA.Seq = "FALSE",

          RNA.Seq.RTN = "FALSE",

          microRNA.Seq = "FALSE",

          microarray_with_mRNA_data = "FALSE",

          microarray_with_microRNA_data = "FALSE",

          methylation = "FALSE")

      }


    }, cbio ) #! mycgds


    # looking for supported techniques at level 2

    list_of_available_data_L2 <- sapply(

      list_of_studies$cancer_study_id, function(cs, cbio_api) {

        # Obtain available techniques

        mols <- molecularProfiles(cbio_api, cs)

        # Converting new format to the old format

        colnames(mols)[colnames(mols) == "name"] <- "genetic_profile_name"

        colnames(mols)[colnames(mols) == "molecularProfileId"] <-

                                                    "genetic_profile_id"

        available_options_2 <- mols

        #! available_options_2 <- getGeneticProfiles(cgds, cs)



        # Update progressbar

        i <<- i + 1

        setTxtProgressBar(pb, i)


        if(length(available_options_2) > 1){

          if(any(colnames(available_options_2) == "genetic_profile_name")){


            description <- available_options_2$genetic_profile_name


            c(RNA.Seq = as.character(

              any(RNA.Seq_L2.terms %in% description)

            ),


            RNA.Seq.RTN = as.character(

              any(RNA.Seq_rtn_L2.terms %in% description)

            ),


            microRNA.Seq = as.character(

              any(microRNA.Seq_L2.terms %in% description)

            ),


            microarray_with_mRNA_data = as.character(

              any(microarray.with.mRNA_L2.terms %in% description)

            ),


            microarray_with_microRNA_data = as.character(

              any(microarray.with.microRNA_L2.terms %in% description)

            ),


            methylation = as.character(

              any(methylation_L2.terms %in% description))

            )


          } else{


            c(RNA.Seq = "FALSE",

              RNA.Seq.RTN = "FALSE",

              microRNA.Seq = "FALSE",

              microarray_with_mRNA_data = "FALSE",

              microarray_with_microRNA_data = "FALSE",

              methylation = "FALSE")


          }

        }else{

          c(RNA.Seq = "FALSE",

            RNA.Seq.RTN = "FALSE",

            microRNA.Seq = "FALSE",

            microarray_with_mRNA_data = "FALSE",

            microarray_with_microRNA_data = "FALSE",

            methylation = "FALSE")

        }


      }, cbio) #! mycgds


    # close progressbar

    close(pb)



    # Find double positives

    for(col_number in seq_len(ncol(list_of_available_data_L1))){

      for(row_number in seq_len(nrow(list_of_available_data_L1))){

        current_L1 <- list_of_available_data_L1[row_number, col_number]

        if(current_L1 == "TRUE"){

          current_L2 <- list_of_available_data_L2[row_number, col_number]

          if(current_L2 == "FALSE"){

            current_L1 <- FALSE

            list_of_available_data_L1[row_number, col_number] <- current_L1
          }

        }

      }

    }




    # Replacing True and False with available and ""

    list_of_available_data_L1[list_of_available_data_L1=="TRUE"] <- "available"

    list_of_available_data_L1[list_of_available_data_L1=="FALSE"] <- "-"





    # joining list.of.available.data to list.of.studies

    combined_list <- cbind(list_of_studies[,"cancer_study_id"],

                           list_of_studies[,"name"],

                           t(list_of_available_data_L1),

                           list_of_studies[,"allSampleCount"],

                           list_of_studies[,"pmid"],

                           list_of_studies[,"description"])

    colnames(combined_list) <-

      c("Cancer_Study_ID", "Cancer_Study_Name", "RNA.Seq", "RNA.Seq.RTN",

        "microRNA.Seq", "microarray_with_mRNA_data" ,

        "microarray_with_microRNA_data", "methylation", "Number of All Samples",

        "PMID", "Description")

    rownames(combined_list) <- seq_len(nrow(combined_list))



    ############################################################################
    # Exporting results

    # Converting matrix to data.frame, and the store as an excel file

    combined_list_dataframe <- data.frame(combined_list)

    colnames(combined_list_dataframe) <-

      gsub("_", " ", colnames(combined_list_dataframe))

    colnames(combined_list_dataframe) <-

      gsub("\\.", "-", colnames(combined_list_dataframe))

    colnames(combined_list_dataframe) <-

      gsub("\\-RTN", " (RTN)", colnames(combined_list_dataframe))

    rownames(combined_list_dataframe) <- seq_len(nrow(combined_list))

    # Store Xlsx file

    ad <- createWorkbook()

    addWorksheet(ad, sheetName = "Available Data")

    writeData(ad,

              sheet = "Available Data",

              x = combined_list_dataframe,

              rowNames = TRUE)

    saveWorkbook(ad, file=paste(excelFileName, ".xlsx", sep = ""))

    # message("[availableData] Finished.")

    message(c("[availableData] The output was stored as '",  excelFileName, ".xlsx","'."))

  }else{

    message("[availableData] Function was haulted!")

  }

}
