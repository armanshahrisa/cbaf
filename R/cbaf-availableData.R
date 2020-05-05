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
#' Version: \tab 1.10.1 \cr
#' Date: \tab 2020-05-05 \cr
#' License: \tab Artistic-2.0 \cr
#' }
#'
#'
#'
#' @importFrom cgdsr CGDS getCancerStudies getCaseLists getGeneticProfiles getProfileData
#'
#' @importFrom xlsx write.xlsx
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

    stop("'excelFileName' must be entered as a character string.")

  }



  if(file.exists(paste(excelFileName, ".xlsx", sep = ""))){

    choiceYesNo <- readline(prompt = "An excel file with the given name already exists. Rewrite? (yes/no):      ")

    if(choiceYesNo == "yes"){

      continue <- TRUE

    }else if(choiceYesNo == "no"){

      continue <- FALSE

    }else{

      stop("please type 'yes' or 'no'!")

    }

  } else{

    continue <- TRUE

  }



  if(continue){

    ############################################################################
    # Prerequisites

    # Prerequisites for cBioportal

    mycgds = CGDS("http://www.cbioportal.org/")

    list_of_studies <- getCancerStudies(mycgds)



    message("Cheching the available data for every cancer study")


    # create progress bar

    pb <- txtProgressBar(min = 0, max = nrow(list_of_studies)*2, style = 3)

    i <- 0



    ############################################################################
    ## core segment

    # looking for supported techniques at level 1

    list_of_available_data_L1 <- sapply(

      list_of_studies[, "cancer_study_id"], function(cs, cgds) {

      # Obtain available techniques

      available_options_1 <- getCaseLists(cgds, cs)



      # Update progressbar

      i <<- i + 1

      setTxtProgressBar(pb, i)


      if(length(available_options_1) > 1){

        if(any(colnames(available_options_1) == "case_list_name")){


          description <- available_options_1[, "case_list_name"]


          c(RNA.Seq = as.character(

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

            microRNA.Seq = "FALSE",

            microarray_with_mRNA_data = "FALSE",

            microarray_with_microRNA_data = "FALSE",

            methylation = "FALSE")


        }

      }else{

        c(RNA.Seq = "FALSE",

          microRNA.Seq = "FALSE",

          microarray_with_mRNA_data = "FALSE",

          microarray_with_microRNA_data = "FALSE",

          methylation = "FALSE")

      }


    }, mycgds)


    # looking for supported techniques at level 2

    list_of_available_data_L2 <- sapply(

      list_of_studies[, "cancer_study_id"], function(cs, cgds) {

        # Obtain available techniques

        available_options_2 <- getGeneticProfiles(cgds, cs)



        # Update progressbar

        i <<- i + 1

        setTxtProgressBar(pb, i)


        if(length(available_options_2) > 1){

          if(any(colnames(available_options_2) == "genetic_profile_name")){


            description <- available_options_2[, "genetic_profile_name"]


            c(RNA.Seq = as.character(

              any(RNA.Seq_L2.terms %in% description)

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

              microRNA.Seq = "FALSE",

              microarray_with_mRNA_data = "FALSE",

              microarray_with_microRNA_data = "FALSE",

              methylation = "FALSE")


          }

        }else{

          c(RNA.Seq = "FALSE",

            microRNA.Seq = "FALSE",

            microarray_with_mRNA_data = "FALSE",

            microarray_with_microRNA_data = "FALSE",

            methylation = "FALSE")

        }


      }, mycgds)


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

                           list_of_studies[,"name"], t(list_of_available_data_L1),

                           list_of_studies[,"description"])

    colnames(combined_list) <-

      c("cancer_study_id", "cancer_study_name", "RNA.Seq", "microRNA.Seq",

        "microarray_with_mRNA_data" , "microarray_with_microRNA_data", "methylation",

        "description")

    rownames(combined_list) <- seq_len(nrow(combined_list))



    ############################################################################
    # Exporting results

    # Converting matrix to data.frame, and the store as an excel file

    combined_list_dataframe <- data.frame(combined_list)

    colnames(combined_list_dataframe) <-

      gsub("_", " ", colnames(combined_list_dataframe))

    colnames(combined_list_dataframe) <-

      gsub("\\.", "-", colnames(combined_list_dataframe))

    rownames(combined_list_dataframe) <- seq_len(nrow(combined_list))

    write.xlsx(

      combined_list_dataframe, file=paste(excelFileName, ".xlsx", sep = "")

      )

  }else{

    message("--- Function 'availableData()' was skipped. ---")

  }

}
