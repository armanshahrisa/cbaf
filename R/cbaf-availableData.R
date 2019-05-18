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
#' Version: \tab 1.7.2 \cr
#' Date: \tab 2019-05-18 \cr
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

    choiceYesNo <- readline(prompt = "An excel file with the given name already exists in the working directory. Proceed anyway and overwrite the file? (yes/no):      ")

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



    # Terms associated with different techniques

    RNA.Seq_terms_L1 <- c(

      "Tumor Samples with mRNA data (RNA Seq V2)",

      "Tumors with mRNA data (RNA Seq V2)",

      "Samples with mRNA data (RNA Seq V2)",

      "Tumor Samples with mRNA data (RNA Seq)",

      "Tumors with mRNA data (RNA Seq)",

      "Samples with mRNA data (RNA Seq)"
    )

    RNA.Seq_terms_L2 <- c(

      "mRNA Expression z-Scores (RNA Seq V2 RSEM)",

      "mRNA expression z-Scores (RNA Seq V2 RSEM)",

      "mRNA Expression Zscores, RSEM (Batch normalized from Illumina HiSeq_RNASeqV2)",

      "mRNA Expression z-Scores (RNA Seq RPKM)",

      "mRNA expression z-Scores"
    )

    microRNA.Seq_terms_L1 <- c("Tumors with microRNA data (microRNA-Seq)",

                               "Tumor Samples with microRNA data (microRNA-Seq)"

                               )

    microRNA.Seq_terms_L2 <- c("microRNA expression Z-scores",

                               "mRNA Expression Z-Scores vs Normals"

                               )

    microarray_for_mRNA_terms_L1 <- c(

      "Tumor Samples with mRNA data (Agilent microarray)",

      "Tumors with mRNA data (Agilent microarray)",

      "Samples with mRNA data (Agilent microarray)",

      "Tumor Samples with mRNA data (U133 microarray only)",

      "Samples with mRNA data (U133 microarray)",

      "Tumors with mRNA data",

      "Tumors with mRNA"

    )

    microarray_for_mRNA_terms_L2 <- c(

      "mRNA Expression z-Scores (microarray)",

      "mRNA expression (microarray) z-scores",

      "mRNA Expression z-Scores (U133 microarray only)",

      "mRNA expression z-scores (Illumina)",

      "mRNA expression Z-scores (all genes)",

      "mRNA Expression Z-Scores vs Normals",

      "mRNA Expression z-Scores (combined microarray)",

      "mRNA Z-scores vs normal fat"

    )

    microarray_for_miRNA_terms_L1 <- "Tumors with microRNA"

    microarray_for_miRNA_terms_L2 <-

      c("mRNA Expression Z-Scores vs Normals",

        "mRNA Expression z-Scores (combined microarray)")

    methylation_terms_L1 <- c(

      "Tumor Samples with methylation data (HM450)",

      "Tumors with methylation data (HM450)",

      "Samples with methylation data (HM450)",

      "Tumor Samples with methylation data (HM27)",

      "Tumors with methylation data (HM27)",

      "Samples with methylation data (HM27)",

      "Tumors with methylation data",

      "Samples with methylation data"

    )

    methylation_terms_L2 <- c(

      "Methylation (HM450)",

      "Methylation (HM27)",

      "Methylation (hm27)",

      "Methylation"

    )



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

            any(RNA.Seq_terms_L1 %in% description)

          ),


          microRNA.Seq = as.character(

            any(microRNA.Seq_terms_L1 %in% description)

          ),


          microarray_of_mRNA = as.character(

            any(microarray_for_mRNA_terms_L1 %in% description)

          ),


          microarray_of_miRNA = as.character(

            any(microarray_for_miRNA_terms_L1 %in% description)

          ),


          methylation = as.character(

            any(methylation_terms_L1 %in% description))

          )


        } else{


          c(RNA.Seq = "FALSE",

            microRNA.Seq = "FALSE",

            microarray_of_mRNA = "FALSE",

            microarray_of_miRNA = "FALSE",

            methylation = "FALSE")


        }

      }else{

        c(RNA.Seq = "FALSE",

          microRNA.Seq = "FALSE",

          microarray_of_mRNA = "FALSE",

          microarray_of_miRNA = "FALSE",

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

              any(RNA.Seq_terms_L2 %in% description)

            ),


            microRNA.Seq = as.character(

              any(microRNA.Seq_terms_L2 %in% description)

            ),


            microarray_of_mRNA = as.character(

              any(microarray_for_mRNA_terms_L2 %in% description)

            ),


            microarray_of_miRNA = as.character(

              any(microarray_for_miRNA_terms_L2 %in% description)

            ),


            methylation = as.character(

              any(methylation_terms_L2 %in% description))

            )


          } else{


            c(RNA.Seq = "FALSE",

              microRNA.Seq = "FALSE",

              microarray_of_mRNA = "FALSE",

              microarray_of_miRNA = "FALSE",

              methylation = "FALSE")


          }

        }else{

          c(RNA.Seq = "FALSE",

            microRNA.Seq = "FALSE",

            microarray_of_mRNA = "FALSE",

            microarray_of_miRNA = "FALSE",

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

        "microarray_of_mRNAs" , "microarray_of_microRNAs", "methylation",

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
