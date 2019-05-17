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
#' Version: \tab 1.7.1 \cr
#' Date: \tab 2019-05-17 \cr
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

    list.of.studies <- getCancerStudies(mycgds)



    # Terms associated with different techniques

    RNA_Seq.terms <- c(

      "Tumor Samples with mRNA data (RNA Seq V2)",

      "Tumors with mRNA data (RNA Seq V2)",

      "Samples with mRNA data (RNA Seq V2)",

      "Tumor Samples with mRNA data (RNA Seq)",

      "Tumors with mRNA data (RNA Seq)",

      "Samples with mRNA data (RNA Seq)"
    )

    microRNA_Seq.terms <- c("Tumors with microRNA data (microRNA-Seq)",

                            "Tumor Samples with microRNA data (microRNA-Seq)")

    microarray.for.mRNA.term <- c(

      "Tumor Samples with mRNA data (Agilent microarray)",

      "Tumors with mRNA data (Agilent microarray)",

      "Samples with mRNA data (Agilent microarray)",

      "Tumor Samples with mRNA data (U133 microarray only)",

      "Samples with mRNA data (U133 microarray)",

      "Tumors with mRNA data",

      "Tumors with mRNA"

    )

    microarray.for.miRNA.term <- "Tumors with microRNA"

    methylation.term <- c(

      "Tumor Samples with methylation data (HM450)",

      "Tumors with methylation data (HM450)",

      "Samples with methylation data (HM450)",

      "Tumor Samples with methylation data (HM27)",

      "Tumors with methylation data (HM27)",

      "Samples with methylation data (HM27)",

      "Tumors with methylation data",

      "Samples with methylation data"

    )



    message("Cheching the available data for every cancer study")


    # create progress bar

    pb <- txtProgressBar(min = 0, max = nrow(list.of.studies), style = 3)

    i <- 0



    ############################################################################
    ## core segment

    # looking for supported technique data

    list.of.available.data <- sapply(

      list.of.studies[, "cancer_study_id"], function(cs, cgds) {

      # Obtain available techniques

      available.options <- getCaseLists(cgds, cs)



      # Update progressbar

      i <<- i + 1

      setTxtProgressBar(pb, i)


      if(length(available.options) > 1){

        if(any(colnames(available.options) == "case_list_name")){


          description <- available.options[, "case_list_name"]


          c(RNA.Seq = as.character(

            any(RNA_Seq.terms %in% description)

          ),


          microRNA.Seq = as.character(

            any(microRNA_Seq.terms %in% description)

          ),


          microarray_of_mRNA = as.character(

            any(microarray.for.mRNA.term %in% description)

          ),


          microarray_of_miRNA = as.character(

            any(microarray.for.miRNA.term %in% description)

          ),


          methylation = as.character(

            any(methylation.term %in% description))

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



    # Replacing True and False with available and ""

    list.of.available.data[list.of.available.data=="TRUE"] <- "available"

    list.of.available.data[list.of.available.data=="FALSE"] <- "-"



    # Overrule studies lacking z-score

    # RNA-Seq

    RNA_Seq_index <- which(colnames(list.of.available.data)
                           %in% c("paad_qcmg_uq_2016",
                                  "aml_target_2018_pub",
                                  "prad_mskcc_cheny1_organoids_2014",
                                  "utuc_cornell_baylor_mdacc_2019"))

    list.of.available.data[1, RNA_Seq_index] <- "-"


    # Microarray

    Microarray_index <- which(colnames(list.of.available.data)
                           %in% c("prad_broad_2013"))

    list.of.available.data[3, Microarray_index] <- "-"



    # joining list.of.available.data to list.of.studies

    combined.list <- cbind(list.of.studies[,"cancer_study_id"],

                           list.of.studies[,"name"], t(list.of.available.data),

                           list.of.studies[,"description"])

    colnames(combined.list) <-

      c("cancer_study_id", "cancer_study_name", "RNA.Seq", "microRNA.Seq",

        "microarray_of_mRNAs" , "microarray_of_microRNAs", "methylation",

        "description")

    rownames(combined.list) <- seq_len(nrow(combined.list))



    ############################################################################
    # Exporting results

    # Converting matrix to data.frame, and the store as an excel file

    combined.list.dataframe <- data.frame(combined.list)

    colnames(combined.list.dataframe) <-

      gsub("_", " ", colnames(combined.list.dataframe))

    colnames(combined.list.dataframe) <-

      gsub("\\.", "-", colnames(combined.list.dataframe))

    rownames(combined.list.dataframe) <- seq_len(nrow(combined.list))

    write.xlsx(

      combined.list.dataframe, file=paste(excelFileName, ".xlsx", sep = "")

      )

  }else{

    message("--- Function 'availableData()' was skipped. ---")

  }

}
