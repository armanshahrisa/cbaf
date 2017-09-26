#' @title Check which Data types are available for each cancer study.
#'
#' @description This function checks all the cancer studies that are registered
#' in 'cbioportal.org' to examine whether or not they contain RNA-seq,
#' microRNA-seq, microarray(mRNA), microarray(miRNA) and methylation data.
#'
#' @details
#' \tabular{lllll}{
#' Package: \tab cbaf \cr
#' Type: \tab Package \cr
#' Version: \tab 1.0.0 \cr
#' Date: \tab 2017-09-26 \cr
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
#' @return A matrix that contains all cancer studies versus available data types
#' . It is available as an excel file with the given name.
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

    RNA_seq.terms <- c(

      "Tumor Samples with mRNA data (RNA Seq V2)",

      "Tumors with mRNA data (RNA Seq V2)",

      "Tumor Samples with mRNA data (RNA Seq)",

      "Tumors with mRNA data (RNA Seq)"
    )

    microRNA_seq.terms <- "Tumors with microRNA data (microRNA-Seq)"

    microarray.for.mRNA.term <- c(

      "Tumor Samples with mRNA data (Agilent microarray)",

      "Tumors with mRNA data (Agilent microarray)",

      "Tumor Samples with mRNA data (U133 microarray only)",

      "Tumors with mRNA data"

    )

    microarray.for.miRNA.term <- "Tumors with microRNA"

    methylation.term <- c(

      "Tumor Samples with methylation data (HM450)",

      "Tumors with methylation data (HM450)",

      "Tumor Samples with methylation data (HM27)",

      "Tumors with methylation data (HM27)",

      "Tumors with methylation data"

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



      if(any(colnames(available.options) == "case_list_name")){


        description <- available.options[, "case_list_name"]


        c(RNA.seq = as.character(

          any(RNA_seq.terms %in% description)

          ),


          microRNA.seq = as.character(

            any(microRNA_seq.terms %in% description)

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


        c(RNA.seq = "FALSE",

          microRNA.seq = "FALSE",

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



    # joining list.of.available.data to list.of.studies

    combined.list <- cbind(list.of.studies[,"cancer_study_id"],

                           list.of.studies[,"name"], t(list.of.available.data),

                           list.of.studies[,"description"])

    colnames(combined.list) <-

      c("cancer_study_id", "cancer_study_name", "RNA.seq", "microRNA.seq",

        "microarray_of_mRNA" , "microarray_of_microRNA", "methylation",

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
