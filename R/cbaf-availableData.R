#' @title Checking which Data types are available for Each Cancer study.
#'
#' @description This function checks all the cancer studies that are registered in 'cbioportal.org' to
#' examine whether or not they contain RNA-seq, microRNA-seq, microarray(mRNA), microarray(miRNA)
#' and methylation data.
#'
#' @details
#' \tabular{lllll}{
#' Package: \tab cBioAutomatedTools \cr
#' Type: \tab Package \cr
#' Version: \tab 0.99.0 \cr
#' Date: \tab 2017-06-22 \cr
#' License: \tab Artistic-2.0 \cr
#' }
#'
#' @import cgdsr xlsxjars xlsx Biobase
#'
#' @param outputName a character string that is required to name the output and, if requested, excel file.
#'
#' @param excelFile a logical value that tells the function whether or not export the results as an excel file. Default vlue is TRUE
#'
#' @return A matrix that contain all cancers versus available data types. It is stored in the
#' global enviroment (user's workspace) as 'available.data.types.output'. For convenience, an
#' excel file - 'Available Data Types Output.xlsx' - is also generated in the working directory.
#'
#' @usage availableData()
#'
#' @author Arman Shahrisa, \email{shahrisa.arman@hotmail.com} [maintainer, copyright holder]
#' @author Maryam Tahmasebi Birgani, \email{tahmasebi-ma@ajums.ac.ir}
#'
#' @export



###################################################################################################
###################################################################################################
########################## Dataset availability in all cBioportal Cancers #########################
###################################################################################################
###################################################################################################

availableData <- function(outputName, excelFile = TRUE){


  ############################################
  # Check input parameters

  if(exists("outputName")){

    if(!is.character(outputName)){

      stop("'outputName' must be entered as a character string")

    }

  } else{

    stop("'outputName' must be entered as a character string")

  }



  if(!is.logical(excelFile)){

    stop("'excelFile' must be either TRUE or FALSE")

  }



  ############################################
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



  print("Cheching the available data for every cancer study")


  # create progress bar

  pb <- txtProgressBar(min = 0, max = nrow(list.of.studies), style = 3)

  i <- 0



  ############################################
  ## core segment

  # looking for supported technique data

  list.of.available.data <- sapply(list.of.studies[, "cancer_study_id"], function(cs, cgds) {

    description <- getCaseLists(cgds, cs)[, "case_list_name"]

    i <<- i + 1

    setTxtProgressBar(pb, i)

    c(RNA.seq = any(RNA_seq.terms %in% description),

      microRNA.seq = any(microRNA_seq.terms %in% description),

      microarray_of_mRNA = any(microarray.for.mRNA.term %in% description),

      microarray_of_miRNA = any(microarray.for.miRNA.term %in% description),

      methylation = any(methylation.term %in% description))



  }, mycgds)


  # close progressbar

  close(pb)



  # Replacing True and False with available and ""

  list.of.available.data[list.of.available.data=="TRUE"] <- "available"

  list.of.available.data[list.of.available.data=="FALSE"] <- "-"



  # joining list.of.available.data to list.of.studies

  combined.list <- cbind(list.of.studies[,"cancer_study_id"], list.of.studies[,"name"], t(list.of.available.data), list.of.studies[,"description"])

  colnames(combined.list) <- c("cancer_study_id", "cancer_study_name", "RNA.seq", "microRNA.seq", "microarray_of_mRNA",

                               "microarray_of_microRNA", "methylation", "description")

  rownames(combined.list) <- 1:nrow(combined.list)



  ################################################
  # Exporting results

  # Storing the output as a variable

  assign(outputName, combined.list, envir = globalenv())



  # Converting matrix to data.frame

  if(excelFile == TRUE){

    combined.list.dataframe <- data.frame(combined.list)

    colnames(combined.list.dataframe) <- gsub("_", " ", colnames(combined.list.dataframe))

    colnames(combined.list.dataframe) <- gsub("\\.", "-", colnames(combined.list.dataframe))

    rownames(combined.list.dataframe) <- 1:nrow(combined.list)

    write.xlsx(combined.list.dataframe, file=paste(outputName, ".xlsx", sep = ""))

  }

}
