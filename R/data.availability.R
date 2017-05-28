#' @title Checking which Datasets are available for Each Cancer
#'
#' @description This function checks all the cancers that are registered in 'cbioportal.org' to
#' examine whether or not they contain RNA-seq, microRNA-seq, microarray(mRNA),
#' microarray(miRNA) and methylation datasets.
#'
#' @details
#' \tabular{lllll}{
#' Package: \tab cBioAutomatedTools \cr
#' Type: \tab Package \cr
#' Version: \tab 0.99.0 \cr
#' Date: \tab 2017-05-30 \cr
#' License: \tab Artistic-2.0 \cr
#' }
#'
#' @return A matrix that contain all cancers and their available datasets. It is
#' available in the global enviroment (user's workspace). For convenience, an excel
#' file will also be generated in the working directory.
#'
#' @details This function checks all the cancers that are registered in 'cbioportal.org' to
#' examine whether or not they contain RNA-seq, microRNA-seq, microarray(mRNA),
#' microarray(miRNA) and methylation datasets.
#'
#' @usage data.availability()
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

data.availability <- function(){


  ############################################
  # Prerequisites

  # Prerequisites for cBioportal

  mycgds = CGDS("http://www.cbioportal.org/")

  all.cancers <- getCancerStudies(mycgds)[,2]


  # Creating Empty matrix in order to fill with availibility data

  availability.matrix <- matrix(rep(" ", length(all.cancers)), nrow = length(all.cancers), ncol = 6)

  rownames(availability.matrix) <- 1:length(all.cancers)

  availability.matrix[,1] <- all.cancers

  colnames(availability.matrix) <- c("Cancers Study" ,"RNA-seq", "MicroRNA-seq", "Microarray (mRNA)", "Microarray (miRNA)", "Methylation")




  ############################################
  # Availability of four data sets

  print("Cheching the availability of requested data sets")


  # create progress bar

  pb <- txtProgressBar(min = 0, max = length(all.cancers), style = 3)


  # Getting information

  for(se in 1:length(all.cancers)){

    mycancerstudy = getCancerStudies(mycgds)[se,1]

    All.types <- getCaseLists(mycgds,mycancerstudy)[,2]



    # Availability of RNA-seq data

    availability.matrix[se,2] <- any(All.types %in% c("Tumor Samples with mRNA data (RNA Seq V2)", "Tumors with mRNA data (RNA Seq V2)", "Tumor Samples with mRNA data (RNA Seq)", "Tumors with mRNA data (RNA Seq)"))



    # Availability of microRNA-seq data

    availability.matrix[se,3] <- any(All.types %in% c("Tumors with microRNA data (microRNA-Seq)"))



    # Availability of Microarray data (mRNA)

    availability.matrix[se,4] <- any(All.types %in% c("Tumor Samples with mRNA data (Agilent microarray)", "Tumors with mRNA data (Agilent microarray)", "Tumor Samples with mRNA data (U133 microarray only)", "Tumors with mRNA data"))



    # Availability of RNA-seq data (miRNA)

    availability.matrix[se,5] <- any(All.types %in% c("Tumors with microRNA"))



    # Availability of methylation data

    availability.matrix[se,6] <- any(All.types %in% c("Tumor Samples with methylation data (HM450)", "Tumors with methylation data (HM450)", "Tumor Samples with methylation data (HM27)", "Tumors with methylation data (HM27)", "Tumors with methylation data"))



    # update progress bar

    setTxtProgressBar(pb, se)

  }

  close(pb)




  ################################################
  # Replacing True and False with available and ""

  availability.matrix[availability.matrix=="TRUE"] <- "Available"

  availability.matrix[availability.matrix=="FALSE"] <- "-"




  ################################################
  # Exporting results

  # Storing data matrix in the global enviornment

  availability.matrix <<- availability.matrix

  # Converting matrix to data.frame

  availability.dataframe <- data.frame(availability.matrix)

  colnames(availability.dataframe) <- colnames(availability.matrix)

  write.xlsx(availability.dataframe, file="Available Datasets.xlsx")

  print(paste("An .xlsx file entitled 'AvailableDataSets.xlsx' that shows the availability of specific datasets in all cancers was saved in", getwd(), sep=" "))

}
