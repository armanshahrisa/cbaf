#' @title Checking which Datasets are available for Each Cancer study.
#'
#' @description This function checks all the cancer studies that are registered in 'cbioportal.org' to
#' examine whether or not they contain RNA-seq, microRNA-seq, microarray(mRNA),
#' microarray(miRNA) and methylation data.
#'
#' @details
#' \tabular{lllll}{
#' Package: \tab cBioAutomatedTools \cr
#' Type: \tab Package \cr
#' Version: \tab 0.99.0 \cr
#' Date: \tab 2017-05-31 \cr
#' License: \tab Artistic-2.0 \cr
#' }
#'
#' @return A matrix that contain all cancers versus available data types. It is
#' stored in the global enviroment (user's workspace). For convenience, an excel
#' file is also generated in the working directory.
#'
#' @usage available.data.types()
#' @example available.data.types()
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

available.data.types <- function(){


  ############################################
  # Prerequisites

  # Prerequisites for cBioportal

  mycgds = cgdsr::CGDS("http://www.cbioportal.org/")

  all.cancers <- getCancerStudies(mycgds)[,2]


  # Creating Empty matrix in order to fill with availibility data

  available.data.matrix <- matrix(rep(" ", length(all.cancers)), nrow = length(all.cancers), ncol = 6)

  rownames(available.data.matrix) <- 1:length(all.cancers)

  available.data.matrix[,1] <- all.cancers

  colnames(available.data.matrix) <- c("Cancers Study" ,"RNA-seq", "MicroRNA-seq", "Microarray (mRNA)", "Microarray (miRNA)", "Methylation")




  ############################################
  # Availability of four data sets

  print("Cheching which data types are available for each cancer")


  # create progress bar

  pb <- txtProgressBar(min = 0, max = length(all.cancers), style = 3)


  # Getting information

  for(se in 1:length(all.cancers)){

    mycancerstudy = getCancerStudies(mycgds)[se,1]

    All.types <- getCaseLists(mycgds,mycancerstudy)[,2]



    # Availability of RNA-seq data

    available.data.matrix[se,2] <- any(All.types %in% c("Tumor Samples with mRNA data (RNA Seq V2)", "Tumors with mRNA data (RNA Seq V2)",

                                                        "Tumor Samples with mRNA data (RNA Seq)", "Tumors with mRNA data (RNA Seq)"))



    # Availability of microRNA-seq data

    available.data.matrix[se,3] <- any(All.types %in% c("Tumors with microRNA data (microRNA-Seq)"))



    # Availability of Microarray data (mRNA)

    available.data.matrix[se,4] <- any(All.types %in% c("Tumor Samples with mRNA data (Agilent microarray)",

                                                        "Tumors with mRNA data (Agilent microarray)", "Tumor Samples with mRNA data (U133 microarray only)",

                                                        "Tumors with mRNA data"))



    # Availability of RNA-seq data (miRNA)

    available.data.matrix[se,5] <- any(All.types %in% c("Tumors with microRNA"))



    # Availability of methylation data

    available.data.matrix[se,6] <- any(All.types %in% c("Tumor Samples with methylation data (HM450)", "Tumors with methylation data (HM450)",

                                                        "Tumor Samples with methylation data (HM27)", "Tumors with methylation data (HM27)",

                                                        "Tumors with methylation data"))



    # update progress bar

    setTxtProgressBar(pb, se)

  }

  close(pb)




  ################################################
  # Replacing True and False with available and ""

  available.data.matrix[available.data.matrix=="TRUE"] <- "Available"

  available.data.matrix[available.data.matrix=="FALSE"] <- "-"




  ################################################
  # Exporting results

  # Storing data matrix in the global enviornment

  available.data.types.output <<- available.data.matrix

  # Converting matrix to data.frame

  available.data.dataframe <- data.frame(available.data.matrix)

  colnames(available.data.dataframe) <- colnames(available.data.matrix)

  write.xlsx(available.data.dataframe, file="Available Data Types Output.xlsx")

  print(paste("An excel file entitled 'Available Data Types Output.xlsx' ,which shows the availabile data types for each cancer study, was saved in",

              getwd(), ".", "In adition, the output is available in the global enviroment as 'available.data.types.output'", sep=" "))

}
