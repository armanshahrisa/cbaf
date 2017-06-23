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
#' @return A matrix that contain all cancers versus available data types. It is stored in the
#' global enviroment (user's workspace) as 'available.data.types.output'. For convenience, an
#' excel file - 'Available Data Types Output.xlsx' - is also generated in the working directory.
#'
#' @usage available.data.types()
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










# The function has three parts -- find all studies, get case lists, output
# to Excel. I revised it into three separate functions, which might be
# documented on the same page. The first gets all studies

available.studies <- function() {
  mycgds <- cgdsr::CGDS("http://www.cbioportal.org/")
  getCancerStudies(mycgds)
}

# I revised available.data.types as

available.data.types <- function(studies = available.studies()) {
  ## relevant descriptions
  RNASeq <- c(
    "Tumor Samples with mRNA data (RNA Seq V2)",
    "Tumors with mRNA data (RNA Seq V2)",
    "Tumor Samples with mRNA data (RNA Seq)",
    "Tumors with mRNA data (RNA Seq)"
  )
  microRNA <- "Tumors with microRNA data (microRNA-Seq)"
  mRNA <- c(
    "Tumor Samples with mRNA data (Agilent microarray)",
    "Tumors with mRNA data (Agilent microarray)",
    "Tumor Samples with mRNA data (U133 microarray only)",
    "Tumors with mRNA data"
  )
  miRNA <- "Tumors with microRNA"
  methylation <- c(
    "Tumor Samples with methylation data (HM450)",
    "Tumors with methylation data (HM450)",
    "Tumor Samples with methylation data (HM27)",
    "Tumors with methylation data (HM27)",
    "Tumors with methylation data"
  )

  ## query each study
  mycgds <- CGDS("http://www.cbioportal.org/")
  pb <- txtProgressBar(min = 0, max = nrow(studies), style = 3)
  i <- 0
  available <- sapply(studies[, "name"], function(se, cgds) {
    description <- getCaseLists(cgds, se)[, "case_list_description"]
    i <<- i + 1
    setTxtProgressBar(pb, i)
    c(RNASeq = any(description %in% RNASeq),
      microRNA = any(description %in% microRNA),
      microarray_mRNA = any(description %in% mRNA),
      microarray_miRNA = any(description %in% miRNA),
      methylation = any(description %in% methylation))
  }, mycgds)
  close(pb)

  cbind(studies, t(available))
}

# The main changes are

# 1. provide an argument with default value, preserving original behavior

# 2. collect constant strings to an 'initialization' section of the function

# 3. use sapply() to manage memory, rather than allocating and filling a
# matrix

# 4. respect R's notion of data types, e.g., logical() values rather than
# "TRUE" / "FALSE" character strings; using a data.frame (because columns
# are of different type) rather than matrix.

#5. avoid querying the web service for getCancerStudies() on each iteration

#6. not printing the result to an Excel file

#7. returning the result to the user for subsequent processing, without
#writing to the global environment.

#I created a third function to write to Excel. It expects a file name
#from the user and will not overwrite an existing file. It returns the
#file name where the output can be found

write.xls <- function(data.types, file) {
stopifnot(!file.exists(file))
## your implementation here
file
}


# On the help page I might illustrate this as

# studies <- available.studies()
# head(studies)

# then use a subset of the data for the next function

# data.types <- available.data.types(head(studies))
# data.types

# and finally write the output to a temporary file

xls.file <- write.xls(data.types, tempfile(fileext=".xls"))
xls.file    # open in Excel

# You might wish to write a final function that integrates these steps

all.available.data.types <- function(file) {
stopifnot(!file.exists(file))
studies = available.studies()
data.types = available.data.types(studies)
write.xls(data.types, file)
}

# and include this as Michael suggests

# \dontrun{
# all.available.data.types("my.xls")
# }

#Bioconductor style would strongly discourage '.'-separated function
# names, instead preferring camelCase() or snake_case().'






















genesList <- list(b=c("HOTAIR", "GAS5", "SNHG6"), hhc = c("PVT1", "GAS5", "DBBBGTZZXCC", "NSD2"))

genesList <- list(b=c("HOTAIR", "GAS5"))

cancers <- c("Acute Myeloid Leukemia (TCGA, NEJM 2013)", "Acute Myeloid Leukemia (TCGA, Provisional)")



# obtainMultipleStudies(genesList, "test", cancers, desiredTechnique = "RNA-seq")








setwd("Z:/2017.04.25-Breast-LysineX3in1-Dr. Tahmasebi")

cancername <- "Breast Invasive Carcinoma (TCGA, Cell 2015)"

library(xlsx)
data.x3 <- read.xlsx("all the gene.xlsx", sheetIndex = 1)

lysine.X3 <- list(Lysine.acetyltransferases = data.x3[1:17,2], Lysine.demethylases = data.x3[18:41,2], Lysine.methyltransferases = data.x3[42:75,2])
lysine1 <- list(Lysinex3 = sort(data.x3[1:75,2])[1:7], r=c("SETDB1", "SUV39H2", "PVT1"))




# obtainOneStudy(lysine1, "test2", "Breast Invasive Carcinoma (TCGA, Cell 2015)", "RNA-seq", desiredCaseList = c(2,3,4,5), validateGenes = TRUE)


# automatedStatistics(obtainedDataType = "single study", "test2", "RNA-seq", calculate = c("frequencyPercentage", "frequencyRatio", "meanValue", "medianValue"), cutoff=NULL, round=TRUE, topGenes = TRUE)








# processMultipleStudies(genesList, "test2", cancers, "RNA-seq", cutoff = 2)

# processOneStudy(lysine1, "test2", cancername, "RNA-seq", c(2,3,4,5))





