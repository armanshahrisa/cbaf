#' @title Generate excel output for various studies/subgroups of a study.
#'
#' @description This function generates excel files containing gene validation
#' and all selected statistical methods. It uses outputs of
#' obtainOneStudy()/obtainMultipleStudies() and automatedStatistics() functions.
#'
#' @details
#' \tabular{lllll}{
#' Package: \tab cbaf \cr
#' Type: \tab Package \cr
#' Version: \tab 1.1.1 \cr
#' Date: \tab 2017-11-11 \cr
#' License: \tab Artistic-2.0 \cr
#' }
#'
#'
#'
#' @importFrom xlsx write.xlsx
#'
#' @importFrom BiocFileCache bfcnew bfcquery bfcpath
#'
#' @importFrom utils head setTxtProgressBar txtProgressBar
#'
#'
#'
#' @include cbaf-obtainOneStudy.R cbaf-obtainMultipleStudies.R
#' cbaf-automatedStatistics.R
#'
#'
#'
#' @usage xlsxOutput(submissionName)
#'
#'
#'
#' @param submissionName a character string containing name of interest. It is
#' used for naming the process.
#'
#'
#'
#' @return It generates one excel file for each gene group. This excel file
#' contains output of automatedStatistics() and validation result from output of
#' either obtainOneStudy() or obtainMultipleStudies().
#'
#'
#'
#' @examples
#' genes <- list(K.demethylases = c("KDM1A", "KDM1B", "KDM2A", "KDM2B", "KDM3A",
#'  "KDM3B", "JMJD1C", "KDM4A"), K.methyltransferases = c("SUV39H1", "SUV39H2",
#'  "EHMT1", "EHMT2", "SETDB1", "SETDB2", "KMT2A", "KMT2A"))
#'
#' obtainOneStudy(genes, "test", "Breast Invasive Carcinoma (TCGA, Cell 2015)",
#' "RNA-seq", desiredCaseList = c(3,4))
#'
#' automatedStatistics("test", obtainedDataType = "single study", calculate =
#' c("frequencyPercentage", "frequencyRatio"))
#'
#' xlsxOutput("test")
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
########## Generating excel file for the processed and validation data #########
################################################################################
################################################################################

xlsxOutput <- function(submissionName){

  ##############################################################################
  ########## Prerequisites

  # Check submissionName

  if(exists("submissionName")){

    if(!is.character(submissionName)){

      stop("'submissionName' must be entered as a character string for naming the process")

    }

  } else{

    stop("'submissionName' must be entered as a character string for naming the process")

  }





  ##############################################################################
  ########## Decide whether functions should stop now!

  # Check wheather the requested data exists

  database <-

    system.file("extdata", submissionName, package="cbaf")

  if(!dir.exists(database)){

    stop("Please run one of the obtainSingleStudy() or obtainMultipleStudies() functions and then the automatedStatistics() function first")

  } else if(dir.exists(database)){

    bfc <- BiocFileCache(

      file.path(system.file("extdata", package = "cbaf"), submissionName)

      )

    if(!nrow(bfcquery(bfc, c("Parameters for automatedStatistics()"))) == 1){

      stop("Please run the automatedStatistics() function first")

    }

  }



  # obtain parameters for prevous function

  previousFunctionParam <-

    readRDS(

      bfcpath(bfc, bfcquery(bfc, c("Parameters for automatedStatistics()"))$rid)

      )




  # fetch an old parameter from the previous function

  desiredTechnique <- previousFunctionParam$desiredTechnique

  cutoff <- previousFunctionParam$cutoff

  obtainedDataType <- previousFunctionParam$obtainedDataType



  # setting the value for cutoff

  if(desiredTechnique == "methylation"){

    cutoff.phrase <- "obs/exp cutoff"

  } else{

    cutoff.phrase <- "z-score cutoff"

  }





  # Store the new parameteres

  newParameters <-list()

  newParameters$submissionName <- submissionName





  # Check wheather the requested data exists

  number.of.rows.parameters <-

    nrow(bfcquery(bfc, "Parameters for xlsxOutput()"))


  if(number.of.rows.parameters == 1){

    oldParameters <-

      readRDS(bfcpath(bfc, bfcquery(bfc, c("Parameters for xlsxOutput()"))$rid))

    # Check whether the previous function is skipped

    if(previousFunctionParam$lastRunStatus == "skipped"){

      if(identical(oldParameters, newParameters)){

        continue <- FALSE

      } else{

        continue <- TRUE

      }

    } else{

      continue <- TRUE

    }

  } else{

    continue <- TRUE

  }







  # Getting the source data

  statisticsData <-

    readRDS(bfcpath(bfc, bfcquery(bfc, c("Calculated statistics"))$rid))

  if(!is.list(statisticsData)){

    stop("Input database must be a list.")

  }





  # Check gene validation dat

  if(obtainedDataType == "multiple studies"){

    validationName <- "Validation data for multiple studie"

  } else if(obtainedDataType == "single study"){

    validationName <- "Validation data for single study"

  }

  if(nrow(bfcquery(bfc, validationName)) == 1){

    validation.data <- readRDS(bfcpath(bfc, bfcquery(bfc, validationName)$rid))

  }





  # get the working directory

  parent.directory <- getwd()



  ##############################################################################
  ########## Set the function ready to work

  # Report

  message("***", " Preparing the requested excel file(s) for ", submissionName, " ***")



  # Count number of skipped excel files

  skipped <- 0



  # Create progressbar

  total.number <- length(statisticsData)

  heatmapOutputProgressBar <-

    txtProgressBar(min = 0, max = total.number, style = 3)

  ExtH <- 0





  ##############################################################################
  ########## Core segment

  # Save heatmaps in separate folder

  for(gr in seq_along(statisticsData)){

    # Subset data that can be presented as heatmap

    subset.name <- names(statisticsData)[gr]

    subset.data <- statisticsData[[gr]]

    sheet.names <- names(subset.data)



    # merge with validation data

    if(exists("validation.data")){

      v.subset.name <- names(validation.data)[gr]

      v.subset.data <- validation.data[[gr]]



      if(identical(subset.name, v.subset.name)){

        subset.data$Validation.Result <- v.subset.data

        reorder <-

          match(c("Validation.Result", sheet.names), names(subset.data))

        subset.data <- subset.data[reorder]

      }

    }







    # Create a directory and set it as desired folder

    child.directory <-

      paste0(gr, ". ", sub(x = subset.name, pattern = "\\.", replacement = "-"))

    dir.create(

      paste(parent.directory, child.directory, sep = "/"),

      showWarnings = FALSE

      )


    setwd(paste(parent.directory, child.directory, sep = "/"))



    # determine ourput file name

    name.of.excel.file <-

      paste0(

        gsub(x = subset.name, pattern = "\\.", replacement = "-"),

        " (",

        cutoff.phrase,

        "=",

        cutoff,

        ")" ,

        ".xlsx"

        )

    # Check if excel file already exists

    if(continue & file.exists(name.of.excel.file)){

      file.remove(name.of.excel.file)

      hault <- FALSE

    }else if(continue| !continue & !file.exists(name.of.excel.file)){

      hault <- FALSE

    }else{

      skipped <- skipped + 1

      hault <- TRUE

    }




    for(possible in seq_along(subset.data)){

      if(!hault){

        # subset statistics

        statistics.data <- subset.data[[possible]]

        name.statistics.data <- names(subset.data)[possible]



        if(name.statistics.data %in% c("Frequency.Percentage",

                                       "Mean.Value",

                                       "Median.Value")){

          # Replace 'NaN' values with character string

          statistics.data[is.nan(statistics.data)] <- "NaN"

          # Replace 'NA' values with character string

          statistics.data[is.na(statistics.data)] <- "NA"

        }



        # Saving the expression profile

        write.xlsx(

          statistics.data,

          file = name.of.excel.file,

          sheetName = gsub(

            x = name.statistics.data,

            pattern = "\\.",

            replacement = " "

            ),

          append = TRUE)

      }

      # Update progressbar

      ExtH <- ExtH + 1

      setTxtProgressBar(heatmapOutputProgressBar, ExtH)

    }

  }

  # Close progressbar

  close(heatmapOutputProgressBar)

  # report number of skipped heatmaps

  if(skipped > 0 & skipped != 1){

    message("--- ", as.character(skipped), " out of ", as.character(total.number), " excel files were skipped: They already exist. ---")

  } else if(skipped > 0 & skipped == 1){

    message("--- ", as.character(skipped), " out of ", as.character(total.number), " excel file was skipped: It already exists. ---")

  }



  # Store the last parameter

  oldParamXlsxOutput <- newParameters


  # Store the parameters for this run

  if(number.of.rows.parameters == 0){

    saveRDS(

      oldParamXlsxOutput,

      file=bfcnew(bfc, "Parameters for xlsxOutput()", ext="RDS")

      )

  } else if(number.of.rows.parameters == 1){

    saveRDS(

      oldParamXlsxOutput,

      file=bfc[[bfcquery(bfc, "Parameters for xlsxOutput()")$rid]]

      )

  }



  # change directory to parent directory

  setwd(parent.directory)

}
