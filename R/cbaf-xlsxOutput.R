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
#' Version: \tab 1.19.5 \cr
#' Date: \tab 2022-07-14 \cr
#' License: \tab Artistic-2.0 \cr
#' }
#'
#'
#'
#' @importFrom openxlsx createWorkbook addWorksheet writeData saveWorkbook
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
#' @usage xlsxOutput(submissionName, transposeResults = FALSE)
#'
#'
#'
#' @param submissionName a character string containing name of interest. It is
#' used for naming the process.
#'
#' @param transposeResults a logical value that enables the function to replace
#' the columns and rows of data.
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
#' "RNA-Seq", desiredCaseList = c(3,4))
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

xlsxOutput <- function(submissionName, transposeResults = FALSE){

  ##############################################################################
  ########## Prerequisites

  # Check submissionName

  if(!missing(submissionName)){

    if(!is.character(submissionName)){

      stop("[xlsxOutput] 'submissionName' must be a character string!")

    }

  } else{

    stop("[xlsxOutput] 'submissionName' must be a character string!")

  }

  # Check transposeResults

  if(!is.logical(transposeResults)){

    stop("[xlsxOutput] 'transposeResults' must be either TRUE or FALSE!")

  }





  ##############################################################################
  ########## Decide whether functions should stop now!

  # Check wheather the requested data exists

  database <-

    system.file("extdata", submissionName, package="cbaf")

  if(!dir.exists(database)){

    stop("[xlsxOutput] Please run one of the obtainSingleStudy() or obtainMultipleStudies() functions first, and then the automatedStatistics() function!")

  } else if(dir.exists(database)){

    bfc <- BiocFileCache(

      file.path(system.file("extdata", package = "cbaf"), submissionName),

      ask = FALSE

      )

    if(!nrow(bfcquery(bfc, c("Parameters for automatedStatistics()"))) == 1){

      stop("[xlsxOutput] Please run the automatedStatistics() function first!")

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

    cutoff.phrase <- "Mean methylation cutoff"

  } else{

    cutoff.phrase <- "log z-score cutoff"

  }





  # Store the new parameteres

  newParameters <-list()

  newParameters$submissionName <- submissionName

  newParameters$transposeResults <- transposeResults





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

    stop("[xlsxOutput] Input database must be a list!")

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

  message("[xlsxOutput] Preparing excel file(s).")



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

    if(is.null(subset.name)){

      subset.name <- " "

    }

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


    # Initiating excel file creation

    xo <- createWorkbook()

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

    if(continue | !continue & !file.exists(name.of.excel.file)){

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

        ## Replacing dots with space in sheet name

        currentSheetName = gsub(

          x = name.statistics.data,

          pattern = "\\.",

          replacement = " "

        )

        ## Fixing the Max length of 31 characters for openxlsx package

        currentSheetName = gsub(

          x = currentSheetName,

          pattern = "Top Genes of Frequency Percentage",

          replacement = "Top Genes of Frequency Percenta"

        )

        ## Save xlsx file

        if(!transposeResults){

          addWorksheet(xo, sheetName = currentSheetName)

          writeData(xo, sheet = currentSheetName, x = statistics.data,

                    rowNames = TRUE)


        } else if(transposeResults){

          addWorksheet(xo, sheetName = currentSheetName)

          writeData(xo, sheet = currentSheetName, x = t(statistics.data),

                    rowNames = TRUE)

        }


      }

      # Update progressbar

      ExtH <- ExtH + 1

      setTxtProgressBar(heatmapOutputProgressBar, ExtH)

    }

    # Avoid removing previos xlsx file

    if(file.exists(name.of.excel.file) & length(names(xo)) != 0){

      file.remove(name.of.excel.file)

    }

    # Avoid storing empty workbook

    if(length(names(xo)) != 0){

      saveWorkbook(xo, name.of.excel.file)

    }


  }

  # Close progressbar

  close(heatmapOutputProgressBar)

  # report number of skipped heatmaps

  if(skipped > 0 & skipped != 1){

    message("[xlsxOutput] ", as.character(skipped), " out of ", as.character(total.number), " excel files were skipped, because they already exist!")

  } else if(skipped > 0 & skipped == 1){

    message("[xlsxOutput] ", as.character(skipped), " out of ", as.character(total.number), " excel file was skipped, because it already exists!")

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

  # message("[xlsxOutput] Finished.")

}
