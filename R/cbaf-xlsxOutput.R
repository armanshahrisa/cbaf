#' @title Generating excel output for various studies/subgroups of a study.
#'
#' @description This function generates excel files containing gene validation and all selected statistical methods. It
#' uses outputs of obtainOneStudy()/obtainMultipleStudies() and automatedStatistics() functions.
#'
#' @details
#' \tabular{lllll}{
#' Package: \tab cbaf \cr
#' Type: \tab Package \cr
#' Version: \tab 0.99.0 \cr
#' Date: \tab 2017-06-22 \cr
#' License: \tab Artistic-2.0 \cr
#' }
#'
#' @import xlsx Biobase
#'
#' @include cbaf-obtainOneStudy.R cbaf-obtainMultipleStudies.R cbaf-automatedStatistics.R
#'
#' @usage xlsxOutput(submissionName)
#'
#' @param submissionName a character string containing name of interest. It is used for naming the process.
#'
#' @return It generates one excel file for each genegroup.
#'
#' @examples
#' # xlsxOutput("test")
#'
#'
#' @author Arman Shahrisa, \email{shahrisa.arman@hotmail.com} [maintainer, copyright holder]
#' @author Maryam Tahmasebi Birgani, \email{tahmasebi-ma@ajums.ac.ir}
#'
#' @export



#########################################################################
#########################################################################
############# Obtain the requested data for multiple Cancer #############
#########################################################################
#########################################################################

xlsxOutput <- function(submissionName){

  ##########################################################################
  ########## Prerequisites

  # Check submissionName

  if(exists("submissionName")){

    if(!is.character(submissionName)){

      stop("'submissionName' must be entered as a character string for naming the process")

    }

  } else{

    stop("'submissionName' must be entered as a character string for naming the process")

  }





  ##########################################################################
  ########## Decide whether functions should stops now!

  # Store the new parameteres

  newParameters <-list()

  newParameters$submissionName <- submissionName





  # Check wheather the requested data exists

  if(!exists(paste("Pa.PrData.", submissionName, sep = ""))){

    Stop("Please run automatedStatistics() function first")

  } else{

    oldParam <- get(paste("Pa.PrData.", submissionName, sep = ""))

    desiredTechnique <- oldParam$desiredTechnique

    cutoff <- oldParam$cutoff

    obtainedDataType <- oldParam$obtainedDataType

  }



  # setting the value for cutoff

  if(desiredTechnique == "methylation"){

    cutoff.phrase <- "obs/exp cutoff"

  } else{

    cutoff.phrase <- "z-score cutoff"

  }



  if(exists(paste("Pa.Excel.", submissionName, sep = ""))){

    if(oldParam$HaultOrder == TRUE){

      if(identical(get(paste("Pa.Excel.", submissionName, sep = "")), newParameters)){

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





  # Get the previous function's result

  whole.data <- get(paste("PrData.", submissionName, sep = ""))

  # Check gene validation dat

  if(obtainedDataType == "multiple studies"){

    validationType <- "Va.Mu."

  } else if(obtainedDataType == "single study"){

    validationType <- "Va.Si."

  }

  if(exists(paste(validationType, submissionName, sep = ""))){

    validation.data <- get(paste(validationType, submissionName, sep = ""))

  }





  # get the working directory

  parent.directory <- getwd()



  ##########################################################################
  ########## Set the function ready to work

  # Report

  print(paste("***", "Preparing the requested excel file(s) for", submissionName, "***", sep = " "))



  # Count number of skipped excel files

  skipped <- 0





  # Create progressbar

  total.number <- length(whole.data)

  heatmapOutputProgressBar <- txtProgressBar(min = 0, max = total.number, style = 3)

  ExtH <- 0





  ##########################################################################
  ########## Core segment

  # Save heatmaps in separate folder

  for(gr in 1:length(whole.data)){

    # Subset data that can be presented as heatmap

    subset.name <- names(whole.data)[gr]

    subset.data <- whole.data[[gr]]

    sheet.names <- names(subset.data)



    # merge with validation data

    if(exists("validation.data")){

      v.subset.name <- names(validation.data)[gr]

      v.subset.data <- validation.data[[gr]]



      if(identical(subset.name, v.subset.name)){

        subset.data$Validation.Result <- v.subset.data

        subset.data <- subset.data[match(c("Validation.Result", sheet.names), names(subset.data))]

      }

    }







    # Create a directory and set it as desired folder

    child.directory <- paste(gr, ". ", sub(x = subset.name, pattern = "\\.", replacement = "-"), sep="")

    dir.create(paste(parent.directory, child.directory, sep = "/"), showWarnings = FALSE)

    setwd(paste(parent.directory, child.directory, sep = "/"))



    # determine ourput file name

    name.of.excel.file <- paste(gsub(x = subset.name, pattern = "\\.", replacement = "-"), " (",cutoff.phrase, "=", cutoff, ")" , ".xlsx", sep = "")

    # Check if excel file already exists

    if(continue == TRUE & file.exists(name.of.excel.file)){

      file.remove(name.of.excel.file)

      hault <- FALSE

    }else if(continue == TRUE | continue == FALSE & !file.exists(name.of.excel.file)){

      hault <- FALSE

    }else{

      skipped <- skipped + 1

      hault <- TRUE

    }




    for(possible in 1:length(subset.data)){

      if(hault == FALSE){

        # subset statistics

        statistics.data <- subset.data[[possible]]

        name.statistics.data <- names(subset.data)[possible]



        if(name.statistics.data %in% c("Frequency.Percentage", "Mean.Value", "Median.Value")){

          # Replace 'NaN' values with character string

          statistics.data[is.nan(statistics.data)] <- "NaN"

          # Replace 'NA' values with character string

          statistics.data[is.na(statistics.data)] <- "NA"

        }



        if(name.statistics.data %in% c("Top.Genes.of.Frequency.Percentage", "Top.Genes.of.Mean.Value", "Top.Genes.of.Median.Value")){

          # Removing index as rownames

          statistics.data <- statistics.data[, 2:length(ncol(statistics.data))]

        }



        # Saving the expression profile

        write.xlsx(statistics.data, file = name.of.excel.file, sheetName = gsub(x = name.statistics.data, pattern = "\\.", replacement = " "), append = TRUE)

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

    print(paste("--- ", as.character(skipped), " out of ", as.character(total.number)," excel files were skipped: They already exist. ---", sep = ""))

  } else if(skipped > 0 & skipped == 1){

    print(paste("--- ", as.character(skipped), " out of ", as.character(total.number)," excel file was skipped: It already exist. ---", sep = ""))

  }

  # Store the parameters for this run

  assign(paste("Pa.Excel.", submissionName, sep = ""), newParameters, envir = globalenv())

  # change directory to parent directory

  setwd(parent.directory)

}
