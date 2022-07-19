#' @title Clean the created database(s)
#'
#' @description This function removes the created databases in the cbaf package
#' directory. This helps users to obtain the fresh data from cbioportal.org.
#'
#'
#' @details
#' \tabular{lllll}{
#' Package: \tab cbaf \cr
#' Type: \tab Package \cr
#' Version: \tab 1.19.5 \cr
#' Date: \tab 2022-07-19 \cr
#' License: \tab Artistic-2.0 \cr
#' }
#'
#'
#'
#' @usage cleanDatabase(databaseNames = NULL)
#'
#'
#'
#' @param databaseNames a character vector that contains name of databases that
#' will be removed. The default value in \code{null}.
#'
#'
#'
#' @return prints the number of removed databases.
#'
#'
#'
#' @examples
#' cleanDatabase(databaseNames = "Whole")
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
############## Obtain the requested data for Subtypes of a Cancer ##############
################################################################################
################################################################################

cleanDatabase <- function(databaseNames = NULL){

  ##############################################################################
  ########## Prerequisites

  # List all folders in extdata folder of cbaf package

  message("[cleanDatabase] Checking for the local data sets.")

  package.path <- system.file("extdata", package = "cbaf")

  just.names <- sapply(

    strsplit(

      as.character(list.dirs(path = package.path)),

      split="extdata/",

      fixed=TRUE

      ) ,

    function(x) (x[2])
    )

  removable.directories <- just.names[!(just.names %in% c(NA, "test", "test2"))]



  ##############################################################################
  ########## Core part

  # Print and ask the folder names

  if(is.null(databaseNames) & length(removable.directories) != 0){

    message("[cleanDatabase] List of removable data set(s):")

    print(removable.directories)

    writeLines("")

    message("[cleanDatabase] Please enter the name of data set(s) you want to remove. Example: test, test2")

    input <- readline(prompt = "Your choice(s): ")

    input <- as.character(unlist(strsplit(input, ", ")))

  } else if(is.character(databaseNames) & length(removable.directories) != 0){

    input <- databaseNames

  }else if(length(removable.directories) == 0){

    message("[cleanDatabase] No removable data set was found.")

  }else{

    stop("[cleanDatabase] Incorrect name of data set(s)!")

  }



  # Check which entered names are correct

  correct.directories <- input[input %in% removable.directories]

  number.of.correct.directories <- length(correct.directories)

  full.path.of.correct.directories <-

    paste0(package.path, "/", correct.directories)



  # Remove directories

  if(number.of.correct.directories  != 0){

    unlink(full.path.of.correct.directories, recursive = TRUE)

  }



  # Inform the user

  if(number.of.correct.directories  != 0){

    if(number.of.correct.directories  == 1){

      message("[cleanDatabase] 1 data set was removed.")

    }else if((number.of.correct.directories  > 1)){

      message("[cleanDatabase] " ,sum(!dir.exists(full.path.of.correct.directories)), " data set(s) were removed.")

    }

  }else if(length(removable.directories) != 0){

    message("[cleanDatabase] 0 data set was removed.")

  }

  # message("[cleanDatabase] Finished.")

}
