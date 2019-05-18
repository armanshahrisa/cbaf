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
#' Version: \tab 1.7.2 \cr
#' Date: \tab 2019-05-18 \cr
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

    print(removable.directories)

    writeLines("")

    message("Please enter the folder names you wish to remove, seperated by comma. For instance, 'test' and 'test2' must be enterd as: test, test2")

    input <- readline(prompt = "Enter the name of folders: ")

    input <- as.character(unlist(strsplit(input, ", ")))

  } else if(is.character(databaseNames) & length(removable.directories) != 0){

    input <- databaseNames

  }else if(length(removable.directories) == 0){

    message("No removable directory was found.")

  }else{

    stop("Incorrect database name!")

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

      message("1 database was removed")

    }else if((number.of.correct.directories  > 1)){

      message(sum(!dir.exists(full.path.of.correct.directories)), " databases were removed")

    }

  }else if(length(removable.directories) != 0){

    message("0 database was removed")

  }

}
