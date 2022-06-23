#' @title Generate heatmaps for various studies/subgroups of a study.
#'
#' @description This function can prepare heatmap for 'frequency percentage',
#' 'mean value' and 'median value' data provided by
#' automatedStatistics() function.
#'
#' @details
#' \tabular{lllll}{
#' Package: \tab cbaf \cr
#' Type: \tab Package \cr
#' Version: \tab 1.19.2 \cr
#' Date: \tab 2022-06-23 \cr
#' License: \tab Artistic-2.0 \cr
#' }
#'
#'
#'
#' @importFrom genefilter rowVars
#'
#' @importFrom RColorBrewer brewer.pal
#'
#' @importFrom gplots heatmap.2 redgreen
#'
#' @importFrom BiocFileCache bfcnew bfcquery bfcpath
#'
#' @importFrom grDevices colorRampPalette dev.off tiff png bmp jpeg pdf
#'
#' @importFrom utils head setTxtProgressBar txtProgressBar
#'
# @importFrom magicK image_read image_crop image_write
#'
#'
#'
#' @include cbaf-obtainOneStudy.R cbaf-obtainMultipleStudies.R
#' cbaf-automatedStatistics.R
#'
#'
#'
#' @usage heatmapOutput(submissionName, shortenStudyNames = TRUE,
#'   geneLimit = 50, rankingMethod = "variation", heatmapFileFormat = "TIFF",
#'   resolution = 600, RowCex = "auto", ColCex = "auto",
#'   heatmapMargines = "auto", rowLabelsAngle = 0, columnLabelsAngle = 45,
#'   heatmapColor = "RdBu", reverseColor = TRUE, transposedHeatmap = FALSE,
#'   simplifyBy = FALSE, genesToDrop = FALSE)
#'
#'
#'
#' @param submissionName a character string containing name of interest. It is
#' used for naming the process.
#'
#' @param shortenStudyNames a logical vector. If the value is set as TRUE,
#' function will try to remove the last part of the cancer names aiming to
#' shorten them. The removed segment usually contains the name of scientific
#' group that has conducted the experiment.
#'
#' @param geneLimit if large number of genes exist in at least one gene group,
#' this option can be used to limit the number of genes that are shown on
#' heatmap. For instance, \code{geneLimit=50} will limit the heatmap to 50 genes
#' that show the most variation across multiple study / study subgroups. The
#' default value is \code{50}.
#'
#' @param rankingMethod a character value that determines how genes will be
#' ranked prior to drawing heatmap. \code{"variation"} orders the genes based on
#' unique values in one or few cancer studies while \code{"highValue"} ranks the
#'  genes when they contain high values in multiple / many cancer studies. This
#'  option is useful when number of genes are too much so that user has to limit
#'  the number of genes on heatmap by \code{geneLimit}.
#'
#' @param heatmapFileFormat This option enables the user to select the desired
#' image file format of the heatmaps. The default value is \code{"TIFF"}. Other
#' supported formats include \code{"JPG"}, \code{"BMP"}, \code{"PNG"}, and
#' \code{"PDF"}.
#'
#' @param resolution a number. This option can be used to adjust the resolution
#' of the output heatmaps as 'dot per inch'. The defalut value is 600.
#'
#' @param RowCex a number that specifies letter size in heatmap row names,
#' which ranges from 0 to 2. If \code{RowCex = "auto"}, the function will
#' automatically determine the best RowCex.
#'
#' @param ColCex a number that specifies letter size in heatmap column names,
#' which ranges from 0 to 2. If \code{ColCex = "auto"}, the function will
#' automatically determine the best ColCex.
#'
#' @param heatmapMargines a numeric vector that is used to set heatmap margins.
#'  If \code{heatmapMargines = "auto"}, the function will automatically
#'  determine the best possible margines. Otherwise, enter the desired margine as
#'  e.g. c(10,10.)
#'
#' @param rowLabelsAngle a number that determines the angle with which the
#' gene names are shown in heatmaps. The default value is 0 degree.
#'
#' @param columnLabelsAngle a number that determines the angle with which the
#' studies/study subgroups names are shown in heatmaps. The default value is 45
#' degree.
#'
#' @param heatmapColor a character string that defines heatmap color. The
#' default value is \code{'RdBu'}. \code{'RdGr'} is also a popular color in
#' genomic studies. To see the rest of colors, please type
#' \code{library(RColorBrewer)} and then \code{display.brewer.all()}.
#'
#' @param reverseColor a logical value that reverses the color gradiant for
#' heatmap(s).
#'
#' @param transposedHeatmap a logical value that transposes heatmap rows to
#' columns and vice versa.
#'
#' @param simplifyBy a number that tells the function to change the values
#' smaller than that to zero. The purpose behind this option is to facilitate
#' recognizing candidate genes. Therefore, it is not suited for publications. It
#' has the same unit as \code{cutoff}.
#'
#' @param genesToDrop a character vector. Gene names within this vector will be
#' omitted from heatmap.The default value is \code{FALSE}.
#'
#'
#'
#' @return Based on preference, three heatmaps for \code{"Frequency.Percentage"}
#' , \code{"Mean.Value"} and \code{"Median.value"} can be generated. If more
#' than one group of genes are entered, output for each group will be strored in
#'  a separate sub-directory.
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
#' heatmapOutput(submissionName = "test")
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
################## Generating heatmap for the processed data ###################
################################################################################
################################################################################

heatmapOutput <- function(

  submissionName,

  shortenStudyNames = TRUE,

  geneLimit = 50,

  rankingMethod = "variation",

  heatmapFileFormat = "TIFF",

  resolution = 600,

  RowCex = "auto",

  ColCex = "auto",

  heatmapMargines = "auto",

  rowLabelsAngle = 0,

  columnLabelsAngle = 45,

  heatmapColor = "RdBu",

  reverseColor = TRUE,

  transposedHeatmap = FALSE,

  simplifyBy = FALSE,

  genesToDrop = FALSE

  ){

  ##############################################################################
  ########## Prerequisites

  # Check submissionName

  if(!missing(submissionName)){

    if(!is.character(submissionName)){

      stop("[heatmapOutput] 'submissionName' must be a character string!")

    }

  } else{

    stop("[heatmapOutput] 'submissionName' is mandatory!")

  }



  # Check shortenStudyNames

  if(!is.logical(shortenStudyNames)){

    stop("[heatmapOutput] 'shortenStudyNames' must be either TRUE or FALSE!")

  }



  # Check simplifyBy

  if(!geneLimit == FALSE & !is.numeric(geneLimit)){

    stop("[heatmapOutput] 'geneLimit' must either specify the maximum number of genes on the heatmap(s) or be FALSE.")

  }



  # Check rankingMethod

  if(!rankingMethod %in% c("variation", "highValue")){

    stop("[heatmapOutput] 'rankingMethod' must be either 'variation' or 'highValue'!")

  }



  # Check heatmap image file format

  if(!(heatmapFileFormat %in% c("TIFF", "JPG", "BMP", "PNG", "PDF"))){

    stop("[heatmapOutput] 'heatmapFileFormat' must be one of these formats: 'TIFF', 'JPG', 'BMP' ,'PNG' or 'PDF'")

  } else if(heatmapFileFormat == "PDF"){

    message("[heatmapOutput] 'resolution' is not applicable for PDF files.")

  }



  # Check resolution

  if(!is.numeric(resolution)){

    stop("[heatmapOutput] 'resolution' must be a number!")

  }



  # Check RowCex

  if(!RowCex == "auto" & !is.numeric(RowCex) |

     is.numeric(RowCex) & ! (RowCex >= 0 & RowCex <= 2)){

    stop("[heatmapOutput] 'RowCex' must be a number between 0 and 2!")

  }

  if(RowCex == "auto"){

    message("[heatmapOutput] Automatically determining 'RowCex'.")

  }



  # Check ColCex

  if(!ColCex == "auto" & !is.numeric(ColCex) |

     is.numeric(ColCex) & ! (ColCex >= 0 & ColCex <= 2)){

    stop("[heatmapOutput] 'ColCex' must be a number between 0 and 2!")

  }

  if(ColCex == "auto"){

    message("[heatmapOutput] Automatically determining 'ColCex'.")

  }



  # Check heatmapMargines

  if(is.character(heatmapMargines)){

    if(length(heatmapMargines) == 1){

      if(!heatmapMargines == "auto"){

        stop("[heatmapOutput] 'heatmapMargines' must be either a numerical vector of two numbers or be 'auto'!")

      } else{

        heatMapMode <- "algorithm"

        message("[heatmapOutput] Automatically determining 'heatmapMargines'.")

      }

    }else{

      stop("[heatmapOutput] 'heatmapMargines' must be either a numerical vector of two numbers or be 'auto'!")

    }

  }else if(is.numeric(heatmapMargines)){

    if(! length(heatmapMargines) == 2){

      stop("[heatmapOutput] 'heatmapMargines' must be either a numerical vector of two numbers or be 'auto'!")

    } else{

      heatMapMode <- "manual"

    }

  }



  # Check rowLabelsAngle

  if(! rowLabelsAngle == "auto" & ! is.numeric(rowLabelsAngle) |

     is.numeric(rowLabelsAngle) &

     ! (rowLabelsAngle >= 0 & rowLabelsAngle <= 360)){

    stop("[heatmapOutput] 'rowLabelsAngle' must be either 'auto' or a number ranging from 0 to 360.")

  }



  # Check columnLabelsAngle

  if(! columnLabelsAngle == "auto" & ! is.numeric(columnLabelsAngle) |

     is.numeric(columnLabelsAngle) &

     ! (columnLabelsAngle >= 0 & columnLabelsAngle <= 360)){

    stop("[heatmapOutput] 'columnLabelsAngle' must be either 'auto' or a number ranging from 0 to 360.")

  }



  # Check heatmapColor

  if(!heatmapColor %in% c("RdGr", "YlOrRd", "YlOrBl", "YlGnBu", "YlGn", "Reds",

                          "RdPu", "Purples", "PuRd", "PuBuGn", "PuBU", "OrRd",

                          "Oranges", "Greys", "Greens", "GnBu", "BuPu", "BuGn",

                          "Blues", "Set3", "Set2", "Set1", "Pastel2", "Pastel1",

                          "Paired", "Dark2", "Accent", "Spectral", "RdYlGn",

                          "RdYlBu", "RdGy", "RdBu", "PuOr", "PRGn", "PiYG",

                          "BrBG")){

    stop("[heatmapOutput] The entered 'heatmapColor' is not supported!")

  }



  # Check reverseColor

  if(!is.logical(reverseColor)){

    stop("[heatmapOutput] 'reverseColor' must be either TRUE or FALSE!")

  }



  # Check transposedHeatmap

  if(!is.logical(transposedHeatmap)){

    stop("[heatmapOutput] 'transposedHeatmap' must be either TRUE or FALSE!")

  }



  # Check transposedHeatmap

  # The FALSE argument is not removable, unfortunately.

  if(!simplifyBy == FALSE & !is.numeric(simplifyBy)){

    stop("[heatmapOutput] 'simplify' must be either FALSE or a numerical value!")

  }



  # Check genesToDrop

  # The FALSE argument is not removable, unfortunately.

  if(!genesToDrop == FALSE & !is.character(genesToDrop)){

    stop("[heatmapOutput] 'genesToDrop' must be a character vector of desired gene names!")

  }





  ##############################################################################
  ########## Decide whether function should stops now!

  # Check wheather the requested data exists

  database <- system.file("extdata", submissionName, package="cbaf")

  if(!dir.exists(database)){

    stop("[heatmapOutput] Please run one of the obtainSingleStudy() or obtainMultipleStudies() functions first, and then the automatedStatistics() function!")

  } else if(dir.exists(database)){

    bfc <- BiocFileCache(

      file.path(system.file("extdata", package = "cbaf"), submissionName),

      ask = FALSE

      )

    if(!nrow(bfcquery(bfc, c("Parameters for automatedStatistics()"))) == 1){

      stop("[heatmapOutput] Please run the automatedStatistics() function first!")

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



  # setting the value for cutoff

  if(desiredTechnique == "methylation"){

    cutoff.phrase <- "Mean methylation cutoff"

  } else{

    cutoff.phrase <- "log z-score cutoff"

  }





  # Store the new parameteres

  newParameters <-list()

  newParameters$submissionName <- submissionName

  newParameters$shortenStudyNames <- shortenStudyNames

  newParameters$geneLimit <- geneLimit

  newParameters$heatmapFileFormat <- heatmapFileFormat

  newParameters$resolution <- resolution

  newParameters$RowCex <- RowCex

  newParameters$ColCex <- ColCex

  newParameters$heatmapMargines <- heatmapMargines

  newParameters$rowLabelsAngle <- rowLabelsAngle

  newParameters$columnLabelsAngle <- columnLabelsAngle

  newParameters$heatmapColor <- heatmapColor

  newParameters$reverseColor <- reverseColor

  newParameters$transposedHeatmap <- transposedHeatmap

  newParameters$simplifyBy <- simplifyBy

  newParameters$genesToDrop <- genesToDrop





  # Check wheather the requested data exists

  number.of.rows.parameters <-

    nrow(bfcquery(bfc, "Parameters for heatmapOutput()"))


  if(number.of.rows.parameters == 1){

    oldParameters <-

      readRDS(

        bfcpath(bfc, bfcquery(bfc, c("Parameters for heatmapOutput()"))$rid)

        )

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

    stop("[heatmapOutput] Input database must be a list!")

  }





  # get the working directory

  parent.directory <- getwd()



  ##############################################################################
  ########## Set the function ready to work

  # Report

  message("[heatmapOutput] Preparing heatmap(s).")

  if(is.numeric(simplifyBy)){

    message("[heatmapOutput] Only significant results will be shown on heatmap(s)!")

  }



  # Count number of skipped heatmaps

  skipped <- 0



  # Create progressbar

  possible.subgroups <- c("Frequency.Percentage", "Mean.Value", "Median.Value")


  idx <- names(statisticsData[[1]]) %in% possible.subgroups

  total.number <- length(statisticsData)*length((statisticsData[[1]])[idx])


  heatmapOutputProgressBar <-

    txtProgressBar(min = 0, max = total.number, style = 3)

  ExtH <- 0





  ##############################################################################
  ########## Core segment

  # Save heatmaps in separate folder

  for(gr in seq_along(statisticsData)){

    # Subset data that can be presented as heatmap

    subset.name <- names(statisticsData)[gr]


    possible.subgroups.idx <-

      names(statisticsData[[gr]]) %in% possible.subgroups


    subset.data <- (statisticsData[[gr]])[possible.subgroups.idx]



    # Create a directory and set it as desired folder

    child.directory <- paste(

      gr, ". ", sub(x = subset.name, pattern = "\\.", replacement = "-"), sep=""

      )

    new.directory <- paste(parent.directory, child.directory, sep = "/")

    dir.create(new.directory, showWarnings = FALSE)

    setwd(new.directory)



    for(possible in seq_along(subset.data)){

      # subset statistics

      statistics.data <- subset.data[[possible]]



      # Remove desired genes

      if(!is.logical(genesToDrop)){

        if(!is.character(genesToDrop)){

          stop("[heatmapOutput] 'genesToDrop' must be a character vector of desired genes!")

        } else{

          filtered.colnames <- !(colnames(statistics.data) %in% genesToDrop)

          statistics.data <- statistics.data[, filtered.colnames, drop = FALSE]

        }

      }



      name.statistics.data <- names(subset.data)[possible]


      # determine ourput file name

      output.file.name <- paste0(

        gsub(x = name.statistics.data, pattern = "\\.", replacement = "-"),

        ", ",

        gsub(x = subset.name, pattern = "\\.", replacement = " "),

        " (",

        cutoff.phrase,

        "=",

        cutoff,

        ")",

        if(heatmapFileFormat == "TIFF"){

          ".tiff"

        }else if(heatmapFileFormat == "PNG"){

          ".png"

        }else if(heatmapFileFormat == "JPG"){

          ".jpg"

        }else if(heatmapFileFormat == "BMP"){

          ".bmp"

        }else if(heatmapFileFormat == "PDF"){

          ".pdf"

        }

      )






      # Check continue permission

      if(continue | !continue & !file.exists(output.file.name)){




        # Check whether study names should be shorted

        if(shortenStudyNames){

          rownames(statistics.data) <-

            sapply(

              strsplit(

                as.character(rownames(statistics.data)),

                split=" (",

                fixed=TRUE

                ),

              function(x) (x[1])

              )

        }



        heatmap.data <- t(statistics.data)

        # Removing NA

        not.just.na <- apply(heatmap.data, 1, function(x) any(!is.na(x)==TRUE))

        heatmap.data <- heatmap.data[not.just.na,]

        # Removing NaN

        heatmap.data[is.nan(heatmap.data)] <- 0

        # Removing rows that contain only 0

        if(is.matrix(heatmap.data)){

          heatmap.data <- heatmap.data[rowSums(heatmap.data, na.rm = TRUE)!=0,]

        }

        if(is.matrix(heatmap.data)){

          if(nrow(heatmap.data) > 1 & ncol(heatmap.data) > 1){




            # Limiting the number of genes in heatmap to get better resolution

            if(geneLimit==FALSE | is.numeric(geneLimit) &

               geneLimit > nrow(heatmap.data)){

              heatmap.data <- heatmap.data

            } else if(is.numeric(geneLimit) & geneLimit <= nrow(heatmap.data) &

                      rankingMethod == "variation"){

              ordering <- order(abs(rowVars(heatmap.data)), decreasing=TRUE)

              heatmap.data <- heatmap.data[ordering[seq_len(geneLimit)],]

            } else if(is.numeric(geneLimit) & geneLimit <= nrow(heatmap.data) &

                      rankingMethod == "highValue"){

              ordering <- order(abs(rowSums(heatmap.data)), decreasing=TRUE)

              heatmap.data <- heatmap.data[ordering[seq_len(geneLimit)],]

            } else{

              stop("[heatmapOutput] 'geneLimit' must be either a numerical value or FALSE!")

            }




            if(is.numeric(simplifyBy)){

              heatmap.data[heatmap.data < simplifyBy] <- 0

            }




            # Heatmap color

            if(reverseColor){

              if(heatmapColor == "RdGr"){

                hmcol <- rev(redgreen(75))

              } else {

                hmcol <- rev(colorRampPalette(brewer.pal(9, heatmapColor))(100))

              }


            } else if (!reverseColor){

              if(heatmapColor == "RdGr"){

                hmcol <- redgreen(75)

              } else {

                hmcol <- colorRampPalette(brewer.pal(9, heatmapColor))(100)

              }

            }





            ######  automatic parameters determination   ######

            if(ColCex == "auto"){

              if(ncol(heatmap.data) <= 18){

                d.ColCex <- 1.8 - ncol(heatmap.data) * 0.0333333333

              }else{

                d.ColCex <- 1 - ((ncol(heatmap.data) - 18)) * 0.0166666666

              }

            }else{

              d.ColCex <- ColCex

            }




            if(RowCex == "auto"){

              if(nrow(heatmap.data) <= 18){

                d.RowCex <- 1.8 - nrow(heatmap.data) * 0.0333333333

              }else{

                d.RowCex <- 1 - ((nrow(heatmap.data) - 18)) * 0.0166666666

              }

            }else{

              d.RowCex <- RowCex

            }


            # Equalizing RowCex and ColCex

            minCex <- min(c(d.ColCex, d.RowCex))

            d.ColCex <- minCex

            d.RowCex <- minCex


            # Check whether heatmapMargines is "auto"

            if(heatMapMode == "algorithm"){

              unitSize <- 1.19


              lengthDeterminant <- function(vector){

                relativeLengthVector <- vector("numeric", length = length(vector))


                for(vNames in seq_along(vector)){

                  currentName <- vector[vNames]

                  vectorLetters <- unlist(strsplit(currentName, ""))

                  startingSize <- 0


                  for(vLetters in seq_along(vectorLetters)){

                    currentLetter <- vectorLetters[vLetters]

                    if(currentLetter %in% c("a", "b", "c", "d", "e", "g", "h",
                                            "k", "n", "o", "p", "q", "s", "u",
                                            "v", "x", "y", "z", "2", "3", "4",
                                            "5", "6", "8", "9", "0", " ")){

                      startingSize <- startingSize + 0.855

                    } else if(currentLetter %in% c("f", "r", "j", "t", "7", "1")
                              ){

                      startingSize <- startingSize + 0.73

                    } else if(currentLetter %in% c("i", "l", "I")){

                      startingSize <- startingSize + 0.25

                    } else if(currentLetter %in% c("m", "w")){

                      startingSize <- startingSize + 1.20

                    } else if(currentLetter %in% c("B", "C", "D", "E", "F", "H",
                                                   "J", "K", "L", "N", "O", "P",
                                                   "Q", "R", "S", "T", "U", "V",
                                                   "X", "Y", "Z")){

                      startingSize <- startingSize + 0.98

                    } else if(currentLetter %in% c("A", "G", "M")){

                      startingSize <- startingSize + 1.10

                    } else if(currentLetter %in% c("W")){

                      startingSize <- startingSize + 1.25

                    }

                  }

                  relativeLengthVector[vNames] <- startingSize

                }

                max(relativeLengthVector)

              }


              # determining the best margin for column names

              # y = ax + b

              longest.study <-

                  lengthDeterminant(colnames(heatmap.data))*unitSize + 7.0

              longest.study.effect <-

                longest.study*abs(sin(columnLabelsAngle*0.0174532925))


              colMargin <- longest.study.effect * d.ColCex * 0.4278074866



              # determining the best margin for row names

              # y = ax + b

              longest.gene <-

                lengthDeterminant(rownames(heatmap.data))*unitSize + 7.0

              longest.gene.effect <-

                longest.gene*abs(cos(rowLabelsAngle*0.0174532925))


              rowMargin <- longest.gene.effect * d.RowCex * 0.4278074866



              # determining which margine influence the final vector

              largestMargine <- max(colMargin, rowMargin)

              d.heatmapMargines <- c(largestMargine, largestMargine)



            }else if(heatMapMode == "manual"){

              d.heatmapMargines <- heatmapMargines

            }



            ###################################################





            # Drawing heatmap

            if(heatmapFileFormat == "TIFF"){


              tiff(

                filename = paste(getwd(), output.file.name, sep="/"),

                width = 11,

                height = 11,

                units = "in",

                res = resolution,

                compression = "lzw"

              )


            }else if(heatmapFileFormat == "PNG"){


              png(

                filename = paste(getwd(), output.file.name, sep="/"),

                width = 11,

                height = 11,

                units = "in",

                res = resolution

              )


            }else if(heatmapFileFormat == "BMP"){


              bmp(

                filename = paste(getwd(), output.file.name, sep="/"),

                width = 11,

                height = 11,

                units = "in",

                res = resolution

              )


            }else if(heatmapFileFormat == "JPG"){


              jpeg(

                filename=paste(getwd(), output.file.name, sep="/"),

                width = 11,

                height = 11,

                units = "in",

                res = resolution

              )


            }else if(heatmapFileFormat == "PDF"){


              pdf(

                file=paste(getwd(), output.file.name, sep="/"),

                width = 11,

                height = 11

              )


            }





            # detemining the oriantation of heatmap

            if(!transposedHeatmap){

                heatmap.input.matrix <- heatmap.data

                labCol <- colnames(heatmap.data)


            } else if(transposedHeatmap){

                heatmap.input.matrix <- t(heatmap.data)

                labCol <- rownames(heatmap.data)


            }



            # Draw heatmap

            heatmap.2(

              heatmap.input.matrix,

              labCol=labCol,

              na.color= "light gray",

              trace="none",

              symbreaks = TRUE,

              col= hmcol,

              cexRow = d.RowCex,

              cexCol= d.ColCex,

              margins = d.heatmapMargines,

              srtRow = rowLabelsAngle,

              srtCol = columnLabelsAngle

            )



            dev.off()



            # Crop margines of the stored image

            # cropped.image <- image_read(output.file.name)

            # cropped.image <- image_crop(cropped.image, "1000x1500+500")

            # image_write(cropped.image,

            #             path = output.file.name,

            #             format = if(heatmapFileFormat == "TIFF"){

            #               "tiff"

            #             }else if(heatmapFileFormat == "PNG"){

            #               "png"

            #             }else if(heatmapFileFormat == "JPG"){

            #               "jpg"

            #             }else if(heatmapFileFormat == "BMP"){

            #               "bmp"

            #            })


          }

        }

      }else{

        skipped <- skipped + 1

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

    message("[heatmapOutput] ", as.character(skipped), " out of ", as.character(total.number)," heatmaps were skipped, because they already exist!")

  } else if(skipped > 0 & skipped == 1){

    message("[heatmapOutput] ", as.character(skipped), " out of ", as.character(total.number)," heatmaps was skipped, because it already exists!")

  }



  # Store the last parameter

  oldParamHeatmapOutput <- newParameters


  # Store the parameters for this run

  if(number.of.rows.parameters == 0){

    saveRDS(

      oldParamHeatmapOutput,

      file=bfcnew(bfc, "Parameters for heatmapOutput()", ext="RDS")

      )

  } else if(number.of.rows.parameters == 1){

    saveRDS(

      oldParamHeatmapOutput,

      file=bfc[[bfcquery(bfc, "Parameters for heatmapOutput()")$rid]]

      )

  }



  # change directory to parent directory

  setwd(parent.directory)

  # message("[heatmapOutput] Finished.")

}
