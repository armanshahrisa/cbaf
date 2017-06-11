#' @title Checking Expression/methylation Profile for various subgroups of a cancer study.
#'
#' @description This function Obtains the requested data for the given genes across multiple subgroups of a cancer. It can
#' check whether or not all genes are included in subgroups of a cancer study and, if not, looks for the alternative gene names.
#' The main part of function calculates frequency percentage, frequency ratio, mean expression and median of samples greather than
#' specific value in the selected subgroups of the cancer. Furthermore, it looks for the genes that comprise the highest values
#' in each cancer subgroup.
#'
#' @details
#' \tabular{lllll}{
#' Package: \tab cBioAutomatedTools \cr
#' Type: \tab Package \cr
#' Version: \tab 0.99.0 \cr
#' Date: \tab 2017-06-01 \cr
#' License: \tab Artistic-2.0 \cr
#' }
#'
#' @import cgdsr Biobase
#'
#' @usage process.one.study(genes, cancername, high.throughput.data.type,
#' data.presented.as = c("Frequency.Percentage", "Frequency.Ratio", "Mean.Value", "Median"),
#' transposedHeatmap=FALSE, desired.case.list="None", genelimit="none",
#' resolution=600, RowCex=0.8, ColCex=0.8, heatmapMargines=c(10,10), cutoff=NULL,
#' angle.for.heatmap.cancernames=45, heatmap.color = "RdBu", reverse.heatmap.color = TRUE,
#' rewrite.output.list = TRUE, round=TRUE, top.genes = TRUE, validate.genes = TRUE,
#' simplify.visulization=FALSE, simplifiction.cuttoff=FALSE)
#'
#' @param genes a list that contains at least one gene group
#'
#' @param cancername a character string showing the desired cancer name. It is an standard cancer study name that can be found on
#' cbioportal.org, such as \code{Acute Myeloid Leukemia (TCGA, NEJM 2013)}.
#'
#' @param high.throughput.data.type a character string that is one of the following techniques: 'RNA-seq', 'microRNA-Seq',
#' 'microarray.mRNA', 'microarray.microRNA' or 'methylation'.
#'
#' @param data.presented.as a character vector that containes the statistical precedures users prefer the function to compute.
#' default is \code{c("Frequency.Percentage", "Frequency.Ratio", "Mean.Value", "Median")}. This will tell the function to
#' compute the following:
#' frequency precentage, which is the number of samples having value greather than specific cutoff divided by the total sample
#' size for each cancer;
#' frequency ratio, which shows the number of selected samples divided by the total number of samples that give the frequency
#' percentage
#' for each cancer-to know selecected and total sample sizes only;
#' Mean Expression, that contains mean value of selected samples for each cancer;
#' Median, which shows the median value of selected samples for each cancer.
#'
#' @param transposedHeatmap a logical value that can transpose the rows and columns of heatmaps.
#'
#' @param desired.case.list a numeric vector that contains the index of rdesired cancer subgroups, if user knows index of desired
#' subgroups. If set to be \code{none}, function will ask the user to enter them during the process. The default value is
#' \code{none}.
#'
#' @param genelimit if large number of genes exists in at least one gene group, this option can be use to limit the number of
#' genes to be shown on hitmap. For instance, \code{genelimit=50} will limit the heatmap to 50 genes showing the most variation.
#' The default value is \code{none}.
#'
#' @param resolution a number. This option can be used to adjust the resolution of the output heatmaps as 'dot per inch'. The
#' defalut value is 600.
#'
#' @param RowCex a number that specifies letter size in heatmap row names.
#'
#' @param ColCex a number that specifies letter size in heatmap column names.
#'
#' @param heatmapMargines a numeric vectors that can be use to set heatmap margins. The default value is
#' \code{heatmapMargines=c(10,10)}.
#'
#' @param cutoff a number used to limit samples to those that are greather than specific number (cutoff). The default value for
#' methylation data is 0.6 while gene expression studies use default value of 2. For methylation studies, it is
#' \code{observed/expected ratio}, for the rest, it is \code{z-score} - Function can understand the correct unit. TO change the
#' cutoff to any desired number, change the option to \code{cutoff = desiredNumber} in which \code{desiredNumber} is the number
#' of interest.
#'
#' @param angle.for.heatmap.cancernames a number that determines the angle with which the cancer names are shown in heatmaps.
#' The default value is 45 degree.
#'
#' @param heatmap.color a character string matches standard color names. The default value is "RdBu". "redgreen" is also a popular
#' color in genomic studies. To see the rest of colors, please type \code{display.brewer.all()}.
#'
#' @param reverse.heatmap.color a logical value that reverses the color gradiant for heatmap.
#'
#' @param rewrite.output.list a logical value. This option can be used to modify heatmap, for instance change margin,
#' without obtaining data from internet again. If set to false, function will use the values previously stored in the
#' global environment to draw heatmaps and save excel file. The default value is \code{TRUE}.
#'
#' @param round a logical value that, if set to be \code{TRUE}, will force the function to round all the calculated values
#' to two decimal places. The default value is \code{TRUE}.
#'
#' @param top.genes a logical value that, if set as \code{TRUE}, cause the function to create three data.frame that contain the
#' gene names with the highest values for each cancer. To get all the three data.frames, `Frequency.Percentage`, `Mean.Value` and
#' `Median` must have been included for \code{data.presented.as}.
#'
#' @param validate.genes a logical value that, if set to be \code{TRUE}, function will checks each cancer study to finds whether
#' or not each gene has a record. If the given cancer doesn't have a record for specific gene, it checks for alternative gene
#' names that cbioportal might use instead of the given gene name.
#'
#' @param simplify.visulization a logical value that tells the function whether or not to change values under
#' \code{simplifiction.cuttoff} to zero. It only affects heatmaps to assist finding the candidate genes faster. Therefore, it is
#' not suited for publications.
#'
#' @param simplifiction.cuttoff a logical value that, if \code{simplify.visulization = TRUE}, needs to be set as a desired cuttoff
#' for \code{simplify.visulization}. It has the same unit as \code{cutoff}.
#'
#' @return a list that containes some or all of the following groups, based on what user has chosen: \code{Validation.Results},
#' \code{Frequency.Percentage}, \code{Top.Genes.of.Frequency.Percentage}, \code{Frequency.Ratio}, \code{Mean.Value},
#' \code{Top.Genes.of.Mean.Value}, \code{Median}, \code{Top.Genes.of.Median}. It also saves these groups in one excel
#' file for convenience. Based on preference, three heatmaps for \code{Frequency.Percentage}, \code{Mean.Value} and
#' \code{Median} can be generated. If more than one gene group is entered, output for each group will be strored in a separate sub-directory.
#'
#' @examples
#' # Creating a list that contains one gene group: 'K.demethylases'
#' genes <- list(K.demethylases = c("KDM1A", "KDM1B", "KDM2A"))
#'
#' # Running the function to obtain and process the selected data
#' process.one.study(genes, "Breast Invasive Carcinoma (TCGA, Cell 2015)", "RNA-seq", desired.case.list = c(3,4,5),
#' data.presented.as = c("Frequency.Percentage", "Frequency.Ratio", "Mean.Value"), heatmap.color = "redgreen")
#'
#' @author Arman Shahrisa, \email{shahrisa.arman@hotmail.com} [maintainer, copyright holder]
#' @author Maryam Tahmasebi Birgani, \email{tahmasebi-ma@ajums.ac.ir}
#'



#########################################################################
#########################################################################
########### Automatically calculate statistical measurements ############
#########################################################################
#########################################################################

automatedStatistics<- function(obtainedDataType = "multiple studies", databaseType, submissionName, calculate = c("FrequencyPercentage", "Frequency.Ratio", "MeanValue", "Median"),

                              cutoff=NULL, round=TRUE, topGenes = TRUE, validateGenes = TRUE){

  ##########################################################################
  ########## Prerequisites

  # Obtain the unprocessed data list

  if(obtainedDataType == "multiple studies"){

    databaseSymbol <- "obM"

    validationSymbol <- "vaM"

    haultType <- "ohaultM"

  } else if(obtainedDataType == "single study"){

    databaseSymbol <- "obS"

    validationSymbol <- "vaS"

    haultType <- "ohaultS"

  }



  # high-throughput data type

  if(is.character(desiredTechnique)){

    if(!(desiredTechnique %in% c("RNA-seq", "microRNA-Seq", "Microarray.mRNA", "Microarray.microRNA", "methylation")) | length(desiredTechnique)!= 1){

      stop("'desiredTechnique' must contain one of the following techniques: 'RNA-seq', 'microRNA-Seq', 'Microarray.mRNA', 'Microarray.microRNA' or

           'methylation'")

    }

  } else {

    stop("'desiredTechnique' must be entered as a character string describing a technique name")

  }



  # setting the calue for cutoff

  if(high.throughput.data.type == "methylation"){

    cutoff.phrase <- "obs/exp cutoff"

    if(is.null(cutoff)){

      cutoff <- 0.6

    }

  } else{

    cutoff.phrase <- "zscore cutoff"

    if(is.null(cutoff)){

      cutoff <- 2

    }

  }



  # Check hault older

  if(exists(paste(paste(haultType, ".", submissionName, sep = "")))){


    if(exists(paste(paste(cut, ".", submissionName, sep = "")))){

      oldCutoff <- get(paste(paste(cut, ".", submissionName, sep = "")))

      if(identical(oldCutoff, cutoff)){

        haultDecision1 <- NULL

      }

    }



    if(exists(paste(paste(cal, ".", submissionName, sep = "")))){

      oldCalculate <- sort(get(paste(paste(cal, ".", submissionName, sep = ""))))

      newCalculate <- sort(calculate)


      if(identical(oldCalculate, newCalculate)){

        haultDecision2 <- NULL

      }

    }



    # Halt the function

    if(all(exists("haultDecision1"), exists("haultDecision2"))){

      return("--- Function 'automatedStatistics()' was skipped: The requested statistical parameters are already calculated ---")

    }

  }















  sourceDataList <- get(paste(databaseSymbol, ".", submissionName, sep = ""))


  *************exists(validation)***********************


  # Create a progressbar

  automatedStatisticsProgressBar <- txtProgressBar(min = 0, max = length(studiesNames), style = 3)



  # Prepairing for a for loop

  if(is.list(sourceDataList)){

    numberDataSets <- length(sourceDataListt)

  } else{

    DataSets <- 1

  }



  ##########################################################################
  ########## Core segment

  # calculating the reuested statistical measurements as a for loop

  for(cal in 1:numberDataSets){

    # Subsetting matrixes from List by 'for' control structure

    if(is.list(sourceDataList)){

      oneDataMatrix <- sourceDataList[[cal]]

      oneDataMatrix <- names(sourceDataList)[[cal]]

    } else{

      oneDataMatrix <- sourceDataList

    }





    setTxtProgressBar(automatedStatisticsProgressBar, cal)

  }

  close(automatedStatisticsProgressBar)

  # Export the obtained data as list

  assign(paste("PrData", ":", submissionName, sep = ""), obtainedData, envir = globalenv())

}





# Get the order of genes in obtained expression data

ordered.Genes <- colnames(ProfileData)

Expected.Genes <- c.genes[order(c.genes)]






# Create a matrix contained genes that are not included in required cancer

if(validate.genes == TRUE){

  # Alter c.genes to be compatible with gene names in cBioPortal output

  Altered.sorted.gene.names <- sort(gsub("-", ".", c.genes))

  # Obtain name of genes that are absent in requested cancer

  Absent.genes <- Altered.sorted.gene.names[!Altered.sorted.gene.names %in% colnames(ProfileData)]

  # For loop for determining changed genes

  Alternative.gene.names <- vector("character", length = length(Absent.genes))

  for(ab in 1:length(Absent.genes)){

    Absent.gene.name.ProfileData <- colnames(data.matrix(getProfileData(mycgds, Absent.genes[ab], mygeneticprofile,mycaselist)))

    # Check wheter gene has an alternative name or missed in the database

    if(length(Absent.gene.name.ProfileData) == 1){

      Alternative.gene.names[ab] <- Absent.gene.name.ProfileData

    } else if(length(Absent.gene.name.ProfileData) != 1){

      Alternative.gene.names[ab] <- "-"

    }

  }

  # Naming Alternative.gene.names

  names(Alternative.gene.names) <- Absent.genes

  # Seperating genes with alternative names from those that are absent

  Gene.with.no.data <- Alternative.gene.names[Alternative.gene.names == "-"]

  Gene.with.data <- Alternative.gene.names[Alternative.gene.names != "-"]




  # Switch for Absent.genes

  if(length(Absent.genes) != 0 & any(Alternative.gene.names == "-")){



    ## Create profileData which containes absent genes

    # Find sample numbers in ProfileData

    numberOfRows <- nrow(ProfileData)

    # Create NA variable to be bind to Profile data with same sample size

    NotFound <- matrix(rep(NA, length(Gene.with.no.data)*numberOfRows), ncol = length(Gene.with.no.data), nrow = numberOfRows)

    # make sample names identical to sample names of ProfileData

    colnames(NotFound) <- names(Gene.with.no.data)

    # Merge two ProfileData

    CompletedProfileData <- cbind(ProfileData, NotFound)

    # Sort the merged ProfileData according to alphabetic name of genes (In order to mix to ProfileData)

    CompletedProfileData <- CompletedProfileData[,sort(colnames(CompletedProfileData))]

  } else if(length(Absent.genes) == 0 | (length(Absent.genes) != 0 & any(Alternative.gene.names != "-"))){

    CompletedProfileData <- ProfileData

  }

}






# Creating Empty matrixes to be filled with 'for' control structure

if(validate.genes == TRUE){

  Complete.ordered.Genes <- colnames(CompletedProfileData)

  Serial.ProfileData.for.Validation <- vector("character", length=length(Expected.Genes))

  names(Serial.ProfileData.for.Validation) <- Complete.ordered.Genes

}

if("Frequency.Percentage" %in% data.presented.as){

  Serial.ProfileData.with.Frequency.Percentage <- vector("numeric", length=length(ordered.Genes))

}

if("Frequency.Ratio" %in% data.presented.as){

  Serial.ProfileData.with.Frequency.Ratio <- vector("character", length=length(ordered.Genes))

}

if("Mean.Value" %in% data.presented.as){

  Serial.ProfileData.with.Mean.Value <- vector("numeric", length=length(ordered.Genes))

}

if("Median" %in% data.presented.as){

  Serial.ProfileData.with.Median <- vector("numeric", length=length(ordered.Genes))

}






### Calculating validity, frequency, Mean.Value and Median values for each requested gene

if(validate.genes == TRUE){

  for(expect in 1:length(Expected.Genes)){

    ## Validating Genes

    # Correct those that are not found

    if(length((CompletedProfileData[,expect])[!is.nan(CompletedProfileData[,expect])]) > 0 & all(!is.finite(CompletedProfileData[,expect])) &

       is.nan(mean(as.vector(CompletedProfileData[,expect])[abs(CompletedProfileData[,expect])], na.rm=TRUE))){

      Serial.ProfileData.for.Validation[expect] <- "-"

    } else {

      Serial.ProfileData.for.Validation[expect] <- "Found"

    }

  }

  # Replacing name and alternative names

  for(re in 1:length(Gene.with.data)){

    replace.index <- which(names(Serial.ProfileData.for.Validation) == Gene.with.data[re])

    Serial.ProfileData.for.Validation[replace.index] <- paste("Alternative name: ", Gene.with.data[re], sep="")

    names(Serial.ProfileData.for.Validation)[replace.index] <- names(Gene.with.data)[re]

  }

  # Sorting genes alphabetically

  Serial.ProfileData.for.Validation <- Serial.ProfileData.for.Validation[sort(names(Serial.ProfileData.for.Validation))]

}



for(j in 1:length(ordered.Genes)){

  # Frequency.Percentage

  if("Frequency.Percentage" %in% data.presented.as){

    # Check all members are under cutoff

    if(!is.na(mean(as.vector(ProfileData[,j])[abs(ProfileData[,j])])) &

       is.nan(mean(as.vector(ProfileData[,j])[abs(ProfileData[,j]) >= cutoff], na.rm=TRUE))){

      Serial.ProfileData.with.Frequency.Percentage[j] <- 0

      # Check all members are NaN

    } else if (length((ProfileData[,j])[!is.nan(ProfileData[,j])]) == 0 & all(!is.finite(ProfileData[,j])) &

               is.nan(mean(as.vector(ProfileData[,j])[abs(ProfileData[,j]) >= cutoff], na.rm=TRUE))){

      Serial.ProfileData.with.Frequency.Percentage[j] <- NaN

      # Check all members are NA

    } else if (length((ProfileData[,j])[!is.nan(ProfileData[,j])]) > 0 & all(!is.finite(ProfileData[,j])) &

               is.nan(mean(as.vector(ProfileData[,j])[abs(ProfileData[,j]) >= cutoff], na.rm=TRUE))){

      Serial.ProfileData.with.Frequency.Percentage[j] <- NA

      # Mean is bigger than 0

    } else if (mean(as.vector(ProfileData[,j])[abs(ProfileData[,j]) >= cutoff], na.rm=TRUE) > 0 &

               !is.nan(mean(as.vector(ProfileData[,j])[abs(ProfileData[,j]) >= cutoff], na.rm=TRUE))){

      Serial.ProfileData.with.Frequency.Percentage[j] <- 100*mean(as.vector(abs(ProfileData[,j]) >= cutoff), na.rm=TRUE)

      # Mean is smaller than 0

    } else if (mean(as.vector(ProfileData[,j])[abs(ProfileData[,j]) >= cutoff], na.rm=TRUE) < 0 &

               !is.nan(mean(as.vector(ProfileData[,j])[abs(ProfileData[,j]) >= cutoff], na.rm=TRUE))){

      Serial.ProfileData.with.Frequency.Percentage[j] <- -100*(mean(as.vector(abs(ProfileData[,j]) >= cutoff), na.rm=TRUE))

    }

  }


  # Frequency Ratio

  if("Frequency.Ratio" %in% data.presented.as){

    # Check all members are under cutoff

    if(!is.na(mean(as.vector(ProfileData[,j])[abs(ProfileData[,j])])) &

       is.nan(mean(as.vector(ProfileData[,j])[abs(ProfileData[,j]) >= cutoff], na.rm=TRUE))){

      Serial.ProfileData.with.Frequency.Ratio[j] <- paste(as.character(0), "/", as.character(length(as.vector(ProfileData[,j]))), sep="")

      # Check all members are NaN

    } else if (length((ProfileData[,j])[!is.nan(ProfileData[,j])]) == 0 & all(!is.finite(ProfileData[,j])) &

               is.nan(mean(as.vector(ProfileData[,j])[abs(ProfileData[,j]) >= cutoff], na.rm=TRUE))){

      Serial.ProfileData.with.Frequency.Ratio[j] <- NaN

      # Check all members are NA

    } else if (length((ProfileData[,j])[!is.nan(ProfileData[,j])]) > 0 & all(!is.finite(ProfileData[,j])) &

               is.nan(mean(as.vector(ProfileData[,j])[abs(ProfileData[,j]) >= cutoff], na.rm=TRUE))){

      Serial.ProfileData.with.Frequency.Ratio[j] <- NA

      # Mean is number

    } else if (!is.nan(mean(as.vector(ProfileData[,j])[abs(ProfileData[,j]) >= cutoff], na.rm=TRUE))){

      Serial.ProfileData.with.Frequency.Ratio[j] <- paste(as.character(length(na.omit(as.vector(ProfileData[,j])[abs(ProfileData[,j]) >= cutoff]))),

                                                          "/", as.character(length(as.vector(ProfileData[,j]))), sep="")

    }

  }


  # Mean.Value

  if("Mean.Value" %in% data.presented.as){

    # Check all members are under cutoff

    if(!is.na(mean(as.vector(ProfileData[,j])[abs(ProfileData[,j])])) &

       is.nan(mean(as.vector(ProfileData[,j])[abs(ProfileData[,j]) >= cutoff], na.rm=TRUE))){

      Serial.ProfileData.with.Mean.Value[j] <- 0

      # Check all members are NaN

    } else if (length((ProfileData[,j])[!is.nan(ProfileData[,j])]) == 0 & all(!is.finite(ProfileData[,j])) &

               is.nan(mean(as.vector(ProfileData[,j])[abs(ProfileData[,j]) >= cutoff], na.rm=TRUE))){

      Serial.ProfileData.with.Mean.Value[j] <- NaN

      # Check all members are NA

    } else if (length((ProfileData[,j])[!is.nan(ProfileData[,j])]) > 0 & all(!is.finite(ProfileData[,j])) &

               is.nan(mean(as.vector(ProfileData[,j])[abs(ProfileData[,j]) >= cutoff], na.rm=TRUE))){

      Serial.ProfileData.with.Mean.Value[j] <- NA

      # Mean is number

    } else if (!is.nan(mean(as.vector(ProfileData[,j])[abs(ProfileData[,j]) >= cutoff], na.rm=TRUE))){

      Serial.ProfileData.with.Mean.Value[j] <- mean(as.vector(ProfileData[,j])[abs(ProfileData[,j]) >= cutoff], na.rm=TRUE)

    }

  }


  # Median

  if("Median" %in% data.presented.as){

    # Check all members are under cutoff

    if(!is.na(mean(as.vector(ProfileData[,j])[abs(ProfileData[,j])])) &

       is.nan(mean(as.vector(ProfileData[,j])[abs(ProfileData[,j]) >= cutoff], na.rm=TRUE))){

      Serial.ProfileData.with.Median[j] <- 0

      # Check all members are NaN

    } else if (length((ProfileData[,j])[!is.nan(ProfileData[,j])]) == 0 & all(!is.finite(ProfileData[,j])) &

               is.nan(mean(as.vector(ProfileData[,j])[abs(ProfileData[,j]) >= cutoff], na.rm=TRUE))){

      Serial.ProfileData.with.Median[j] <- NaN

      # Check all members are NA

    } else if (length((ProfileData[,j])[!is.nan(ProfileData[,j])]) > 0 & all(!is.finite(ProfileData[,j])) &

               is.nan(mean(as.vector(ProfileData[,j])[abs(ProfileData[,j]) >= cutoff], na.rm=TRUE))){

      Serial.ProfileData.with.Median[j] <- NA

      # Mean is number

    } else if (!is.nan(mean(as.vector(ProfileData[,j])[abs(ProfileData[,j]) >= cutoff], na.rm=TRUE))){

      Serial.ProfileData.with.Median[j] <- median(as.vector(ProfileData[,j])[abs(ProfileData[,j]) >= cutoff], na.rm=TRUE)

    }

  }

}


# Converting each ProfileData to independant cartical matrix with one column

if(validate.genes == TRUE){

  Original.validate.gene.names <- names(Serial.ProfileData.for.Validation)

  dim(Serial.ProfileData.for.Validation) <- c(length(Expected.Genes),1)

  rownames(Serial.ProfileData.for.Validation) <- Original.validate.gene.names

}

if("Frequency.Percentage" %in% data.presented.as){

  dim(Serial.ProfileData.with.Frequency.Percentage) <- c(length(ordered.Genes),1)

}

if("Frequency.Ratio" %in% data.presented.as){

  dim(Serial.ProfileData.with.Frequency.Ratio) <- c(length(ordered.Genes),1)

}

if("Mean.Value" %in% data.presented.as){

  dim(Serial.ProfileData.with.Mean.Value) <- c(length(ordered.Genes),1)

}

if("Median" %in% data.presented.as){

  dim(Serial.ProfileData.with.Median) <- c(length(ordered.Genes),1)

}




# Combining serries of columns to obtain four matrixes presenting 'Frequency.Percentage', 'Frequency in ratio', 'Mean.Value', 'Median' in a specific cancer

if(i==1){

  if(validate.genes == TRUE){

    Validation.matrix <- Serial.ProfileData.for.Validation

  }

  if("Frequency.Percentage" %in% data.presented.as){

    Frequency.Percentage.Matrix <- Serial.ProfileData.with.Frequency.Percentage

  }

  if("Frequency.Ratio" %in% data.presented.as){

    Frequency.Ratio.Matrix <- Serial.ProfileData.with.Frequency.Ratio

  }

  if("Mean.Value" %in% data.presented.as){

    Mean.Value.Matrix <- Serial.ProfileData.with.Mean.Value

  }

  if("Median" %in% data.presented.as){

    Median.Matrix <- Serial.ProfileData.with.Median

  }

}else if(i!=1){

  if(validate.genes == TRUE){

    Validation.matrix <- cbind(Validation.matrix, Serial.ProfileData.for.Validation)

  }

  if("Frequency.Percentage" %in% data.presented.as){

    Frequency.Percentage.Matrix <- cbind(Frequency.Percentage.Matrix, Serial.ProfileData.with.Frequency.Percentage)

  }

  if("Frequency.Ratio" %in% data.presented.as){

    Frequency.Ratio.Matrix <- cbind(Frequency.Ratio.Matrix, Serial.ProfileData.with.Frequency.Ratio)

  }

  if("Mean.Value" %in% data.presented.as){

    Mean.Value.Matrix <- cbind(Mean.Value.Matrix, Serial.ProfileData.with.Mean.Value)

  }

  if("Median" %in% data.presented.as){

    Median.Matrix <- cbind(Median.Matrix, Serial.ProfileData.with.Median)

  }

}






# update progress bar

setTxtProgressBar(obtainOneStudyProgressBar, i)

}

# Close progress bar

close(obtainOneStudyProgressBar)






### Adding dimention names to four matrixes

print("Preparing a 'list' to fill with the requested data")






## Creating the required matrixes and Overal list

# Create empty list

Data.list <- list()

# filling the list

if(validate.genes == TRUE){

  colnames(Validation.matrix) <- inputCases.names

  Data.list$Validation.Results <- t(Validation.matrix)

}

if("Frequency.Percentage" %in% data.presented.as){

  dimnames(Frequency.Percentage.Matrix) <- list(ordered.Genes, inputCases.names)

  Data.list$Frequency.Percentage <- t(Frequency.Percentage.Matrix)

  if(top.genes==TRUE){

    # Check if manual naming is requested

    Exp.matrix <- t(Frequency.Percentage.Matrix)

    # Removing NaN

    Exp.matrix[is.nan(Exp.matrix)] <- 0

    # Removing NA

    Exp.matrix[is.na(Exp.matrix)] <- 0

    # Removing rows that contain only 0

    Exp.matrix <- Exp.matrix[rowSums(Exp.matrix)!=0,]

    # Empty data.frame to be filled with data

    Frequency.Percentage.top.genes.dataframe <- data.frame(Gene.with.highest.value=vector("character", length = nrow(Exp.matrix)), Frequency.Percentage=vector("numeric", length = nrow(Exp.matrix)), stringsAsFactors = FALSE)

    # Finding genes with highest value

    for (m in 1:nrow(Exp.matrix)){

      Frequency.Percentage.top.genes.dataframe[m,1] <- names(which.max(Exp.matrix[m,]))

    }

    # Finding values corresponding to these genes

    for (n in 1:nrow(Exp.matrix)){

      Frequency.Percentage.top.genes.dataframe[n,2] <- max(Exp.matrix[n,])

    }

    # Creating data.frame for desired profile

    rownames(Frequency.Percentage.top.genes.dataframe) <- rownames(Exp.matrix)

    Data.list$Top.Genes.of.Frequency.Percentage <- Frequency.Percentage.top.genes.dataframe

  }

}

if("Frequency.Ratio" %in% data.presented.as){

  dimnames(Frequency.Ratio.Matrix) <- list(ordered.Genes, inputCases.names)

  Data.list$Frequency.Ratio <- t(Frequency.Ratio.Matrix)

}

if("Mean.Value" %in% data.presented.as){

  dimnames(Mean.Value.Matrix) <- list(ordered.Genes, inputCases.names)

  Data.list$Mean.Value <- t(Mean.Value.Matrix)

  if(top.genes==TRUE){

    # Check if manual naming is requested

    Exp.matrix <- t(Mean.Value.Matrix)

    # Removing NaN

    Exp.matrix[is.nan(Exp.matrix)] <- 0

    # Removing NA

    Exp.matrix[is.na(Exp.matrix)] <- 0

    # Removing rows that contain only 0

    Exp.matrix <- Exp.matrix[rowSums(Exp.matrix)!=0,]

    # Empty data.frame to be filled with data

    Mean.Value.top.genes.dataframe <- data.frame(Gene.with.highest.value=vector("character", length = nrow(Exp.matrix)), Mean.Value=vector("numeric", length = nrow(Exp.matrix)), stringsAsFactors = FALSE)

    # Finding genes with highest value

    for (m in 1:nrow(Exp.matrix)){

      Mean.Value.top.genes.dataframe[m,1] <- names(which.max(Exp.matrix[m,]))

    }

    # Finding values corresponding to these genes

    for (n in 1:nrow(Exp.matrix)){

      Mean.Value.top.genes.dataframe[n,2] <- max(Exp.matrix[n,])

    }

    # Creating data.frame for desired profile

    rownames(Mean.Value.top.genes.dataframe) <- rownames(Exp.matrix)

    Data.list$Top.Genes.of.Mean.Value <- Mean.Value.top.genes.dataframe

  }

}

if("Median" %in% data.presented.as){

  dimnames(Median.Matrix) <- list(ordered.Genes, inputCases.names)

  Data.list$Median <- t(Median.Matrix)

  if(top.genes==TRUE){

    # Check if manual naming is requested

    Exp.matrix <- t(Median.Matrix)

    # Removing NaN

    Exp.matrix[is.nan(Exp.matrix)] <- 0

    # Removing NA

    Exp.matrix[is.na(Exp.matrix)] <- 0

    # Removing rows that contain only 0

    Exp.matrix <- Exp.matrix[rowSums(Exp.matrix)!=0,]

    # Empty data.frame to be filled with data

    Median.top.genes.dataframe <- data.frame(Gene.with.highest.value=vector("character", length = nrow(Exp.matrix)), Median=vector("numeric", length = nrow(Exp.matrix)), stringsAsFactors = FALSE)

    # Finding genes with highest value

    for (m in 1:nrow(Exp.matrix)){

      Median.top.genes.dataframe[m,1] <- names(which.max(Exp.matrix[m,]))

    }

    # Finding values corresponding to these genes

    for (n in 1:nrow(Exp.matrix)){

      Median.top.genes.dataframe[n,2] <- max(Exp.matrix[n,])

    }

    # Creating data.frame for desired profile

    rownames(Median.top.genes.dataframe) <- rownames(Exp.matrix)

    Data.list$Top.Genes.of.Median <- Median.top.genes.dataframe

  }

}


assign(paste(paste(cut, ".", submissionName, sep = "")), cutoff, envir = globalenv())

assign(paste(paste(cal, ".", submissionName, sep = "")), calculate, envir = globalenv())


# Saving data list to enviroment with an informative name



assign(paste(names(genes)[[g]], ".", "Data.List",  "_", "for", "_", Shortened.cancername, sep=""), Data.list, envir = globalenv())

}
