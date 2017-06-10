#' @title Checking Expression/methylation Profile for various cancer studies.
#'
#' @description This function Obtains the requested data for the given genes across multiple cancer studies. It can check
#' whether or not all genes are included in cancer studies and, if not, looks for the alternative gene names. Tha main part
#' of function calculates frequency percentage, frequency ratio, mean expression and median of samples greather than specific
#' value in the selected cancer studies. Furthermore, it looks for the genes that comprise the highest values in each cancer study.
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
#' @usage process.multiple.studies(genes, cancers, high.throughput.data.type,
#' data.presented.as = c("Frequency.Percentage", "Frequency.Ratio", "Mean.Value", "Median"),
#' shorteded.cancer.names = TRUE, genelimit="none", resolution=600, RowCex=0.8, ColCex=0.8,
#' heatmapMargines=c(15,07), cutoff=NULL, angle.for.heatmap.cancernames=45, heatmap.color = "RdBu",
#' reverse.heatmap.color = TRUE, rewrite.output.list = TRUE, round=TRUE, top.genes = TRUE,
#' validate.genes = TRUE, Use.CancerCode.as.Name = FALSE, simplify.visulization=FALSE,
#' simplifiction.cuttoff=FALSE)
#'
#' @param genes a list that contains at least one gene group
#'
#' @param cancers a character vector or a matrix that containes desired cancer names. The character vector containes standard
#' cancer names that can be found on cbioportal.org, such as \code{Acute Myeloid Leukemia (TCGA, NEJM 2013)}. Alternatively,
#' a matrix can be used if users prefer user-defined cancer names, in which the first column of matrix comprises the standard
#' cancer names while the second column must contain the desired cancer names.
#'
#' @param high.throughput.data.type a character string that is one of the following techniques: 'RNA-seq', 'microRNA-Seq',
#' 'microarray.mRNA', 'microarray.microRNA' or 'methylation'.
#'
#' @param data.presented.as a character vector that containes the statistical precedures users prefer the function to compute.
#' default is \code{c("Frequency.Percentage", "Frequency.Ratio", "Mean.Value", "Median")}. This will tell the function to
#' compute the following:
#' frequency precentage, which is the number of samples having value greather than specific cutoff
#' divided by the total sample size for each cancer study;
#' frequency ratio, which shows the number of selected samples divided by the total number of samples that give the frequency
#' percentage for each cancer-to know selecected and total sample sizes only;
#' Mean Expression, that contains mean value of selected samples for each cancer;
#' Median, which shows the median value of selected samples for each cancer.
#'
#' @param shorteded.cancer.names a logical vector. If the value is set as true, function will try to remove the end part of
#' cancer names aiming to shorten them. The removed segment usually contains the name of scientific group that has conducted
#' the experiment.
#'
#' @param genelimit if large number of genes exists in at least one gene group, this option can be use to limit the number of
#' genes to be shown on hitmap. For instance, \code{genelimit=50} will limit the heatmap to 50 genes showing the most variation.
#' The default value is \code{none}.
#'
#' @param resolution a number. This option can be used to adjust the resolution of the output heatmaps as 'dot per inch'.
#' The defalut value is 600.
#'
#' @param RowCex a number that specifies letter size in heatmap row names.
#'
#' @param ColCex a number that specifies letter size in heatmap column names.
#'
#' @param heatmapMargines a numeric vectors that can be use to set heatmap margins. The default value is
#' \code{heatmapMargines=c(15,07)}.
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
#' or not each gene has a record. If a cancer study doesn't have a record for specific gene, it checks for alternative gene names
#' that cbioportal might use instead of the given gene name.
#'
#' @param Use.CancerCode.as.Name a logical value that tells the function to use sandard abbreviated cancer names instead of complete cancer
#' names, if set to be \code{TRUE}. For example, \code{laml_tcga_pub} is the shortened name for \code{Acute Myeloid Leukemia (TCGA, NEJM 2013)}.
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
#' \code{Median} can be generated. If more than one group of genes is entered, output for each group will be strored in a separate sub-directory.
#'
#' @examples
#' # Creating a list that contains one gene group: 'K.acetyltransferases'
#' genes <- list(K.acetyltransferases = c("CLOCK", "CREBBP", "ELP3", "EP300"))
#'
#' # creating a character vector of study names.
#' cancernames <- c("Acute Myeloid Leukemia (TCGA, Provisional)", "Adrenocortical Carcinoma (TCGA, Provisional)",
#' "Bladder Urothelial Carcinoma (TCGA, Provisional)", "Brain Lower Grade Glioma (TCGA, Provisional)",
#' "Breast Invasive Carcinoma (TCGA, Provisional)")
#'
#' # Running the function to obtain and process the selected data
#' process.multiple.studies(genes, cancernames, "RNA-seq",
#' data.presented.as = c("Frequency.Percentage", "Frequency.Ratio", "Mean.Value"),
#' heatmapMargines=c(15,5), heatmap.color = "redgreen")
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

obtainMultipleStudies <- function(genesList, submissionName, studiesNames, desiredTechnique, cancerCode = FALSE){

  ##########################################################################
  ########## Prerequisites

  # Check genes

  if(!is.list(genesList)){

    stop("'genes' must be entered as a list containing at list one group of genes with descriptive group name for logistical purposes")

  } else if(is.list(genesList) & exists("formerGeneList")){

    oldList <- unname(sort(unlist(formerGeneList)))

    newList <- unname(sort(unlist(genesList)))

    if(identical(newList, oldList)){

       return("--- Function 'obtainMultipleStudies()' was skipped: Data for the requested genes already exist ---")

    }

  }





  # studiesNames names

  if(!is.character(studiesNames)){

    stop("'studiesNames' must be entered a character vector containing at least one cancer name")

  }



  # high.throughput.data.type

  if(is.character(desiredTechnique)){

    if(!(desiredTechnique %in% c("RNA-seq", "microRNA-Seq", "Microarray.mRNA", "Microarray.microRNA", "methylation")) | length(desiredTechnique)!= 1){

      stop("'desiredTechnique' must contain one of the following techniques: 'RNA-seq', 'microRNA-Seq', 'Microarray.mRNA', 'Microarray.microRNA' or

           'methylation'")

    }

  } else {

    stop("'desiredTechnique' must be entered as a character string describing a technique name")

  }




  # Choice of high-throughput data type

  if(desiredTechnique == "RNA-seq"){

    L1.characteristics <- c("Tumor Samples with mRNA data (RNA Seq V2)", "Tumors with mRNA data (RNA Seq V2)", "Tumor Samples with mRNA data (RNA Seq)", "Tumors with mRNA data (RNA Seq)")

    L2.characteristics <- c("mRNA Expression z-Scores (RNA Seq V2 RSEM)", "mRNA Expression z-Scores (RNA Seq RPKM)")

  } else if(desiredTechnique == "microRNA-Seq"){

    L1.characteristics <- c("Tumors with microRNA data (microRNA-Seq)")

    L2.characteristics <- c("microRNA expression Z-scores")

  } else if(desiredTechnique == "microarray.mRNA"){

    L1.characteristics <- c("Tumor Samples with mRNA data (Agilent microarray)", "Tumors with mRNA data (Agilent microarray)", "Tumor Samples with mRNA data (U133 microarray only)", "Tumors with mRNA data", "Tumors with mRNA")

    L2.characteristics <- c("mRNA Expression z-Scores (microarray)", "mRNA Expression z-Scores (U133 microarray only)", "mRNA expression z-scores (Illumina)", "mRNA expression Z-scores (all genes)", "mRNA Expression Z-Scores vs Normals", "mRNA Expression z-Scores (combined microarray)")

  } else if(desiredTechnique == "microarray.microRNA"){

    L1.characteristics <- c("Tumors with microRNA")

    L2.characteristics <- c("mRNA Expression Z-Scores vs Normals", "mRNA Expression z-Scores (combined microarray)")

  } else if(desiredTechnique == "methylation"){

    L1.characteristics <- c("Tumor Samples with methylation data (HM450)", "Tumors with methylation data (HM450)", "Tumor Samples with methylation data (HM27)", "Tumors with methylation data (HM27)", "Tumors with methylation data")

    L2.characteristics <- c("Methylation (HM450)", "Methylation (HM27)", "Methylation")

  } else{

    stop("desiredTechnique field can not be left empety. It should be chosen as 'RNA-seq', 'microRNA-Seq', 'microarray.mRNA', 'microarray.microRNA'or 'methylation'")

  }



  ##########################################################################
  ########## Set the function ready to work

  # Set cgdsr

  mycgds = CGDS("http://www.cbioportal.org/")



  # Creating a vector for cancer names and subsequent list sebset name

  if(!is.matrix(studiesNames)){

    studiesNames <- studiesNames[order(studiesNames)]

    studiesNamesMatrix <- studiesNames



    if(cancerCode == TRUE){

      groupNames <- getCancerStudies(mycgds)[which(getCancerStudies(mycgds)[,2] %in% as.character(studiesNames)),1]

    } else if(cancerCode == FALSE){

      groupNames <- as.character(studiesNames)

    }


  } else if(is.matrix(studiesNames)){

    studiesNamesMatrix <- studiesNames[order(studiesNames[,1]),]

    studiesNames <- studiesNamesMatrix[,2]

    groupNames <- studiesNames

  }


  ##########################################################################
  ########## Core segment

  # Report

  print(paste("***", "Obtaining the requested data for", submissionName, "***", sep = " "))



  # create progress bar

  obtainMultipleStudiesProgressBar <- txtProgressBar(min = 0, max = length(studiesNames), style = 3)



  # Create parent list for storing final results in the global environment

  rawList <- list()

  # Creating child lists

  for(nname in 1:length(genesList)){

       rawList[[nname]] <- list()

       names(rawList)[nname] <- names(genesList)[nname]

  }



  ## Getting the required gene expresssion profile ...

  # 'for' control structure for obtaining data and calculating the requested parameters

  for(c in 1:length(studiesNames)){

    # Determining name for list subset of study name

    groupName <- groupNames[c]

    # Correcting possible errors of list names

    groupName <- gsub(groupName, pattern = "\\+ ", replacement = " possitive ", ignore.case = TRUE)

    groupName <- gsub(groupName, pattern = "\\- ", replacement = " negative ", ignore.case = TRUE)

    groupName <- gsub(groupName, pattern = " ", replacement = "_", ignore.case = TRUE)





    ##--## print(as.character(studiesNames[c]))

    mycancerstudy = getCancerStudies(mycgds)[which(getCancerStudies(mycgds)[,2]==as.character(studiesNames[c])),1]



    # Finding the first characteristics of data in the cancer

    f.condition <- (getCaseLists(mycgds,mycancerstudy)[,2])[getCaseLists(mycgds,mycancerstudy)[,2] %in% L1.characteristics]

    mycaselist = getCaseLists(mycgds,mycancerstudy)[which(getCaseLists(mycgds,mycancerstudy)[,2] == if(length(f.condition) >= 1){

      f.condition[1]

    } else if(length(f.condition) == 0){

      stop(paste(studiesNames[c], "lacks", desiredTechnique, "data!", sep=" "))

    }) ,1]



    # Finding the second characteristics of data in the cancer

    s.condition <- (getGeneticProfiles(mycgds,mycancerstudy)[,2])[getGeneticProfiles(mycgds,mycancerstudy)[,2] %in% L2.characteristics]

    mygeneticprofile = getGeneticProfiles(mycgds,mycancerstudy)[which(getGeneticProfiles(mycgds,mycancerstudy)[,2] == if(length(s.condition) >= 1){

      s.condition[1]

    } else if (length(s.condition) == 0){

      stop(paste(studiesNames[i], "doesn't have an appropriate 'level 2' condition for", desiredTechnique, "data!", sep=" "))

    }) ,1]



    # obtaining data for every genegroup

    for(group in 1:length(genesList)){

      # Chose one group of genes

      genesNames <- genesList[[group]]



      # Obtaining Expression x-scores fore the requested genes

      # Assaign data to specific list member

      rawList[[group]][[c]] <- data.matrix(getProfileData(mycgds,genesNames[order(genesNames)],mygeneticprofile,mycaselist))

      names(rawList[[group]])[c] <- groupName

    }

    # Update progressbar

    setTxtProgressBar(obtainMultipleStudiesProgressBar, c)


  }

  # Closing progress bar

  close(obtainMultipleStudiesProgressBar)

  # Store name of the genes

  formerGeneList <<-  genesList

  # Export the obtained data as list

  assign(paste("ObM", ":", submissionName, sep = ""), rawList, envir = globalenv())

}
