








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



# Saving data list to enviroment with an informative name

assign(paste(names(genes)[[g]], ".", "Data.List",  "_", "for", "_", Shortened.cancername, sep=""), Data.list, envir = globalenv())

}
