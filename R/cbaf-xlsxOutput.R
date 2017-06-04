
as <- function(x){
  x+3
}
ad <- function(y, as){
  as+y
}

 ghj <- function(x=4, y=5){
   c <- as(x=x)
   ad (c, y = y)
 }





##########################################################################
##########################################################################
#################### Exporting Data as a excel file ######################
##########################################################################
##########################################################################




xlsxOutput <- function(pocessedDatabase){

  # Break

  List.to.Go <- get(paste(names(genes)[[g]], "Data", "List", sep = "."))

  if(is.list(pocessedDatabase)){

    tabsNumber <- length(pocessedDatabase)

  } else{

    tabsNumber <- 1

  }


  # Create progress bar

  xlsxProgressBar <- txtProgressBar(min = 0, max = tabsNumber, style = 3)

  ## 'for' to prepare ouputs

  for(segment in 1:tabsNumber){

    ## Subsetting matrixes from List by 'for' control structure

    if(is.list(pocessedDatabase)){

      Temporary.source <- pocessedDatabase[[segment]]

      Temporary.source.name <- names(pocessedDatabase)[[segment]]

    }

  ## Exporting the expression profile

  Exp.matrix2 <- Temporary.source

  if(is.matrix(cancers.matrix)){

    rownames(Exp.matrix2) <- cancers.matrix[,2]

  }

  # Rounding values

  if(round==TRUE & Temporary.source.name %in% c("Frequency.Percentage","Mean.Value" ,"Median")){

    Exp.matrix2 <- round(Exp.matrix2, digits=2)

  }

  # Copy data to a variable to modify

  Save.matrix <- Exp.matrix2

  if(!(Temporary.source.name %in% c("Top.Genes.of.Frequency.Percentage","Top.Genes.of.Mean.Value" ,"Top.Genes.of.Median"))){

    # Replace 'NaN' values with character string

    Save.matrix[is.nan(Save.matrix)] <- "NaN"

    # Replace 'NA' values with character string

    Save.matrix[is.na(Save.matrix)] <- "NA"

  }

  if(Temporary.source.name %in% c("Top.Genes.of.Frequency.Percentage","Top.Genes.of.Mean.Value" ,"Top.Genes.of.Median")){

    colnames(Save.matrix)[1:2] <- c("Gene with highest value", if(Temporary.source.name == "Top.Genes.of.Frequency.Percentage"){"Frequency Percentage"}

                                    else if(Temporary.source.name == "Top.Genes.of.Mean.Value"){"Mean Expression"}

                                    else if(Temporary.source.name == "Top.Genes.of.Median"){"Median"})

  }



  # Saving the expression profile

  write.xlsx(Save.matrix, file = paste(gsub(x = names(genes)[g], pattern = "\\.", replacement = "-"), " (",cutoff.phrase, "=", cutoff, ")" , ".xlsx", sep = ""), sheetName = gsub(x = Temporary.source.name, pattern = "\\.", replacement = " "), append = TRUE)




  # Update progressbar

  setTxtProgressBar(xlsxProgressBar, segment)

}

# Close progress bar

close(xlsxProgressBar)

# Inform the user about the directory

# Print information

print(paste("Requested excel files were saved in", getwd(), sep = " "))

}
