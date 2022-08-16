
##MAKING loopformula function to run all combinations of covariates in regressions (REML)


loopformula <- function(data,columns, folderpath, idealvector,start,chunks) {
  
  
  
  
  
  
  
  
  
  #make df without changing but only including specified columns.
  
  df1 <- data[,c(columns)]
  

  #make a vector of variable names to use in each regression formula
  #drop variables that are in every formula (standardised T, precision, and study ID)
  X <- names(df1[,-c(1,2,3)])
  
  
  
  #make list of formulas
  
  out <- unlist(lapply(1:length(X), function(n) combn(X, n, FUN=function(row) paste0("Standardised_T ~  Precision + (1 | Study_ID) + ", paste0(row, collapse = "+")))))
  
  #add in FAT-PET with no covariates 
  out <- c("Standardised_T ~  Precision + (1 | Study_ID)", out)
  
  

  #=============================================================
  #THIS IS MAIN LOOP
  
  #Set up parallel computing
  no_cores <- detectCores() - 2 
  registerDoParallel(cores=no_cores) 
  cl <- makeCluster(no_cores) 
  
  
  #splitting function for chunks, you may not need
  
  n <- split(out, cut(seq(length(out)), chunks))
  
  mods <- list()
  
  for (N in start:length(n)) {
    
    mods <- parLapply(cl, n[[N]], fun =  lme4::lmer, data=df1, REML = TRUE, start = 0,
                      control = lme4:: lmerControl(calc.derivs = FALSE, optimizer = "nloptwrap", check.scaleX = "ignore"))
    
    
    
    p <- length(mods)
    
    
    #Make empty data frame to fill out with coefficients and standard errors
    a <- as.data.frame(matrix(nrow = length(mods), ncol = (length(X) + 2)*2 + 1))
    names(a) <- c("PET", "FAT", substring(X, first=2), 
                  "PET_SE", "FAT_SE",  paste0(substring(X, first=2),"_SE"), "frml")
    
    
    for (i in 1:p) {
      
      #Start off by getting the FAT and PET with standard errors
      #These are present in all models so it is easy
      
      a$FAT[i] <- mods[[i]] @beta[1]
      a$FAT_SE[i] <- sqrt(vcov(mods[[i]]) [1,1])
      a$PET[i] <- mods[[i]] @beta[2]
      a$PET_SE[i] <- sqrt(vcov(mods[[i]]) [2,2])
      
      
      
      
      
      #Now call a loop to find in which order a variable appears in the regression equation (or not)
      #And then put the reference coefficient and standard error into the data frame
      
      for (j in 1:length(X)) {
        
        
        #These lines find the coefficient order in the formula (if it is in this specification) 
        #And then put coefficient into data frame
        
        varname <- paste0(substring(X[j],first = 2))
        
        
        a[[varname]][i] <- ifelse(identical(mods[[i]]@beta[grep(X[j], x = labels(terms(mods[[i]]))) +1], numeric(0)) == FALSE, 
                                  mods[[i]]@beta[grep(X[j], x = labels(terms(mods[[i]]))) +1],NA)
        
        
        #now we do the same for standard error
        varname <- paste0(substring(X[j],first = 2), "_SE") 
        
        #this part checks if coefficient standard error exists in each specification
        #if not we skip to next 
        
        skip_to_next <- FALSE
        
        tryCatch(sqrt(vcov(mods[[i]]) [X[j], X[j]]), error = function(e) { skip_to_next <<- TRUE})
        
        if(skip_to_next) { next }  
        
        
        a[[varname]][i] <-  sqrt(vcov(mods[[i]]) [X[j], X[j]])
        
        
        
      }
      
      #end of loops+
      
    }
    
  #write feather for this chunk before moving to next
    path <- paste0(folderpath, N)
    write_feather(a, path)
  }
  #stop parallel cluster
  stopCluster(cl)  
  
  
  
  #================================================================
  #bind chunks back together
  
  tmp1 <- data.frame()
  
  for (N in 1:length(n)) {
    
    path <- paste0(folderpath, N)
    tmp1 <- rbind(tmp1, read_feather(path))
    
  }
  
  #add in the formula for each specification
  
  tmp1$frml <- out
  
  
  #add in names
  tmp1 <- cbind(tmp1, sapply(X, FUN = grepl, x = tmp1$frml))
  
  
  #=================================================
  #This gets the sample average and then combines them with coefficients to get estimates of average effect size
  
  
  
  #get sample averages
  mn <- sapply(data[,substring(X, first=2)],mean, na.rm = TRUE)
  
  #Now combine with 0 vector that we fill out after
  a <- as.data.frame(rbind(mn,0), stringsAsFactors = F)
  
  #NOW WE PUT IN "IDEAL" SPECIFICATION (max or mean sample and covariates as you like)
  


  a[2,] <- c(idealvector,mean(data$PubYear), mean(data$Zsample), mean(data$ZCovariates))

  
  
  #add in a column of 1's for PET
  a <- cbind(PET = 1,a)
  #names for each
  rownames(a) <- c( "average","idealNAAv")
  
  
  #Now we make them into matrices and just use columns that we need
  #also turn NA to 0 for matrix multiplication
  amat <- as.matrix(a[,c("PET", substring(X, first=2))])
  
  amat[is.na(amat)] <- 0
  
  bmat <- as.matrix(tmp1[, c("PET", substring(X, first=2))])
  
  
  bmat[is.na(bmat)] <- 0
  
  #combine with previous data frame
  tmp1 <- cbind(tmp1, bmat %*% t(amat))
  
  
  

  
  

  
  
}