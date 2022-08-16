#This script makes a table of random effects meta analysis estimates 
#and heterogeneity estimates for full and subsamples


#First make a subset for each subsample
hettablenames <- c("Control_gender" ,          "Control_race"   ,          "Control_income" ,          "Control_education",       
                   "Homicide"      ,           "Violent"       ,          "Non_Violent"    ,          "Both"             ,      
                   "Area_dummy"    ,           "OLS"           ,           "ML"            ,           "Odds_Ratio"      ,        
                   "Panel_dummy"   ,           "Endogeneity"   ,          "North_America"  ,         "Europe"            ,      
                         "Direct_Lead_Measure" ,     "Rep_estimate" )
  
  
  
  
p <- lapply(hettablenames, function(x) subset(df, df[,x] == 1))


#make empty data frame for table
a <- data.frame(Variable =  rep(0, length(p)),RE_Estimate = rep(0, length(p)), SE = rep(0, length(p)), 
                Tau2 = rep(0, length(p)) , I2 = rep(0, length(p)), H2 = rep(0, length(p),
Studies = rep(0, length(p)), Estimates = rep(0, length(p)) ))


#empty list for loop
x <- list()         


#Loop which gets RE and heterogeneity estimates   
for (i in 1:length(p)) {

x <- rma.uni(p[[i]] $Standardised_effect_size, p[[i]]$var, method = "DL")

a$Studies[i] <- length(unique(p[[i]] $Study_ID))
a$Estimates[i] <- length(p[[i]]$UID)
a$Variable[i] <- paste(colnames(df[(20+i)]), " = TRUE")
a$RE_Estimate[i] <- conv(x$beta [1])
a$SE[i] <- x$se
a$Tau2[i] <- x$tau2
a$I2[i] <- x$I2
a$H2[i] <- x$H2 
}



#Now we do same thing but for opposite subsamples.  
#i.e instead of a subsample which controls for income, we have a subsample which does not control for income

d <- lapply(hettablenames, function(x) subset(df, df[,x] == 0))



b <- data.frame(Variable =  rep(0, length(d)),RE_Estimate = rep(0, length(d)), SE = rep(0, length(d)),
                Tau2 = rep(0, length(d)) ,
                I2 = rep(0, length(d)), H2 = rep(0, length(d)
                                  , Studies = rep(0, length(d)), Estimates = rep(0, length(d)) ))
x <- list()         

for (i in 1:length(d)) {
  
  x <- rma.uni(d[[i]] $Standardised_effect_size, d[[i]]$var, method = "DL")
  
  
 b$Studies[i] <- length(unique(d[[i]] $Study_ID))
 b$Estimates[i] <- length(d[[i]]$UID)
 b$Variable[i] <- paste(colnames(df[(20+i)]), " = FALSE")
 b$RE_Estimate[i] <-conv( x$beta [1])
 b$SE[i] <- x$se
 b$Tau2[i] <- x$tau2
 b$I2[i] <- x$I2
 b$H2[i] <- x$H2 

}


#==================================================
#finally we get random effects and heterogenity estimates for full sample


#bind them together

a  <- rbind(a,b)




a[nrow(a) + 1, ] <- list("Full Sample", conv(MET2$beta)[1], MET2$se, MET2$tau2, MET2$I2, 
                      MET2$H2, length(unique(df$Study_ID)), nrow(df))





#============================================================
#now get the elasticity estimates

a[nrow(a) + 1, ] <- list("Elasticity Sample*", ELASMET2$beta[1], ELASMET2$se,ELASMET2$tau2, ELASMET2$I2, 
                         ELASMET2$H2, dim(table(dfE$Short_name)), nrow(dfE))





#and now the elasticity and endogeneity subsample
a[nrow(a) + 1, ] <- list("Elasticity Sample (Addressing Endogeneity)*",ELASENDOMET2$beta[1], ELASENDOMET2$se ,ELASENDOMET2$tau2, ELASENDOMET2$I2, 
                         ELASENDOMET2$H2, dim(table(dfE$Short_name[dfE$Endogeneity==1])), nrow(dfE[dfE$Endogeneity==1,]))


#make the table
write.csv(a, "heterogeneitytable1.csv", row.names = FALSE)

