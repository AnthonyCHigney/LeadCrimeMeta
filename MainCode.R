
#Meta-analysis code



#WARNING - The full code takes about 40 hours on my machine 

#You need Rtools for some packages https://cran.r-project.org/bin/windows/Rtools/ (for windows)
#for Bayesian model averaging packages you need JAGS https://sourceforge.net/projects/mcmc-jags/


#=================================================================================
#Preliminaries

#Import the dataset
df <- read.csv("Lead_Meta_Data14.csv", stringsAsFactors=FALSE)

#make vector of needed packages to be used
pack <- c( "Cairo",  "tidyverse", "ggplot2", "gridExtra", "reshape2" , #general packages
           "metafor","forestplot","metasens", "meta" , #standard meta analysis packages
           "lme4", "mice", "miceadds",   "fixest",  "sandwich",  "AER" ,    "corrplot",           #additional estimation packages
           "foreach", "parallel", "doParallel", "feather", #for parallel computing and feather for quick saving
           "BMS"    , "RoBMA"                                     ) #Bayesian Model averaging packages  
                                                         
     



#Install or just use packages as required
#lapply(pack, FUN = install.packages)
lapply(pack, FUN = require, character.only = TRUE)

#set number of cores to use for parallel computing. I use total - 2
#cores variable is "no_cores"

no_cores <- detectCores() - 2


#================================================================================
#functions

#make a function to convert standardized effects back to PCCs
#this is used with the common and random effects averages

conv <- function(x) {
  x <- (exp(2*x)-1)/(exp(2*x)+1)
  return(x)}



#make a function to weigh covariates for WLS meta-regressions


weigh <- function(x) {
  x/df$Standarised_se
}


#The function we use to run all the regressions
source("Formulaloop.R")

#separate one for elasticities
source("Loopelast.R")

#Functions for Andrews-Kasy method
source("AndrewsKasyFunctions.R")

#==================================================================================
#Data prep

#First Normalise the three continuous variables, otherwise REML has trouble converging in some estimations
scaled <- as.data.frame(scale(df$Covariates))
df$ZCovariates <- scaled$V1

scaled <- as.data.frame(scale(df$Sample_size))
df$Zsample <- scaled$V1

#make a variable for publication year (normalised to aid REML convergence)
df$PubYear <- scale(df$Publication_Year) 


#Make new columns for the weighted covariates
#first make a name vector for the new columns
orignames <- c("Control_gender" ,          "Control_race"   ,          "Control_income" ,          "Control_education",       
                "Homicide"      ,           "Violent"       ,          "Non_Violent"    ,          "Both"             ,      
                "Area_dummy"    ,           "OLS"           ,           "ML"            ,           "Odds_Ratio"      ,        
                "Panel_dummy"   ,           "Endogeneity"   ,          "North_America"  ,         "Europe"            ,      
                         "Direct_Lead_Measure" ,     "Rep_estimate" ,
               "ZCovariates"    ,          "Zsample" , "PubYear")

namevector <- paste0("W", orignames)


#Weigh the covariates by the PCC standard errors using namevector to name the columns
df[,namevector] <- lapply(df[orignames] , FUN = weigh)






#==================================================================================
#Now run the code
#Elasticity and PCC code is mostly separate. This is mainly a legacy thing as elasticity was added after, but also because they use different main variables
#So order of script running does matter.

#Do the elasticity subsample workings 
source("elast.R")

#This makes a forest plot of the studies and a table of their mean, median effect sizes etc
source("StudySummary.R")

#this makes a table to show the heterogeneity measures by different sub samples
source("HeterogeneityTable.R")

#here is the main estimates of pub bias and effect after accounting for pub bias
source("pubbias.R")

#makes a table of the moderator variables (means and SDs)
source("metaregression.R")



#do Bayesian model averaging over the meta-regression covariates
source("BMA.R")
#==================================================================================
##MAKE DENSITY PLOTS
#These are the main density plots
#how many chunks to split code up. Depends on RAM
chunks = 300
#WARNING THIS WILL TAKE A LONG TIME!


#data <- df is nonsense but i can't be bothered changing code now
data <- df

#the columns used for the regressions. The order is important in the function!
columns <- c("Standardised_T" ,           "Precision",                "Study_ID",                  "WControl_gender" ,        
              "WControl_race"  ,          "WControl_income" ,         "WControl_education",        "WHomicide" ,              
              "WViolent" ,                "WNon_Violent",             "WArea_dummy" ,            
              "WOLS"   ,                  "WML"    ,                  "WOdds_Ratio"  ,            "WPanel_dummy" ,           
              "WEndogeneity"  ,           "WNorth_America" ,          "WEurope"   ,                        
              "WDirect_Lead_Measure",     "WPubYear",                  "WZCovariates",             "WZsample" )

folderpath = "./TMP/tmp"


#this is the "ideal" specification as we use in the paper, this is why order important
idealvector = c(1,1,1,1,0,0,0,0,0,0,0,1,1,1,0,1)

#now function made in FormulaLoop.R is called
tmp <- loopformula(data,columns, folderpath, idealvector,295,chunks)


#write output as a csv
write.csv(tmp, "tmp.csv")




#Make dataset for endogeneity subsample
data = df[df$Endogeneity ==1,]

#use only columns with some variation as covariates
#also have to be aware that if base case has no variation, you need to omit the new dummy base case
#otherwise collinearity induced by only selecting a subsample
columns = c("Standardised_T" ,           "Precision",                "Study_ID",                  "WControl_gender" ,        
            "WControl_race"  ,          "WControl_income" ,        "WHomicide" ,              
            "WViolent" ,                "WNon_Violent",            "WArea_dummy" ,          
            "WOLS"   ,                  "WPanel_dummy" ,           "WPubYear"  ,
            "WZCovariates",             "WZsample" )

#where you store feather files before they get bundled together (you may not need to use chunks)
folderpath = "./TMPQUAS/tmp"

#the vector of dummy variables you use to construct the effect size at end for the "ideal" case
idealvector = c(1,1,1,0,0,0,0,0,1)

#Run Endogeneity loop
tmpquas <- loopformula(data,columns, folderpath, idealvector,1,chunks)

#write as a csv
write.csv(tmpquas, "tmpquas.csv")


#now for non-endogeneity

data = df[df$Endogeneity ==0,]

columns = c("Standardised_T" ,           "Precision",                "Study_ID",                  "WControl_gender" ,        
            "WControl_race"  ,          "WControl_income" ,         "WControl_education",        "WHomicide" ,              
            "WViolent" ,                "WNon_Violent",             "WArea_dummy" ,            
            "WOLS"   ,                  "WML"    ,                  "WOdds_Ratio"  ,            "WPanel_dummy" ,           
            "WNorth_America" ,                
            "WDirect_Lead_Measure",     "WPubYear",   "WZCovariates",             "WZsample" )

folderpath = "./TMPNONQUAS/tmp"

idealvector = c(1,1,1,1,0,0,0,0,0,0,0,1,1,1)


tmpnonquas <- loopformula(data,columns, folderpath, idealvector,1,chunks)


write.csv(tmpnonquas, "tmpnonquas.csv")


#now for homicide

data = df[df$Homicide ==1,]

columns = c("Standardised_T" ,           "Precision",                "Study_ID",       
            "WControl_race"  ,          "WControl_income" ,                    
            "WOLS"   ,                   "WPanel_dummy" ,           
            "WEndogeneity"  ,        
            "WPubYear",                 "WZCovariates",             "WZsample" )

folderpath = "./TMPHOMICIDE/tmp"

idealvector = c(1,1,0,1,1)


tmphomicide <- loopformula(data,columns, folderpath, idealvector, 1, chunks)


write.csv(tmphomicide, "tmphomicide.csv")




#now for violent crime

data = df[df$Violent ==1,]

columns = c("Standardised_T" ,           "Precision",                "Study_ID",                "WControl_gender" ,        
            "WControl_race"  ,          "WControl_income" ,         "WControl_education",       "WArea_dummy" ,            
            "WOLS"   ,                  "WML"    ,                   "WPanel_dummy" ,           
            "WEndogeneity"  ,           "WNorth_America" ,            
            "WDirect_Lead_Measure",     "WPubYear",                "WZCovariates",             "WZsample" )

folderpath = "./TMPVIOLENT/tmp"

idealvector = c(1,1,1,1,0,0,0,1,1,1,1)


tmpviolent <- loopformula(data,columns, folderpath, idealvector,1, chunks)


write.csv(tmpviolent, "tmpviolent.csv")





#now for non-violent crime

data = df[df$Non_Violent ==1,]


columns = c("Standardised_T" ,           "Precision",                "Study_ID",                "WControl_gender" ,        
            "WControl_race"  ,          "WControl_income" ,         "WControl_education",       "WArea_dummy" ,            
            "WOLS"   ,                  "WML"    ,                  "WOdds_Ratio"  ,            "WPanel_dummy" ,           
            "WEndogeneity"  ,           "WNorth_America"   ,        
            "WDirect_Lead_Measure",     "WPubYear",                 "WZCovariates",             "WZsample" )

folderpath = "./TMPNONVIOLENT/tmp"

idealvector = c(1,1,1,1,0,0,0,0,1,1,1,1)


tmpnonviolent <- loopformula(data,columns, folderpath, idealvector,1,chunks)


write.csv(tmpnonviolent, "tmpnonviolent.csv")



#now do areas

data = df[df$Area_dummy ==1,]

columns <- c("Standardised_T" ,           "Precision",                "Study_ID",
              "WControl_race"  ,          "WControl_income" ,         "WControl_education",        "WHomicide" ,              
              "WViolent" ,                "WNon_Violent",                  
              "WOLS"   ,                  "WML"    ,                  "WPanel_dummy" ,           
              "WEndogeneity"  ,           "WNorth_America" ,                 
                 "WPubYear",                 "WZCovariates",             "WZsample" )

folderpath = "./TMPAREA/tmp"

idealvector = c(1,1,1,0,0,0,0,0,1,1,1)


tmparea <- loopformula(data,columns, folderpath, idealvector,1,chunks)


write.csv(tmparea, "tmparea.csv")




#now do individuals

data = df[df$Area_dummy ==0,]

 columns <- c("Standardised_T" ,           "Precision",                "Study_ID",                  "WControl_gender" ,        
                       "WControl_race"  ,          "WControl_income" ,         "WControl_education",              
                       "WViolent" ,                "WNon_Violent",                      
                       "WOLS"   ,                  "WML"    ,                  "WOdds_Ratio"  ,            "WPanel_dummy" ,           
                       "WEndogeneity"  ,           "WNorth_America" ,                 
                       "WDirect_Lead_Measure",     "WPubYear",                 "WZCovariates",             "WZsample" )

folderpath = "./TMPIND/tmp"

idealvector = c(1,1,1,1,0,0,0,0,0,1,1,1,1)


tmpind <- loopformula(data,columns, folderpath, idealvector,1,chunks)


write.csv(tmpind, "tmpind.csv")


#elasticities
data = dfE

columns <- c("T_Elast" ,           "Elas_Precision",                "Study_ID",                  "WControl_gender" ,        
              "WControl_race"  ,          "WControl_income" ,         "WControl_education",        "WHomicide" ,              
              "WViolent" ,                "WNon_Violent",             "WArea_dummy" ,            
              "WOLS"   ,                  "WML"    ,                  "WPanel_dummy" ,           
              "WEndogeneity"  ,           "WNorth_America" ,                   
              "WDirect_Lead_Measure",     "WPubYear",                 "WZCovariates",             "WZsample" )

folderpath = "./TMPELAS/tmp"

idealvector = c(1,1,1,1,0,0,0,0,0,0,1,1,1,1)


tmpElas <- loopformulaELAS(data,columns, folderpath, idealvector,1,chunks)


write.csv(tmpElas, "tmpELAS.csv")




#elasticities with addressing endogeneity

data = dfE1

columns = c("T_Elast" ,           "Elas_Precision",                  "Study_ID",                  "WControl_gender" ,        
            "WControl_race"  ,          "WControl_income" ,        "WHomicide" ,              
            "WViolent" ,                "WNon_Violent",            "WArea_dummy" ,          
                   "WPubYear"  ,
            "WZCovariates",             "WZsample" )

folderpath = "./TMPELASENDO/tmp"

idealvector = c(1,1,1,0,0,0,0)


tmpElasEndo <- loopformulaELAS(data,columns, folderpath, idealvector,1,chunks)


write.csv(tmpElasEndo, "tmpElasEndo.csv")





#===========================================================
#Now plot the densities

source("histcharts.R")

