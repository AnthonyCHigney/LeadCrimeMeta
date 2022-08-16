#--------------------------------------------------------------------------------------------------------------------------------
#Meta-regression with moderators code

#make a name vector for sunning desc stats

orignames2 <- c("Control_gender" ,          "Control_race"   ,          "Control_income" ,          "Control_education",       
               "Homicide"      ,           "Violent"       ,          "Non_Violent"    ,          "Both"             ,      
               "Area_dummy"    ,           "OLS"           ,           "ML"            ,           "Odds_Ratio"      ,        
               "Panel_dummy"   ,           "Endogeneity"   ,          "North_America"  ,         "Europe"            ,      
               "Direct_Lead_Measure" ,     "Rep_estimate" ,
               "Covariates"    ,          "Sample_size" , "Publication_Year")

namevector <- paste0("W", orignames)

#Get means and standard deviations of unweighted covariates for summary table
means <- as.data.frame(lapply(df[orignames2], FUN = mean))

SD <- as.data.frame(lapply(df[orignames2], FUN = sd))

medians <- as.data.frame(lapply(df[orignames2], FUN = median))


tableX <- cbind(names(df[orignames2]), t(means), t(medians), t(SD))


colnames(tableX) <- c("Name", "Mean","Median", "Standard Deviation")
write.csv(tableX, "Moderators.csv", row.names = FALSE)

#-------------------------------------------------------------------------------------
#Carry out meta-regression with all covariates


FullModel <- lmer(data = df, formula =  Standardised_T ~ Precision + (1 | Study_ID)
                                   
                            + WHomicide + WViolent + WNon_Violent + 
                              WOdds_Ratio + WML + WPanel_dummy   + 
                              WEndogeneity + WArea_dummy + WZCovariates + 
                              WZsample   + WNorth_America + WEurope + 
                              WControl_gender + WControl_race + WControl_income + 
                              WControl_education  + WOLS  + WDirect_Lead_Measure + WPubYear, REML = TRUE,
                            start = 0, control = lme4:: lmerControl(calc.derivs = FALSE, optimizer ="Nelder_Mead"))

write.csv(summary(FullModel)$coefficients, "fullModel.csv")


#Now just with representative estimates




RepModel <- lm(data = df[df$Rep_estimate == 1,], formula =  Standardised_T ~ Precision +
                 WHomicide + WViolent + WNon_Violent + 
                 WOdds_Ratio + WML + WPanel_dummy   + 
                 WEndogeneity + WArea_dummy + WZCovariates + 
                 WZsample   + WNorth_America + WEurope + 
                 WControl_gender + WControl_race + WControl_income + 
                 WControl_education  + WOLS  + WDirect_Lead_Measure + WPubYear )




