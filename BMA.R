
#========================================================================
#Bayesian model averaging over many different pub bias models
#This ends up putting a lot of weight on andrews-kasy selection type models, so is also negative.


# Parameters for in-sample estimation 

burn_       = 0.5e5
iter_       = 1e5
gprior      = "UIP"
modelprior  = "uniform"
order_PIP   = F
separator   = ";"



fit <- RoBMA( y = df$PCC , se = df$Standarised_se, study_names = df$Short_name
            , seed = 56189 #from Random.org  
        ,parallel = TRUE, burnin = burn_  )
#



#BMA model for representative estimates

REPfit <- RoBMA( y = df$PCC[df$Rep_estimate ==1] , se = df$Standarised_se[df$Rep_estimate ==1], study_names = df$Short_name[df$Rep_estimate ==1]
            , seed = 1737,  burnin = burn_ )



#BMA model for elasticity sample

ELASfit <- RoBMA( y = dfE$Elasticity , se = dfE$Elasticity_SE, study_names = dfE$Short_name
           , seed = 20053,  burnin = burn_  )



ELASEndofit <- RoBMA( y = dfE1$Elasticity , se = dfE1$Elasticity_SE, study_names = dfE1$Short_name
                      , seed = 20053 , burnin = burn_ )


#================================================================

#Bayesian meta regression model averaging code





# Estimation #

BMA_nw  = bms(data = df, X.data =  Standardised_T ~ Precision +
                WControl_gender + WControl_race + WControl_income + 
                WControl_education +
              WHomicide + WViolent + WNon_Violent + 
               WArea_dummy  + WOLS + WML + WOdds_Ratio + WPanel_dummy   + 
                WEndogeneity  + WNorth_America + WEurope +
                + WDirect_Lead_Measure +  WPubYear +  
                WZCovariates +   WZsample ,
              burn=burn_,iter=iter_, g=gprior, mprior=modelprior, nmodel=5000, mcmc="bd", user.int=FALSE)



BMA_nw_res<-coef(BMA_nw,std.coefs=F, order.by.pip = order_PIP, exact=T,include.constant = T)
write.csv(BMA_nw_res,"BMAPCC.csv")

image(BMA_nw, cex=0.5)
summary(BMA_nw)
plot(BMA_nw)


#get average effect using sample averages
X <- c(         "WControl_gender" ,        
               "WControl_race"  ,          "WControl_income" ,         "WControl_education",        "WHomicide" ,              
               "WViolent" ,                "WNon_Violent",             "WArea_dummy" ,            
               "WOLS"   ,                  "WML"    ,                  "WOdds_Ratio"  ,            "WPanel_dummy" ,           
               "WEndogeneity"  ,           "WNorth_America" ,          "WEurope"   ,                      
               "WDirect_Lead_Measure",     "WPubYear",                  "WZCovariates",             "WZsample" )

mn <- sapply(df[,substring(X, first=2)],mean, na.rm = TRUE)


mn <- matrix(c(1,mn))

BMAAverage <- t(mn) %*% matrix(BMA_nw_res[-nrow(BMA_nw_res),2])


# Save 

save.image("results_excess.Rdata")

#==============================================
#Now same for elasticity


data = dfE


BMA_nw2  = bms(data = dfE, X.data =  T_Elast ~ Elas_Precision +
                WControl_gender + WControl_race + WControl_income + 
                WControl_education +
                WHomicide + WViolent + WNon_Violent + 
                WArea_dummy  + WOLS + WML  + WPanel_dummy   + 
                WEndogeneity  + WNorth_America + 
                + WDirect_Lead_Measure +  WPubYear +  
                WZCovariates +   WZsample ,
              burn=burn_,iter=iter_, g=gprior, mprior=modelprior, nmodel=5000, mcmc="bd", user.int=FALSE)



BMA_nw_res2<-coef(BMA_nw2,std.coefs=F, order.by.pip = order_PIP, exact=T,include.constant = T)
write.csv(BMA_nw_res2,"BMAelas.csv")

X <- c(              "WControl_gender" ,        
             "WControl_race"  ,          "WControl_income" ,         "WControl_education",        "WHomicide" ,              
             "WViolent" ,                "WNon_Violent",             "WArea_dummy" ,            
             "WOLS"   ,                  "WML"    ,                  "WPanel_dummy" ,           
             "WEndogeneity"  ,           "WNorth_America" ,                   
             "WDirect_Lead_Measure",     "WPubYear",                 "WZCovariates",             "WZsample" )



mn <- sapply(dfE[,substring(X, first=2)],mean, na.rm = TRUE)


mn <- matrix(c(1,mn))

BMAAverage2 <- t(mn) %*% matrix(BMA_nw_res2[-nrow(BMA_nw_res2),2])


