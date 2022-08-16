
#--------------------------------------------------------------------------------------------------------------------------------
#Testing for publication bias


#=================================================================

#Funnel plot

#first plot the funnel plot to examine for bias

#we need to make a column showing significance for the plot (at 95% level of 2-tail)
df$sig <- ifelse(abs(df$PCC/df$Standarised_se) > 1.96,"Significant","Not Significant")


funnelPCC <- ggplot(df, aes(x = PCC, y = Precision, group = sig)) + 
  geom_point(cex = 4, aes(shape = sig, col = sig, fill = sig, group = sig), alpha=min(.5,max(40/5000,.3))) +
  
  theme_bw()   + theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) + 
  scale_shape_manual( values = c(25,21)) + scale_color_manual(values = c("#B95465", "#235D96")) + 
  scale_fill_manual(values = c("#B95465", "#235D96")) +
  ylab("Precision \n1/SE")  +  
  theme(panel.grid = element_blank()) +
  theme(text = element_text(size = 15), legend.title = element_blank(), legend.position = "top") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  scale_x_continuous(expand = c(0,0),limits = c(-1,1)) + 
  geom_vline(xintercept = 0) 


tiff(file="funnel1.tiff",width=1024,height=768,res=100, type = "cairo")


  
funnelPCC
  dev.off()
  
  


#=============================================================


  

#Visual inspection done, but now we test more formally
#First with linear methods
  
  
  
#First is standard FAT-PET Egger test. regression of Effect = B1 + B2*SE + e.  
#To account for heteroskedasticity we do WLS, with weights of 1/SE, so it becomes Effect/SE = B1*1/SE + B2 + v.
#The coefficient B2 is the Funnel Asymmetry Test. If it is not equal to 0, and the t-test is significant,
#then we reject the null of no publication bias. The sign indicates the direction of bias. 
#B1 is the Precision Effect Test.  This is the estimated "true" effect size once we have adjusted for the bias. 


FATPET <- lm.cluster(formula = Standardised_T ~ Precision, cluster = df$Study_ID, data = df)




#Next the FAT-PEESE which regresses on the variance instead of the standard error, 
#Effect = B1 + B2*SE^2 + e, which becomes Effect/SE = B1*1/SE + B2*SE + v when we weight.  
##The interpretation of coefficients is the same.

#now with clustered standard errors
FATPEESE <- lm.cluster(formula = Standardised_T ~ 0 + Precision + Standarised_se, cluster = df$Study_ID, data = df)


#Now we use a  multi-level model to adjust for the heteroskedasticity 
#and dependence between estimates
#This model is fitted with restricted maximum likelihood rather than OLS.

lmer1 <- lmer(df$Standardised_T ~ df$Precision + (1 | df$Study_ID), REML = TRUE,
                      start = 0, control = lme4:: lmerControl(calc.derivs = FALSE, optimizer ="Nelder_Mead"))




#Now do IV FAT-PET. This uses sqrt of sample size as IV instead of standard error (or rather the precision as we use WLS)
IV = as.numeric(sqrt(df$Sample_size))

IVEstimate <- ivreg(df$Standardised_T ~ df$Precision | IV)





#==============================================================================

#Non-linear methods


#trim and fill - adds in studies on left hand side of funnel
g <- metagen(TE = df$PCC, seTE = df$Standarised_se, studlab = df$Study_ID)

tf <- meta::trimfill(g)

#plot funnel of trim and fill
tf$sig <- ifelse(abs(tf$TE/tf$seTE) > 1.96,"Significant","Not Significant")
tf$added <- as.factor(ifelse(index(tf$TE) > nrow(df),"Added","Not Added"))

funnelTrim <- ggplot(as.data.frame(tf), aes(x = TE, y = 1/tf$seTE, group = added)) + 
  geom_point(cex = 4, aes(col = added),alpha=min(.7,max(40/5000,.9))) +
             theme_bw()   + theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) + 
               scale_shape_manual( values = c(25,21)) + scale_color_manual(values = c("Black","#235D96")) + 
  scale_fill_manual(values = c( "Black", "#235D96"))  +
  ylab("Precision \n1/SE")  +  
  theme(panel.grid = element_blank()) +
  theme(text = element_text(size = 15), legend.title = element_blank(), legend.position = "top") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  scale_x_continuous(expand = c(0,0),limits = c(-1,1)) + 
  geom_vline(xintercept = 0) 


tiff(file="funneltrim.tiff",width=1024,height=768,res=100, type = "cairo")



funnelTrim
dev.off()

#WAAP method only uses estimates deemed to have adequate power. Need some estimate of the" "true" effect to estimate power,
#We use the RE estimate (FE estimate is more common but this leaves only 1 study deemed to have adequate power and results in smaller WAAP estimate)

#only get adequate power estimates
WAAPTRUE <- df[df$Standarised_se <= as.numeric(abs(conv(MET2$beta))/2.8),]

#estimate WAAP with unrestricted WLS
WAAP <- lm.cluster(formula = Standardised_T ~ Precision + 0, cluster = WAAPTRUE$Short_name, data = WAAPTRUE)


#table with included studies
write.csv(table(WAAPTRUE$Short_name), "WAAPstudies.csv")


#===============================================================================

#This code performs the Andrews-Kasy (2019) method for detecting and adjusting for publication bias
# I use the code they have generously provided at https://github.com/maxkasy/MetaStudiesApp with only slight modifications


#Set parameters and vaklues for estimation
X = df$PCC
sigma = df$Standarised_se
cluster_ID = as.numeric(df$Study_ID)
symmetric = FALSE
critval=1.96
cutoffs <- c(-1.96,0,1.96)


#get estimates from method
##warning will be produced due to the low standard errors of study 7 (Grönqvist,Nilsson, Robling)####
#you can omit study 7 and get an estimate of about -0.55)
#We suppress warnings for these lines
suppressWarnings(

estimates <- metastudies_estimation(X,sigma, cutoffs,symmetric, model="normal", cluster_ID)

)

AKtable <- estimatestable(estimates$Psihat, estimates$SE, cutoffs, symmetric = FALSE, model = "normal")


suppressWarnings(
AK2 <- estimates_plot2(X, sigma, cutoffs, symmetric, estimates, model="normal"))


tiff(file="AKPCC.tiff",width=1024,height=768,res=100, type = "cairo")

suppressWarnings(
AK2 )
dev.off()
  
#======================================================================== 
#Now write up into a table


pubbias <- data.frame(NAME = rep(NA,7),Effect_beyond_bias = rep(NA,7),
                      Effect_beyond_bias_SE = rep(NA,7), Pub_Bias = rep(NA,7),  Pub_Bias_SE = rep(NA,7) )



pubbias[1,1] <- "FATPET"
pubbias[1,2] <- FATPET$lm_res$coefficients[[2]]
pubbias[1,3] <- sqrt(vcov(FATPET)[2,2])
pubbias[1,4] <- FATPET$lm_res$coefficients[[1]]
pubbias[1,5] <- sqrt(vcov(FATPET)[1,1]) 
  
pubbias[2,1] <- "FATPEESE"
pubbias[2,2] <- FATPEESE$lm_res$coefficients[[1]]
pubbias[2,3] <- sqrt(vcov(FATPEESE)[1,1])
pubbias[2,4] <- FATPEESE$lm_res$coefficients[[2]]
pubbias[2,5] <- sqrt(vcov(FATPEESE)[2,2]) 

pubbias[3,1] <- "Multi-Level FATPET"
pubbias[3,2] <- lmer1@beta[2]
pubbias[3,3] <- sqrt(vcov(lmer1)[2,2])
pubbias[3,4] <- lmer1@beta[1]
pubbias[3,5] <- sqrt(vcov(lmer1)[1,1])

pubbias[4,1] <- "IV"
pubbias[4,2] <- IVEstimate$coefficients[[2]]
pubbias[4,3] <- sqrt(vcovCL(IVEstimate, cluster = df$Study_ID)[2,2])
pubbias[4,4] <- IVEstimate$coefficients[[1]]
pubbias[4,5] <- sqrt(vcovCL(IVEstimate, cluster = df$Study_ID)[1,1])
  
pubbias[5,1] <- "Trim and Fill"
pubbias[5,2] <- tf$TE.random
pubbias[5,3] <- tf$seTE.random


pubbias[6,1] <- "WAAP"
pubbias[6,2] <- WAAP$lm_res$coefficients[[1]]
pubbias[6,3] <- sqrt(vcov(WAAP)[1.1])
  


pubbias[7,1] <- "Andrews-Kasy"
pubbias[7,2] <- estimates$Psihat[1]
pubbias[7,3] <- estimates$SE[1]


pubbias <- t(pubbias)


#write the table
write.csv(pubbias, "pubbias.csv")


  

#========================================================================
#Now with representative estimates only (so cannot do multi-level model)
#this is in the online appendix



REPFATPET  <- lm.cluster(formula = Standardised_T ~ Precision, data = df[df$Rep_estimate ==1,],
                         cluster = df$Study_ID[df$Rep_estimate==1])



REPFATPEESE  <- lm.cluster(formula = Standardised_T ~ 0 + Precision + Standarised_se, 
                           data = df[df$Rep_estimate ==1,], cluster = df$Study_ID[df$Rep_estimate==1])



REPIV = as.numeric(sqrt(df$Sample_size[df$Rep_estimate ==1]))

REPIVEstimate <- ivreg(df$Standardised_T[df$Rep_estimate ==1] ~ df$Precision[df$Rep_estimate ==1]| REPIV)




REPg <- metagen(TE = df$PCC[df$Rep_estimate ==1], seTE = df$Standarised_se[df$Rep_estimate ==1], studlab = df$Study_ID[df$Rep_estimate ==1])

REPtf <- meta::trimfill(REPg)


#make an RE estimate for the WAAP and only rep estimates

REPMET2 <- rma.uni(df$Standardised_effect_size[df$Rep_estimate==1], vi = df$var[df$Rep_estimate==1], method = "DL")

#only get adequate power estimates
REPWAAPTRUE <- df[df$Standarised_se <= as.numeric(abs(conv(REPMET2$beta))/2.8) & df$Rep_estimate ==1,]

REPWAAP <- lm.cluster(formula = Standardised_T ~ Precision + 0, 
                      cluster = REPWAAPTRUE$Short_name,
                      data = REPWAAPTRUE)



#andrews-kasy did not converge#####


#now put into table


REPpubbias <- data.frame(NAME = rep(NA,7),Effect_beyond_bias = rep(NA,7),
                         Effect_beyond_bias_SE = rep(NA,7), Pub_Bias = rep(NA,7),  Pub_Bias_SE = rep(NA,7) )



REPpubbias[1,1] <- "FATPET"
REPpubbias[1,2] <- REPFATPET$lm_res$coefficients[[2]]
REPpubbias[1,3] <- sqrt(vcov(REPFATPET)[2,2])
REPpubbias[1,4] <- REPFATPET$lm_res$coefficients[[1]]
REPpubbias[1,5] <- sqrt(vcov(REPFATPET)[1,1]) 

REPpubbias[2,1] <- "FATPEESE"
REPpubbias[2,2] <- REPFATPEESE$lm_res$coefficients[[1]]
REPpubbias[2,3] <- sqrt(vcov(REPFATPEESE)[1,1])
REPpubbias[2,4] <- REPFATPEESE$lm_res$coefficients[[2]]
REPpubbias[2,5] <- sqrt(vcov(REPFATPEESE)[2,2]) 


REPpubbias[4,1] <- "IV"
REPpubbias[4,2] <- REPIVEstimate$coefficients[[2]]
REPpubbias[4,3] <- sqrt(vcovCL(REPIVEstimate, cluster = df$Study_ID[df$Rep_estimate ==1] )[2,2])
REPpubbias[4,4] <- REPIVEstimate$coefficients[[1]]
REPpubbias[4,5] <- sqrt(vcovCL(REPIVEstimate, cluster = df$Study_ID[df$Rep_estimate ==1] )[1,1])

REPpubbias[5,1] <- "Trim and Fill"
REPpubbias[5,2] <- REPtf$TE.random
REPpubbias[5,3] <- REPtf$seTE.random


REPpubbias[6,1] <- "WAAP"
REPpubbias[6,2] <- REPWAAP$lm_res$coefficients[[1]]
REPpubbias[6,3] <- sqrt(vcov(REPWAAP)[1.1])



REPpubbias[7,1] <- "Andrews-Kasy"



REPpubbias <- t(REPpubbias)


#write the table
write.csv(REPpubbias, "REPpubbias.csv")


#========================================================================
#Now do pub bias tests just for addressing endogeneity sub sample

dfEndo <- df[df$Endogeneity==1 ,]



FATPET <- lm.cluster(formula = Standardised_T ~ Precision, cluster = dfEndo$Study_ID, data = dfEndo)




#Next the FAT-PEESE which regresses on the variance instead of the standard error, 
#Effect = B1 + B2*SE^2 + e, which becomes Effect/SE = B1*1/SE + B2*SE + v when we weight.  
##The interpretation of coefficients is the same.

#now with clustered standard errors
FATPEESE <- lm.cluster(formula = Standardised_T ~ 0 + Precision + Standarised_se, cluster = dfEndo$Study_ID, data = dfEndo)


#Now we use a mixed-effects multi-level model to adjust for the heteroskedasticity 
#and dependence between estimates
#This model is fitted with restricted maximum likelihood rather than OLS.

lmer1 <- lmer(dfEndo$Standardised_T ~ dfEndo$Precision + (1 | dfEndo$Study_ID), REML = TRUE,
              start = 0, control = lme4:: lmerControl(calc.derivs = FALSE, optimizer ="Nelder_Mead"))




#Now do IV FAT-PET. This uses sqrt of sample size as IV instead of standard error (or rather the precision as we use WLS)
IV = as.numeric(sqrt(dfEndo$Sample_size))

IVEstimate <- ivreg(dfEndo$Standardised_T ~ dfEndo$Precision | IV)





#==============================================================================

#Non-linear methods


#trim and fill - adds in studies on left hand side of funnel
g <- metagen(TE = dfEndo$PCC, seTE = dfEndo$Standarised_se, studlab = dfEndo$Study_ID)

tf <- meta::trimfill(g)

#plot funnel of trim and fill
tf$sig <- ifelse(abs(tf$TE/tf$seTE) > 1.96,"Significant","Not Significant")
tf$added <- as.factor(ifelse(index(tf$TE) > nrow(dfEndo),"Added","Not Added"))

funnelTrimEndog <- ggplot(as.data.frame(tf), aes(x = TE, y = 1/tf$seTE, group = added)) + 
  geom_point(cex = 4, aes(col = added),alpha=min(.7,max(40/5000,.9))) +
  theme_bw()   + theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) + 
  scale_shape_manual( values = c(25,21)) + scale_color_manual(values = c("Black","#235D96")) + 
  scale_fill_manual(values = c( "Black", "#235D96"))  +
  ylab("Precision \n1/SE")  +  
  theme(panel.grid = element_blank()) +
  theme(text = element_text(size = 15), legend.title = element_blank(), legend.position = "top") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  scale_x_continuous(expand = c(0,0),limits = c(-1,1)) + 
  geom_vline(xintercept = 0) 


tiff(file="funneltrimEndog.tiff",width=1024,height=768,res=100, type = "cairo")



funnelTrimEndog
dev.off()


#make an RE estimate for the WAAP and only addressing endog estimates

EndogMET2 <- rma.uni(dfEndo$Standardised_effect_size, vi = dfEndo$var, method = "DL")
#only get adequate power estimates
WAAPTRUE <- dfEndo[dfEndo$Standarised_se <= as.numeric(abs(conv(EndogMET2$beta))/2.8),]

#estimate WAAP with unrestricted WLS
WAAP <- lm.cluster(formula = Standardised_T ~ Precision + 0, cluster = WAAPTRUE$Short_name, data = WAAPTRUE)




#===============================================================================


#Set parameters and vaklues for estimation
X = dfEndo$PCC
sigma = dfEndo$Standarised_se
cluster_ID = as.numeric(dfEndo$Study_ID)
symmetric = FALSE
critval=1.96
cutoffs <- c(-1.96,0,1.96)


suppressWarnings(
  
  estimates <- metastudies_estimation(X,sigma, cutoffs,symmetric, model="normal", cluster_ID)
  
)

AKtable <- estimatestable(estimates$Psihat, estimates$SE, cutoffs, symmetric = FALSE, model = "normal")


suppressWarnings(
  AK2 <- estimates_plot2(X, sigma, cutoffs, symmetric, estimates, model="normal"))


tiff(file="AKPCCEndog.tiff",width=1024,height=768,res=100, type = "cairo")

suppressWarnings(
  AK2 )
dev.off()


#now make file

pubbias <- data.frame(NAME = rep(NA,7),Effect_beyond_bias = rep(NA,7),
                      Effect_beyond_bias_SE = rep(NA,7), Pub_Bias = rep(NA,7),  Pub_Bias_SE = rep(NA,7) )



pubbias[1,1] <- "FATPET"
pubbias[1,2] <- FATPET$lm_res$coefficients[[2]]
pubbias[1,3] <- sqrt(vcov(FATPET)[2,2])
pubbias[1,4] <- FATPET$lm_res$coefficients[[1]]
pubbias[1,5] <- sqrt(vcov(FATPET)[1,1]) 

pubbias[2,1] <- "FATPEESE"
pubbias[2,2] <- FATPEESE$lm_res$coefficients[[1]]
pubbias[2,3] <- sqrt(vcov(FATPEESE)[1,1])
pubbias[2,4] <- FATPEESE$lm_res$coefficients[[2]]
pubbias[2,5] <- sqrt(vcov(FATPEESE)[2,2]) 

pubbias[3,1] <- "Multi-Level FATPET"
pubbias[3,2] <- lmer1@beta[2]
pubbias[3,3] <- sqrt(vcov(lmer1)[2,2])
pubbias[3,4] <- lmer1@beta[1]
pubbias[3,5] <- sqrt(vcov(lmer1)[1,1])

pubbias[4,1] <- "IV"
pubbias[4,2] <- IVEstimate$coefficients[[2]]
pubbias[4,3] <- sqrt(vcovCL(IVEstimate, cluster = dfEndo$Study_ID)[2,2])
pubbias[4,4] <- IVEstimate$coefficients[[1]]
pubbias[4,5] <- sqrt(vcovCL(IVEstimate, cluster = dfEndo$Study_ID)[1,1])

pubbias[5,1] <- "Trim and Fill"
pubbias[5,2] <- tf$TE.random
pubbias[5,3] <- tf$seTE.random


pubbias[6,1] <- "WAAP"
pubbias[6,2] <- WAAP$lm_res$coefficients[[1]]
pubbias[6,3] <- sqrt(vcov(WAAP)[1.1])



pubbias[7,1] <- "Andrews-Kasy"
pubbias[7,2] <- estimates$Psihat[1]
pubbias[7,3] <- estimates$SE[1]


pubbias <- t(pubbias)


#write the table
write.csv(pubbias, "pubbiasEndog.csv")

#========================================================================






#

  