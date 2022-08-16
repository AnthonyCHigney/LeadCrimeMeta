
#============================================
#Data frame subset of just the elasticities
dfE <- df[!is.na(df$Elasticity),]

#new elasticity variables for computations
dfE$Elas_Precision <- 1/dfE$Elasticity_SE
dfE$T_Elast <- dfE$Elasticity/dfE$Elasticity_SE


weigh1 <- function(x) {
  x/dfE$Elasticity_SE
}

#Weigh the covariates by the PCC standard errors using namevector to name the columns
dfE[,namevector] <- lapply(dfE[orignames] , FUN = weigh1)


#Now make funnel plot
dfE$sig <- ifelse(abs(dfE$Elasticity/dfE$Elasticity_SE) > 1.96,"Significant","Not Significant")


funnelELAS <- ggplot(dfE, aes(x = Elasticity, y = Elas_Precision, group = sig)) + 
  geom_point(cex = 4, aes(shape = sig, col = sig, fill = sig, group = sig), alpha=min(.5,max(40/5000,.3))) +
  
  theme_bw()   + theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) + 
  scale_shape_manual( values = c(25,21)) + scale_color_manual(values = c("#B95465", "#235D96")) + 
  scale_fill_manual(values = c("#B95465", "#235D96")) +
  ylab("Precision \n1/SE")  +  
  theme(panel.grid = element_blank()) +
  theme(text = element_text(size = 15), legend.title = element_blank(), legend.position = "top") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
 scale_x_continuous(expand = c(0,0),limits = c(-4,4)) + 
  geom_vline(xintercept = 0) 





tiff(file="funnelelasticty.tiff",width=1024,height=768,res=100, type = "cairo")



funnelELAS
dev.off()


tiff(file="funnelbothELAS.tiff",width=800,height=1024,res=100, type = "cairo")
grid.arrange(funnelPCC, funnelELAS, nrow=2, ncol=1)
dev.off()




#==============================================================================

#now forest plot


#Fixed effects estimate for full sample
ELASMET1 <- rma.uni(dfE$Elasticity, sei = dfE$Elasticity_SE, method = "FE")


#random effects estimate for full sample
ELASMET2 <- rma.uni(dfE$Elasticity, sei = dfE$Elasticity_SE, method = "DL")





fr <- by(dfE, dfE$Study_ID,
         function(x) rma.uni(Elasticity, sei = Elasticity_SE, data = x, method = "FE"))

#study ID no longer conecutive so chnage

names(fr) <- 1:length(unique(dfE$Short_name))

#make a dataframe to hold FE averages, standard errors and confidence bounds

f <- data.frame(rep(0,length(unique(dfE$Short_name))),
                rep(0,length(unique(dfE$Short_name))),
                rep(0,length(unique(dfE$Short_name))),
                rep(0,length(unique(dfE$Short_name))), 
                rep(0,length(unique(dfE$Short_name))))
names(f) <- c("Average", "SE", "CI_L", "CI_U", "names")

#fill columns with PCCs, standard errors and lower and upper bounds
nn <- length(unique(dfE$Short_name))


for (i in 1:nn){
  betas <- eval(parse(text = paste0("fr$`", i, "`$beta")))
  SEs <- eval(parse(text = paste0("fr$`", i, "`$se")))
  CI_L <- eval(parse(text = paste0("fr$`", i, "`$ci.lb")))
  CI_U <- eval(parse(text = paste0("fr$`", i, "`$ci.ub")))
  id <- i
  f[i,1] <- conv(betas)
  f[i,2] <- SEs
  f[i,3] <- CI_L
  f[i,4] <- CI_U
  f[i,5] <- id
  
}



#get study names for forest plot

f$name <- unique(dfE$Short_name)


#order data frame by smallest effect size to largest
f <- f[order(f$Average),]


#make forest plot and export
tiff(file="forestplotElas.tiff",width=2400, height=1800,res=200, type = "cairo")
metafor::forest(f$Average, sei = f$SE, slab = f$name, xlab = "Elasticity", cex = 1.5, ylim = c(-3,15))
abline(0,0)
addpoly(x = ELASMET2$beta, sei = ELASMET2$se, efac = 2, mlab = c("Random Effects Model"), row = -2, cex = 1.5)
addpoly(x = ELASMET1$beta, sei = ELASMET1$se, efac = 2, mlab = c("Common Effects Model"), row = -1, cex = 1.5)
par(font=2)
text(-1.5, 15, "Author(s) and Year",  pos=1, cex = 1.5)

dev.off()



#================================================================
#Pub bias tests


#FAT-PET
ELASFATPET  <- lm.cluster(formula = T_Elast ~ Elas_Precision, data = dfE,
                         cluster = dfE$Study_ID)


#PEESE
ELASFATPEESE  <- lm.cluster(formula = T_Elast ~ 0 + Elas_Precision + Elasticity_SE, 
                           data = dfE, cluster = dfE$Study_ID)

#mixed model


ELASlmer1 <- lmer(dfE$T_Elast ~ dfE$Elas_Precision + (1 | dfE$Study_ID), REML = TRUE,
              start = 0, control = lme4:: lmerControl(calc.derivs = FALSE, optimizer ="Nelder_Mead"))


#IV
ELASIV = as.numeric(sqrt(dfE$Sample_size))

ELASIVEstimate <- ivreg(dfE$T_Elast ~ dfE$Elas_Precision| ELASIV)



#Trim and fill
ELASg <- metagen(TE = dfE$Elasticity, seTE = dfE$Elasticity_SE, studlab = dfE$Study_ID)

ELAStf <- meta::trimfill(ELASg)

#WAAP

#need to get some estimate for effect for use in WAAP



WAAPTRUE <- dfE[dfE$Elasticity_SE <= as.numeric(abs(ELASMET2$beta)/2.8),]

ELASWAAP <- lm.cluster(formula = T_Elast ~ Elas_Precision + 0, cluster = WAAPTRUE$Short_name,
                       data = WAAPTRUE)


#Andrews -kasy

X = dfE$Elasticity
sigma = dfE$Elasticity_SE
cluster_ID = as.numeric(dfE$Study_ID)
symmetric = FALSE
critval=1.96
cutoffs <- c(-1.96,0,1.96)


suppressWarnings(
  ELASestimates <- metastudies_estimation(X,sigma, cutoffs,symmetric, model="normal", cluster_ID)
  
)


suppressWarnings(
elastAK2 <- estimates_plot2(X, sigma, cutoffs, symmetric, ELASestimates, model="normal"))

tiff(file="ELASAK.tiff",width=1024,height=768,res=100, type = "cairo")
elastAK2
dev.off()




#now put into table


ELASpubbias <- data.frame(NAME = rep(NA,7),Effect_beyond_bias = rep(NA,7),
                         Effect_beyond_bias_SE = rep(NA,7), Pub_Bias = rep(NA,7),  Pub_Bias_SE = rep(NA,7) )



ELASpubbias[1,1] <- "FATPET"
ELASpubbias[1,2] <- ELASFATPET$lm_res$coefficients[[2]]
ELASpubbias[1,3] <- sqrt(vcov(ELASFATPET)[2,2])
ELASpubbias[1,4] <- ELASFATPET$lm_res$coefficients[[1]]
ELASpubbias[1,5] <- sqrt(vcov(ELASFATPET)[1,1]) 

ELASpubbias[2,1] <- "FATPEESE"
ELASpubbias[2,2] <- ELASFATPEESE$lm_res$coefficients[[1]]
ELASpubbias[2,3] <- sqrt(vcov(ELASFATPEESE)[1,1])
ELASpubbias[2,4] <- ELASFATPEESE$lm_res$coefficients[[2]]
ELASpubbias[2,5] <- sqrt(vcov(ELASFATPEESE)[2,2]) 

ELASpubbias[3,1] <- "Multi-Level FATPET"
ELASpubbias[3,2] <- ELASlmer1@beta[2]
ELASpubbias[3,3] <- sqrt(vcov(ELASlmer1)[2,2])
ELASpubbias[3,4] <- ELASlmer1@beta[1]
ELASpubbias[3,5] <- sqrt(vcov(ELASlmer1)[1,1])

ELASpubbias[4,1] <- "IV"
ELASpubbias[4,2] <- ELASIVEstimate$coefficients[[2]]
ELASpubbias[4,3] <- sqrt(vcovCL(ELASIVEstimate, cluster = dfE$Study_ID)[2,2])
ELASpubbias[4,4] <- ELASIVEstimate$coefficients[[1]]
ELASpubbias[4,5] <- sqrt(vcovCL(ELASIVEstimate, cluster = dfE$Study_ID)[1,1])

ELASpubbias[5,1] <- "Trim and Fill"
ELASpubbias[5,2] <- ELAStf$TE.random
ELASpubbias[5,3] <- ELAStf$seTE.random


ELASpubbias[6,1] <- "WAAP"
ELASpubbias[6,2] <- ELASWAAP$lm_res$coefficients[[1]]
ELASpubbias[6,3] <- sqrt(vcov(ELASWAAP)[1.1])



ELASpubbias[7,1] <- "Andrews-Kasy"
ELASpubbias[7,2] <- ELASestimates$Psihat[1]
ELASpubbias[7,3] <- ELASestimates$SE[1]


ELASpubbias <- t(ELASpubbias)





write.csv(ELASpubbias, "Elaspubbias.csv")



#==================================
studies <- data.frame( rep(0,length(unique(dfE$Short_name))), 0, 0, 0 ,0 ,0 ,0 , 0, 0,0,0)
names(studies) <-  c("Study & Year"	,"Median",	"Mean",	"Fixed Effects Average" , 
                     "Type of Crime" , "Individual or Area","Addressing Endogeneity", "elasticity")

studies[,1]  <- (dfE %>% group_by(Study_ID) %>% summarize(name = unique(Short_name)))[2]
studies[,2]  <-  (dfE %>% group_by(Study_ID) %>% summarize(median = median(Elasticity)))[2]
studies[,3]  <-  (dfE %>% group_by(Study_ID) %>% summarize(mean = mean(Elasticity)))[2]


for (i in 1:length(unique(dfE$Short_name))){
  betas <- eval(parse(text = paste0("fr$`", i, "`$beta")))
  studies[i,4] <- conv(betas)
  next
}


write.csv(studies, "StudiesRElas.csv", row.names = FALSE)






#====================================================================


#Now just do Elasticity and addressing endogeneity



dfE1 <- dfE[dfE$Endogeneity==1,]






#FAT-PET
ELASFATPET  <- lm.cluster(formula = T_Elast ~ Elas_Precision, data = dfE1,
                          cluster = dfE1$Study_ID)


#PEESE
ELASFATPEESE  <- lm.cluster(formula = T_Elast ~ 0 + Elas_Precision + Elasticity_SE, 
                            data = dfE1, cluster = dfE1$Study_ID)

#mixed model


ELASlmer1 <- lmer(dfE1$T_Elast ~ dfE1$Elas_Precision + (1 | dfE1$Study_ID), REML = TRUE,
                  start = 0, control = lme4:: lmerControl(calc.derivs = FALSE, optimizer ="Nelder_Mead"))


#IV
ELASIV = as.numeric(sqrt(dfE1$Sample_size))

ELASIVEstimate <- ivreg(dfE1$T_Elast ~ dfE1$Elas_Precision| ELASIV)



#Trim and fill
ELASg <- metagen(TE = dfE1$Elasticity, seTE = dfE1$Elasticity_SE, studlab = dfE1$Study_ID)

ELAStf <- meta::trimfill(ELASg)

#WAAP

#need to get some estimate for effect for use in WAAP

#random effects estimate for full sample
ELASENDOMET2 <- rma.uni(dfE1$Elasticity, sei = dfE1$Elasticity_SE, method = "DL")



WAAPTRUE <- dfE1[dfE1$Elasticity_SE <= as.numeric(abs(ELASENDOMET2$beta)/2.8),]

ELASWAAP <- lm.cluster(formula = T_Elast ~ Elas_Precision + 0, cluster = WAAPTRUE$Short_name,
                       data = WAAPTRUE)


#Andrews -kasy

X = dfE1$Elasticity
sigma = dfE1$Elasticity_SE
cluster_ID = as.numeric(dfE1$Study_ID)
symmetric = FALSE
critval=1.96
cutoffs <- c(-1.96,0,1.96)


suppressWarnings(
  ELASestimates <- metastudies_estimation(X,sigma, cutoffs,symmetric, model="normal", cluster_ID)
  
)


suppressWarnings(
  elastAK2 <- estimates_plot2(X, sigma, cutoffs, symmetric, ELASestimates, model="normal"))

tiff(file="ELASAKEndog.tiff",width=1024,height=768,res=100, type = "cairo")
elastAK2
dev.off()




#now put into table


ELASpubbias <- data.frame(NAME = rep(NA,7),Effect_beyond_bias = rep(NA,7),
                          Effect_beyond_bias_SE = rep(NA,7), Pub_Bias = rep(NA,7),  Pub_Bias_SE = rep(NA,7) )



ELASpubbias[1,1] <- "FATPET"
ELASpubbias[1,2] <- ELASFATPET$lm_res$coefficients[[2]]
ELASpubbias[1,3] <- sqrt(vcov(ELASFATPET)[2,2])
ELASpubbias[1,4] <- ELASFATPET$lm_res$coefficients[[1]]
ELASpubbias[1,5] <- sqrt(vcov(ELASFATPET)[1,1]) 

ELASpubbias[2,1] <- "FATPEESE"
ELASpubbias[2,2] <- ELASFATPEESE$lm_res$coefficients[[1]]
ELASpubbias[2,3] <- sqrt(vcov(ELASFATPEESE)[1,1])
ELASpubbias[2,4] <- ELASFATPEESE$lm_res$coefficients[[2]]
ELASpubbias[2,5] <- sqrt(vcov(ELASFATPEESE)[2,2]) 

ELASpubbias[3,1] <- "Multi-Level FATPET"
ELASpubbias[3,2] <- ELASlmer1@beta[2]
ELASpubbias[3,3] <- sqrt(vcov(ELASlmer1)[2,2])
ELASpubbias[3,4] <- ELASlmer1@beta[1]
ELASpubbias[3,5] <- sqrt(vcov(ELASlmer1)[1,1])

ELASpubbias[4,1] <- "IV"
ELASpubbias[4,2] <- ELASIVEstimate$coefficients[[2]]
ELASpubbias[4,3] <- sqrt(vcovCL(ELASIVEstimate, cluster = dfE1$Study_ID)[2,2])
ELASpubbias[4,4] <- ELASIVEstimate$coefficients[[1]]
ELASpubbias[4,5] <- sqrt(vcovCL(ELASIVEstimate, cluster = dfE1$Study_ID)[1,1])

ELASpubbias[5,1] <- "Trim and Fill"
ELASpubbias[5,2] <- ELAStf$TE.random
ELASpubbias[5,3] <- ELAStf$seTE.random


ELASpubbias[6,1] <- "WAAP"
ELASpubbias[6,2] <- ELASWAAP$lm_res$coefficients[[1]]
ELASpubbias[6,3] <- sqrt(vcov(ELASWAAP)[1.1])



ELASpubbias[7,1] <- "Andrews-Kasy"
ELASpubbias[7,2] <- ELASestimates$Psihat[1]
ELASpubbias[7,3] <- ELASestimates$SE[1]


ELASpubbias <- t(ELASpubbias)



#write the table
write.csv(ELASpubbias, "pubbiasEndogElas.csv")



#==========================================================================
#what studies included

studies <- data.frame( rep(0,length(unique(dfE1$Short_name))), 0, 0, 0 ,0 ,0 ,0 , 0, 0,0,0)
names(studies) <-  c("Study & Year"	,"Median",	"Mean",	"Fixed Effects Average" , 
                     "Type of Crime" , "Individual or Area","Addressing Endogeneity", "elasticity")

studies[,1]  <- (dfE1 %>% group_by(Study_ID) %>% summarize(name = unique(Short_name)))[2]
studies[,2]  <-  (dfE1 %>% group_by(Study_ID) %>% summarize(median = median(Elasticity)))[2]
studies[,3]  <-  (dfE1 %>% group_by(Study_ID) %>% summarize(mean = mean(Elasticity)))[2]


for (i in 1:length(unique(dfE1$Short_name))){
  betas <- eval(parse(text = paste0("fr$`", i, "`$beta")))
  studies[i,4] <- conv(betas)
  next
}


write.csv(studies, "StudiesRElasEndog.csv", row.names = FALSE)




