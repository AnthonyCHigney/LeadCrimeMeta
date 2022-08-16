#forest plot code

#pool the estimates from each study with weighted average for each individual study

fr <- by(df, df$Study_ID,
         function(x) rma.uni(Standardised_effect_size, var, data = x, method = "FE"))

#Get Common Effects and random effects estimates for all estimates combined from every study


#This is common effects estimate for full sample
MET1 <- rma.uni(df$Standardised_effect_size, df$var, method = "FE")



#random effects estimate for full sample
MET2 <- rma.uni(df$Standardised_effect_size, df$var, method = "DL")







#make a dataframe to hold weighted averages, standard errors and confidence bounds

f <- data.frame(rep(0,24),rep(0,24),rep(0,24),rep(0,24), rep(0,24))
names(f) <- c("Average", "SE", "CI_L", "CI_U", "names")

#fill columns with PCCs, standard errors and lower and upper bounds

for (i in 1:24){
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

f$name <- unique(df$Short_name)
  

#order data frame by smallest effect size to largest
f <- f[order(f$Average),]



#make forest plot and export
tiff(file="forestplot.tiff",width=2400, height=2400,res=200, type = "cairo")
metafor::forest(f$Average, sei = f$SE, slab = f$name, xlab = "PCC", ylim=c(-6, 27), cex = 1.5)
abline(0,0)
addpoly(x = conv(MET1$beta), sei = MET1$se, efac = 2, mlab = c("Common Effects Model"), row = -2, cex = 1.5)
addpoly(x = conv(MET2$beta), sei = MET2$se, efac = 2, mlab = c("Random Effects Model"), row = -4, cex = 1.5)
par(font=2)
text(-1.3, 27, "Author(s) and Year",  pos=1, cex = 1.5)

dev.off()



#===========================================================
#Make a table for studies


studies <- data.frame( rep(0,24), 0, 0, 0 ,0 ,0 ,0 , 0, 0,0,0)
names(studies) <-  c("Study & Year"	,"Median",	"Mean",	"Fixed Effects Average" , 
                     "Type of Crime" , "Individual or Area","Addressing Endogeneity", "elasticity")

studies[,1]  <- (df %>% group_by(Study_ID) %>% summarize(name = unique(Short_name)))[2]
studies[,2]  <-  (df %>% group_by(Study_ID) %>% summarize(pcc = median(PCC)))[2]
studies[,3]  <-  (df %>% group_by(Study_ID) %>% summarize(pcc = mean(PCC)))[2]
studies[,8]  <-  (df %>% group_by(Study_ID) %>% summarize(Elasticity = mean(Elasticity, na.rm = TRUE)))[2]

for (i in 1:24){
  betas <- eval(parse(text = paste0("fr$`", i, "`$beta")))
  studies[i,4] <- conv(betas)
  next
}


write.csv(studies, "StudiesR.csv", row.names = FALSE)