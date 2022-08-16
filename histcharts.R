
tmp <- read.csv("tmp.csv")

tmpquas <- read.csv("tmpquas.csv")
tmpnonquas <- read.csv("tmpnonquas.csv")
tmparea <- read.csv("tmparea.csv")
tmpind <- read.csv("tmpind.csv")
tmphomicide <- read.csv("tmphomicide.csv")
tmpviolent <- read.csv("tmpviolent.csv")
tmpnonviolent <- read.csv("tmpnonviolent.csv")
tmpELAS <- read.csv("tmpELAS.csv")
tmpELASENDO <- read.csv("tmpElasEndo.csv")





#===================================================================
#make a big dataframe with all estimates in together.

#first need to add a name in each dataframe
tmp$name <- "Full Sample"
tmpquas$name <- "Addressing Endogeneity Sample"
tmpnonquas$name <- "Correlational Sample"
tmparea$name <- "Area Sample"
tmpind$name <- "Individual Sample"
tmphomicide$name <- "Homicide Sample"
tmpviolent$name <- "Violent Crime Sample"
tmpnonviolent$name <- "Non-violent Crime Sample"
tmpELAS$name <- "Elasticity Sample"
tmpELASENDO$name <- "Elasticity and Addressing Endogeneity Sample"
#function to help with binding

rb <- function(x) {
  x[, c("average", "idealNAAv", "FAT", "name")]
}

bigtmp <- rbind(rb(tmp), rb(tmpquas), rb(tmpnonquas), rb(tmparea), 
                rb(tmpind), rb(tmphomicide), rb(tmpviolent), rb(tmpnonviolent), rb(tmpELAS), rb(tmpELASENDO))




#======================================================================
#make density plots for Average and FAT

#full sample first




#Sample Averages
gg1 <- ggplot(data = tmp, aes(x = average)) + geom_density(color = "#004488", lwd = 2)



gg1 <- gg1 + ylab("") + xlab(label ="PCC") + theme_bw() + scale_x_continuous(limits = c(-1,1)) + theme(text = element_text(size = 15))

gg1 <- gg1 +  theme(panel.grid = element_blank(), text = element_text(size=20))

gg1

tiff(file="densfullsampAVG.tiff",width=1024,height=600,res=100, type = "cairo")
gg1
dev.off()




#"ideal" specification
gg2 <- ggplot(data = tmp, aes(x = idealNAAv)) + geom_density(color = "#004488", lwd = 2)



gg2 <- gg2 + ylab("") + xlab(label ="PCC") + theme_bw() + scale_x_continuous(limits = c(-1,1)) + theme(text = element_text(size = 15))

gg2 <- gg2 +  theme(panel.grid = element_blank(), text = element_text(size=20))

gg2

tiff(file="densfullideal.tiff",width=1024,height=600,res=100, type = "cairo")
gg2
dev.off()

tiff(file="densfullBoth.tiff",width=1024,height=600,res=100, type = "cairo")
grid.arrange(gg1, gg2, nrow=1, ncol=2)
dev.off()



#FAT
gg <- ggplot(data = tmp, aes(x = FAT)) + geom_density(color = "#004488", lwd = 2)



gg <- gg + ylab("") + xlab("FAT") + theme_bw()

gg <- gg +  theme(panel.grid = element_blank()) + theme(text = element_text(size = 15))


gg

tiff(file="densfullsampFAT.tiff",width=1024,height=600,res=100, type = "cairo")
gg
dev.off()
#===============================================================================
#same for elasticity
#Sample Averages
gg1 <- ggplot(data = tmpELAS, aes(x = average)) + geom_density(color = "#004488", lwd = 2)



gg1 <- gg1 + ylab("") + xlab(label ="Elasticity") + theme_bw() + scale_x_continuous(limits = c(-1,1)) + theme(text = element_text(size = 15))

gg1 <- gg1 +  theme(panel.grid = element_blank(), text = element_text(size=20))

gg1

tiff(file="densELASsampAVG.tiff",width=1024,height=600,res=100, type = "cairo")
gg1
dev.off()




#"ideal" specification
gg2 <- ggplot(data = tmpELAS, aes(x = idealNAAv)) + geom_density(color = "#004488", lwd = 2)



gg2 <- gg2 + ylab("") + xlab(label ="Elasticity") + theme_bw() + scale_x_continuous(limits = c(-1,1)) + theme(text = element_text(size = 15))

gg2 <- gg2 +  theme(panel.grid = element_blank(), text = element_text(size=20))

gg2

tiff(file="densELASideal.tiff",width=1024,height=600,res=100, type = "cairo")
gg2
dev.off()

tiff(file="densfullBothELAS.tiff",width=1024,height=600,res=100, type = "cairo")
grid.arrange(gg1, gg2, nrow=1, ncol=2)
dev.off()

#=============================================================

#now "Elasticity and Addressing Endogeneity Sample"
#Sample Averages
gg1 <- ggplot(data = tmpELASENDO, aes(x = average)) + geom_density(color = "#004488", lwd = 2)



gg1 <- gg1 + ylab("") + xlab(label ="Elasticity") + theme_bw() + scale_x_continuous(limits = c(-1,1)) + theme(text = element_text(size = 15))

gg1 <- gg1 +  theme(panel.grid = element_blank(), text = element_text(size=20))

gg1

tiff(file="densELASENDOsampAVG.tiff",width=1024,height=600,res=100, type = "cairo")
gg1
dev.off()




#"ideal" specification
gg2 <- ggplot(data = tmpELASENDO, aes(x = idealNAAv)) + geom_density(color = "#004488", lwd = 2)



gg2 <- gg2 + ylab("") + xlab(label ="Elasticity") + theme_bw() + scale_x_continuous(limits = c(-1,1)) + theme(text = element_text(size = 15))

gg2 <- gg2 +  theme(panel.grid = element_blank(), text = element_text(size=20))

gg2

tiff(file="densELASENDOideal.tiff",width=1024,height=600,res=100, type = "cairo")
gg2
dev.off()

tiff(file="densfullBothELASENDO.tiff",width=1024,height=600,res=100, type = "cairo")
grid.arrange(gg1, gg2, nrow=1, ncol=2)
dev.off()


























#==========================================================
#Now do same for addressing endogeneity sample
#Sample Averages
gg1 <- ggplot(data = tmpquas, aes(x = average)) + geom_density(color = "#004488", lwd = 2)



gg1 <- gg1 + ylab("") + xlab(label ="PCC") + theme_bw() + scale_x_continuous(limits = c(-1,1)) + theme(text = element_text(size = 15))

gg1 <- gg1 +  theme(panel.grid = element_blank(), text = element_text(size=20))

gg1

tiff(file="densEndosampAVG.tiff",width=1024,height=600,res=100, type = "cairo")
gg1
dev.off()




#"ideal" specification
gg2 <- ggplot(data = tmpquas, aes(x = idealNAAv)) + geom_density(color = "#004488", lwd = 2)



gg2 <- gg2 + ylab("") + xlab(label ="PCC") + theme_bw() + scale_x_continuous(limits = c(-1,1)) + theme(text = element_text(size = 15))

gg2 <- gg2 +  theme(panel.grid = element_blank(), text = element_text(size=20))

gg2

tiff(file="densEndoideal.tiff",width=1024,height=600,res=100, type = "cairo")
gg2
dev.off()

tiff(file="densEndoBoth.tiff",width=1024,height=600,res=100, type = "cairo")
grid.arrange(gg1, gg2, nrow=1, ncol=2)
dev.off()



#FAT
gg <- ggplot(data = tmpquas, aes(x = FAT)) + geom_density(color = "#004488", lwd = 2)



gg <- gg + ylab("") + xlab("FAT") + theme_bw()

gg <- gg +  theme(panel.grid = element_blank()) + theme(text = element_text(size = 15))


gg

tiff(file="densEndosampFAT.tiff",width=1024,height=600,res=100, type = "cairo")
gg
dev.off()
#==========================================================

#subsample
#non-scaled version

gg <- ggplot(data = bigtmp[bigtmp$name != "Full Sample" & bigtmp$name != "Addressing Endogeneity Sample" & 
                             bigtmp$name != "Elasticity Sample" ,], 
             aes(x = average)) + geom_density(color = "#004488", lwd = 2) + xlab(label =expression(hat(beta)[M])) + 
  facet_wrap(~ name, scales = "free", nrow = 4, ncol = 3 ) 

gg <- gg + ylab("")  + theme_bw() + theme(text = element_text(size = 15))

gg <- gg +  theme(panel.grid = element_blank()) + theme(text = element_text(size = 15))


gg


tiff(file="densnonscaleAVG.tiff",width=1024,height=600,res=100, type = "cairo")
gg
dev.off()


#scaled version

gg <- ggplot(data = bigtmp[bigtmp$name != "Full Sample" & bigtmp$name != "Addressing Endogeneity Sample"
                           & bigtmp$name != "Elasticity Sample" ,], 
             aes(x = average)) + geom_density(color = "#004488", lwd = 2) + 
  facet_wrap(~ name, scales = "free", nrow = 4, ncol = 3 ) + scale_x_continuous(limits = c(-1,1))

gg <- gg + ylab("") + xlab("") + theme_bw() + theme(text = element_text(size = 15))


gg <- gg +  theme(panel.grid = element_blank())

gg

tiff(file="densscaleAVG.tiff",width=1024,height=600,res=100, type = "cairo")
gg
dev.off()

#=======================================
#same with "ideal"

gg <- ggplot(data = bigtmp[bigtmp$name != "Full Sample" & bigtmp$name != "Addressing Endogeneity Sample"
                           & 
                             bigtmp$name != "Elasticity Sample" ,], aes(x = idealNAAv)) + geom_density(color = "#004488", lwd = 2) + xlab(label =expression(hat(beta)[M])) + facet_wrap(~ name, scales = "free", nrow = 4, ncol = 3 ) 

gg <- gg + ylab("")  + theme_bw() + theme(text = element_text(size = 15))

gg <- gg +  theme(panel.grid = element_blank()) + theme(text = element_text(size = 15))


gg


tiff(file="densnonscaleIDEAL.tiff",width=1024,height=600,res=100, type = "cairo")
gg
dev.off()


#scaled version

gg <- ggplot(data = bigtmp[bigtmp$name != "Full Sample" & bigtmp$name != "Addressing Endogeneity Sample" & 
                             bigtmp$name != "Elasticity Sample" ,], aes(x = idealNAAv)) + geom_density(color = "#004488", lwd = 2) + facet_wrap(~ name, scales = "free", nrow = 4, ncol = 3 ) + scale_x_continuous(limits = c(-1,1))

gg <- gg + ylab("") + xlab("") + theme_bw() + theme(text = element_text(size = 15))


gg <- gg +  theme(panel.grid = element_blank())

gg

tiff(file="densscaleIDEAL.tiff",width=1024,height=600,res=100, type = "cairo")
gg
dev.off()














#==================================================

gg <- ggplot(data = bigtmp, aes(x = FAT)) + geom_density(color = "#004488", lwd = 2) + facet_wrap(~ name, scales = "free", nrow = 4, ncol = 3 ) 

gg <- gg + ylab("") + xlab("FAT") + theme_bw() + theme(text = element_text(size = 15))


gg <- gg +  theme(panel.grid = element_blank())

gg


tiff(file="densFAT.tiff",width=1024,height=600,res=100, type = "cairo")
gg
dev.off()



#=====================================================================
#make a table for all means and medians




#add in the N for each table

Avgtab <- bigtmp %>%  group_by(name) %>% summarize(mean = mean(average), median = median(average),
                                        N = length(average), sd =sd(average), pos = length(average[average>1 | average < -1])/length(average))




#now for FAT

FATtab <- bigtmp %>%  group_by(name) %>% summarize(mean = mean(FAT), median = median(FAT), 
                                        N = length(FAT))

Idealtab <- bigtmp %>%  group_by(name) %>% summarize(mean = mean(idealNAAv), median = median(idealNAAv),
                                                  N = length(idealNAAv), sd =sd(idealNAAv), pos = length(idealNAAv[idealNAAv>1 | idealNAAv < -1])/length(idealNAAv))




write.csv(FATtab, "FATtab.csv")



write.csv(Avgtab, "Avgtab.csv")





write.csv(Idealtab, "Idealtab.csv")


#==========================================================================================
#chart all histograms of covariate coefficients


X <- c(         "Control_gender" ,        
                "Control_race"  ,          "Control_income" ,         "Control_education",        "Homicide" ,              
                "Violent" ,                "Non_Violent",             "Area_dummy" ,            
                "OLS"   ,                  "ML"    ,                  "Odds_Ratio"  ,            "Panel_dummy" ,           
                "Endogeneity"  ,           "North_America" ,          "Europe"   ,                      
                "Direct_Lead_Measure",     "PubYear",                  "ZCovariates",             "Zsample" )



toq <- reshape2::melt(tmp, measure.vars = X)


toq <- toq[c("variable", "value")]
toq$variable <- as.character(toq$variable)

toq$variable[toq$variable == "PET"] <- "Precision"

toq$variable[toq$variable == "Endogeneity"] <- "Address Endogeneity"

ord <- c("FAT", "Precision", X)

toq <- toq[order(match(toq$variable, ord)),]

gg <- ggplot(toq[!is.na(toq$value),], aes(x = value)) + geom_density(color = "#004488", lwd = 1) + 
  facet_wrap(~variable, scales = "free", nrow = 6, ncol = 4)  +
  scale_x_continuous(expand = c(0,0),limits = c(-1,1)) + theme_bw() 
gg <- gg +  theme(panel.grid = element_blank(), text = element_text(size=10))

gg <- gg + ylab("") + xlab(label ="") + theme(text = element_text(size = 10))

tiff(file="denscovary.tiff",width=600,height=700,res=100, type = "cairo")
gg
dev.off()




#make table of means and medians


x <- sapply(tmp[X], mean, na.rm = T)
y <- sapply(tmp[X], median, na.rm = T)
z <- sapply(tmp[X], sd, na.rm = T)

covarytable <- cbind(x,y,z)

write.csv(covarytable, "covarytable.csv")



#==========================================

#table of covariates used

usedtable <- c(paste0(tmp$frml[length(tmp$frml)]),
               paste0(tmpquas$frml[length(tmpquas$frml)]),
               paste0(   tmpnonquas$frml[length(tmpnonquas$frml)]),
               paste0(   tmparea$frml[length(tmparea$frml)]),
               paste0(  tmpind$frml[length(tmpind$frml)]),
               paste0(  tmphomicide$frml[length(tmphomicide$frml)]),
               paste0(  tmpviolent$frml[length(tmpviolent$frml)]),
               paste0( tmpnonviolent$frml[length(tmpnonviolent$frml)]),
               paste0( tmpELAS$frml[length(tmpELAS$frml)]),
               paste0( tmpELASENDO$frml[length(tmpELASENDO$frml)])
               )

name <- unique(bigtmp$name)

usedtable <- cbind(usedtable, name)

write.csv(usedtable, "used.table.csv")




