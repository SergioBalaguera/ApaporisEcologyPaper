##=========================================##
##Population ecology paper Apaporis caiman+##
##=========================================##
##StepbyStep##

## Clean----
rm(list=ls())
graphics.off()

#Libraries----
library(dplyr) ## A Grammar of Data Manipulation
library(readxl)
library(lattice)
library(unmarked)
library(AICcmodavg)
library(AHMbook)
library(RColorBrewer) #color <- brewer.pal(5, "Set1")
library(MuMIn)

#Barplots----
AllData <- read_xlsx("Workbook.xlsx", sheet = "All Caiman")

tiff("Model_Figs/fig2.tif", width=250, height=200, units="mm", res=500, compression="lzw")
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
par(mar=c(4, 4, 2, 1), cex = 0.8)

#Barpot Observations per Site per Month
color2 = c("lemonchiffon", "khaki1", "lemonchiffon3", "khaki3", "khaki4")
BP2 <- table(AllData$`Especific locality`, AllData$Month)
barplot(BP2, ylab = "Number of observations", ylim = c(0,250), xlab = "Months", 
        col = color2)
legend("topright", legend = c("Inana Lagoon", "Churuco Lagoon", "Arriba Lagoon", 
                              "Apaporis River T3", "Apaporis River T1"), fill = color2, 
       title = "Transects",  box.lty=0, cex = 1)

#Barplot Observations per SizeClass
color1 = c("lemonchiffon", "khaki1", "lemonchiffon3", "khaki3")
BP3 <- table(AllData$`Size class`)
barplot(BP3, ylab = "Number of observations", ylim = c(0,600), xlab = "Size Class", 
        col = color1)

#Barplot Observations per Size Class per Month
color1 = c("lemonchiffon", "khaki1", "lemonchiffon3", "khaki3")
BP1 <- table(AllData$`Size class`, AllData$Month)
barplot(BP1, ylab = "Number of observations", ylim = c(0,250), xlab = "Months", 
        col = color1)
legend("topright", legend = c("I", "II", "III", "IV"), fill = color1, 
       title = "Size class",  box.lty=0, cex = 1)

dev.off()

## Analysis fisical variables----
(FV = read_xlsx("Workbook.xlsx", sheet = "Fisical variables"))
mean(FV$`Air Temp S`)
sd(FV$`Air Temp S`)
mean(FV$`Air Temp E`)
sd(FV$`Air Temp E`)

mean(FV$`Water Temp S`)
sd(FV$`Water Temp S`)
mean(FV$`Water Temp E`)
sd(FV$`Water Temp E`)

mean(FV$`Relative Humidity S`)
sd(FV$`Relative Humidity S`)
mean(FV$`Relative Humidity E`)
sd(FV$`Relative Humidity E`)

TempA = FV[,1:2]
TempA = data.frame(TempA)
t.test(TempA$Air.Temp.S, TempA$Air.Temp.E, var.equal = T)

TempW = FV[,3:4]
TempW = data.frame(TempW)
t.test(TempW$Water.Temp.S, TempW$Water.Temp.E, var.equal = T)

ReHum = FV[,5:6]
ReHum = data.frame(ReHum)
t.test(ReHum$Relative.Humidity.S, ReHum$Relative.Humidity.E, var.equal = T)

##Colinearity----
Co.co<- read_xlsx("Workbook.xlsx", sheet = "Abundance_R")
(Co.co <- as.data.frame(Co.co))
C.lm <- lm(Co.co$`Air Temp Mean` ~  Co.co$`Water Temp Mean`)
C.lm1 <- lm(Co.co$`Air Temp Mean` ~  Co.co$`Relative Humidity Mean`)
C.lm2 <- lm(Co.co$`Water Temp Mean` ~  Co.co$`Relative Humidity Mean`)

summary(C.lm)$r.squared
summary(C.lm1)$r.squared
summary(C.lm2)$r.squared

##Spatial autocorrelation analysis----
library(ape)
#Loading data Segments
SurveysSA <- read.csv("Abundance_S.csv", header = T, sep = ",")
SurveysSA

SurveysSA.dists <- as.matrix(dist(cbind(SurveysSA$Lat, SurveysSA$Long)))
SurveysSA.dists.inv <- 1/SurveysSA.dists
diag(SurveysSA.dists.inv) <- 0
SurveysSA.dists.inv[is.infinite(SurveysSA.dists.inv)] <- 0
SurveysSA.dists.inv[1:5, 1:5]
set.seed(4587)
Moran.I(SurveysSA$Crocs, SurveysSA.dists.inv)

SurveysSA.dists.bin <- (SurveysSA.dists > 0 & SurveysSA.dists <= 0.75)
Moran.I(SurveysSA$Crocs, SurveysSA.dists.bin)

#Loading data Transcts
SurveysSAT <- read.csv("Abundance_T.csv", header = T, sep = ",")
SurveysSAT

SurveysSAT.dists <- as.matrix(dist(cbind(SurveysSAT$Lat, SurveysSAT$Long)))
SurveysSAT.dists.inv <- 1/SurveysSAT.dists
diag(SurveysSAT.dists.inv) <- 0
SurveysSAT.dists.inv[is.infinite(SurveysSAT.dists.inv)] <- 0
SurveysSAT.dists.inv

set.seed(4587)
Moran.I(SurveysSAT$Crocs, SurveysSAT.dists.inv)

SurveysSAT.dists.bin <- (SurveysSAT.dists > 0 & SurveysSAT.dists <= 1000)
Moran.I(SurveysSAT$Crocs, SurveysSAT.dists.bin)

## N-mixed models----
##By Segments----
#Loading and setting data by segments
Surveys <- read.csv("Abundance_S.csv", header = T, sep = ",")
Surveys$Mm.tr <- paste(Surveys$Month, Surveys$Transect, sep="_")
Surveys <- Surveys %>%
  mutate(Route_name = factor(Route_name), 
         Month_Text = factor(Month_Text))
#Center and scale environmental variables
at.mn <- mean(Surveys$Air.Temp.Mean)
at.sd <- sd(Surveys$Air.Temp.Mean)
wt.mn <- mean(Surveys$Water.Temp.Mean)
wt.sd <- sd(Surveys$Water.Temp.Mean)
hum.mn <- mean(Surveys$Relative.Humidity.Mean)
hum.sd <- sd(Surveys$Relative.Humidity.Mean)
mt.mn <- mean(Surveys$Month)
mt.sd <- sd(Surveys$Month)

Surveys.sc <- Surveys %>%
  mutate(Air.Temp.Mean = (Air.Temp.Mean-at.mn)/at.sd,
         Water.Temp.Mean = (Water.Temp.Mean-wt.mn)/wt.sd,
         Relative.Humidity.Mean = (Relative.Humidity.Mean-hum.mn)/hum.sd,
         Month = (Month-mt.mn)/mt.sd)
#Model data preparation
#Count data
y <- tapply(Surveys.sc$Crocs, list(Surveys.sc$Mm.tr, Surveys.sc$Segment), max, 
            na.rm=T)
#Site Covariates
site.cov <- Surveys.sc %>%
  select(Mm.tr, Month, Transect, Hab) %>% 
  group_by(Mm.tr) %>%
  summarise(transect = getElement(Transect, 1),
            hab = getElement(Hab, 1),
            month = getElement(Month, 1))

(site.cov <- data.frame(site.cov))
#Observation Covariates
obs.cov <- list()
obs.cov.nms <- c("Month", "Air.Temp.Mean", "Water.Temp.Mean", 
                 "Relative.Humidity.Mean")
obs.cov.nms <- obs.cov.nms[order(obs.cov.nms)]
for (i in 1:length(obs.cov.nms)){
  obs.cov[[i]] <- tapply(Surveys.sc[,obs.cov.nms[i]], list(Surveys.sc$Mm.tr, 
                                                           Surveys.sc$Segment), getElement, 1)
  names(obs.cov)[i] <- obs.cov.nms[i]
}

obs.cov

pc.data <- unmarkedFramePCount(y=y, siteCovs = site.cov,  obsCovs = obs.cov)
summary(pc.data)

pc.data <- as(pc.data, "data.frame")
write.csv(pc.data, "Model_Tables/pc.data.csv", row.names=F)
str(pc.data)
#Model fitting
starts.P <- pcount(~ Month + Air.Temp.Mean + Water.Temp.Mean + Relative.Humidity.Mean
                   ~ transect + month + hab, pc.data, mixture = "P", 
                   engine = "C", se = F); beepr::beep(11)
m.full.P <- pcount(~ Month + Air.Temp.Mean + Water.Temp.Mean + Relative.Humidity.Mean
                   ~ transect + month + hab, pc.data, mixture="P", engine="C", 
                   starts=coef(starts.P)); beepr::beep(11)

starts.ZIP <- pcount(~ Month + Air.Temp.Mean + Water.Temp.Mean  + Relative.Humidity.Mean
                     ~ transect + month + hab, pc.data, mixture="ZIP", se=F); beepr::beep(11)
m.full.ZIP <- pcount(~ Month + Air.Temp.Mean + Water.Temp.Mean + Relative.Humidity.Mean
                     ~ transect + month + hab, pc.data, mixture="ZIP", starts=coef(starts.ZIP)); beepr::beep(11)

starts.NB <- pcount(~ Month + Air.Temp.Mean + Water.Temp.Mean + Relative.Humidity.Mean
                    ~ transect + month + hab, pc.data, mixture="NB", se=F); beepr::beep(11)
m.full.NB <- pcount(~ Month + Air.Temp.Mean + Water.Temp.Mean + Relative.Humidity.Mean
                    ~ transect + month + hab, pc.data, mixture="NB", starts=coef(starts.NB)); beepr::beep(11)
#Evaluate
(full.mods <- fitList("Pois" = m.full.P,
                      "ZIP" = m.full.ZIP,
                      "NB" = m.full.NB))
#AIC table
(aic.models.table.S <- modSel(full.mods))
(aic.models.table.S <- as(aic.models.table.S, "data.frame"))
write.csv(aic.models.table.S, "Model_Tables/aic.models.table.S.csv", row.names=F)
#Goodness of Fit
system.time(gof.P <- Nmix.gof.test(m.full.P, engine = "C", nsim=10000, #run 100 to test, 10000 for final model
                                   parallel = F)); beepr::beep(11)

system.time(gof.ZIP <- Nmix.gof.test(m.full.ZIP, engine = "C", nsim=10000, #run 100 to test, 10000 for final model
                                     parallel = F)); beepr::beep(11)

system.time(gof.NB <- Nmix.gof.test(m.full.NB, engine = "C", nsim=10000, #run 100 to test, 10000 for final model
                                    parallel = F)); beepr::beep(11)

gof.P ; gof.ZIP ; gof.NB
#Examine residuals
plot_Nmix_resi(fmP=m.full.P, fmNB=m.full.NB, fmZIP=m.full.ZIP)
#RMSE
(RMSE.P <- sqrt(mean((y-fitted(m.full.P))^2, na.rm=T)))
(RMSE.ZIP <- sqrt(mean((y-fitted(m.full.ZIP))^2, na.rm=T)))
(RMSE.NB <- sqrt(mean((y-fitted(m.full.NB))^2, na.rm=T)))
##Models (Note: formula is for detection, then abundance)-
m0 <- pcount(~ 1 
             ~ 1, pc.data,  mixture="NB", engine="C"); beepr::beep(11)
m1 <- pcount(~ Month 
             ~ 1, pc.data,  mixture="NB", engine="C"); beepr::beep(11)
m2 <- pcount(~ Air.Temp.Mean 
             ~ 1, pc.data,  mixture="NB", engine="C"); beepr::beep(11)
m3 <- pcount(~ Water.Temp.Mean 
             ~ 1, pc.data,  mixture="NB", engine="C"); beepr::beep(11)
m4 <- pcount(~ Relative.Humidity.Mean 
             ~ 1, pc.data,  mixture="NB", engine="C"); beepr::beep(11)
m5 <- pcount(~ Month + Air.Temp.Mean 
             ~ 1, pc.data,  mixture="NB", engine="C"); beepr::beep(11)
m6 <- pcount(~ Month + Water.Temp.Mean 
             ~ 1, pc.data,  mixture="NB", engine="C"); beepr::beep(11)
m7 <- pcount(~ Month + Relative.Humidity.Mean 
             ~ 1, pc.data,  mixture="NB", engine="C"); beepr::beep(11)
m8 <- pcount(~ Air.Temp.Mean + Water.Temp.Mean 
             ~ 1, pc.data,  mixture="NB", engine="C"); beepr::beep(11)
m9 <- pcount(~ Air.Temp.Mean + Relative.Humidity.Mean 
             ~ 1, pc.data,  mixture="NB", engine="C"); beepr::beep(11)
m10 <- pcount(~ Water.Temp.Mean + Relative.Humidity.Mean 
              ~ 1, pc.data,  mixture="NB", engine="C"); beepr::beep(11)
m11 <- pcount(~ Month + Air.Temp.Mean + Water.Temp.Mean 
              ~ 1, pc.data,  mixture="NB", engine="C"); beepr::beep(11)
m12 <- pcount(~ Month + Air.Temp.Mean + Relative.Humidity.Mean 
              ~ 1, pc.data, mixture="NB", engine="C"); beepr::beep(11)
m13 <- pcount(~ Air.Temp.Mean + Water.Temp.Mean + Relative.Humidity.Mean 
              ~ 1, pc.data, mixture="NB", engine="C"); beepr::beep(11)
m14 <- pcount(~ Month + Air.Temp.Mean + Water.Temp.Mean + Relative.Humidity.Mean
              ~ 1, pc.data,  mixture="NB", engine="C"); beepr::beep(11)
m15 <- pcount(~ Month + Air.Temp.Mean + Water.Temp.Mean + Relative.Humidity.Mean
              ~ route_name, pc.data,  mixture="NB", engine="C"); beepr::beep(11)
m16 <- pcount(~ Month + Air.Temp.Mean + Water.Temp.Mean + Relative.Humidity.Mean
              ~ hab, pc.data,  mixture="NB", engine="C"); beepr::beep(11)
m17 <- pcount(~ Month + Air.Temp.Mean + Water.Temp.Mean + Relative.Humidity.Mean
              ~ month, pc.data,  mixture="NB", engine="C"); beepr::beep(11)
m18 <- pcount(~ Month + Air.Temp.Mean + Water.Temp.Mean + Relative.Humidity.Mean
              ~ route_name + hab, pc.data,  mixture="NB", engine="C"); beepr::beep(11)
m19 <- pcount(~ Month + Air.Temp.Mean + Water.Temp.Mean + Relative.Humidity.Mean
              ~ route_name + month, pc.data,  mixture="NB", engine="C"); beepr::beep(11)
m20 <- pcount(~ Month + Air.Temp.Mean + Water.Temp.Mean + Relative.Humidity.Mean
              ~ hab + month, pc.data,  mixture="NB", engine="C"); beepr::beep(11)
m21 <- pcount(~ Month + Air.Temp.Mean + Water.Temp.Mean + Relative.Humidity.Mean     #### Full
              ~ route_name + hab + month, pc.data,  mixture="NB", engine="C"); beepr::beep(11)
m22 <- pcount(~ Month + Air.Temp.Mean + Water.Temp.Mean  
              ~ route_name + hab + month, pc.data,  mixture="NB", engine="C"); beepr::beep(11)
m23 <- pcount(~ Month + Air.Temp.Mean + Relative.Humidity.Mean 
              ~ route_name + hab + month, pc.data,  mixture="NB", engine="C"); beepr::beep(11)
m24 <- pcount(~ Month + Water.Temp.Mean + Relative.Humidity.Mean 
              ~ route_name + hab + month, pc.data,  mixture="NB", engine="C"); beepr::beep(11)
m25 <- pcount(~ Air.Temp.Mean + Water.Temp.Mean + Relative.Humidity.Mean 
              ~ route_name + hab + month, pc.data,  mixture="NB", engine="C"); beepr::beep(11)
m26 <- pcount(~ Month + Air.Temp.Mean 
              ~ route_name + hab + month, pc.data, mixture="NB", engine="C"); beepr::beep(11)
m27 <- pcount(~ Month + Water.Temp.Mean 
              ~ route_name + hab + month, pc.data,  mixture="NB", engine="C"); beepr::beep(11)
m28 <- pcount(~ Month + Relative.Humidity.Mean 
              ~ route_name + hab + month, pc.data,  mixture="NB", engine="C"); beepr::beep(11)
m29 <- pcount(~ Air.Temp.Mean + Water.Temp.Mean 
              ~ route_name + hab + month, pc.data, mixture="NB", engine="C"); beepr::beep(11)
m30 <- pcount(~ Air.Temp.Mean + Relative.Humidity.Mean 
              ~ route_name + hab + month, pc.data, mixture="NB", engine="C"); beepr::beep(11)
m31 <- pcount(~ Water.Temp.Mean + Relative.Humidity.Mean 
              ~ route_name + hab + month, pc.data, mixture="NB", engine="C"); beepr::beep(11)
m32 <- pcount(~ Month 
              ~ route_name + hab + month, pc.data,  mixture="NB", engine="C"); beepr::beep(11)
m33 <- pcount(~ Air.Temp.Mean 
              ~ route_name + hab + month, pc.data,  mixture="NB", engine="C"); beepr::beep(11)
m34 <- pcount(~ Water.Temp.Mean 
              ~ route_name + hab + month, pc.data,  mixture="NB", engine="C"); beepr::beep(11)
m35 <- pcount(~ Relative.Humidity.Mean 
              ~ route_name + hab + month, pc.data,  mixture="NB", engine="C"); beepr::beep(11)
m36 <- pcount(~ 1 
              ~ route_name + hab + month, pc.data,  mixture="NB", engine="C"); beepr::beep(11)
m37 <- pcount(~ 1 
              ~ route_name + hab, pc.data, mixture = "NB", engine="C"); beepr::beep(11)
m38 <- pcount(~ 1 
              ~ route_name + month, pc.data, mixture = "NB", engine="C"); beepr::beep(11)
m39 <- pcount(~ 1 
              ~ hab + month, pc.data, mixture = "NB", engine="C"); beepr::beep(11)
m40 <- pcount(~ 1 
              ~ route_name, pc.data, mixture = "NB", engine="C"); beepr::beep(11)
m41 <- pcount(~ 1 
              ~ hab, pc.data, mixture = "NB", engine="C"); beepr::beep(11)
m42 <- pcount(~ 1 
              ~ month, pc.data, mixture = "NB", engine="C"); beepr::beep(11)

mods <- fitList(
  `p(.)lambda(.)` = m0,
  `p(mt)lambda(.)` = m1,
  `p(airt)lambda(.)` = m2,
  `p(wt)lambda(.)` = m3,
  `p(relhum)lambda(.)` = m4,
  `p(mt, airt)lambda(.)` = m5,
  `p(mt, wt)lambda(.)` = m6,
  `p(mt, relhum)lambda(.)` = m7,
  `p(airt, wt)lambda(.)` = m8,
  `p(airt, relhum)lambda(.)` = m9,
  `p(wt, relhum)lambda(.)` = m10,
  `p(mt, airt, wt)lambda(.)` = m11,
  `p(mt, airt, relhum)lambda(.)` = m12,
  `p(airt, wt, relhum)lambda(.)` = m13,
  `p(full)lambda(.)` = m14,
  `p(full)lambda(route_name)` = m15,
  `p(full)lambda(hab)` = m16,
  `p(full)lambda(mt)` = m17,
  `p(full)lambda(route_name, hab)` = m18,
  `p(full)lambda(route_name, mt)` = m19,
  `p(full)lambda(hab, mt)` = m20,
  `p(full)lambda(full)` = m21,
  `p(mt, airt, wt)lambda(full)` = m22,
  `p(mt, airt, relhum)lambda(full)` = m23,
  `p(mt, wt, relhum)lambda(full)` = m24,
  `p(airt, wt, relhum)lambda(full)` = m25,
  `p(mt, airt)lambda(full)` = m26,
  `p(mt, wt)lambda(full)` = m27,
  `p(mt, relhum)lambda(full)` = m28,
  `p(airt, wt)lambda(full)` = m29,
  `p(airt, relhum)lambda(full)` = m30,
  `p(wt, relhum)lambda(full)` = m31,
  `p(mt)lambda(full)` = m32,
  `p(airt)lambda(full)` = m33,
  `p(wt)lambda(full)` = m34,
  `p(relhum)lambda(full)` = m35,
  `p(.)lambda(full)` = m36,
  `p(.)lambda(route_name, hab)` = m37,
  `p(.)lambda(route_name, mt)` = m38,
  `p(.)lambda(hab, mt)` = m39,
  `p(.)lambda(route_name)` = m40,
  `p(.)lambda(hab)` = m41,
  `p(.)lambda(mt)` = m42
)

aic.table.all.S <- modSel(mods, nullmod="p(.)lambda(.)")
(aic.table.all.S <- as(aic.table.all.S, "data.frame"))
write.csv(aic.table.all.S, "Model_Tables/aic_table.all.S.csv", row.names=F)

#Best model and general model-GoF figure
palette(brewer.pal(n=3, name="Dark2"))
system.time(gof.NB <- Nmix.gof.test(m.full.NB, engine = "C", nsim=100, #run 100 to test, 10000 for final model
                                    parallel = F)); beepr::beep(11)
system.time(gof.m20 <- Nmix.gof.test(m20, engine = "C", nsim=100, #run 100 to test, 10000 for final model
                                    parallel = F)); beepr::beep(11)
gof.m20

#Figures N-mixture model
#Detection probabilities
m20
(a <- predict(m20, type = "det"))
mean(a$Predicted)
mean(a$lower)
mean(a$upper)

#Color palette
palette(brewer.pal(n=3, name="Dark2"))
tiff("Model_Figs/p_full.tif", width=190, height=100, units="mm", res=500, compression="lzw")
par(mfrow=c(1,4), mar=c(4, 4, 2, 1), cex = 0.7)

#Month
pred.df <- data.frame(mt.un = seq(1, 12, length.out = 58), airt.un=at.mn, wt.un=wt.mn, 
                       rh.un=hum.mn)  %>%
    mutate(Month = (mt.un-mt.mn)/mt.sd,
           Air.Temp.Mean = (airt.un-at.mn)/at.sd,
           Water.Temp.Mean = (wt.un-wt.mn)/wt.sd,
           Relative.Humidity.Mean = (rh.un - hum.mn)/hum.sd)

preds <- predict(m20, pred.df, type="det")
(pred.df <- cbind(pred.df, preds))

pred.df.m20.mt <- as(pred.df, "data.frame")
write.csv(pred.df.m20.mt, "Model_Tables/pred.df.m20.mt.csv", row.names=F)

plot(Predicted~mt.un, pred.df, type="l", lwd=3, ylim=c(0,1), 
     xlab=expression("Month"), 
     ylab="Detection Probability")

lines(lower~mt.un, pred.df, lty=2)
lines(upper~mt.un, pred.df, lty=2)

#Air temp
pred.df <- data.frame(mt.un=mt.mn, airt.un=seq(20, 30, length.out = 58),  wt.un=wt.mn, 
                       rh.un=hum.mn) %>%
  mutate(Month = (mt.un-mt.mn)/mt.sd,
         Air.Temp.Mean = (airt.un-at.mn)/at.sd,
         Water.Temp.Mean = (wt.un-wt.mn)/wt.sd,
         Relative.Humidity.Mean = (rh.un - hum.mn)/hum.sd,
  )

preds <- predict(m20, pred.df, type="det")
pred.df <- cbind(pred.df, preds)

pred.df.m20.airt <- as(pred.df, "data.frame")
write.csv(pred.df.m20.airt, "Model_Tables/pred.df.m20.airt.csv", row.names=F)

plot(Predicted~airt.un, pred.df, type="l", lwd=3, ylim=c(0,1), 
    xlab=expression("Air Temperature"~(degree*C)),
    ylab="")

lines(lower~airt.un, pred.df, lty=2)
lines(upper~airt.un, pred.df, lty=2)

#Water temp
pred.df <- data.frame(mt.un=mt.mn, airt.un=at.mn,  wt.un=seq(20, 30, length.out = 58), 
                       rh.un=hum.mn) %>%
    mutate(Month = (mt.un-mt.mn)/mt.sd,
           Air.Temp.Mean = (airt.un-at.mn)/at.sd,
           Water.Temp.Mean = (wt.un-wt.mn)/wt.sd,
           Relative.Humidity.Mean = (rh.un - hum.mn)/hum.sd,
    )

preds <- predict(m20, pred.df, type="det")
pred.df <- cbind(pred.df, preds)

pred.df.m20.wt <- as(pred.df, "data.frame")
write.csv(pred.df.m20.wt, "Model_Tables/pred.df.m20.wt.csv", row.names=F)

plot(Predicted~wt.un, pred.df, type="l", lwd=3, ylim=c(0,1), 
     xlab=expression("Water Temperature"~(degree*C)),
     ylab="")

lines(lower~wt.un, pred.df, lty=2)
lines(upper~wt.un, pred.df, lty=2)

#Relative humidity
pred.df <- data.frame(mt.un=mt.mn, airt.un=at.mn,  wt.un=wt.mn, 
                       rh.un=seq(60, 100, length.out = 58)) %>%
    mutate(Month = (mt.un-mt.mn)/mt.sd,
           Air.Temp.Mean = (airt.un-at.mn)/at.sd,
           Water.Temp.Mean = (wt.un-wt.mn)/wt.sd,
           Relative.Humidity.Mean = (rh.un - hum.mn)/hum.sd,
    )

preds <- predict(m20, pred.df, type="det")
pred.df <- cbind(pred.df, preds)

pred.df.m20.rh <- as(pred.df, "data.frame")
write.csv(pred.df.m20.rh, "Model_Tables/pred.df.m20.rh.csv", row.names=F)

plot(Predicted~rh.un, pred.df, type="l", lwd=3, ylim=c(0,1), 
     xlab=expression("Relative humidity (%)"),
     ylab="")

lines(lower~rh.un, pred.df, lty=2)
lines(upper~rh.un, pred.df, lty=2)

dev.off()

#Abundance Model predictions
#lambda
m20
(Nhat.pred <- predict(m20, type = "state"))
(Nhat <- sum(Nhat.pred$Predicted))
(Nhat <- sum(Nhat.pred$lower))
(Nhat <- sum(Nhat.pred$upper))

average.density.m16 <- predict(m16, type="state")
ave.dens.m16 <- mean(average.density.m16$Predicted)

average.density.m16 <- as(average.density.m16, "data.frame")
write.csv(average.density.m16, "Model_Tables/average.density.m16.csv", row.names=F)

rts <- factor(Surveys.sc$Hab)
pred.dat <- data.frame(Hab=factor("1", levels=levels(rts))) %>%
  mutate(hab = Hab)
(predict(m16, pred.dat, type="state"))

rts <- factor(Surveys.sc$Hab)
pred.dat <- data.frame(Hab=factor("2", levels=levels(rts))) %>%
  mutate(hab = Hab)
(predict(m16, pred.dat, type="state"))

#lambda (route_name) m15
average.density.m15 <- predict(m15, type="state")
mean(average.density.m15$Predicted)

rts <- factor(Surveys.sc$Route_name)
pred.dat <- data.frame(Route_name=factor("Apaporis River T1", levels=levels(rts))) %>%
  mutate(route_name = Route_name)
(predict(m15, pred.dat, type="state"))

  rts <- factor(Surveys.sc$Route_name)
pred.dat <- data.frame(Route_name=factor("Inana Lagoon", levels=levels(rts))) %>%
  mutate(route_name = Route_name)
(predict(m15, pred.dat, type="state"))

rts <- factor(Surveys.sc$Route_name)
pred.dat <- data.frame(Route_name=factor("Apaporis River T3", levels=levels(rts))) %>%
  mutate(route_name = Route_name)
(predict(m15, pred.dat, type="state"))

rts <- factor(Surveys.sc$Route_name)
pred.dat <- data.frame(Route_name=factor("Churuco Lagoon", levels=levels(rts))) %>%
  mutate(route_name = Route_name)
(predict(m15, pred.dat, type="state"))

rts <- factor(Surveys.sc$Route_name)
pred.dat <- data.frame(Route_name=factor("Arriba Lagoon", levels=levels(rts))) %>%
  mutate(route_name = Route_name)
(predict(m15, pred.dat, type="state"))

## By Transects----
set.seed(1985)
#Loading and setting data by transects
Data <- read.csv("Abundance_T.csv", header = T, sep = ",")
as.factor(Data$Hab)
Data$Mt.Tr <- paste(Data$Month, Data$Transect, sep = "_")
Data$ilength <- paste(1/Data$Transect.length)
Data$ilength <- as.numeric(as.character(Data$ilength))
str(Data)

mean(Data$Crocs)
var(Data$Crocs)

#Center and scale environmental variables
atmn <- mean(Data$Air.Temp.Mean)
atsd <- sd(Data$Air.Temp.Mean)
wtmn <- mean(Data$Water.Temp.Mean)
wtsd <- sd(Data$Water.Temp.Mean)
hummn <- mean(Data$Relative.Humidity.Mean)
humsd <- sd(Data$Relative.Humidity.Mean)
mtmn <- mean(Data$Month)
mtsd <- sd(Data$Month)
at.wtmn <- mean(Data$AT.WT)
at.wtsd <- sd(Data$AT.WT)

Data.sc <- Data %>%
  mutate(Air.Temp.Mean = (Air.Temp.Mean-atmn)/atsd,
         Water.Temp.Mean = (Water.Temp.Mean-wtmn)/wtsd,
         Relative.Humidity.Mean = (Relative.Humidity.Mean-hummn)/humsd,
         Month = (Month-mtmn)/mtsd,
         AT.WT = (AT.WT-at.wtmn)/at.wtsd)

#Count data
Counts <- tapply(Data.sc$Crocs, list(Data.sc$Mt.Tr), max, 
                 na.rm=T)
Counts <- data.frame(Counts)

#Site Covariates
SC <- list()
SC.nms <- c("Month", "Transect", "Hab", "ilength")
SC.nms <- SC.nms[order(SC.nms)]
for (i in 1:length(SC.nms)){
  SC[[i]] <- tapply(Data.sc[,SC.nms[i]], list(Data.sc$Mt.Tr), getElement, 1)
  names(SC)[i] <- SC.nms[i]
}

SC <- data.frame(SC)

#Observation Covariates
OC <- list()
OC.nms <- c("Month", "Air.Temp.Mean", "Water.Temp.Mean", 
            "Relative.Humidity.Mean", "AT.WT")
OC.nms <- OC.nms[order(OC.nms)]
for (i in 1:length(OC.nms)){
  OC[[i]] <- tapply(Data.sc[,OC.nms[i]], list(Data.sc$Mt.Tr), getElement, 1)
  names(OC)[i] <- OC.nms[i]
}

OC <- data.frame(OC)

pcountdata <- unmarkedFramePCount(y=Counts, siteCovs = SC, obsCovs = OC)
summary(pcountdata)

pcountdata <- as(pcountdata, "data.frame")
write.csv(pcountdata, "Model_Tables/Transects/pcountdata.csv", row.names=F)

#Model fitting
#Likelihood evaluation
fm <- pcount(~ Month + Air.Temp.Mean + Water.Temp.Mean + Relative.Humidity.Mean
             ~ Month + ilength + Hab, pcountdata, control=list(trace=TRUE, REPORT=5), 
             engine = "C", se = F); beepr::beep(11)
summary(fm) ; fm@AIC
fm.k500 <- pcount(~ Month + Air.Temp.Mean + Water.Temp.Mean + Relative.Humidity.Mean
                  ~ Hab + Month + ilength, pcountdata, control=list(trace=T, REPORT=5), K = 500, 
                  engine = "C", se = F); beepr::beep(11)
summary(fm.k500) ; fm.k500@AIC

#Defining covariate structure
fm <- pcount(~ Month * Air.Temp.Mean * Water.Temp.Mean * Relative.Humidity.Mean
             ~ Month * ilength * Hab, pcountdata,  mixture = "NB", control=list(trace=TRUE, REPORT=5),
             engine = "C", se = F); beepr::beep(11)
fmf <- pcount(~ Month * Air.Temp.Mean * Water.Temp.Mean * Relative.Humidity.Mean
              ~ Month * ilength * Hab, pcountdata, mixture="NB", engine="C", control=list(trace=TRUE, REPORT=5),
              starts=coef(fm)); beepr::beep(11)
fmf
fm1 <- pcount(~ Month * Air.Temp.Mean * Water.Temp.Mean * Relative.Humidity.Mean
              ~ Month + Month:ilength + Month:Hab + Month:ilength:Hab, pcountdata, mixture = "NB", control=list(trace=TRUE, REPORT=5),
              engine = "C", se = F); beepr::beep(11)
fmf1 <- pcount(~ Month * Air.Temp.Mean * Water.Temp.Mean * Relative.Humidity.Mean
               ~ Month + Month:ilength + Month:Hab + Month:ilength:Hab, pcountdata, mixture="NB", engine="C", control=list(trace=TRUE, REPORT=5),
               starts=coef(fm1)); beepr::beep(11)
fmf1
fm2 <- pcount(~ Month + Water.Temp.Mean + Month:Relative.Humidity.Mean + Water.Temp.Mean:Relative.Humidity.Mean
              ~ Month + Month:ilength + Month:Hab + Month:ilength:Hab, pcountdata, mixture = "NB", control=list(trace=TRUE, REPORT=5),
              engine = "C", se = F); beepr::beep(11)
fmf2 <- pcount(~ Month + Water.Temp.Mean + Month:Relative.Humidity.Mean + Water.Temp.Mean:Relative.Humidity.Mean
               ~ Month + Month:ilength + Month:Hab + Month:ilength:Hab, pcountdata, mixture="NB", engine="C", control=list(trace=TRUE, REPORT=5),
               starts=coef(fm2)); beepr::beep(11)
fmf2
fm3 <- pcount(~ Month + Water.Temp.Mean + Month:Relative.Humidity.Mean + Water.Temp.Mean:Relative.Humidity.Mean
              ~ Month + ilength + Month:ilength + Month:Hab + Month:ilength:Hab, pcountdata, mixture = "NB", control=list(trace=TRUE, REPORT=5),
              engine = "C", se = F); beepr::beep(11)
fmf3 <- pcount(~ Month + Water.Temp.Mean + Month:Relative.Humidity.Mean + Water.Temp.Mean:Relative.Humidity.Mean
               ~ Month + ilength +  Month:ilength + Month:Hab + Month:ilength:Hab, pcountdata, mixture="NB", engine="C", control=list(trace=TRUE, REPORT=5),
               starts=coef(fm3)); beepr::beep(11)
fmf3

#Defining variance structure
Ps <- pcount(fmf3@formula, pcountdata, mixture = "P", engine = "C", control=list(trace=TRUE, REPORT=5), 
             se = F); beepr::beep(11)
Pfull <- pcount(fmf3@formula, pcountdata, mixture= "P", engine= "C", starts=coef(Ps), 
                control=list(trace=TRUE, REPORT=5)); beepr::beep(11)

ZIPs <- pcount(fmf3@formula, pcountdata, mixture= "ZIP", engine = "C", control=list(trace=TRUE, REPORT=5), 
               se=F); beepr::beep(11)
ZIPfull <- pcount(fmf3@formula, pcountdata, mixture= "ZIP", engine= "C", starts=coef(ZIPs), 
                  control=list(trace=TRUE, REPORT=5)); beepr::beep(11)

NBs <- pcount(fmf3@formula, pcountdata, mixture= "NB", engine= "C", control=list(trace=TRUE, REPORT=5), 
              se=F); beepr::beep(11)
NBfull <- pcount(fmf3@formula, pcountdata, mixture= "NB", engine= "C", starts=coef(NBs), 
                 control=list(trace=TRUE, REPORT=5)); beepr::beep(11)
#Evaluate
(full.models <- fitList("Pois" = Pfull,
                        "ZIP" = ZIPfull,
                        "NB" = NBfull))

#AIC table
(aic.models.table.T <- modSel(full.models))
aic.models.table.T <- as(aic.models.table.T, "data.frame")
write.csv(aic.models.table.T, "Model_Tables/Transects/aic.models.table.T.csv", row.names=F)

##General p and lamda values
Pfulllambda <- predict(Pfull, type = "state")
mean(Pfulllambda$Predicted)
PfullP <- predict(Pfull, type = "det")
mean(PfullP$Predicted)

ZIPfulllambda <- predict(ZIPfull, type = "state")
mean(ZIPfulllambda$Predicted)
ZIPfullP <- predict(ZIPfull, type = "det")
mean(ZIPfullP$Predicted)

NBfulllambda <- predict(NBfull, type = "state")
mean(NBfulllambda$Predicted)
NBfullP <- predict(NBfull, type = "det")
mean(NBfullP$Predicted)

#Goodness of Fit
system.time(Pgof <- Nmix.gof.test(Pfull, engine = "C", nsim=10000, #run 100 to test, 10000 for final model
                                  parallel = F)); beepr::beep(11)

system.time(ZIPgof <- Nmix.gof.test(ZIPfull, engine = "C", nsim=10000, #run 100 to test, 10000 for final model
                                    parallel = F)); beepr::beep(11)

system.time(NBgof <- Nmix.gof.test(NBfull, engine = "C", nsim=10000, #run 100 to test, 10000 for final model
                                   parallel = F)); beepr::beep(11)

Pgof ; ZIPgof ; NBgof

#Examine residuals
plot_Nmix_resi(fmP=Pfull, fmZIP=ZIPfull, fmNB=NBfull)

##Models (Note: formula is for detection, then abundance)
m1 <- pcount(~ 1 
             ~ 1, pcountdata,  mixture="NB", engine="C"); beepr::beep(11)
m2 <- pcount(~ Month 
             ~ 1, pcountdata,  mixture="NB", engine="C"); beepr::beep(11)
m3 <- pcount(~ Water.Temp.Mean 
             ~ 1, pcountdata,  mixture="NB", engine="C"); beepr::beep(11)
m4 <- pcount(~ Month:Relative.Humidity.Mean  
             ~ 1, pcountdata,  mixture="NB", engine="C"); beepr::beep(11)
m5 <- pcount(~ Water.Temp.Mean:Relative.Humidity.Mean  
             ~ 1, pcountdata,  mixture="NB", engine="C"); beepr::beep(11)
m6 <- pcount(~ Month + Water.Temp.Mean 
             ~ 1, pcountdata,  mixture="NB", engine="C"); beepr::beep(11)
m7 <- pcount(~ Month + Month:Relative.Humidity.Mean 
             ~ 1, pcountdata,  mixture="NB", engine="C"); beepr::beep(11)
m8 <- pcount(~ Month + Water.Temp.Mean:Relative.Humidity.Mean  
             ~ 1, pcountdata,  mixture="NB", engine="C"); beepr::beep(11)
m9 <- pcount(~ Water.Temp.Mean + Month:Relative.Humidity.Mean 
             ~ 1, pcountdata,  mixture="NB", engine="C"); beepr::beep(11)
m10 <- pcount(~ Water.Temp.Mean + Water.Temp.Mean:Relative.Humidity.Mean  
              ~ 1, pcountdata,  mixture="NB", engine="C"); beepr::beep(11)
m11 <- pcount(~ Month:Relative.Humidity.Mean + Water.Temp.Mean:Relative.Humidity.Mean  
              ~ 1, pcountdata,  mixture="NB", engine="C"); beepr::beep(11)
m12 <- pcount(~ Month + Water.Temp.Mean + Month:Relative.Humidity.Mean  
              ~ 1, pcountdata,  mixture="NB", engine="C"); beepr::beep(11)
m13 <- pcount(~ Month + Water.Temp.Mean + Water.Temp.Mean:Relative.Humidity.Mean  
              ~ 1, pcountdata,  mixture="NB", engine="C"); beepr::beep(11)
m14 <- pcount(~ Water.Temp.Mean + Month:Relative.Humidity.Mean + Water.Temp.Mean:Relative.Humidity.Mean     
              ~ 1, pcountdata,  mixture="NB", engine="C"); beepr::beep(11)
m15 <- pcount(~ Month + Water.Temp.Mean + Month:Relative.Humidity.Mean + Water.Temp.Mean:Relative.Humidity.Mean  ## Full p
              ~ 1, pcountdata,  mixture="NB", engine="C"); beepr::beep(11)
m16 <- pcount(~ Month + Water.Temp.Mean + Month:Relative.Humidity.Mean + Water.Temp.Mean:Relative.Humidity.Mean
              ~ Month, pcountdata,  mixture="NB", engine="C"); beepr::beep(11)
m17 <- pcount(~ Month + Water.Temp.Mean + Month:Relative.Humidity.Mean + Water.Temp.Mean:Relative.Humidity.Mean
              ~ ilength, pcountdata,  mixture="NB", engine="C"); beepr::beep(11)
m18 <- pcount(~ Month + Water.Temp.Mean + Month:Relative.Humidity.Mean + Water.Temp.Mean:Relative.Humidity.Mean
              ~ Month:ilength, pcountdata,  mixture="NB", engine="C"); beepr::beep(11)
m19 <- pcount(~ Month + Water.Temp.Mean + Month:Relative.Humidity.Mean + Water.Temp.Mean:Relative.Humidity.Mean
              ~ Month:Hab, pcountdata,  mixture="NB", engine="C"); beepr::beep(11)
m20 <- pcount(~ Month + Water.Temp.Mean + Month:Relative.Humidity.Mean + Water.Temp.Mean:Relative.Humidity.Mean
              ~ Month:ilength:Hab, pcountdata,  mixture="NB", engine="C"); beepr::beep(11)
m21 <- pcount(~ Month + Water.Temp.Mean + Month:Relative.Humidity.Mean + Water.Temp.Mean:Relative.Humidity.Mean
              ~ Month + ilength, pcountdata,  mixture="NB", engine="C"); beepr::beep(11)
m22 <- pcount(~ Month + Water.Temp.Mean + Month:Relative.Humidity.Mean + Water.Temp.Mean:Relative.Humidity.Mean
              ~ Month + Month:ilength, pcountdata, mixture="NB", engine="C"); beepr::beep(11)
m23 <- pcount(~ Month + Water.Temp.Mean + Month:Relative.Humidity.Mean + Water.Temp.Mean:Relative.Humidity.Mean
              ~ Month + Month:Hab, pcountdata, mixture="NB", engine="C"); beepr::beep(11)
m24 <- pcount(~ Month + Water.Temp.Mean + Month:Relative.Humidity.Mean + Water.Temp.Mean:Relative.Humidity.Mean
              ~ Month + Month:ilength:Hab, pcountdata, mixture="NB", engine="C"); beepr::beep(11)
m25 <- pcount(~ Month + Water.Temp.Mean + Month:Relative.Humidity.Mean + Water.Temp.Mean:Relative.Humidity.Mean
              ~ ilength + Month:ilength, pcountdata,  mixture="NB", engine="C"); beepr::beep(11)
m26 <- pcount(~ Month + Water.Temp.Mean + Month:Relative.Humidity.Mean + Water.Temp.Mean:Relative.Humidity.Mean
              ~ ilength + Month:Hab, pcountdata,  mixture="NB", engine="C"); beepr::beep(11)
m27 <- pcount(~ Month + Water.Temp.Mean + Month:Relative.Humidity.Mean + Water.Temp.Mean:Relative.Humidity.Mean
              ~ ilength + Month:ilength:Hab, pcountdata,  mixture="NB", engine="C"); beepr::beep(11)
m28 <- pcount(~ Month + Water.Temp.Mean + Month:Relative.Humidity.Mean + Water.Temp.Mean:Relative.Humidity.Mean
              ~ Month:ilength + Month:Hab, pcountdata,  mixture="NB", engine="C"); beepr::beep(11)
m29 <- pcount(~ Month + Water.Temp.Mean + Month:Relative.Humidity.Mean + Water.Temp.Mean:Relative.Humidity.Mean
              ~ Month:ilength + Month:ilength:Hab, pcountdata,  mixture="NB", engine="C"); beepr::beep(11)
m30 <- pcount(~ Month + Water.Temp.Mean + Month:Relative.Humidity.Mean + Water.Temp.Mean:Relative.Humidity.Mean
              ~ Month:Hab + Month:ilength:Hab, pcountdata,  mixture="NB", engine="C"); beepr::beep(11)
m31 <- pcount(~ Month + Water.Temp.Mean + Month:Relative.Humidity.Mean + Water.Temp.Mean:Relative.Humidity.Mean
              ~ Month + ilength + Month:ilength, pcountdata,  mixture="NB", engine="C"); beepr::beep(11)
m32 <- pcount(~ Month + Water.Temp.Mean + Month:Relative.Humidity.Mean + Water.Temp.Mean:Relative.Humidity.Mean
              ~ Month + ilength + Month:Hab, pcountdata,  mixture="NB", engine="C"); beepr::beep(11)
m33 <- pcount(~ Month + Water.Temp.Mean + Month:Relative.Humidity.Mean + Water.Temp.Mean:Relative.Humidity.Mean
              ~ Month + ilength + Month:ilength:Hab, pcountdata,  mixture="NB", engine="C"); beepr::beep(11)
m34 <- pcount(~ Month + Water.Temp.Mean + Month:Relative.Humidity.Mean + Water.Temp.Mean:Relative.Humidity.Mean
              ~ ilength + Month:ilength + Month:Hab, pcountdata,  mixture="NB", engine="C"); beepr::beep(11)
m35 <- pcount(~ Month + Water.Temp.Mean + Month:Relative.Humidity.Mean + Water.Temp.Mean:Relative.Humidity.Mean
              ~ ilength + Month:ilength + Month:ilength:Hab, pcountdata,  mixture="NB", engine="C"); beepr::beep(11)
m36 <- pcount(~ Month + Water.Temp.Mean + Month:Relative.Humidity.Mean + Water.Temp.Mean:Relative.Humidity.Mean
              ~ Month:ilength + Month:Hab + Month:ilength:Hab, pcountdata,  mixture="NB", engine="C"); beepr::beep(11)
m37 <- pcount(~ Month + Water.Temp.Mean + Month:Relative.Humidity.Mean + Water.Temp.Mean:Relative.Humidity.Mean
              ~ Month + ilength + Month:ilength + Month:Hab, pcountdata,  mixture="NB", engine="C"); beepr::beep(11)
m38 <- pcount(~ Month + Water.Temp.Mean + Month:Relative.Humidity.Mean + Water.Temp.Mean:Relative.Humidity.Mean
              ~ Month + ilength + Month:ilength + Month:ilength:Hab, pcountdata,  mixture="NB", engine="C"); beepr::beep(11)
m39 <- pcount(~ Month + Water.Temp.Mean + Month:Relative.Humidity.Mean + Water.Temp.Mean:Relative.Humidity.Mean
              ~ ilength + Month:ilength + Month:Hab + Month:ilength:Hab, pcountdata,  mixture="NB", engine="C"); beepr::beep(11)
m40 <- pcount(~ Month + Water.Temp.Mean + Month:Relative.Humidity.Mean + Water.Temp.Mean:Relative.Humidity.Mean ## Full p and lambda
              ~ Month + ilength + Month:ilength + Month:Hab + Month:ilength:Hab, pcountdata,  mixture="NB", engine="C"); beepr::beep(11)  
m41 <- pcount(~ Water.Temp.Mean + Month:Relative.Humidity.Mean + Water.Temp.Mean:Relative.Humidity.Mean  
              ~ Month + ilength + Month:ilength + Month:Hab + Month:ilength:Hab, pcountdata, mixture="NB", engine="C"); beepr::beep(11)
m42 <- pcount(~ Month + Water.Temp.Mean + Water.Temp.Mean:Relative.Humidity.Mean
              ~ Month + ilength + Month:ilength + Month:Hab + Month:ilength:Hab, pcountdata, mixture="NB", engine="C"); beepr::beep(11)
m43 <- pcount(~ Month + Water.Temp.Mean + Month:Relative.Humidity.Mean 
              ~ Month + ilength + Month:ilength + Month:Hab + Month:ilength:Hab, pcountdata, mixture="NB", engine="C"); beepr::beep(11)
m44 <- pcount(~ Month:Relative.Humidity.Mean + Water.Temp.Mean:Relative.Humidity.Mean
              ~ Month + ilength + Month:ilength + Month:Hab + Month:ilength:Hab, pcountdata, mixture="NB", engine="C"); beepr::beep(11)
m45 <- pcount(~ Water.Temp.Mean + Water.Temp.Mean:Relative.Humidity.Mean 
              ~ Month + ilength + Month:ilength + Month:Hab + Month:ilength:Hab, pcountdata, mixture="NB", engine="C"); beepr::beep(11)
m46 <- pcount(~ Water.Temp.Mean + Month:Relative.Humidity.Mean 
              ~ Month + ilength + Month:ilength + Month:Hab + Month:ilength:Hab, pcountdata, mixture="NB", engine="C"); beepr::beep(11)
m47 <- pcount(~ Month + Water.Temp.Mean:Relative.Humidity.Mean 
              ~ Month + ilength + Month:ilength + Month:Hab + Month:ilength:Hab, pcountdata, mixture="NB", engine="C"); beepr::beep(11)
m48 <- pcount(~ Month + Month:Relative.Humidity.Mean 
              ~ Month + ilength + Month:ilength + Month:Hab + Month:ilength:Hab, pcountdata, mixture="NB", engine="C"); beepr::beep(11)
m49 <- pcount(~ Month + Water.Temp.Mean 
              ~ Month + ilength + Month:ilength + Month:Hab + Month:ilength:Hab, pcountdata, mixture="NB", engine="C"); beepr::beep(11)
m50 <- pcount(~ Water.Temp.Mean:Relative.Humidity.Mean
              ~ Month + ilength + Month:ilength + Month:Hab + Month:ilength:Hab, pcountdata, mixture="NB", engine="C"); beepr::beep(11)
m51 <- pcount(~ Month:Relative.Humidity.Mean 
              ~ Month + ilength + Month:ilength + Month:Hab + Month:ilength:Hab, pcountdata, mixture="NB", engine="C"); beepr::beep(11)
m52 <- pcount(~ Water.Temp.Mean 
              ~ Month + ilength + Month:ilength + Month:Hab + Month:ilength:Hab, pcountdata, mixture="NB", engine="C"); beepr::beep(11)
m53 <- pcount(~ Month 
              ~ Month + ilength + Month:ilength + Month:Hab + Month:ilength:Hab, pcountdata, mixture="NB", engine="C"); beepr::beep(11)
m54 <- pcount(~ 1 
              ~ Month + ilength + Month:ilength + Month:Hab + Month:ilength:Hab, pcountdata, mixture="NB",  engine="C"); beepr::beep(11) ## Full lamda

models <- fitList(
  `p(.)lambda(.)` = m1,
  `p(mt)lambda(.)` = m2,
  `p(wt)lambda(.)` = m3,
  `p(mt.rh)lambda(.)` = m4,
  `p(wt.rh)lambda(.)` = m5,
  `p(mt, wt)lambda(.)` = m6,
  `p(mt, mt.rh)lambda(.)` = m7,
  `p(mt, wt.rh)lambda(.)` = m8,
  `p(wt, mt.rh)lambda(.)` = m9,
  `p(wt, wt.rh)lambda(.)` = m10,
  `p(mt.rh, wt.rh)lambda(.)` = m11,
  `p(mt, wt, mt.rh)lambda(.)` = m12,
  `p(mt, wt, wt.rh)lambda(.)` = m13,
  `p(wt, mt.rh, wt.rh)lambda(.)` = m14,
  `p(full)lambda(.)` = m15,
  `p(full)lambda(mt)` = m16,
  `p(full)lambda(ilength)` = m17,
  `p(full)lambda(mt.ilength)` = m18,
  `p(full)lambda(mt.hab)` = m19,
  `p(full)lambda(mt.ilength.hab)` = m20,
  `p(full)lambda(mt, ilength)` = m21,
  `p(full)lambda(mt, mt.ilength)` = m22,
  `p(full)lambda(mt, mt.hab)` = m23,
  `p(full)lambda(mt, mt.ilength.hab)` = m24,
  `p(full)lambda(ilength, mt.ilength)` = m25,
  `p(full)lambda(ilength, mt.hab)` = m26,
  `p(full)lambda(ilength, mt.ilength.hab)` = m27,
  `p(full)lambda(mt.ilength, mt.hab)` = m28,
  `p(full)lambda(mt.ilength, mt.ilength.hab)` = m29,
  `p(full)lambda(mt.hab, mt.ilength.hab)` = m30,
  `p(full)lambda(mt, ilength, mt.ilength)` = m31,
  `p(full)lambda(mt, ilength, mt.hab)` = m32,
  `p(full)lambda(mt, ilength, mt.ilength.hab)` = m33,
  `p(full)lambda(ilength, mt.ilength, mt.hab)` = m34,
  `p(full)lambda(ilength, mt.ilength, mt.ilength.hab)` = m35,
  `p(full)lambda(mt.ilength, mt.hab, mt.ilength.hab)` = m36,
  `p(full)lambda(mt, ilength, mt.ilength, mt.hab)` = m37,
  `p(full)lambda(mt, ilength, mt.ilength, mt.ilength.hab)` = m38,
  `p(full)lambda(ilength, mt.ilength, mt.hab, mt.ilength.hab)` = m39,
  `p(full)lambda(full)` = m40,
  `p(wt, mt.rh, wt.rh)lambda(full)` = m41,
  `p(mt, wt, wt.rh)lambda(full)` = m42,
  `p(mt, wt, mt.rh)lambda(full)` = m43,
  `p(mt.rh, wt.rh)lambda(full)` = m44,
  `p(wt, wt.rh)lambda(full)` = m45,
  `p(wt, mt.rh)lambda(full)` = m46,
  `p(mt, wt.rh)lambda(full)` = m47,
  `p(mt, mt.rh)lambda(full)` = m48,
  `p(mt, wt)lambda(full)` = m49,
  `p(wt.rh)lambda(full)` = m50,
  `p(mt.rh)lambda(full)` = m51,
  `p(wt)lambda(full)` = m52,
  `p(mt)lambda(full)` = m53,
  `p(.)lambda(full)` = m54)

aic.table.all.T <- modSel(models, nullmod="p(.)lambda(.)")
(aic.table.all.T <- as(aic.table.all.T, "data.frame"))
write.csv(aic.table.all.T, "Model_Tables/Transects/aic_table.all.T.csv", row.names=F)

#Best model and general model-GoF figure
palette(brewer.pal(n=3, name="Dark2"))
system.time(gof.NB <- Nmix.gof.test(NBfull, engine = "C", nsim=100, #run 100 to test, 10000 for final model
                                    parallel = F)); beepr::beep(11)
system.time(gof.m19 <- Nmix.gof.test(m19, engine = "C", nsim=100, #run 100 to test, 10000 for final model
                                     parallel = F)); beepr::beep(11)
gof.NB; gof.m19

#Model averaging (AIC <2)
Cond.models <- list( )
Cond.models[[1]] <- m19
Cond.models[[2]] <- m26
Cond.models[[3]] <- m6
Cond.models[[4]] <- m16
Cond.models[[5]] <- m13
Cond.models[[6]] <- m30
Cond.models[[7]] <- m21
Cond.models[[8]] <- m28
Cond.models[[9]] <- m23

Modnames <- paste("mod", 1:length(Cond.models), sep = "")

summary(model.avg(Cond.models, subset = delta < 2, revised.var = TRUE))
ma <- model.avg(Cond.models, subset = delta < 2, revised.var = TRUE)
print(ma)

(mapP <- predict(ma, type = "det"))
mean(mapP$fit)
mean(mapP$se.fit)

maplambda <- predict(ma, type = "state")
mean(maplambda$fit)
mean(maplambda$se.fit)

sum(maplambda$fit)
sum(maplambda$se.fit)

Cond.models1 <- list( )
Cond.models1[[1]] <- m19
Cond.models1[[2]] <- m26
Cond.models1[[3]] <- m6
Cond.models1[[4]] <- m16
Cond.models1[[5]] <- m13
Cond.models1[[6]] <- m30
Cond.models1[[7]] <- m21
Cond.models1[[8]] <- m28
Cond.models1[[9]] <- m23
Cond.models1[[10]] <- m18
Cond.models1[[11]] <- m20
Cond.models1[[12]] <- m34
Cond.models1[[13]] <- m32
Cond.models1[[14]] <- m40
Cond.models1[[15]] <- m12
Cond.models1[[16]] <- m24
Cond.models1[[17]] <- m22
Cond.models1[[18]] <- m27
Cond.models1[[19]] <- m15
Cond.models1[[20]] <- m36
Cond.models1[[21]] <- m33
Cond.models1[[22]] <- m31

Modnames <- paste("mod", 1:length(Cond.models1), sep = "")

summary(model.avg(Cond.models1, subset = delta < 2, revised.var = TRUE))
ma1 <- model.avg(Cond.models1)
print(ma1)
coef(ma1)

(map1P <- predict(ma1, type = "det"))
mean(map1P$fit)
mean(map1P$se.fit)

map1 <- predict(ma1, type = "state")
sum(map1$fit)
sum(map1$se.fit)

avgmod.95p <- model.avg(Cond.models1, cumsum(weight) <= .95)
confint(avgmod.95p)
model.avg(ma1, cumsum(weight) <= .95, fit = TRUE)

##Model average prediction lambda and p
palette(brewer.pal(n=3, name="Dark2"))
tiff("Model_Figs/Transects/mod.ave.pre.tif", width=210, height=100, units="mm", res=500, compression="lzw")
par(mfrow=c(1,5), mar=c(4, 4, 2, 1), cex = 0.7)
##layout(matrix(c(1,1,2,3,4,5), 3, 2, byrow = TRUE))

#ilength lambda
rlength <- seq(3.5, 8.1, 0.135)
pred.df <- data.frame(Month = 0, 
                      ilength = 1/rlength, 
                      Hab = 0)

(pred.il.lambda <- modavgPred(cand.set = Cond.models, subset = delta < 2, modnames = Modnames, newdata = pred.df, 
                              parm.type = "lambda", type = "response"))

plot(rlength, pred.il.lambda$mod.avg.pred, type="l", lwd=3, ylim=c(10, 45), 
     xlab=expression("ilength"), 
     ylab="λ")

matlines(rlength, cbind(pred.il.lambda$mod.avg.pred - pred.il.lambda$uncond.se, 
                        pred.il.lambda$mod.avg.pred + pred.il.lambda$uncond.se), type = "l",
         lty = 1, lwd = 2, col = "gray")

pred.il.lambda
mean(pred.il.lambda$mod.avg.pred)
mean(pred.il.lambda$uncond.se)

#Month lambda
pred.df <- data.frame(Month = seq(1, 12, length.out = 60), 
                      ilength = 0, 
                      Hab = 0)

(pred.mt.lambda <- modavgPred(cand.set = Cond.models, modnames = Modnames, newdata = pred.df, 
           parm.type = "lambda", type = "response"))

plot(pred.df$Month, pred.mt.lambda$mod.avg.pred, type="l", lwd=3, ylim=c(0, 60), 
     xlab=expression("Month"), 
     ylab="λ")

matlines(pred.df$Month, cbind(pred.mt.lambda$mod.avg.pred - pred.mt.lambda$uncond.se, 
                        pred.mt.lambda$mod.avg.pred + pred.mt.lambda$uncond.se), type = "l",
         lty = 1, lwd = 2, col = "gray")

pred.mt.lambda
mean(pred.mt.lambda$mod.avg.pred)
mean(pred.mt.lambda$uncond.se)

##Month p
pred.df <- data.frame(Month = seq(1, 12, length.out = 60), Water.Temp.Mean=0, Relative.Humidity.Mean=0)

(pred.mt.p <- modavgPred(cand.set = Cond.models, modnames = Modnames, newdata = pred.df, 
                         parm.type = "detect", type = "response"))

plot(pred.df$Month, pred.mt.p$mod.avg.pred, type="l", lwd=3, ylim=c(0, 1), 
     xlab=expression("Month"), 
     ylab="p")

matlines(pred.df$Month, cbind(pred.mt.p$mod.avg.pred - pred.mt.p$uncond.se, 
                              pred.mt.p$mod.avg.pred + pred.mt.p$uncond.se), type = "l",
         lty = 1, lwd = 2, col = "gray")

pred.mt.p
mean(pred.mt.p$mod.avg.pred)
mean(pred.mt.p$uncond.se)

#Water temp p
min(Data$Water.Temp.Mean)
max(Data$Water.Temp.Mean)

pred.df <- data.frame(Month = 0, Water.Temp.Mean=seq(24, 30, length.out = 60), Relative.Humidity.Mean=0)

(pred.wt.p <- modavgPred(cand.set = Cond.models, modnames = Modnames, newdata = pred.df, 
                         parm.type = "detect", type = "response"))

plot(pred.df$Water.Temp.Mean, pred.wt.p$mod.avg.pred, type="l", lwd=3, ylim=c(0, 1), 
     xlab=expression("Water Temperature ºC"), 
     ylab="p")

matlines(pred.df$Water.Temp.Mean, cbind(pred.wt.p$mod.avg.pred - pred.wt.p$uncond.se, 
                              pred.wt.p$mod.avg.pred + pred.wt.p$uncond.se), type = "l",
         lty = 1, lwd = 2, col = "gray")

mean(pred.wt.p$mod.avg.pred)
mean(pred.wt.p$uncond.se)

#Relative humidity p
min(Data$Relative.Humidity.Mean)
max(Data$Relative.Humidity.Mean)

pred.df <- data.frame(Month = 0, Water.Temp.Mean=0, Relative.Humidity.Mean=seq(68, 91, length.out = 60))

(pred.rh.p <- modavgPred(cand.set = Cond.models, modnames = Modnames, newdata = pred.df, 
                         parm.type = "detect", type = "response"))

plot(pred.df$Relative.Humidity.Mean, pred.rh.p$mod.avg.pred, type="l", lwd=3, ylim=c(0, 1), 
     xlab=expression("Relative humidity (%)"), 
     ylab="p")

matlines(pred.df$Relative.Humidity.Mean, cbind(pred.rh.p$mod.avg.pred - pred.rh.p$uncond.se, 
                              pred.rh.p$mod.avg.pred + pred.rh.p$uncond.se), type = "l",
         lty = 1, lwd = 2, col = "gray")

mean(pred.rh.p$mod.avg.pred)
mean(pred.rh.p$uncond.se)

dev.off(); beepr::beep(11)

#Other variables lambda
Data.sc$ilength
pred.df <- data.frame(Month = 0, 
                      ilength = 0.2631579, # Transect 1
                      Hab = 0)

(pred.il.lambda <- modavgPred(cand.set = Cond.models, subset = delta < 2, modnames = Modnames, newdata = pred.df, 
                              parm.type = "lambda", type = "response"))

pred.df <- data.frame(Month = 0, 
                      ilength = 0.2857143, #Transect 2 
                      Hab = 0)

(pred.il.lambda <- modavgPred(cand.set = Cond.models, subset = delta < 2, modnames = Modnames, newdata = pred.df, 
                              parm.type = "lambda", type = "response"))

pred.df <- data.frame(Month = 0, 
                      ilength = 0.123458, #Transect 3 
                      Hab = 0)

(pred.il.lambda <- modavgPred(cand.set = Cond.models, subset = delta < 2, modnames = Modnames, newdata = pred.df, 
                              parm.type = "lambda", type = "response"))

pred.df <- data.frame(Month = 0, 
                      ilength = 0.2564103, #Transect 4 
                      Hab = 0)

(pred.il.lambda <- modavgPred(cand.set = Cond.models, subset = delta < 2, modnames = Modnames, newdata = pred.df, 
                              parm.type = "lambda", type = "response"))

pred.df <- data.frame(Month = 0, 
                      ilength = 0.1492537, #Transect 5
                      Hab = 0)

(pred.il.lambda <- modavgPred(cand.set = Cond.models, subset = delta < 2, modnames = Modnames, newdata = pred.df, 
                              parm.type = "lambda", type = "response"))


#month
pred.df <- data.frame(Month = 1, 
                      ilength = 0,
                      Hab = 0)

(pred.il.lambda <- modavgPred(cand.set = Cond.models, subset = delta < 2, modnames = Modnames, newdata = pred.df, 
                              parm.type = "lambda", type = "response"))

pred.df <- data.frame(Month = 2, 
                      ilength = 0,
                      Hab = 0)

(pred.il.lambda <- modavgPred(cand.set = Cond.models, subset = delta < 2, modnames = Modnames, newdata = pred.df, 
                              parm.type = "lambda", type = "response"))

pred.df <- data.frame(Month = 3, 
                      ilength = 0,
                      Hab = 0)

(pred.il.lambda <- modavgPred(cand.set = Cond.models, subset = delta < 2, modnames = Modnames, newdata = pred.df, 
                              parm.type = "lambda", type = "response"))

pred.df <- data.frame(Month = 4, 
                      ilength = 0,
                      Hab = 0)

(pred.il.lambda <- modavgPred(cand.set = Cond.models, subset = delta < 2, modnames = Modnames, newdata = pred.df, 
                              parm.type = "lambda", type = "response"))

pred.df <- data.frame(Month = 5, 
                      ilength = 0,
                      Hab = 0)

(pred.il.lambda <- modavgPred(cand.set = Cond.models, subset = delta < 2, modnames = Modnames, newdata = pred.df, 
                              parm.type = "lambda", type = "response"))

pred.df <- data.frame(Month = 6, 
                      ilength = 0,
                      Hab = 0)

(pred.il.lambda <- modavgPred(cand.set = Cond.models, subset = delta < 2, modnames = Modnames, newdata = pred.df, 
                              parm.type = "lambda", type = "response"))

pred.df <- data.frame(Month = 7, 
                      ilength = 0,
                      Hab = 0)

(pred.il.lambda <- modavgPred(cand.set = Cond.models, subset = delta < 2, modnames = Modnames, newdata = pred.df, 
                              parm.type = "lambda", type = "response"))

pred.df <- data.frame(Month = 8, 
                      ilength = 0,
                      Hab = 0)

(pred.il.lambda <- modavgPred(cand.set = Cond.models, subset = delta < 2, modnames = Modnames, newdata = pred.df, 
                              parm.type = "lambda", type = "response"))

pred.df <- data.frame(Month = 9, 
                      ilength = 0,
                      Hab = 0)

(pred.il.lambda <- modavgPred(cand.set = Cond.models, subset = delta < 2, modnames = Modnames, newdata = pred.df, 
                              parm.type = "lambda", type = "response"))

pred.df <- data.frame(Month = 10, 
                      ilength = 0,
                      Hab = 0)

(pred.il.lambda <- modavgPred(cand.set = Cond.models, subset = delta < 2, modnames = Modnames, newdata = pred.df, 
                              parm.type = "lambda", type = "response"))

pred.df <- data.frame(Month = 11, 
                      ilength = 0,
                      Hab = 0)

(pred.il.lambda <- modavgPred(cand.set = Cond.models, subset = delta < 2, modnames = Modnames, newdata = pred.df, 
                              parm.type = "lambda", type = "response"))

pred.df <- data.frame(Month = 12, 
                      ilength = 0,
                      Hab = 0)

(pred.il.lambda <- modavgPred(cand.set = Cond.models, subset = delta < 2, modnames = Modnames, newdata = pred.df, 
                              parm.type = "lambda", type = "response"))

#Hab
pred.df <- data.frame(Month = 0, 
                      ilength = 0,
                      Hab = 1) #River

(pred.il.lambda <- modavgPred(cand.set = Cond.models, subset = delta < 2, modnames = Modnames, newdata = pred.df, 
                              parm.type = "lambda", type = "response"))

pred.df <- data.frame(Month = 0, 
                      ilength = 0,
                      Hab = 2) #Oxbow

(pred.il.lambda <- modavgPred(cand.set = Cond.models, subset = delta < 2, modnames = Modnames, newdata = pred.df, 
                              parm.type = "lambda", type = "response"))

##Generalized linear model----
library(MASS)
##Loading and setting up data
Data <- read.csv("Abundance_T.csv", header = T, sep = ",")
as.factor(Data$Hab)
Data$Mt.Tr <- paste(Data$Month, Data$Transect, sep = "_")
Data$ilength <- paste(1/Data$Transect.length)
Data$ilength <- as.numeric(as.character(Data$ilength))
#Center and scale environmental variables
atmn <- mean(Data$Air.Temp.Mean)
atsd <- sd(Data$Air.Temp.Mean)
wtmn <- mean(Data$Water.Temp.Mean)
wtsd <- sd(Data$Water.Temp.Mean)
hummn <- mean(Data$Relative.Humidity.Mean)
humsd <- sd(Data$Relative.Humidity.Mean)
mtmn <- mean(Data$Month)
mtsd <- sd(Data$Month)
at.wtmn <- mean(Data$AT.WT)
at.wtsd <- sd(Data$AT.WT)

Data.sc <- Data %>%
  mutate(Air.Temp.Mean = (Air.Temp.Mean-atmn)/atsd,
         Water.Temp.Mean = (Water.Temp.Mean-wtmn)/wtsd,
         Relative.Humidity.Mean = (Relative.Humidity.Mean-hummn)/humsd,
         Month = (Month-mtmn)/mtsd,
         AT.WT = (AT.WT-at.wtmn)/at.wtsd)
Data.sc
## Exploring data
Var <- c("Air.Temp.Mean", "Water.Temp.Mean", "Relative.Humidity.Mean", "Month", "Hab", "Route_name", "Crocs")

dotplot(as.matrix(as.matrix(Data.sc[,Var])),
        groups = F, strip = strip.custom(bg = "white", par.strip.text = list(cex = 1.2)),
        scales = list(x = list(relation = "free", draw = T),
                      y = list(relation = "free", draw = F)),
        col = 1, cex = 1.0, pch = 16,
        xlab = list(label = "Data range", cex = 1.2),
        ylab = list(label = "Data order", cex = 1.2))

(mean(Data.sc$Crocs))
(var(Data.sc$Crocs))
sum(Data.sc$Crocs == 0)

Data.fit <- Data.sc[, c(1,4,7,8,10,11,13)]
as.factor(as.factor(Data.fit$Transect))
str(Data.fit)

##model fitting
#Poisson
fit1 <- glm(Crocs ~ ., data = Data.fit, family = poisson(), na.action = "na.fail")
fit1
(ms1 <- dredge(fit1))
summary(model.avg(ms1, subset = delta < 4, fit = T))

(ods <- fit1$deviance/fit1$df.residual)

#Negative Bionmial
fit <- glm.nb(Crocs ~ ., data = Data.fit, na.action = "na.fail")
(ms2 <- dredge(fit))
summary(model.avg(ms2, subset = delta < 2, fit = T))
(ms2.ma <- model.avg(ms2, subset = delta < 2, fit = T))

(ods <- fit$deviance/fit$df.residual)

#Predict
ms2p <- predict(ms2.ma)
mean(ms2p)
sd(ms2p)

pred.df <- data.frame(Air.Temp.Mean  = 0, 
                      Month = seq(1, 12, length.out = 60),
                      Relative.Humidity.Mean = 0,
                      Transect = 0,
                      Water.Temp.Mean = 0)

(pred.glm <- predict(ms2.ma, pred.df))

mean(pred.glm)
