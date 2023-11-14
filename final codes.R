library("SSN")
file.copy(system.file("lsndata/MiddleFork04.ssn", package = "SSN"),
          to = tempdir(), recursive = TRUE, copy.mode = FALSE)
setwd(tempdir())
#importSSN(Path, predpts = NULL, o.write = FALSE)
mf04p <- importSSN("./MiddleFork04.ssn",
                   predpts = "pred1km")
head(mf04p)



a <- mf04p@obspoints
obs <- a@SSNPoints[[1]]

b <- mf04p@predpoints
pred <- b@SSNPoints[[1]]

obs_data <- data.frame(obs@network.point.coords,obs@point.data,
                       obs@point.coords)
pred_data <- data.frame(pred@network.point.coords,pred@point.data
                        ,pred@point.coords)

dim(obs_data)
names(obs_data)
dat_obs <- obs_data[,c(39,40,3,8,9,10,11,33,34,35,19)] #final observed data
dat_pred <- pred_data[,c(30,31,3,6,9,10,11,24,25,26)] #final predicted data


library(sp)
library(fields)
library(gstat)




# Simulation Study



set.seed(100)
raw.ssn <- createSSN(n = c(10, 10, 10),
                     obsDesign = binomialDesign(c(40, 40, 40)),
                     predDesign = systematicDesign(c(0.2, 0.4, 0.8)), 
                     importToR = TRUE,
                     path = "./raw8.ssn")
plot(raw.ssn, lwdLineCol = "addfunccol", lwdLineEx = 8,
     lineCol = "green", cex = 2, xlab = "x-coordinate",
     ylab = "y-coordinate", pch = 1, cex.axis = 1.5, cex.lab=1.5)
plot(raw.ssn, PredPointsID = "preds", add = TRUE, cex = 1, pch = 19,
     col = "red", cex.axis = 1.5, cex.lab=1.5)

rawDFobs <- getSSNdata.frame(raw.ssn, "Obs")
rawDFpred <- getSSNdata.frame(raw.ssn, "preds")

rawDFobs[,"X1"] <- rnorm(length(rawDFobs[,1]))
rawDFpred[,"X1"] <- rnorm(length(rawDFpred[,1]))
rawDFobs[,"X2"] <- rnorm(length(rawDFobs[,1]))
rawDFpred[,"X2"] <- rnorm(length(rawDFpred[,1]))

rawDFobs[,"F1"] <- as.factor(sample.int(4,length(rawDFobs[,1]),
                                        replace = TRUE))
rawDFpred[,"F1"] <- as.factor(sample.int(4,length(rawDFpred[,1]),
                                         replace = TRUE))
rawDFobs[,"RE1"] <- as.factor(sample(1:3,length(rawDFobs[,1]),
                                     replace = TRUE))
rawDFobs[,"RE2"] <- as.factor(sample(1:4,length(rawDFobs[,1]),
                                     replace = TRUE))
rawDFpred[,"RE1"] <- as.factor(sample(1:3,length(rawDFpred[,1]),
                                      replace = TRUE))
rawDFpred[,"RE2"] <- as.factor(sample(1:4,length(rawDFpred[,1]),
                                      replace = TRUE))

createDistMat(raw.ssn, "preds", o.write=TRUE, amongpred = TRUE)

set.seed(42)
sim.out <- SimulateOnSSN(raw.ssn, ObsSimDF = rawDFobs,
                         PredSimDF = rawDFpred, PredID = "preds",
                         formula = ~ X1 + X2 + F1, coefficients = c(10,1,0,-2,0,2),
                         CorModels = c("LinearSill.taildown","Exponential.tailup",
                                       "Exponential.Euclid", "RE1", "RE2"), use.nugget = TRUE,
                         CorParms = c(3, 10, 2, 10, 1, 5, 1, .5, .1),
                         addfunccol = "addfunccol")
sim.ssn <- sim.out$ssn.object
plot(sim.ssn, "Sim_Values",
     xlab = "x-coordinate", ylab = "y-coordinate",
     cex = 1.5, cex.axis=1.5, cex.lab=1.5)

plot(mf04p, "Summer_mn",
     xlab = "x-coordinate", ylab = "y-coordinate",
     cex = 1.5, cex.axis=1.5, cex.lab=1.5)


simDFobs <- getSSNdata.frame(sim.ssn, "Obs")
simDFpred <- getSSNdata.frame(sim.ssn, "preds")

simpreds <- simDFpred[,"Sim_Values"]
simDFpred[,"Sim_Values"] <- NA
sim.ssn <- putSSNdata.frame(simDFpred, sim.ssn, "preds")

glmssn.out <- glmssn(Sim_Values ~ X1 + X2 + F1, sim.ssn,
                     CorModels = c("LinearSill.taildown","Exponential.tailup",
                                   "Exponential.Euclid", "RE1", "RE2"),
                     addfunccol = "addfunccol")

summary(glmssn.out)
glmssn.pred <- predict(glmssn.out,"preds")
predDF <- getSSNdata.frame(glmssn.pred, "preds")
plot(simpreds, predDF[,"Sim_Values"], xlab = "True", col="deeppink",
     ylab = "Predicted", pch = 19, cex= 1.5, cex.axis=1.5,cex.lab=1.5)

sim.out$FixedEffects
sim.out$CorParms









# Data Analysis



plot(mf04p, lwdLineCol = "afvArea", lwdLineEx = 8,
     lineCol = "deeppink", cex = 2, xlab = "x-coordinate",
     ylab = "y-coordinate", pch = 1, cex.axis = 1.5, cex.lab=1.5)
plot(mf04p, PredPointsID = "pred1km", add = TRUE, cex = 1, pch = 19,
     col = "darkblue", cex.axis = 1.5, cex.lab=1.5)

mf04.Torg <- Torgegram(mf04p, "Summer_mn", nlag = 70, maxlag = 100000)
plot(mf04.Torg,col=c("deeppink","darkblue"),pch=c(13,8),cex.axis=1.4,cex.lab=1.4)


rawDFobs <- getSSNdata.frame(mf04p, "Obs")
rawDFpred <- getSSNdata.frame(mf04p, "pred1km")

createDistMat(mf04p, "pred1km", o.write=TRUE, amongpred = TRUE)
plot(mf04p,"Summer_mn", xlab = "x-coordinate",
     ylab = "y-coordinate", pch = 19, cex.axis = 1.5, cex.lab=1.5)

glm <- glmssn(Summer_mn ~ ELEV_DEM + SLOPE+CUMDRAINAG+AREAWTMAP+MAXELEVSMO+SLOPE+ratio+afvArea,mf04p,
              CorModels =c("LinearSill.taildown","Exponential.tailup",
                           "Exponential.Euclid"),addfunccol = "afvArea")
glm.pred <- predict(glm,"pred1km")
predfg <- getSSNdata.frame(glm.pred,"pred1km")

tail_up <- glmssn(Summer_mn ~ ELEV_DEM + SLOPE+CUMDRAINAG+AREAWTMAP+MAXELEVSMO+SLOPE+ratio+afvArea,mf04p,
                  CorModels = "Exponential.tailup",addfunccol = "afvArea")
summary(tail_up)
p1 <- predict(tail_up,predpointsID = "pred1km")
plot(p1)

tail_down <- glmssn(Summer_mn ~ ELEV_DEM + SLOPE+CUMDRAINAG+AREAWTMAP+MAXELEVSMO+SLOPE+ratio+afvArea,mf04p,
                    CorModels = "Exponential.taildown",addfunccol = "afvArea")
summary(tail_down)
p2 <- predict(tail_down,predpointsID = "pred1km")
plot(p2)

tail_el <- glmssn(Summer_mn ~ ELEV_DEM + SLOPE+CUMDRAINAG+AREAWTMAP+MAXELEVSMO+SLOPE+ratio+afvArea,mf04p,
                  CorModels = "Exponential.Euclid",addfunccol = "afvArea")
summary(tail_el)
p3 <- predict(tail_el,predpointsID = "pred1km")
plot(p3)



