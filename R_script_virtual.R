# Code Authors: Lukas Gabor & Vojtech Bartak, 2021

## Lukas Gabor
# https://scholar.google.cz/citations?user=pLQXY5wAAAAJ&hl=cs
# https://www.researchgate.net/profile/Lukas-Gabor

## Vojtech Bartak
# https://www.researchgate.net/profile/Vojtech-Bartak


# Load libraries
library(raster)
library(virtualspecies)
library(dismo)
library(rJava)
library(sdm)
library(usdm)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)

# ----- Load functions-----
# Function - Occ shifting
inside <- function(x, y, mask){
  ins <- ifelse(x < mask@extent@xmin | x > mask@extent@xmax | y < mask@extent@ymin | y > mask@extent@ymax, FALSE, TRUE)
  if (ins == TRUE) {
    val <- raster::extract(mask, matrix(c(x,y), ncol=2))
    ins <- ifelse(is.na(val), FALSE, TRUE)
  }
  return(ins)
}
shift <- function(x, y, dist, mask){
  new.x <- x
  new.y <- y
  inside <- FALSE
  while (!inside){
    angle <- runif(1,0,360)*pi/180
    new.x <- x + dist*sin(angle)
    new.y <- y + dist*cos(angle)
    inside <- inside(new.x, new.y, mask)
  }
  return(c(new.x, new.y))
}

# Function - Evaluation metrics from SDM package
sdm.package.evaluation <- function (fit.model){
  th <- mean(getEvaluation(fit.model, stat= "threshold", opt = 2)[,2])
  
  cm1 <- as.table(sdm:::.cmx(o = as.vector(fit.model@models$occ$maxent$`1`@evaluation$test.dep@observed),
                             p = as.vector(ifelse(fit.model@models$occ$maxent$`1`@evaluation$test.dep@predicted[] >= th, 1, 0))))
  
  cm2 <- as.table(sdm:::.cmx(o = as.vector(fit.model@models$occ$maxent$`2`@evaluation$test.dep@observed),
                             p = as.vector(ifelse(fit.model@models$occ$maxent$`2`@evaluation$test.dep@predicted[] >= th, 1, 0))))
  
  cm3 <- as.table(sdm:::.cmx(o = as.vector(fit.model@models$occ$maxent$`3`@evaluation$test.dep@observed),
                             p = as.vector(ifelse(fit.model@models$occ$maxent$`3`@evaluation$test.dep@predicted[] >= th, 1, 0))))
  
  cm4 <- as.table(sdm:::.cmx(o = as.vector(fit.model@models$occ$maxent$`4`@evaluation$test.dep@observed),
                             p = as.vector(ifelse(fit.model@models$occ$maxent$`4`@evaluation$test.dep@predicted[] >= th, 1, 0))))
  
  cm5 <- as.table(sdm:::.cmx(o = as.vector(fit.model@models$occ$maxent$`5`@evaluation$test.dep@observed),
                             p = as.vector(ifelse(fit.model@models$occ$maxent$`5`@evaluation$test.dep@predicted[] >= th, 1, 0))))
  
  eval <- getEvaluation(fit.model, stat= c("AUC", "Kappa"))
  
  perf <- rbind(data.frame(performance(cm1)), (performance(cm2)), data.frame(performance(cm3)),
                data.frame(performance(cm4)),data.frame(performance(cm5)))
  
  data.frame(TPR = mean(perf$TPR), TNR = mean(perf$TNR), FPR = mean(perf$FPR), FNR = mean(perf$FNR), Sensitivity = mean(perf$Sensitivity),
             Specificity = mean(perf$Specificity), AUC = mean(eval$AUC), Kappa = mean(eval$Kappa),TSS = mean(perf$TSS),
             Jaccard = mean(perf$Jaccard), Sorensen = mean(perf$Sorensen), F_measure = mean(perf$F_measure), OPR = mean(perf$OPR),
             UPR = mean(perf$UPR))
  
}

performance <- function (confusion.matrix) {
  tp <- confusion.matrix[1,1]
  fp <- confusion.matrix[1,2]
  fn <- confusion.matrix[2,1]
  tn <- confusion.matrix[2,2]
  TPR <- tp / (tp+fn)
  TNR <- tn / (tn+fp)
  FPR <- fp / (fp+tn)
  FNR <- fn / (fn+tp)
  Sensitivity <- TPR
  Specificity <- TNR
  TSS = Sensitivity + Specificity - 1
  Jaccard = TPR/(FNR + TPR + FPR) 
  Sorensen = 2*tp/(fn + 2*tp + fp)
  F_measure= 2 * Jaccard
  OPR = fp/(tp+fp) 
  UPR = 1 - Sensitivity
  data.frame(TPR = TPR, TNR = TNR, FPR = FPR, FNR = FNR, Sensitivity = Sensitivity, Specificity = Specificity,
             TSS = TSS, Jaccard = Jaccard, Sorensen = Sorensen, F_measure = F_measure, OPR = OPR, UPR = UPR)
}

# Function - Extract Response Curves from SDM package
extract_curves <- function(rc, nfolds = 5){
  rc %>%
    map(~ mutate(.x,
                 response = rowMeans(.[,2:(nfolds + 1)]),
                 variable = names(.x)[1]) %>%
          rename(value = 1) %>%
          select(-c(2:(nfolds + 1)))) %>%
    bind_rows()
}

# Function for mean om model predictions
my_rowmeans <- function(...) Reduce(`+`, list(...))/length(list(...))

predictions <- function(fit.model, nfolds = 5){
  
  df <- data.frame(cbind(fit.model@models$occ$maxent$`1`@evaluation$test.dep@predicted,
                         fit.model@models$occ$maxent$`2`@evaluation$test.dep@predicted,
                         fit.model@models$occ$maxent$`3`@evaluation$test.dep@predicted,
                         fit.model@models$occ$maxent$`4`@evaluation$test.dep@predicted,
                         fit.model@models$occ$maxent$`5`@evaluation$test.dep@predicted))
  
  df2 <- df %>% mutate(rms = my_rowmeans(X1, X2, X3, X4, X5))
  
  data.frame(df2$rms)
}


## ------ Import envi data -----
twi <- raster("dem.tif") # load mask for shift function

Bioall <- stack("dem.tif")

## ----- Virtual Species -----
my.parametrs <- formatFunctions(dem = c(fun = 'dnorm', mean = 1500, sd = 100))


virtual.species <- generateSpFromFun(raster.stack = Bioall,
                                     parameters = my.parametrs,
                                     species.type = "multiplicative",
                                     plot = F)

# Conversion to presence-absence
PA.raster <- convertToPA(virtual.species, alpha = -0.05, beta = 0.3, plot = F)


response.curves <- numeric()
var.importance <- numeric()
perf <- numeric()


# ----- For Cyclus -----
for (i in 1:100) {
  print(paste('Cycle:', i))
  
  # Sampling of '"presence-absence"'" and " presence - only "occurrences
  PA.sampling <- sampleOccurrences(PA.raster, n = 100, type = "presence only", plot = F)
  
  
  ## ----- Data Preparation and Shifting Occurrences -----
  DATA.Prep <- as.data.frame(PA.sampling$sample.points) # Import original coordinates
  s <- nrow(DATA.Prep)
  ID <- 1:s
  Observed <- DATA.Prep$Observed
  
  # Shifting data
  print(paste("Shifting data"))
  DATA.Prep[,c("S1x", "S1y")] <-  t(apply(DATA.Prep, 1, function(.) shift(.["x"], .["y"], dist=runif(1, 50, 200), mask=twi)))
  DATA.Prep[,c("S2x", "S2y")] <-  t(apply(DATA.Prep, 1, function(.) shift(.["x"], .["y"], dist=runif(1, 200, 500), mask=twi)))
  DATA.Prep[,c("S3x", "S3y")] <-  t(apply(DATA.Prep, 1, function(.) shift(.["x"], .["y"], dist=runif(1, 500, 1000), mask=twi)))
  DATA.Prep[,c("S4x", "S4y")] <-  t(apply(DATA.Prep, 1, function(.) shift(.["x"], .["y"], dist=runif(1, 1000, 5000), mask=twi)))
  DATA.Prep[,c("S5x", "S5y")] <-  t(apply(DATA.Prep, 1, function(.) shift(.["x"], .["y"], dist=runif(1, 5000, 10000), mask=twi)))
  DATA.Prep[,c("S6x", "S6y")] <-  t(apply(DATA.Prep, 1, function(.) shift(.["x"], .["y"], dist=runif(1, 10000, 20000), mask=twi)))
  
  # ------ Prepare of input spatial points data frames -----
  sp0 <- na.omit(data.frame(occ=DATA.Prep$Real, DATA.Prep[,1:2])) 
  coordinates(sp0) <- ~x+y
  sp1 <- na.omit(data.frame(occ=DATA.Prep$Real, x = DATA.Prep[,5], y = DATA.Prep[,6])) 
  coordinates(sp1) <- ~x+y
  sp2 <- na.omit(data.frame(occ=DATA.Prep$Real, x = DATA.Prep[,7], y = DATA.Prep[,8]))
  coordinates(sp2) <- ~x+y
  sp3 <- na.omit(data.frame(occ=DATA.Prep$Real, x = DATA.Prep[,9], y = DATA.Prep[,10])) 
  coordinates(sp3) <- ~x+y
  sp4 <- na.omit(data.frame(occ=DATA.Prep$Real, x = DATA.Prep[,11], y = DATA.Prep[,12])) 
  coordinates(sp4) <- ~x+y
  sp5 <- na.omit(data.frame(occ=DATA.Prep$Real, x = DATA.Prep[,13], y = DATA.Prep[,14])) 
  coordinates(sp5) <- ~x+y
  sp6 <- na.omit(data.frame(occ=DATA.Prep$Real, x = DATA.Prep[,15], y = DATA.Prep[,16]))
  coordinates(sp6) <- ~x+y
  
  # ----- SDM -----
  d0 <- sdmData(train = sp0, predictors = Bioall, bg=list(n=1000,method='gRandom',remove=TRUE))
  d1 <- sdmData(train = sp1, predictors = Bioall, bg=list(n=1000,method='gRandom',remove=TRUE))
  d2 <- sdmData(train = sp2, predictors = Bioall, bg=list(n=1000,method='gRandom',remove=TRUE))
  d3 <- sdmData(train = sp3, predictors = Bioall, bg=list(n=1000,method='gRandom',remove=TRUE))
  d4 <- sdmData(train = sp4, predictors = Bioall, bg=list(n=1000,method='gRandom',remove=TRUE))
  d5 <- sdmData(train = sp5, predictors = Bioall, bg=list(n=1000,method='gRandom',remove=TRUE))
  d6 <- sdmData(train = sp6, predictors = Bioall, bg=list(n=1000,method='gRandom',remove=TRUE))
  
  
  #parallel::detectCores()
  m0 <- sdm (occ~., data = d0, methods = "maxent", replication =  "cv", cv.folds=5, parallelSettings=list(ncore=4))
  m1 <- sdm (occ~., data = d1, methods = "maxent", replication =  "cv", cv.folds=5, parallelSettings=list(ncore=4))
  m2 <- sdm (occ~., data = d2, methods = "maxent", replication =  "cv", cv.folds=5, parallelSettings=list(ncore=4))
  m3 <- sdm (occ~., data = d3, methods = "maxent", replication =  "cv", cv.folds=5, parallelSettings=list(ncore=4))
  m4 <- sdm (occ~., data = d4, methods = "maxent", replication =  "cv", cv.folds=5, parallelSettings=list(ncore=4))
  m5 <- sdm (occ~., data = d5, methods = "maxent", replication =  "cv", cv.folds=5, parallelSettings=list(ncore=4))
  m6 <- sdm (occ~., data = d6, methods = "maxent", replication =  "cv", cv.folds=5, parallelSettings=list(ncore=4))
  
  
  perf <- rbind(perf,
                data.frame(sdm.package.evaluation(m0), data = "unaltered"),
                data.frame(sdm.package.evaluation(m1), data = "S1"),
                data.frame(sdm.package.evaluation(m2), data = "S2"),
                data.frame(sdm.package.evaluation(m3), data = "S3"),
                data.frame(sdm.package.evaluation(m4), data = "S4"),
                data.frame(sdm.package.evaluation(m5), data = "S5"),
                data.frame(sdm.package.evaluation(m6), data = "S6"))
  
  
  response.curves <- rbind(response.curves,
                           data.frame(extract_curves(getResponseCurve(m0)@response, nfolds = 5), data = "unaltered"),
                           data.frame(extract_curves(getResponseCurve(m1)@response, nfolds = 5), data = "S1"),
                           data.frame(extract_curves(getResponseCurve(m2)@response, nfolds = 5), data = "S2"),
                           data.frame(extract_curves(getResponseCurve(m3)@response, nfolds = 5), data = "S3"),
                           data.frame(extract_curves(getResponseCurve(m4)@response, nfolds = 5), data = "S4"),
                           data.frame(extract_curves(getResponseCurve(m5)@response, nfolds = 5), data = "S5"),
                           data.frame(extract_curves(getResponseCurve(m6)@response, nfolds = 5), data = "S6"))
  
  
  var.importance <- rbind(var.importance,
                          data.frame(getVarImp(m0)@varImportanceMean$corTest[,1:2], getVarImp(m0)@varImportanceMean$AUCtest[,1:2], data = "unaltered"),
                          data.frame(getVarImp(m1)@varImportanceMean$corTest[,1:2], getVarImp(m1)@varImportanceMean$AUCtest[,1:2], data = "S1"),
                          data.frame(getVarImp(m2)@varImportanceMean$corTest[,1:2], getVarImp(m2)@varImportanceMean$AUCtest[,1:2], data = "S2"),
                          data.frame(getVarImp(m3)@varImportanceMean$corTest[,1:2], getVarImp(m3)@varImportanceMean$AUCtest[,1:2], data = "S3"),
                          data.frame(getVarImp(m4)@varImportanceMean$corTest[,1:2], getVarImp(m4)@varImportanceMean$AUCtest[,1:2], data = "S4"),
                          data.frame(getVarImp(m5)@varImportanceMean$corTest[,1:2], getVarImp(m5)@varImportanceMean$AUCtest[,1:2], data = "S5"),
                          data.frame(getVarImp(m6)@varImportanceMean$corTest[,1:2], getVarImp(m6)@varImportanceMean$AUCtest[,1:2], data = "S6"))
  
  
}


write.table(x = perf, file = "perfromance.csv", row.names=F, col.names = T, sep = ";")
write.table(x = response.curves, file = "response_curves.csv", row.names=F, col.names = T, sep = ";")
write.table(x = var.importance, file = "var_importance.csv", row.names=F, col.names = T, sep = ";")

# ---- Example of ggplots -----
# Model Performance
virtual <- read.csv("perfromance.csv", sep = ";", dec = ".")
virtual$data <- factor(virtual$data, levels=c('unaltered', 'S1', 'S2', 'S3','S4', 'S5', 'S6'), labels=c('Unaltered', 'S1 (50-200 m)', 'S2 (200-500 m)','S3 (500-1000 m)', 'S4 (1-5 km)', 'S5 (5-10 km)', 'S6 (10-20 km)'))


ggplot(virtual, aes(data, Sorensen, fill = data)) +
  geom_violin(trim=FALSE) +
  scale_fill_brewer(palette="RdYlGn", direction=-1) +
  stat_summary(fun.y=median, geom="point", size=1, color = "black") +
  scale_y_continuous(limits=c(min(virtual$Sorensen) - 0.02, max(virtual$Sorensen) + 0.02), breaks = seq(0.67, 1, by = 0.02)) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size=12, face="bold"),
    axis.title.y=element_text(size=10, face="bold"),
    axis.title = element_text(size=9),
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    legend.position = "",
    legend.title = element_blank()) +
  labs(y="Sorensen index", x="", title = "Virtual species")


# Response Curves
df.virtual <- read.csv("response_curves.csv", sep = ";")
df.virtual$data <- factor(df.virtual$data, levels=c('unaltered', 'S1', 'S2', 'S3','S4', 'S5', 'S6'), labels=c('Unaltered', 'S1 (50-200 m)', 'S2 (200-500 m)','S3 (500-1000 m)', 'S4 (1-5 km)', 'S5 (5-10 km)', 'S6 (10-20 km)'))


ggplot(df.virtual, aes(x=value, y=response, color = data, linetype = data)) + 
  geom_smooth(method="loess", se=F, size = 0.65) +
  scale_color_manual(values=c("black", "#006837", "#66bd63", "#d9ef8b", "#fee08b", "#f46d43", "#d73027")) +
  scale_linetype_manual(values = c("dashed", "solid", "solid", "solid", "solid", "solid", "solid")) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, .2, .4, .6, .8, 1)) +
  scale_x_continuous(limits = c(0, 000), breaks = c(0, 500, 1000, 1500, 2000, 2500, 3000)) +
  theme(
    plot.title = element_text(hjust = 0.5, size=11, face="bold"),
    axis.title.y=element_text(size=10, face="bold"),
    axis.title = element_text(size=9),
    axis.title.x=element_text(size=10, face="bold"),
    axis.text.x = element_text(angle = 270, hjust=1, vjust = 0.5, size = 7),
    legend.title = element_blank()) +
  labs(title="",
       x="Elevation",
       y="Occurrence probability",
       fill="Data Accuracy")

ggsave(filename="rc_virtualSpecies.png", dpi=600, height = 10, width = 10, units = "cm")