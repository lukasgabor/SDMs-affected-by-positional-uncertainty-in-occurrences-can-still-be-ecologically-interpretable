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
library(sf)
library(tidyverse)
library(ggplot2)

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

## ------ Import envi data -----
twi <- raster("twi.tif") # load mask for shoft function

Bioall <- stack("aspect.tif", "dem.tif", "forest.tif", "grassland.tif", "twi.tif")


response.curves <- numeric()
var.importance <- numeric()
perf <- numeric()
## ----- Load bear data -----
DATA.Prep <- read.csv("data/bears/bears_projected.csv", sep = ";", dec = ".")

# ----- For Cyclus -----
for (i in 1:100) {
  print(paste('Cycle:', i))
  
  # Shifting data
  print(paste("Shifting data"))
  
  DATA.Prep[,c("S1x", "S1y")] <-  t(apply(DATA.Prep, 1, function(.) shift(.["x"], .["y"], dist=runif(1, 50, 200), mask=twi)))
  DATA.Prep[,c("S2x", "S2y")] <-  t(apply(DATA.Prep, 1, function(.) shift(.["x"], .["y"], dist=runif(1, 200, 500), mask=twi)))
  DATA.Prep[,c("S3x", "S3y")] <-  t(apply(DATA.Prep, 1, function(.) shift(.["x"], .["y"], dist=runif(1, 500, 1000), mask=twi)))
  DATA.Prep[,c("S4x", "S4y")] <-  t(apply(DATA.Prep, 1, function(.) shift(.["x"], .["y"], dist=runif(1, 1000, 5000), mask=twi)))
  DATA.Prep[,c("S5x", "S5y")] <-  t(apply(DATA.Prep, 1, function(.) shift(.["x"], .["y"], dist=runif(1, 5000, 10000), mask=twi)))
  DATA.Prep[,c("S6x", "S6y")] <-  t(apply(DATA.Prep, 1, function(.) shift(.["x"], .["y"], dist=runif(1, 10000, 20000), mask=twi)))
  
  DATA.Prep$occ <- "1"
  
  # write.csv(DATA.Prep, "data/bears/shifted_data.csv")
  
  # ------ Prepare of input spatial points data frames -----
  sp0 <- na.omit(data.frame(occ=DATA.Prep$occ, DATA.Prep[,1:2])) 
  coordinates(sp0) <- ~x+y
  sp1 <- na.omit(data.frame(occ=DATA.Prep$occ, x = DATA.Prep[,3], y = DATA.Prep[,4])) 
  coordinates(sp1) <- ~x+y
  sp2 <- na.omit(data.frame(occ=DATA.Prep$occ, x = DATA.Prep[,5], y = DATA.Prep[,6]))
  coordinates(sp2) <- ~x+y
  sp3 <- na.omit(data.frame(occ=DATA.Prep$occ, x = DATA.Prep[,7], y = DATA.Prep[,8])) 
  coordinates(sp3) <- ~x+y
  sp4 <- na.omit(data.frame(occ=DATA.Prep$occ, x = DATA.Prep[,9], y = DATA.Prep[,10])) 
  coordinates(sp4) <- ~x+y
  sp5 <- na.omit(data.frame(occ=DATA.Prep$occ, x = DATA.Prep[,11], y = DATA.Prep[,12])) 
  coordinates(sp5) <- ~x+y
  sp6 <- na.omit(data.frame(occ=DATA.Prep$occ, x = DATA.Prep[,13], y = DATA.Prep[,14]))
  coordinates(sp6) <- ~x+y
  
  # ----- SDM -----
  print(paste("SDM"))
  
  d0 <- sdmData(occ~.,train=sp0,predictors = Bioall,bg=list(n=10000,method='gRandom',remove=TRUE))
  d1 <- sdmData(occ~.,train=sp1,predictors = Bioall,bg=list(n=10000,method='gRandom',remove=TRUE))
  d2 <- sdmData(occ~.,train=sp2,predictors = Bioall,bg=list(n=10000,method='gRandom',remove=TRUE))
  d3 <- sdmData(occ~.,train=sp3,predictors = Bioall,bg=list(n=10000,method='gRandom',remove=TRUE))
  d4 <- sdmData(occ~.,train=sp4,predictors = Bioall,bg=list(n=10000,method='gRandom',remove=TRUE))
  d5 <- sdmData(occ~.,train=sp5,predictors = Bioall,bg=list(n=10000,method='gRandom',remove=TRUE))
  d6 <- sdmData(occ~.,train=sp6,predictors = Bioall,bg=list(n=10000,method='gRandom',remove=TRUE))
  
  print(paste("m0"))
  m0 <- sdm (occ~., data = d0, methods = "maxent", replication =  "cv", cv.folds=5, parallelSettings=list(ncore=4))
  m1 <- sdm (occ~., data = d1, methods = "maxent", replication =  "cv", cv.folds=5, parallelSettings=list(ncore=4))
  m2 <- sdm (occ~., data = d2, methods = "maxent", replication =  "cv", cv.folds=5, parallelSettings=list(ncore=4))
  print(paste("m3"))
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


write.table(x = perf, file = "outputs/bears/perfromance9.csv", row.names=F, col.names = T, sep = ";")
write.table(x = response.curves, file = "outputs/bears/response_curves9.csv", row.names=F, col.names = T, sep = ";")
write.table(x = var.importance, file = "outputs/bears/var_importance9.csv", row.names=F, col.names = T, sep = ";")

# ---- Example of ggplots -----
# Model Performance
bears <- read.csv("perfromance.csv", sep = ";", dec = ".")
bears$data <- factor(bears$data, levels=c('unaltered', 'S1', 'S2', 'S3','S4', 'S5', 'S6'), labels=c('Unaltered', 'S1 (50-200 m)', 'S2 (200-500 m)','S3 (500-1000 m)', 'S4 (1-5 km)', 'S5 (5-10 km)', 'S6 (10-20 km)'))


ggplot(bears, aes(data, Sorensen, fill = data)) +
  geom_violin(trim=FALSE, color="black") +
  scale_fill_brewer(palette="RdYlGn", direction=-1) +
  stat_summary(fun.y=median, geom="point", size=1, color = "black") +
  scale_y_continuous(limits=c(min(bears$Sorensen) - 0.02, max(bears$Sorensen) + 0.02), breaks = seq(0, 1, by = 0.02)) +
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
  labs(y="", x="", title = "Cantabrian Brown Bear")
ggsave(filename="si.png", dpi=600, height = 10, width = 10, units = "cm")

# Response Curves
df.bear <- read.csv("outputs/bears/response_curves9.csv", sep = ";")
df.bear$data <- factor(df.bear$data, levels=c('unaltered', 'S1', 'S2', 'S3','S4', 'S5', 'S6'), labels=c('Unaltered', 'S1 (50-200 m)', 'S2 (200-500 m)','S3 (500-1000 m)', 'S4 (1-5 km)', 'S5 (5-10 km)', 'S6 (10-20 km)'))

ggplot(df.bear[ which(df.bear$variable=='aspect'), ], aes(x=value, y=response, color = data, linetype = data)) + 
  geom_smooth(method="loess", se=F, size = 0.65) +
  scale_color_manual(values=c("black", "#006837", "#66bd63", "#d9ef8b", "#fee08b", "#f46d43", "#d73027")) +
  scale_linetype_manual(values = c("dashed", "solid", "solid", "solid", "solid", "solid", "solid")) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size=11, face="bold"),
    axis.title.y=element_text(size=10, face="bold"),
    axis.title = element_text(size=9),
    axis.title.x=element_text(size=10, face="bold"),
    axis.text.x = element_text(angle = 270, hjust=1, vjust = 0.5, size = 7),
    legend.title = element_blank()) +
  labs(title="VI = 0.33 %",
       x="Aspect",
       y="Occurrence probability",
       fill="Data Accuracy")
ggsave(filename="rc_aspect.png", dpi=600, height = 10, width = 10, units = "cm")

#Variable importance
df.bear.varimp <- read.csv("outputs/bears/var_importance9.csv", sep = ";")
df.bear.varimp$data <- factor(df.bear.varimp$data, levels=c('unaltered', 'S1', 'S2', 'S3','S4', 'S5', 'S6'), labels=c('Unaltered', 'S1 (50-200 m)', 'S2 (200-500 m)','S3 (500-1000 m)', 'S4 (1-5 km)', 'S5 (5-10 km)', 'S6 (10-20 km)'))
df.bear.varimp$variables <- factor(df.bear.varimp$variables, levels = c("aspect", "forest", "grassland", "dem", "twi"), labels = c("Aspect", "Amount of forest", "Amount of grassland", "Elevation", "Topography wetness index"))


df.bear.varimp2 <- df.bear.varimp %>% 
  group_by(variables, data) %>%
  summarize(mean = mean(AUCtest))

vi.bear <- ggplot(df.bear.varimp2, aes(x = data, y = mean, fill = variables)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99")) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size=11, face="bold"),
    axis.title.y=element_text(size=10, face="bold"),
    axis.title = element_text(size=9),
    axis.title.x=element_text(size=10, face="bold"),
    axis.text.x = element_text(angle = 270, hjust=1, vjust = 0.5, size = 7),
    legend.title = element_blank()) +
  labs(title="Cantabrian Brown Bear",
       x="Positional uncertainty",
       y="Variable importance")
ggsave(filename="variableImportance.png", plot=plot.varimp, dpi=600, height = 10, width = 25, units = "cm")

