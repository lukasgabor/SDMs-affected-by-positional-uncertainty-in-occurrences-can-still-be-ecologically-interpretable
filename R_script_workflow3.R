########################################################################################
#Code Authors:
# Lukas Gabor - https://scholar.google.cz/citations?user=pLQXY5wAAAAJ&hl=cs
# Vojtech Bartak - https://scholar.google.cz/citations?user=p8WAo8oAAAAJ&hl=cs&oi=ao

# Date: 10/14/2022
########################################################################################

# Worfklow 2 R script - real environmental data + virtual species 

# R script used to develope SDMs ------------------------------------------------------------------

# Load packages
library(raster)
library(virtualspecies)
library(dismo)
library(rJava)
library(sdm)
library(usdm)
library(tidyverse)
library(gstat)
library(climateStability)

# Load functions ----------------------------------------------------------

# Shift species occurrences
inside <- function(x, y, mask) {
  # limits shifting to study are
  ins <-
    ifelse(
      x < mask@extent@xmin |
        x > mask@extent@xmax |
        y < mask@extent@ymin | y > mask@extent@ymax,
      FALSE,
      TRUE
    )
  if (ins == TRUE) {
    val <- raster::extract(mask, matrix(c(x, y), ncol = 2))
    ins <- ifelse(is.na(val), FALSE, TRUE)
  }
  return(ins)
}

shift <- function(x, y, dist, mask) {
  # shift species occurrences
  new.x <- x
  new.y <- y
  inside <- FALSE
  while (!inside) {
    angle <- runif(1, 0, 360) * pi / 180
    new.x <- x + dist * sin(angle)
    new.y <- y + dist * cos(angle)
    inside <- inside(new.x, new.y, mask)
  }
  return(c(new.x, new.y))
}

# Evaluate SDM performance
performance <- function (confusion.matrix) {
  tp <- confusion.matrix[1, 1]
  fp <- confusion.matrix[1, 2]
  fn <- confusion.matrix[2, 1]
  tn <- confusion.matrix[2, 2]
  TPR <- tp / (tp + fn)
  TNR <- tn / (tn + fp)
  FPR <- fp / (fp + tn)
  FNR <- fn / (fn + tp)
  Sensitivity <- TPR
  Specificity <- TNR
  TSS = Sensitivity + Specificity - 1
  Jaccard = TPR / (FNR + TPR + FPR)
  Sorensen = 2 * TPR / (FNR + 2 * TPR + FPR)
  F_measure = 2 * Jaccard
  OPR = FPR / (TPR + FPR)
  UPR = 1 - Sensitivity
  data.frame(
    tp = tp,
    fp = fp,
    fn = fn,
    tn = tn,
    TPR = TPR,
    TNR = TNR,
    FPR = FPR,
    FNR = FNR,
    Sensitivity = Sensitivity,
    Specificity = Specificity,
    TSS = TSS,
    Jaccard = Jaccard,
    Sorensen = Sorensen,
    F_measure = F_measure,
    OPR = OPR,
    UPR = UPR
  )
}

sdm.evaluation <- function (fit.model){
  
  fit.model <- fit.model
  th <- mean(getEvaluation(fit.model, stat= "threshold", opt = 2)[,2])
  
  cm1 <- as.table(sdm:::.cmx(o = as.vector(fit.model@models$occ$maxent$`1`@evaluation$test.dep@observed),
                             p = as.vector(ifelse(fit.model@models$occ$maxent$`1`@evaluation$test.dep@predicted >= th, 1, 0))))
  
  cm2 <- as.table(sdm:::.cmx(o = as.vector(fit.model@models$occ$maxent$`2`@evaluation$test.dep@observed),
                             p = as.vector(ifelse(fit.model@models$occ$maxent$`2`@evaluation$test.dep@predicted>= th, 1, 0))))
  
  cm3 <- as.table(sdm:::.cmx(o = as.vector(fit.model@models$occ$maxent$`3`@evaluation$test.dep@observed),
                             p = as.vector(ifelse(fit.model@models$occ$maxent$`3`@evaluation$test.dep@predicted>= th, 1, 0))))
  
  cm4 <- as.table(sdm:::.cmx(o = as.vector(fit.model@models$occ$maxent$`4`@evaluation$test.dep@observed),
                             p = as.vector(ifelse(fit.model@models$occ$maxent$`4`@evaluation$test.dep@predicted>= th, 1, 0))))
  
  cm5 <- as.table(sdm:::.cmx(o = as.vector(fit.model@models$occ$maxent$`5`@evaluation$test.dep@observed),
                             p = as.vector(ifelse(fit.model@models$occ$maxent$`5`@evaluation$test.dep@predicted>= th, 1, 0))))
  
  df <- cbind(getEvaluation(fit.model, stat= c("AUC", "Kappa"))[, 2:3],
              rbind(performance(cm1), performance(cm2), performance(cm3), performance(cm4), performance(cm5)))
  
  df.perf <- data.frame(tp = mean(df$tp), fp = mean(df$fp), fn = mean(df$fn), tn = mean(df$tn), TPR = mean(df$TPR), TNR = mean(df$TNR),
                        FPR = mean(df$FPR), FNR = mean(df$FNR), Sensitivity = mean(df$Sensitivity), Specificity = mean(df$Specificity),
                        TSS = mean(df$TSS), Jaccard = mean(df$Jaccard), Sorensen = mean(df$Sorensen), F_measure = mean(df$F_measure),
                        OPR = mean(df$OPR), UPR = mean(df$UPR), AUC = mean(df$AUC), Kappa = mean(df$Kappa))
  df.perf
  
}

# Extract Response Curves from SDM package
extract_curves <- function(rc, nfolds = 5){
  rc %>%
    map(~ mutate(.x,
                 response = rowMeans(.[,2:(nfolds + 1)]),
                 variable = names(.x)[1]) %>%
          rename(value = 1) %>%
          dplyr::select(-c(2:(nfolds + 1)))) %>%
    bind_rows()
}

# Load predictors
twi <- raster("PathToYourFile/bio1.tif") # load mask for shift function

rastlist <- list.files(path = "PathToYourFiles", pattern='.tif', all.files=TRUE, full.names=TRUE)
Bioall <- stack(rastlist)

# Read and prepare unaltered species data 
DATA.Prep <- read.csv("PathToYourFile/band_tailed_pigeon.csv", sep = ";", dec = ".")
DATA.Prep <- DATA.Prep[DATA.Prep$occ == "1", ]

sp0 <- as.data.frame(raster::extract(Bioall, DATA.Prep[, 2:3]))
sp0$occ <- "1"

# Get extent for sampling background points
ext <- DATA.Prep[, 2:3]
coordinates(ext) <- ~x+y
ext <- extent(ext)

# Sampling background points
bg0 <- as.data.frame(randomPoints(Bioall[[1]], p = sp0, n=10000, ext=ext, extf = 1.1))
bg0 <- as.data.frame(raster::extract(Bioall, bg0))
bg0$occ <- "0"

# Prepare empty data frames for results
response.curves <- numeric()
var.importance <- numeric()
perf <- numeric()

for (i in 1:50) {
  
  # Shifting data
  print(paste("Shifting data"))
  
  DATA.Prep[,c("S1x", "S1y")] <-  t(apply(DATA.Prep, 1, function(.) shift(.["x"], .["y"], dist=runif(1, 1000, 2000), mask=twi))) # 1-2 kilometers
  DATA.Prep[,c("S2x", "S2y")] <-  t(apply(DATA.Prep, 1, function(.) shift(.["x"], .["y"], dist=runif(1, 2000, 5000), mask=twi))) # 2-5 kilometers
  DATA.Prep[,c("S3x", "S3y")] <-  t(apply(DATA.Prep, 1, function(.) shift(.["x"], .["y"], dist=runif(1, 5000, 10000), mask=twi))) # 5-10 kilometers
  DATA.Prep[,c("S4x", "S4y")] <-  t(apply(DATA.Prep, 1, function(.) shift(.["x"], .["y"], dist=runif(1, 10000, 30000), mask=twi))) # 10-30 kilometers 
  
  # Prepare occurrences for SDMs ----------------------------------------------------------
  sp1 <- as.data.frame(raster::extract(Bioall, DATA.Prep[, 4:5]))
  sp1$occ <- "1"
  
  sp2 <- as.data.frame(raster::extract(Bioall, DATA.Prep[, 6:7]))
  sp2$occ <- "1"
  
  sp3 <- as.data.frame(raster::extract(Bioall, DATA.Prep[, 8:9]))
  sp3$occ <- "1"
  
  sp4 <- as.data.frame(raster::extract(Bioall, DATA.Prep[, 10:11]))
  sp4$occ <- "1"
  
  # Sampling background points
  bg1 <- as.data.frame(randomPoints(Bioall[[1]], p = sp1, n=10000, ext=ext, extf = 1.1))
  bg1 <- as.data.frame(raster::extract(Bioall, bg1))
  bg1$occ <- "1"
  
  bg2 <- as.data.frame(randomPoints(Bioall[[1]], p = sp2, n=10000, ext=ext, extf = 1.1))
  bg2 <- as.data.frame(raster::extract(Bioall, bg2))
  bg2$occ <- "0"
  
  bg3 <- as.data.frame(randomPoints(Bioall[[1]], p = sp3, n=10000, ext=ext, extf = 1.1))
  bg3 <- as.data.frame(raster::extract(Bioall, bg3))
  bg3$occ <- "0"
  
  bg4 <- as.data.frame(randomPoints(Bioall[[1]], p = sp3, n=10000, ext=ext, extf = 1.1))
  bg4 <- as.data.frame(raster::extract(Bioall, bg4))
  bg4$occ <- "0"
  
  # Species distribution modeling ----------------------------------------------------------
  # Create sdmData object
  d0 <- sdmData(occ~.,train=sp0, bg=bg0)
  d1 <- sdmData(occ~.,train=sp1, bg=bg1)
  d2 <- sdmData(occ~.,train=sp2, bg=bg2)
  d3 <- sdmData(occ~.,train=sp3, bg=bg3)
  d4 <- sdmData(occ~.,train=sp4, bg=bg4)
  
  # Fit and evaluate models
  m0 <- sdm (occ~., data = d0, methods = "maxent", replication =  "cv", cv.folds=5, parallelSettings=list(ncore=2),
             modelSettings=list(maxent=list(beta=0.5, args=c('nolinear', 'nothreshold', 'nohinge', 'noproduct'))))
  m1 <- sdm (occ~., data = d1, methods = "maxent", replication =  "cv", cv.folds=5, parallelSettings=list(ncore=2),
             modelSettings=list(maxent=list(beta=0.5, args=c('nolinear', 'nothreshold', 'nohinge', 'noproduct'))))
  m2 <- sdm (occ~., data = d2, methods = "maxent", replication =  "cv", cv.folds=5, parallelSettings=list(ncore=2),
             modelSettings=list(maxent=list(beta=0.5, args=c('nolinear', 'nothreshold', 'nohinge', 'noproduct'))))
  m3 <- sdm (occ~., data = d3, methods = "maxent", replication =  "cv", cv.folds=5, parallelSettings=list(ncore=2),
             modelSettings=list(maxent=list(beta=0.5, args=c('nolinear', 'nothreshold', 'nohinge', 'noproduct'))))
  m4 <- sdm (occ~., data = d4, methods = "maxent", replication =  "cv", cv.folds=5, parallelSettings=list(ncore=2),
             modelSettings=list(maxent=list(beta=0.5, args=c('nolinear', 'nothreshold', 'nohinge', 'noproduct'))))
  
  print(paste('Results'))
  
  # Combine results ----------------------------------------------------------
  # Performance metrics
  perf <- rbind(perf,
                data.frame(sdm.evaluation(m0), data = "unaltered"),
                data.frame(sdm.evaluation(m1), data = "S1"),
                data.frame(sdm.evaluation(m2), data = "S2"),
                data.frame(sdm.evaluation(m3), data = "S3"),
                data.frame(sdm.evaluation(m4), data = "S4"))
  
  # Response curves
  response.curves <- rbind(response.curves,
                           data.frame(extract_curves(getResponseCurve(m0)@response, nfolds = 5), data = "unaltered", id = 1:nrow(extract_curves(getResponseCurve(m0)@response, nfolds = 5))),
                           data.frame(extract_curves(getResponseCurve(m1)@response, nfolds = 5), data = "S1", id = 1:nrow(extract_curves(getResponseCurve(m0)@response, nfolds = 5))),
                           data.frame(extract_curves(getResponseCurve(m2)@response, nfolds = 5), data = "S2", id = 1:nrow(extract_curves(getResponseCurve(m0)@response, nfolds = 5))),
                           data.frame(extract_curves(getResponseCurve(m3)@response, nfolds = 5), data = "S3", id = 1:nrow(extract_curves(getResponseCurve(m0)@response, nfolds = 5))),
                           data.frame(extract_curves(getResponseCurve(m4)@response, nfolds = 5), data = "S4", id = 1:nrow(extract_curves(getResponseCurve(m0)@response, nfolds = 5))))
  
  # Variable importance
  var.importance <- rbind(var.importance,
                          data.frame(getVarImp(m0)@varImportanceMean$corTest[,1:2], getVarImp(m0)@varImportanceMean$AUCtest[,1:2], data = "unaltered"),
                          data.frame(getVarImp(m1)@varImportanceMean$corTest[,1:2], getVarImp(m1)@varImportanceMean$AUCtest[,1:2], data = "S1"),
                          data.frame(getVarImp(m2)@varImportanceMean$corTest[,1:2], getVarImp(m2)@varImportanceMean$AUCtest[,1:2], data = "S2"),
                          data.frame(getVarImp(m3)@varImportanceMean$corTest[,1:2], getVarImp(m3)@varImportanceMean$AUCtest[,1:2], data = "S3"),
                          data.frame(getVarImp(m4)@varImportanceMean$corTest[,1:2], getVarImp(m4)@varImportanceMean$AUCtest[,1:2], data = "S4"))
  
  print(paste("Finished running set", i))
  
} # End of loops

# Save results ------------------------------------------------------------
write.table(x = perf, file = "PathToYourFile/performance.csv", row.names=F, col.names = T, sep = ";")
write.table(x = response.curves, file = "PathToYourFile/response_curves.csv", row.names=F, col.names = T, sep = ";")
write.table(x = var.importance, file = "PathToYourFile/var_importance.csv", row.names=F, col.names = T, sep = ";")


# Process results ------------------------------------------------------------
library(dplyr)

## Response curves
# Load data
response.curves <- read.csv(file = "PathToYourFile/response_curves.csv", sep = ";")

# Get mean of values across all scenarios and loops
response.curves.df.sum <- response.curves %>% dplyr::group_by(data, variable, id)  %>% 
  dplyr::summarise(
    value = mean(value, na.rm=T),
    response = mean(response, na.rm=T))

# Save summarized data
write.table(response.curves.df.sum, file = "PathToYourFile/response_curves_sum.csv", sep = ";")

## Variable importance
df.var_importance <- read.csv(file = "PathToYourFile/var_importance.csv", sep = ";")

# Get mean of values across all scenarios and loops
var_importance.sum <- df.var_importance %>% dplyr::group_by(data, variables)  %>% 
  dplyr::summarise(
    corTest = mean(corTest, na.rm=T),
    AUCtest = mean(AUCtest, na.rm=T))

write.table(var_importance.sum, file = "PathToYourFile/var_importance_sum.csv", sep = ";")


# Plot results ------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(RColorBrewer)

## Model performance
# Load and process data
# Load data
df.sum <- read.csv(file = "PathToYourFile/performance.csv", sep = ";")
df.sum$data <- factor(df.sum$data, levels=c("unaltered", "S1", "S2", "S3", "S4"), labels = c("Unaltered", "1-2km", "2-5km", "5-10km", "10-30km"))

# Sorensen index
ggplot(df.sum, aes(data, Sorensen, fill = data)) +
  geom_violin(trim=FALSE, color="black") +
  scale_fill_brewer(palette="RdYlBu", direction=-1) +
  stat_summary(fun.y=median, geom="point", size=1, color = "black") +
  scale_y_continuous(limits=c(0.78, .90), breaks = c(0.78, 0.80, 0.82, 0.84,  0.86, 0.88, 0.9, .92, .94, .96, .98, 1)) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size=12, face="bold"),
    axis.title.y=element_text(size=10, face="bold"),
    axis.title = element_text(size=9),
    axis.title.x=element_text(size=10, face="bold"),
    axis.text.x=element_text(size=9),
    legend.position = "",
    legend.title = element_blank()) +
  labs(y="Sorensen index", x="Positional uncertainty", title = "")
ggsave(file="PathToYourFile/Sorensen.png", dpi = 600, width = 10, height = 10, units = "cm")

# Overprediction rate (OPR)
ggplot(df.sum, aes(data, OPR, fill = data)) +
  geom_violin(trim=FALSE, color="black") +
  scale_fill_brewer(palette="RdYlBu", direction=-1) +
  stat_summary(fun.y=median, geom="point", size=1, color = "black") +
  scale_y_continuous(limits=c(0.09, 0.2), breaks = seq(0.09, 0.25, 0.05)) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size=12, face="bold"),
    axis.title.y=element_text(size=10, face="bold"),
    axis.title = element_text(size=9),
    axis.title.x=element_text(size=10, face="bold"),
    axis.text.x=element_text(size=9),
    legend.position = "",
    legend.title = element_blank()) +
  labs(y="Overprediction rate", x="Positional uncertainty", title = "")
ggsave(file="PathToYourFile/OPR.png", dpi = 600, width = 10, height = 10, units = "cm")

# Underprediction rate (UPR)
ggplot(df.sum, aes(data, UPR, fill = data)) +
  geom_violin(trim=FALSE, color="black") +
  scale_fill_brewer(palette="RdYlBu", direction=-1) +
  stat_summary(fun.y=median, geom="point", size=1, color = "black") +
  scale_y_continuous(limits=c(0.06, 0.28), breaks = seq(0.06, 0.28, 0.05)) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size=12, face="bold"),
    axis.title.y=element_text(size=10, face="bold"),
    axis.title = element_text(size=9),
    axis.title.x=element_text(size=10, face="bold"),
    axis.text.x=element_text(size=9),
    legend.position = "",
    legend.title = element_blank()) +
  labs(y="Underprediction rate", x="Positional uncertainty", title = "")
ggsave(file="PathToYourFile/UPR.png", dpi = 600, width = 10, height = 10, units = "cm")

# TSS
ggplot(df.sum, aes(data, TSS, fill = data)) +
  geom_violin(trim=FALSE, color="black") +
  scale_fill_brewer(palette="RdYlBu", direction=-1) +
  stat_summary(fun.y=median, geom="point", size=1, color = "black") +
  scale_y_continuous(limits=c(0.60, .78), breaks = seq(0.6, 0.78, 0.02)) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size=12, face="bold"),
    axis.title.y=element_text(size=10, face="bold"),
    axis.title = element_text(size=9),
    axis.title.x=element_text(size=10, face="bold"),
    axis.text.x=element_text(size=9),
    legend.position = "",
    legend.title = element_blank()) +
  labs(y="TSS", x="Positional uncertainty", title = "")
ggsave(file="PathToYourFile/TSS.png", dpi = 600, width = 10, height = 10, units = "cm")


## Response curves
# Load and process data
df.work3 <- read.csv(file = "PathToYourFile/response_curves_sum.csv", sep = ";")
df.work3$data <- factor(df.work3$data, levels=c("unaltered", "S1", "S2", "S3", "S4"), labels = c("Unaltered", "1-2km", "2-5km", "5-10km", "10-30km"))
df.work3$variable <- factor(df.work3$variable, levels=c("modfc", "silt", "bio", 'bio1', 'evi', "bio15", "clay", "clay2", 'tri'),
                            labels = c("Cloud cover", "EVI spatial heterogeneity","Growing season precipitation", "Mean annual temperature",
                                       "Mean winter EVI", "Seasonality of precipitation", "Soil clay content", "Soil silt content", "Terrain ruggedness index"))

# Load plot function
plot.responseCurves_text3 <- function(data, var){
  
  my.title <- paste(var)
  
  df <- filter(data, variable == var)
  
  ggplot(df, aes(x=value, y=response, color = data)) + 
    geom_line(size = .3) +
    scale_color_manual(values=c("#313695", "#4575b4", "#fdae61", "#f46d43", "#a50026"), name = "Data Accuracy") +
    scale_linetype_manual(values = c("solid", "solid", "solid", "solid", "solid"), name = "Data Accuracy") +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 5, margin = margin(0,0,2,0)),
      legend.margin=margin(0,0,0,0),
      legend.box.margin=margin(-5,-5,-5,-5),
      axis.title = element_text(size=5, face = "bold"),
      axis.text.x = element_text(size=4),
      axis.text.y = element_text(size=4),
      legend.title = element_blank(),
      legend.text = element_text(size = 5),
      strip.background = element_rect(colour="transparent", fill="transparent"),
      strip.text.x = element_text(face = "bold.italic", size = 5),
      strip.text.y = element_text(face = "bold.italic", size = 5),
      panel.grid.minor = element_blank(),
      axis.ticks.length=unit(.07, "cm"),
      legend.position = "") +
    labs(title=my.title,
         x="Predictor value",
         y="Occurrence probability")
  
  ggsave(file=paste0("PathToYourFile/plot__rc_", var,".png"), dpi = 600, width = 4, height = 4.3, units = "cm")
}

plot.responseCurves_text3(df.work3.sum, "Growing season precipitation")
plot.responseCurves_text3(df.work3.sum, "Mean annual temperature")
plot.responseCurves_text3(df.work3.sum, "Seasonality of precipitation")
plot.responseCurves_text3(df.work3.sum, "Soil clay content")
plot.responseCurves_text3(df.work3.sum, "Soil silt content")
plot.responseCurves_text3(df.work3.sum, "Mean winter EVI")
plot.responseCurves_text3(df.work3.sum, "Cloud cover")
plot.responseCurves_text3(df.work3.sum, "EVI spatial heterogeneity")
plot.responseCurves_text3(df.work3.sum, "Terrain ruggedness index")

## Variable importance
# Load and process data
df.bird.varimp <- read.csv("PathToYourFile/var_importance_sum.csv", sep = ";")
df.bird.varimp$data <- factor(df.bird.varimp$data, levels=c("unaltered", "S1", "S2", "S3", "S4"), labels = c("Unaltered", "1-2km", "2-5km", "5-10km", "10-30km"))
df.bird.varimp$variables <- factor(df.bird.varimp$variables, levels=c("modfc", "silt", "bio", 'bio1', 'evi', "bio15", "clay", "clay2", 'tri'),
                                   labels = c("Cloud cover", "EVI spatial heterogeneity","Growing season precipitation", "Mean annual temperature",
                                              "Mean winter EVI", "Seasonality of precipitation", "Soil clay content", "Soil silt content", "Terrain ruggedness index"))

ggplot(df.bird.varimp, aes(x = data, y = AUCtest, fill = variables)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a")) +
  theme_bw() +
  theme(
    axis.title = element_text(size=5, face = "bold"),
    axis.text.x = element_text(size=4),
    axis.text.y = element_text(size=4),
    legend.title = element_text(size=5, face="bold"),
    legend.text = element_text(size = 5),
    strip.background = element_rect(colour="transparent", fill="transparent"),
    strip.text.x = element_text(face = "bold.italic", size = 5),
    strip.text.y = element_text(face = "bold.italic", size = 5),
    legend.position = "",) +
  labs(title="",
       x="Positional uncertainty",
       y="Variables’ importance")
ggsave(file="PathToYourFile/var_importance.png", dpi = 600, width = 5, height = 5, units = "cm")