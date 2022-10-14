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
preds <- stack("PathToYourFile/aspect.tif", "PathToYourFile/dem.tif", "PathToYourFile/forest.tif", "PathToYourFile/grassland.tif", "PathToYourFile/twi.tif")

# Create virtual Species ------------------------------------------------------------
# Define species response to aspoect and elevation
my.parametrs <- formatFunctions(dem = c(fun = 'dnorm', mean = 1000, sd = 100), # sd = 300 / 500 for species width medium  respectively wide niche
                                aspect = c(fun = 'dnorm', mean = 100, sd = 10)) # # sd = 30 / 50 for species width medium  respectively wide niche

# Generate virtual species
virtual.species <- generateSpFromFun(raster.stack = preds[[c(1,2)]],
                                     parameters = my.parametrs,
                                     species.type = "multiplicative",
                                     plot = T,
                                     rescale.each.response = FALSE)

# Convert probability of occurrence to presence-absence raster
PA.raster <- convertToPA(virtual.species, alpha = -0.05, beta = 0.3, plot = T)

# Define range of species occurences
occ.range <- c(20, 50, 100, 300, 1000)

# Prepare empty data frames for results
response.curves <- numeric()
var.importance <- numeric()
perf <- numeric()

for (i in 1:50) {
  for (j in occ.range){
    
    # Sample species occurrences
    PA.sampling <- sampleOccurrences(PA.raster, n = j, type = "presence only", plot = F)
    
    # Shift Occurrences ----------------------------------------------------------
    DATA.Prep <- as.data.frame(PA.sampling$sample.points) # Import original coordinates
    s <- nrow(DATA.Prep)
    ID <- 1:s
    Observed <- DATA.Prep$Observed
    
    # Shift data
    print(paste("Shifting data"))
    DATA.Prep[,c("S1x", "S1y")] <-  t(apply(DATA.Prep, 1, function(.) shift(.["x"], .["y"], dist=runif(1, 50, 100), mask=preds[[1]]))) # 50-100 meters
    DATA.Prep[,c("S2x", "S2y")] <-  t(apply(DATA.Prep, 1, function(.) shift(.["x"], .["y"], dist=runif(1, 100, 250), mask=preds[[1]]))) # 100-250 meters
    DATA.Prep[,c("S3x", "S3y")] <-  t(apply(DATA.Prep, 1, function(.) shift(.["x"], .["y"], dist=runif(1, 250, 500), mask=preds[[1]]))) # 250-500 meters
    DATA.Prep[,c("S4x", "S4y")] <-  t(apply(DATA.Prep, 1, function(.) shift(.["x"], .["y"], dist=runif(1, 500, 1500), mask=preds[[1]]))) # 500-1500 meters
    
    # prepare occurrences for SDMs ----------------------------------------------------------
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
    
    # Species distribution modeling ----------------------------------------------------------
    # Create sdmData object
    d0 <- sdmData(train = sp0, predictors = pred, bg=list(n=10000,method='gRandom',remove=TRUE))
    d1 <- sdmData(train = sp1, predictors = pred, bg=list(n=10000,method='gRandom',remove=TRUE))
    d2 <- sdmData(train = sp2, predictors = pred, bg=list(n=10000,method='gRandom',remove=TRUE))
    d3 <- sdmData(train = sp3, predictors = pred, bg=list(n=10000,method='gRandom',remove=TRUE))
    d4 <- sdmData(train = sp4, predictors = pred, bg=list(n=10000,method='gRandom',remove=TRUE))
    
    # Fit and evaluate models
    m0 <- sdm (occ~., data = d0, methods = "maxent", replication =  "cv", cv.folds=5, parallelSettings=list(ncore=1),
               modelSettings=list(maxent=list(beta=0.5,args=c('noquadratic', 'nothreshold', 'nolinear', 'noproduct'))))
    m1 <- sdm (occ~., data = d1, methods = "maxent", replication =  "cv", cv.folds=5, parallelSettings=list(ncore=1),
               modelSettings=list(maxent=list(beta=0.5, args=c('noquadratic', 'nothreshold', 'nolinear', 'noproduct'))))
    m2 <- sdm (occ~., data = d2, methods = "maxent", replication =  "cv", cv.folds=5, parallelSettings=list(ncore=1),
               modelSettings=list(maxent=list(beta=0.5, args=c('noquadratic', 'nothreshold', 'nolinear', 'noproduct'))))
    m3 <- sdm (occ~., data = d3, methods = "maxent", replication =  "cv", cv.folds=5, parallelSettings=list(ncore=1),
               modelSettings=list(maxent=list(beta=0.5, args=c('noquadratic', 'nothreshold', 'nolinear', 'noproduct'))))
    m4 <- sdm (occ~., data = d4, methods = "maxent", replication =  "cv", cv.folds=5, parallelSettings=list(ncore=1),
               modelSettings=list(maxent=list(beta=0.5, args=c('noquadratic', 'nothreshold', 'nolinear', 'noproduct'))))
    
    print(paste('Results'))
    
    # Combine results ----------------------------------------------------------
    # Performance metrics
    perf <- rbind(perf,
                  data.frame(sdm.evaluation(m0), data = "unaltered", n_occ = paste(j)),
                  data.frame(sdm.evaluation(m1), data = "S1", n_occ = paste(j)),
                  data.frame(sdm.evaluation(m2), data = "S2", n_occ = paste(j)),
                  data.frame(sdm.evaluation(m3), data = "S3", n_occ = paste(j)),
                  data.frame(sdm.evaluation(m4), data = "S4", n_occ = paste(j)))
    
    response.curves <- rbind(response.curves,
                             data.frame(extract_curves(getResponseCurve(m0)@response, nfolds = 5), data = "unaltered", n_occ = paste(j), id = rep_len(1:nrow(as.data.frame(getResponseCurve(m0)@response[1])), nrow(extract_curves(getResponseCurve(m0)@response, nfolds = 5)))),
                             data.frame(extract_curves(getResponseCurve(m1)@response, nfolds = 5), data = "S1", n_occ = paste(j), id = rep_len(1:nrow(as.data.frame(getResponseCurve(m0)@response[1])), nrow(extract_curves(getResponseCurve(m0)@response, nfolds = 5)))),
                             data.frame(extract_curves(getResponseCurve(m2)@response, nfolds = 5), data = "S2", n_occ = paste(j), id = rep_len(1:nrow(as.data.frame(getResponseCurve(m0)@response[1])), nrow(extract_curves(getResponseCurve(m0)@response, nfolds = 5)))),
                             data.frame(extract_curves(getResponseCurve(m3)@response, nfolds = 5), data = "S3", n_occ = paste(j), id = rep_len(1:nrow(as.data.frame(getResponseCurve(m0)@response[1])), nrow(extract_curves(getResponseCurve(m0)@response, nfolds = 5)))),
                             data.frame(extract_curves(getResponseCurve(m4)@response, nfolds = 5), data = "S4", n_occ = paste(j), id = rep_len(1:nrow(as.data.frame(getResponseCurve(m0)@response[1])), nrow(extract_curves(getResponseCurve(m0)@response, nfolds = 5)))))
    
    var.importance <- rbind(var.importance,
                            data.frame(getVarImp(m0)@varImportanceMean$corTest[,1:2], getVarImp(m0)@varImportanceMean$AUCtest[,1:2], data = "unaltered", n_occ = paste(j)),
                            data.frame(getVarImp(m1)@varImportanceMean$corTest[,1:2], getVarImp(m1)@varImportanceMean$AUCtest[,1:2], data = "S1", n_occ = paste(j)),
                            data.frame(getVarImp(m2)@varImportanceMean$corTest[,1:2], getVarImp(m2)@varImportanceMean$AUCtest[,1:2], data = "S2", n_occ = paste(j)),
                            data.frame(getVarImp(m3)@varImportanceMean$corTest[,1:2], getVarImp(m3)@varImportanceMean$AUCtest[,1:2], data = "S3", n_occ = paste(j)),
                            data.frame(getVarImp(m4)@varImportanceMean$corTest[,1:2], getVarImp(m4)@varImportanceMean$AUCtest[,1:2], data = "S4", n_occ = paste(j)))
    
    
    print(paste("Finished running set", i))
    
  }
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
response.curves.df.sum <- response.curves %>% dplyr::group_by(niche_width, n_occ, data, id, variable)  %>% 
  dplyr::summarise(
    value = mean(value, na.rm=T),
    response = mean(response, na.rm=T))

# Save summarized data
write.table(response.curves.df.sum, file = "PathToYourFile/response_curves_sum.csv", sep = ";")

## Model performance
df.perf <- read.csv(file = "PathToYourFile/performance_sum.csv", sep = ";")

# Get mean of values across all scenarios and loops
performance.df.sum <- df.perf %>% dplyr::group_by(data, niche_width, n_occ)  %>% 
  dplyr::summarise(
    TPR = mean(TPR, na.rm=T),
    TNR = mean(TNR, na.rm=T),
    FPR = mean(FPR, na.rm=T),
    FNR = mean(FNR, na.rm=T),
    Sensitivity = mean(Sensitivity, na.rm=T),
    Specificity = mean(Specificity, na.rm=T),
    TSS = mean(TSS, na.rm=T),
    Jaccard = mean(Jaccard, na.rm=T),
    Sorensen = mean(Sorensen, na.rm=T),
    F_measure = mean(F_measure, na.rm=T),
    OPR = mean(OPR, na.rm=T),
    UPR = mean(UPR, na.rm=T),
    AUC = mean(AUC, na.rm=T))

write.table(performance.df.sum, file = "PathToYourFile/performance_sum.csv", sep = ";")

## Variable importance
df.var_importance <- read.csv(file = "PathToYourFile/var_importance.csv", sep = ";")

# Get mean of values across all scenarios and loops
var_importance.sum <- df.var_importance %>% dplyr::group_by(data, niche_width, n_occ, variables)  %>% 
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
df.sum <- read.csv(file = "PathToYourFile/performance_sum.csv", sep = ";")
df.sum$data <- factor(df.sum$data, levels=c("unaltered", "S1", "S2", "S3", "S4"), labels = c("Unaltered", "50-100m", "100-250m", "250-500m", "500-1500m"))
df.sum$n_occ <- factor(df.sum$n_occ, levels=c("20", "100", "300", "1000"), labels = c("20", "100", "300", "1000"))
df.sum$niche_width <- factor(df.sum$niche_width, levels=c("Narrow", "Medium", "Wide"), labels=c("Narrow niche", "Medium niche", "Wide niche"))

# Sorensen index
ggplot(df.sum, aes(data, n_occ, fill = Sorensen)) + 
  geom_tile(color = "grey70", lwd = 0.3, linetype = 1)+
  geom_text(aes(label = round(Sorensen, 2)), color = "black", size = 2.7) +
  scale_fill_distiller(palette = "RdYlBu", direction = 1) +
  facet_grid(.~niche_width) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  coord_fixed() +
  guides(fill = guide_colourbar(barwidth = 0.6, barheight = 8)) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 10, hjust = 0.5, margin = margin(0,0,-1,0)),
    axis.title = element_text(size=9, face = "bold"),
    axis.text.x = element_text(size=8, angle = 290, vjust =-1),
    axis.text.y = element_text(size=8),
    legend.title = element_blank(),
    strip.background = element_rect(colour="transparent", fill="transparent"),
    strip.text.x = element_text(face = "bold.italic"),
    strip.text.y = element_text(face = "bold.italic"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()) +
  labs(title = "Sorensen index",
       x = "Positional uncertainty",
       y = "Number of occurrences")
ggsave(file="PathToYourFile/Sorensen.png", dpi = 600, width = 20, height = 20, units = "cm")

# Overprediction rate (OPR)
ggplot(df.sum, aes(data, n_occ, fill = OPR)) + 
  geom_tile(color = "grey70", lwd = .3, linetype = 1)+
  geom_text(aes(label = round(OPR, 2)), color = "black", size = 2.7) +
  scale_fill_distiller(palette = "RdYlBu", direction = 1) +
  facet_grid(.~niche_width) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  coord_fixed() +
  guides(fill = guide_colourbar(barwidth = 0.6, barheight = 8)) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 10, hjust = 0.5, margin = margin(0,0,-1,0)),
    axis.title = element_text(size=9, face = "bold"),
    axis.text.x = element_text(size=8, angle = 290, vjust =-1),
    axis.text.y = element_text(size=8),
    legend.title = element_blank(),
    strip.background = element_rect(colour="transparent", fill="transparent"),
    strip.text.x = element_text(face = "bold.italic"),
    strip.text.y = element_text(face = "bold.italic"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()) +
  labs(title = "Overprediction rate",
       x = "Positional uncertainty",
       y = "Number of occurrences")
ggsave(file="PathToYourFile/_OPR.png", dpi = 600, width = 20, height = 20, units = "cm")


# Underprediction rate (UPR)
ggplot(df.sum, aes(data, n_occ, fill = UPR)) + 
  geom_tile(color = "grey70", lwd = .3, linetype = 1)+
  geom_text(aes(label = round(UPR, 2)), color = "black", size = 2.7) +
  scale_fill_distiller(palette = "RdYlBu", direction = 1) +
  facet_grid(.~niche_width) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  coord_fixed() +
  guides(fill = guide_colourbar(barwidth = 0.6, barheight = 8)) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 10, hjust = 0.5, margin = margin(0,0,-1,0)),
    axis.title = element_text(size=9, face = "bold"),
    axis.text.x = element_text(size=8, angle = 290, vjust =-1),
    axis.text.y = element_text(size=8),
    legend.title = element_blank(),
    strip.background = element_rect(colour="transparent", fill="transparent"),
    strip.text.x = element_text(face = "bold.italic"),
    strip.text.y = element_text(face = "bold.italic"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()) +
  labs(title = "Underprediction rate",
       x = "Positional uncertainty",
       y = "Number of occurrences")
ggsave(file="PathToYourFile/_UPR.png", dpi = 600, width = 20, height = 20, units = "cm")

# TSS
ggplot(df.sum, aes(data, n_occ, fill = TSS)) + 
  geom_tile(color = "grey70", lwd = .3, linetype = 1)+
  geom_text(aes(label = round(TSS, 2)), color = "black", size = 2.7) +
  scale_fill_distiller(palette = "RdYlBu", direction = 1) +
  facet_grid(.~niche_width) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  coord_fixed() +
  guides(fill = guide_colourbar(barwidth = 0.6, barheight = 8)) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 10, hjust = 0.5, margin = margin(0,0,-1,0)),
    axis.title = element_text(size=9, face = "bold"),
    axis.text.x = element_text(size=8, angle = 290, vjust =-1),
    axis.text.y = element_text(size=8),
    legend.title = element_blank(),
    strip.background = element_rect(colour="transparent", fill="transparent"),
    strip.text.x = element_text(face = "bold.italic"),
    strip.text.y = element_text(face = "bold.italic"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()) +
  labs(title = "TSS",
       x = "Positional uncertainty",
       y = "Number of occurrences")
ggsave(file="PathToYourFile/TSS.png", dpi = 600, width = 20, height = 20, units = "cm")

## Response curves
# Load and process data
df.work2 <- read.csv(file = "PathToYourFile/response_curves_sum.csv", sep = ";")
df.work2$data <- factor(df.work2$data, levels=c("unaltered", "S1", "S2", "S3", "S4"), labels = c("Unaltered", "50-100m", "100-250m", "250-500m", "500-1500m"))
df.work2$n_occ <- factor(df.work2$n_occ, levels=c("20", "100", "300", "1000"), labels = c("20 occurrences", "100 occurrences", "300 occurrences", "1000 occurrences"))
df.work2$niche_width <- factor(df.work2$niche_width, levels=c("Narrow", "Medium", "Wide"), labels=c("Narrow niche", "Medium niche", "Wide niche"))

# Load function
plot.responseCurves2 <- function(data, niche){
  
  df <- dplyr::filter(data, niche_width == niche)
  
  my.title <- paste(niche)
  
  ggplot(df, aes(x=value, y=response, color = data)) +
    geom_line(size = .5) +
    scale_color_manual(values=c("#313695", "#4575b4", "#fdae61", "#f46d43", "#a50026"), name = "Positional uncertainty") +
    scale_linetype_manual(values = c("solid", "solid", "solid", "solid", "solid"), name = "Positional uncertainty") +
    facet_grid(n_occ~variable, scales="free") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size=10, face="bold"),
      axis.title.y=element_text(size=8, face="bold"),
      axis.title = element_text(size=9),
      axis.title.x=element_text(size=8, face="bold"),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      legend.title = element_text(size=8, face="bold"),
      strip.background = element_rect(colour="black", fill="grey80"),
      strip.text.x = element_text(face = "bold.italic", margin = margin(.08,0,.08,0, "cm")),
      strip.text.y = element_text(face = "bold.italic"),
      panel.grid.minor = element_blank(),
      legend.position = "bottom") +
    labs(title=my.title,
         x="Predictor value",
         y="Occurrence probability")
  ggsave(file=paste0("PathToYourFile/plot_", paste(niche), ".png"), dpi = 600, width = 18, height = 25, units = "cm")
}

plot.responseCurves2(data = df.work2.sum, niche = "Narrow niche")
plot.responseCurves2(data = df.work2.sum, niche = "Medium niche")
plot.responseCurves2(data = df.work2.sum, niche = "Wide niche")

## Variable importance
# Load and process data
df.sum <- read.csv(file = "PathToYourFile/var_importance_sum.csv", sep = ";")
df.sum$data <- factor(df.sum$data, levels=c("unaltered", "S1", "S2", "S3", "S4"), labels = c("Unaltered", "50-100m", "100-250m", "250-500m", "500-1500m"))
df.sum$n_occ <- factor(df.sum$n_occ, levels=c("30", "100", "300", "1000"), labels = c("20 occurrences", "100 occurrences", "300 occurrences", "1000 occurrences"))
df.sum$niche_width <- factor(df.sum$niche_width, levels=c("Narrow", "Medium", "Wide"), labels=c("Narrow niche", "Medium niche", "Wide niche"))
df.sum$variables <- factor(df.sum$variables, levels=c("aspect", "forest", "grassland", "dem", "twi"), labels=c("Aspect", "Amount of forest", "Amount of grassland", "Elevation", "Topography wetness index"))

# Load function
plot.varImportance <- function(data, niche, n.occurrences){
  
  df <- dplyr::filter(df.sum, niche_width == niche, n_occ == n.occurrences)
  
  ggplot(df, aes(x = data, y = AUCtest, fill = variables)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99")) +
    facet_grid(n_occ~niche_width) +
    theme_bw() +
    theme(
      axis.title = element_text(size=9, face = "bold"),
      axis.text.x = element_text(size=8, angle = 290, vjust =-1),
      axis.text.y = element_text(size=8),
      legend.title = element_blank(),
      legend.position = "",
      strip.background = element_rect(colour="transparent", fill="transparent"),
      strip.text.x = element_text(face = "bold.italic"),
      strip.text.y = element_text(face = "bold.italic")) +
    labs(title="",
         x="Positional uncertainty",
         y="Variables importance")
  ggsave(file=paste0("PathToYourFile/plot_", niche, "_", n.occurrences, ".png"), dpi = 600, width = 9, height = 10, units = "cm")
}

plot.varImportance(df, "Wide niche", "20 occurrences")
plot.varImportance(df, "Wide niche", "100 occurrences")
plot.varImportance(df, "Wide niche", "300 occurrences")
plot.varImportance(df, "Wide niche", "1000 occurrences")

plot.varImportance(df, "Medium niche", "20 occurrences")
plot.varImportance(df, "Medium niche", "100 occurrences")
plot.varImportance(df, "Medium niche", "300 occurrences")
plot.varImportance(df, "Medium niche", "1000 occurrences")

plot.varImportance(df, "Narrow niche", "20 occurrences")
plot.varImportance(df, "Narrow niche", "100 occurrences")
plot.varImportance(df, "Narrow niche", "300 occurrences")
plot.varImportance(df, "Narrow niche", "1000 occurrences")