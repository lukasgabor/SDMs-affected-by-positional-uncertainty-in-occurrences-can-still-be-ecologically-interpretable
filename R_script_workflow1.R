########################################################################################
#Code Authors:
# Lukas Gabor - https://scholar.google.cz/citations?user=pLQXY5wAAAAJ&hl=cs
# Vojtech Bartak - https://scholar.google.cz/citations?user=p8WAo8oAAAAJ&hl=cs&oi=ao
# Matthew S. Rogan - https://scholar.google.cz/citations?user=OTgp4V8AAAAJ&hl=cs&oi=ao

# Date: 10/14/2022
########################################################################################

# Because we had fitted over half milion models here we used batch processing to get results in reasonable amount of time.
# For more details see https://www.statmethods.net/interface/batch.html
# The script below thus contain also processes related with batch processing data preparation

# Worfklow 1 R script - virtual landscape + virtual species 

# Prepare batch job file ------------------------------------------------------------------

# Load packages
library(tibble)
library(emdbook)

# Create ranges for various SAC ranges, species niche width, number of occurrences and number of loops
SAC.ranges <- seq(1, 50, 2)
species.ranges <- lseq(0.005,0.2, 25)
occ.ranges <- c(20, 100, 300, 1000)
run.sets = 1:5 # Number of Loops

combinations <- expand.grid(SAC.ranges,species.ranges, occ.ranges, run.sets)

# Write batch job file
df <- tibble(
  text = paste(
    "module load R/4.2.0-foss-2020b;",
    "Rscript Matt.R",
    # always keep Rscript as it is command!, Matt.R is Rscript file with code
    combinations[, 1],
    combinations[, 2],
    combinations[, 3],
    combinations[, 4],
    sep = " "
  )
)

# Save the file
write.table(
  df,
  file = file.path("Alejandra_job_file8.txt"),
  append = FALSE,
  quote = FALSE,
  sep = " ",
  eol = "\n",
  na = "NA",
  dec = ".",
  row.names = FALSE,
  col.names = FALSE,
  qmethod = c("escape", "double"),
  fileEncoding = ""
)

## Batch job file example
# module load R/4.2.0-foss-2020b; Rscript Matt.R 1 0.005 20 1
# module load R/4.2.0-foss-2020b; Rscript Matt.R 3 0.005 20 1
# module load R/4.2.0-foss-2020b; Rscript Matt.R 5 0.005 20 1
# module load R/4.2.0-foss-2020b; Rscript Matt.R 7 0.005 20 1
# module load R/4.2.0-foss-2020b; Rscript Matt.R 9 0.005 20 1
# module load R/4.2.0-foss-2020b; Rscript Matt.R 11 0.005 20 1
# module load R/4.2.0-foss-2020b; Rscript Matt.R 13 0.005 20 1
# module load R/4.2.0-foss-2020b; Rscript Matt.R 15 0.005 20 1
# module load R/4.2.0-foss-2020b; Rscript Matt.R 17 0.005 20 1
# module load R/4.2.0-foss-2020b; Rscript Matt.R 19 0.005 20 1
# module load R/4.2.0-foss-2020b; Rscript Matt.R 21 0.005 20 1
# module load R/4.2.0-foss-2020b; Rscript Matt.R 23 0.005 20 1
# module load R/4.2.0-foss-2020b; Rscript Matt.R 25 0.005 20 1


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

setwd("/vast/palmer/scratch/jetz/lg763/Alejandra/outputs/virtual") # set my wd

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


# Get information about SAC value, niche width and number of occurrences from batch file
run_vars <- as.numeric(commandArgs(trailingOnly=TRUE))

SAC.range <- run_vars[1]
species.range <- run_vars[2]
occ.range <- run_vars[3]

run.set <- run_vars[4]

# Prepare empty data frames for results
response.curves <- numeric()
var.importance <- numeric()
perf <- numeric()

# Set size of artificial variables (in our case 200 * 200 cells)
xy <- expand.grid(1:200, 1:200)
names(xy) <- c('x','y')


print(paste('SAC:', SAC.range))
print(paste('Species:', species.range))
print(paste('Occ:', occ.range))

for (i in 1:5) {
  
  # Create artificial variables ----------------------------------------------------------
  g.dummy <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=1, model=vgm(psill=0.025, range= SAC.range, model='Exp'), nmax=20)
  
  pred <- predict(g.dummy, newdata=xy, nsim=1)
  gridded(pred) = ~x+y
  pred <- raster(pred) 
  pred <- rescale0to1(pred) # rescale values from 0-1 for purpose of fitting virtual species
  
  # Create virtual species ---------------------------------------------------------- 
  my.parametrs <- formatFunctions(sim1 = c(fun = 'dnorm', mean = 0.5, sd = species.range)) # define species response to variable
  
  # Generate virtual species ----------------------------------------------------------
  virtual.species <- generateSpFromFun(raster.stack = pred, 
                                       parameters = my.parametrs,
                                       species.type = "multiplicative",
                                       plot = F)
  
  # Convert probability of occurrence to presence-absence raster
  PA.raster <- convertToPA(virtual.species, alpha = -0.05, beta = 0.3, plot = F)
  
  # Sample species occurrences
  PA.sampling <- sampleOccurrences(PA.raster, n = occ.range, type = "presence only", plot = F)
  
  # Shift Occurrences ----------------------------------------------------------
  DATA.Prep <- as.data.frame(PA.sampling$sample.points) # Import original coordinates
  s <- nrow(DATA.Prep)
  ID <- 1:s
  Observed <- DATA.Prep$Observed
  
  # Shift data
  print(paste("Shifting data"))
  DATA.Prep[,c("S1x", "S1y")] <-  t(apply(DATA.Prep, 1, function(.) shift(.["x"], .["y"], dist=runif(1, 1, 2), mask=pred))) # 1-2 pixels
  DATA.Prep[,c("S2x", "S2y")] <-  t(apply(DATA.Prep, 1, function(.) shift(.["x"], .["y"], dist=runif(1, 2, 5), mask=pred))) # 2-5 pixels
  DATA.Prep[,c("S3x", "S3y")] <-  t(apply(DATA.Prep, 1, function(.) shift(.["x"], .["y"], dist=runif(1, 5, 10), mask=pred))) # 5-10 pixels
  DATA.Prep[,c("S4x", "S4y")] <-  t(apply(DATA.Prep, 1, function(.) shift(.["x"], .["y"], dist=runif(1, 10, 30), mask=pred))) # 10-30 pixels
  
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
             modelSettings=list(maxent=list(beta=0.5, args=c('noquadratic', 'nothreshold', 'nolinear', 'noproduct'))))
  m1 <- sdm (occ~., data = d1, methods = "maxent", replication =  "cv", cv.folds=5, parallelSettings=list(ncore=1),
             modelSettings=list(maxent=list(beta=0.5, args=c('noquadratic', 'nothreshold', 'nolinear', 'noproduct'))))
  m2 <- sdm (occ~., data = d2, methods = "maxent", replication =  "cv", cv.folds=5, parallelSettings=list(ncore=1),
             modelSettings=list(maxent=list(beta=0.5, args=c('noquadratic', 'nothreshold', 'nolinear', 'noproduct'))))
  m3 <- sdm (occ~., data = d3, methods = "maxent", replication =  "cv", cv.folds=5, parallelSettings=list(ncore=1),
             modelSettings=list(maxent=list(beta=0.5, args=c('noquadratic', 'nothreshold', 'nolinear', 'noproduct'))))
  m4 <- sdm (occ~., data = d4, methods = "maxent", replication =  "cv", cv.folds=5, parallelSettings=list(ncore=1),
             modelSettings=list(maxent=list(beta=0.5, args=c('noquadratic', 'nothreshold', 'nolinear', 'noproduct'))))
  
  print(paste('Results'))
  # get species prevalence
  raster <- PA.raster$pa.raster # virtual species presence / absence raster
  
  raster_data <- as.data.frame(rasterToPoints(raster)) # raster to data frame
  
  n_presences <- raster_data %>% dplyr::filter(layer == 1) # filter presence cells
  
  n_presences_sum <- nrow(n_presences) # get sum of presence cells
  
  n_cells_sum <- nrow(raster_data) # get sum of all cells
  
  spec.prevalence <- n_presences_sum / n_cells_sum # species prevalence (presence cells / all cells)
  
  # get sampling prevalence
  n_sampled_occurrences <- as.numeric(paste(occ.range)) # get number of sampled occurences 
  
  sampl.prevalence1 <- n_sampled_occurrences / n_cells_sum # sampling prevalence against all study area
  
  sampl.prevalence2 <- n_sampled_occurrences / n_presences_sum # sampling prevalence against number of presence cells
  
  # Combine results ----------------------------------------------------------
  # Performance metrics
  perf <- rbind(perf,
                data.frame(sdm.evaluation(m0), data = "unaltered", sac = paste(SAC.range), niche_width = paste(species.range), n_occ = paste(occ.range), species_prevalence = spec.prevalence, sampling_prevalence_studyArea = sampl.prevalence1, sampling_prevalence_presCells = sampl.prevalence2),
                data.frame(sdm.evaluation(m1), data = "S1", sac = paste(SAC.range), niche_width = paste(species.range), n_occ = paste(occ.range), species_prevalence = spec.prevalence, sampling_prevalence_studyArea = sampl.prevalence1, sampling_prevalence_presCells = sampl.prevalence2),
                data.frame(sdm.evaluation(m2), data = "S2", sac = paste(SAC.range), niche_width = paste(species.range), n_occ = paste(occ.range), species_prevalence = spec.prevalence, sampling_prevalence_studyArea = sampl.prevalence1, sampling_prevalence_presCells = sampl.prevalence2),
                data.frame(sdm.evaluation(m3), data = "S3", sac = paste(SAC.range), niche_width = paste(species.range), n_occ = paste(occ.range), species_prevalence = spec.prevalence, sampling_prevalence_studyArea = sampl.prevalence1, sampling_prevalence_presCells = sampl.prevalence2),
                data.frame(sdm.evaluation(m4), data = "S4", sac = paste(SAC.range), niche_width = paste(species.range), n_occ = paste(occ.range), species_prevalence = spec.prevalence, sampling_prevalence_studyArea = sampl.prevalence1, sampling_prevalence_presCells = sampl.prevalence2))
  
  # Response curves
  response.curves <- rbind(response.curves,
                           data.frame(extract_curves(getResponseCurve(m0)@response, nfolds = 5), data = "unaltered", sac = paste(SAC.range), niche_width = paste(species.range), n_occ = paste(occ.range), species_prevalence = spec.prevalence, sampling_prevalence_studyArea = sampl.prevalence1, sampling_prevalence_presCells = sampl.prevalence2, id = as.numeric(1:paste(occ.range))),
                           data.frame(extract_curves(getResponseCurve(m1)@response, nfolds = 5), data = "S1", sac = paste(SAC.range), niche_width = paste(species.range), n_occ = paste(occ.range), species_prevalence = spec.prevalence, sampling_prevalence_studyArea = sampl.prevalence1, sampling_prevalence_presCells = sampl.prevalence2, id = as.numeric(1:paste(occ.range))),
                           data.frame(extract_curves(getResponseCurve(m2)@response, nfolds = 5), data = "S2", sac = paste(SAC.range), niche_width = paste(species.range), n_occ = paste(occ.range), species_prevalence = spec.prevalence, sampling_prevalence_studyArea = sampl.prevalence1, sampling_prevalence_presCells = sampl.prevalence2, id = as.numeric(1:paste(occ.range))),
                           data.frame(extract_curves(getResponseCurve(m3)@response, nfolds = 5), data = "S3", sac = paste(SAC.range), niche_width = paste(species.range), n_occ = paste(occ.range), species_prevalence = spec.prevalence, sampling_prevalence_studyArea = sampl.prevalence1, sampling_prevalence_presCells = sampl.prevalence2, id = as.numeric(1:paste(occ.range))),
                           data.frame(extract_curves(getResponseCurve(m4)@response, nfolds = 5), data = "S4", sac = paste(SAC.range), niche_width = paste(species.range), n_occ = paste(occ.range), species_prevalence = spec.prevalence, sampling_prevalence_studyArea = sampl.prevalence1, sampling_prevalence_presCells = sampl.prevalence2, id = as.numeric(1:paste(occ.range))))
  
  print(paste("Finished running set", i))
  
} # End of loop


# Save results ------------------------------------------------------------
write.table(x = perf, file = paste0("PathToYourFile/performance_", SAC.range, "_", species.range, "_", occ.range, "_", run.set, ".csv"), row.names=F, col.names = T, sep = ";")
write.table(x = response.curves, file = paste0("PathToYourFile/response.curves_", SAC.range, "_", species.range, "_", occ.range, "_", run.set, ".csv"), row.names=F, col.names = T, sep = ";")

# Process results ------------------------------------------------------------
library(dplyr)

## Response curves
setwd("PathToYourFile")# Set working directory
temp <- list.files(pattern="*.csv") # Read all csv files
response.curves <- lapply(temp, read.table,sep = ";", header = T)
response.curves.df <- do.call(rbind.data.frame, response.curves) # Create dataframe with all results

# Get mean of values across all scenarios and loops
response.curves.df.sum <- response.curves.df %>% dplyr::group_by(niche_width, n_occ, data, sac, id)  %>% 
  dplyr::summarise(
    value = mean(value, na.rm=T),
    response = mean(response, na.rm=T),
    species_prevalence = mean(species_prevalence, na.rm=T),
    sampling_prevalence_studyArea = mean(sampling_prevalence_studyArea, na.rm=T),
    sampling_prevalence_presCells = mean(sampling_prevalence_presCells, na.rm=T))

# Save summarized data
write.table(response.curves.df.sum, file = "PathToYourFile/fileName", sep = ";")

## Model performance
setwd("PathToYourFile") # Set working directory
temp <- list.files(pattern="*.csv") # Read all csv files
performance <- lapply(temp, read.table,sep = ";", header = T)
performance.df <- do.call(rbind.data.frame, performance) # Create dataframe with all results

# Get mean of values across all scenarios and loops
performance.df.sum <- performance.df %>% dplyr::group_by(data, sac, niche_width, n_occ)  %>% 
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
    AUC = mean(AUC, na.rm=T),
    species_prevalence = mean(species_prevalence, na.rm=T),
    sampling_prevalence_studyArea = mean(sampling_prevalence_studyArea, na.rm=T),
    sampling_prevalence_presCells = mean(sampling_prevalence_presCells, na.rm=T))

write.table(performance.df.sum, file = "PathToYourFile/fileName", sep = ";")

# Plot results ------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(RColorBrewer)

## Model performance
# load and process data
df.sum <- read.csv(file = "PathToYourFile/fileName", sep = ";")
df.sum$data <- factor(df.sum$data,  levels = c("Unaltered", "PU 1-2 pixels", "PU 2-5 pixels", "PU 5-10 pixels", "PU 10-30 pixels"))
df.sum$n_occ <- factor(df.sum$n_occ,  levels = c("20 occurrences", "100 occurrences", "300 occurrences", "1000 occurrences"))
df.sum$niche_widthRound <- round(df.sum$niche_width, digits = 3)
df.sum$niche_widthRound <- factor(df.sum$niche_widthRound)

# Sorensen index
ggplot(df.sum, aes(sac, niche_widthRound, fill= Sorensen)) + 
  geom_tile()+
  scale_fill_distiller(palette = "RdYlBu", direction = 1) +
  scale_x_continuous(breaks = seq(5, 50, 15)) +
  scale_y_discrete(breaks = c("0.005", "0.017", "0.058", "0.2")) +
  facet_grid(n_occ~data) +
  guides(fill = guide_colourbar(barwidth = 0.6, barheight = 8)) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 8, margin = margin(0,0,5,0)),
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(-5,-5,-5,-5),
    axis.title = element_text(size=6, face = "bold"),
    axis.text.x = element_text(size=5),
    axis.text.y = element_text(size=5),
    legend.title = element_blank(),
    legend.text = element_text(size = 5),
    strip.background = element_rect(colour="transparent", fill="transparent"),
    strip.text.x = element_text(face = "bold.italic", size = 5),
    strip.text.y = element_text(face = "bold.italic", size = 5),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.ticks.length=unit(.07, "cm")) +
  labs(title = "Sorensen index",
       x = "Spatial autocorrelation",
       y = "Niche width")
ggsave(file=paste0("PathToYourFile/fileName.png"), dpi = 600, width = 11, height = 10, units = "cm")

# TSS
ggplot(df.sum, aes(sac, niche_widthRound, fill= TSS)) + 
  geom_tile()+
  scale_fill_distiller(palette = "RdYlBu", direction = 1) +
  scale_x_continuous(breaks = seq(5, 50, 15)) +
  scale_y_discrete(breaks = c("0.005", "0.017", "0.058", "0.2")) +
  facet_grid(n_occ~data) +
  guides(fill = guide_colourbar(barwidth = 0.6, barheight = 8)) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 7, margin = margin(0,0,5,0)),
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(-5,-5,-5,-5),
    axis.title = element_text(size=6, face = "bold"),
    axis.text.x = element_text(size=4),
    axis.text.y = element_text(size=4),
    legend.title = element_blank(),
    legend.text = element_text(size = 5),
    strip.background = element_rect(colour="transparent", fill="transparent"),
    strip.text.x = element_text(face = "bold.italic", size = 5),
    strip.text.y = element_text(face = "bold.italic", size = 5),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.ticks.length=unit(.07, "cm")) +
  labs(title = "TSS",
       x = "Spatial autocorrelation",
       y = "Niche width")
ggsave(file=paste0("PathToYourFile/fileName.png"), dpi = 600, width = 11, height = 10, units = "cm")

# Overprediction rate (OPR)
ggplot(df.sum, aes(sac, niche_widthRound, fill= OPR)) + 
  geom_tile()+
  scale_fill_distiller(palette = "RdYlBu", direction = 1) +
  scale_x_continuous(breaks = seq(5, 50, 15)) +
  scale_y_discrete(breaks = c("0.005", "0.017", "0.058", "0.2")) +
  facet_grid(n_occ~data) +
  guides(fill = guide_colourbar(barwidth = 0.6, barheight = 8)) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 8, margin = margin(0,0,5,0)),
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(-5,-5,-5,-5),
    axis.title = element_text(size=6, face = "bold"),
    axis.text.x = element_text(size=5),
    axis.text.y = element_text(size=5),
    legend.title = element_blank(),
    legend.text = element_text(size = 5),
    strip.background = element_rect(colour="transparent", fill="transparent"),
    strip.text.x = element_text(face = "bold.italic", size = 5),
    strip.text.y = element_text(face = "bold.italic", size = 5),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.ticks.length=unit(.07, "cm")) +
  labs(title = "Overprediction rate",
       x = "Spatial autocorrelation",
       y = "Niche width")
ggsave(file=paste0("PathToYourFile/fileName.png"), dpi = 600, width = 11, height = 10, units = "cm")

# Underprediction rate (UPR)
ggplot(df.sum, aes(sac, niche_widthRound, fill= UPR)) + 
  geom_tile()+
  scale_fill_distiller(palette = "RdYlBu", direction = 1) +
  scale_x_continuous(breaks = seq(5, 50, 15)) +
  scale_y_discrete(breaks = c("0.005", "0.017", "0.058", "0.2")) +
  facet_grid(n_occ~data) +
  guides(fill = guide_colourbar(barwidth = 0.6, barheight = 8)) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 8, margin = margin(0,0,5,0)),
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(-5,-5,-5,-5),
    axis.title = element_text(size=6, face = "bold"),
    axis.text.x = element_text(size=5),
    axis.text.y = element_text(size=5),
    legend.title = element_blank(),
    legend.text = element_text(size = 5),
    strip.background = element_rect(colour="transparent", fill="transparent"),
    strip.text.x = element_text(face = "bold.italic", size = 5),
    strip.text.y = element_text(face = "bold.italic", size = 5),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.ticks.length=unit(.07, "cm")) +
  labs(title = "Underprediction rate",
       x = "Spatial autocorrelation",
       y = "Niche width")
ggsave(file=paste0("PathToYourFile/fileName.png"), dpi = 600, width = 11, height = 10, units = "cm")

## Response curves
# Load and process data
df <- read.csv(file = "PathToYourFile/fileName.csv", sep = ";")
df$data <- factor(df$data, levels = c("Unaltered", "PU 1-2 pixels", "PU 2-5 pixels", "PU 5-10 pixels", "PU 10-30 pixels"))
df$n_occ <- factor(df$n_occ, levels = c("20 occurrences", "100 occurrences", "300 occurrences", "1000 occurrences"))

# Get unique species niche
niche.range <- unique(df$niche_width)

# Load plot function
plot.responseCurves <- function(data, niche, n.occurrences){
  
  df <- dplyr::filter(data, niche_width == niche, n_occ == n.occurrences)
  
  niche_title <- round(niche, digits = 3)
  
  my.title <- paste("Number of occurrences:",paste(n.occurrences),"/ Value of standart deviation:",paste(niche_title))
  
  ggplot(df, aes(x=value, y=response, color = data)) +
    geom_line(size = .5) +
    scale_color_manual(values=c("#313695", "#4575b4", "#fdae61", "#f46d43", "#a50026"), name = "Positional uncertainty") +
    scale_linetype_manual(values = c("solid", "solid", "solid", "solid", "solid"), name = "Positional uncertainty") +
    facet_wrap(.~sac, scales="free", ncol = 3) +
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
  ggsave(file=paste0("PathToYourFile/plot_", niche_title, "_", n.occurrences, ".png"), dpi = 600, width = 22, height = 35, units = "cm")
}

# Plot and save response curves 
plot.responseCurves(data = df, niche = niche.range[1], n.occurrences = "20 occurrences")
plot.responseCurves(data = df, niche = niche.range[1], n.occurrences = "100 occurrences")
plot.responseCurves(data = df, niche = niche.range[1], n.occurrences = "300 occurrences")
plot.responseCurves(data = df, niche = niche.range[1], n.occurrences = "1000 occurrences")

plot.responseCurves(data = df, niche = niche.range[2], n.occurrences = "20 occurrences")
plot.responseCurves(data = df, niche = niche.range[2], n.occurrences = "100 occurrences")
plot.responseCurves(data = df, niche = niche.range[2], n.occurrences = "300 occurrences")
plot.responseCurves(data = df, niche = niche.range[2], n.occurrences = "1000 occurrences")

plot.responseCurves(data = df, niche = niche.range[3], n.occurrences = "20 occurrences")
plot.responseCurves(data = df, niche = niche.range[3], n.occurrences = "100 occurrences")
plot.responseCurves(data = df, niche = niche.range[3], n.occurrences = "300 occurrences")
plot.responseCurves(data = df, niche = niche.range[3], n.occurrences = "1000 occurrences")

plot.responseCurves(data = df, niche = niche.range[4], n.occurrences = "20 occurrences")
plot.responseCurves(data = df, niche = niche.range[4], n.occurrences = "100 occurrences")
plot.responseCurves(data = df, niche = niche.range[4], n.occurrences = "300 occurrences")
plot.responseCurves(data = df, niche = niche.range[4], n.occurrences = "1000 occurrences")

plot.responseCurves(data = df, niche = niche.range[5], n.occurrences = "20 occurrences")
plot.responseCurves(data = df, niche = niche.range[5], n.occurrences = "100 occurrences")
plot.responseCurves(data = df, niche = niche.range[5], n.occurrences = "300 occurrences")
plot.responseCurves(data = df, niche = niche.range[5], n.occurrences = "1000 occurrences")

plot.responseCurves(data = df, niche = niche.range[6], n.occurrences = "20 occurrences")
plot.responseCurves(data = df, niche = niche.range[6], n.occurrences = "100 occurrences")
plot.responseCurves(data = df, niche = niche.range[6], n.occurrences = "300 occurrences")
plot.responseCurves(data = df, niche = niche.range[6], n.occurrences = "1000 occurrences")

plot.responseCurves(data = df, niche = niche.range[7], n.occurrences = "20 occurrences")
plot.responseCurves(data = df, niche = niche.range[7], n.occurrences = "100 occurrences")
plot.responseCurves(data = df, niche = niche.range[7], n.occurrences = "300 occurrences")
plot.responseCurves(data = df, niche = niche.range[7], n.occurrences = "1000 occurrences")

plot.responseCurves(data = df, niche = niche.range[8], n.occurrences = "20 occurrences")
plot.responseCurves(data = df, niche = niche.range[8], n.occurrences = "100 occurrences")
plot.responseCurves(data = df, niche = niche.range[8], n.occurrences = "300 occurrences")
plot.responseCurves(data = df, niche = niche.range[8], n.occurrences = "1000 occurrences")

plot.responseCurves(data = df, niche = niche.range[9], n.occurrences = "20 occurrences")
plot.responseCurves(data = df, niche = niche.range[9], n.occurrences = "100 occurrences")
plot.responseCurves(data = df, niche = niche.range[9], n.occurrences = "300 occurrences")
plot.responseCurves(data = df, niche = niche.range[9], n.occurrences = "1000 occurrences")

plot.responseCurves(data = df, niche = niche.range[10], n.occurrences = "20 occurrences")
plot.responseCurves(data = df, niche = niche.range[10], n.occurrences = "100 occurrences")
plot.responseCurves(data = df, niche = niche.range[10], n.occurrences = "300 occurrences")
plot.responseCurves(data = df, niche = niche.range[10], n.occurrences = "1000 occurrences")

plot.responseCurves(data = df, niche = niche.range[11], n.occurrences = "20 occurrences")
plot.responseCurves(data = df, niche = niche.range[11], n.occurrences = "100 occurrences")
plot.responseCurves(data = df, niche = niche.range[11], n.occurrences = "300 occurrences")
plot.responseCurves(data = df, niche = niche.range[11], n.occurrences = "1000 occurrences")

plot.responseCurves(data = df, niche = niche.range[12], n.occurrences = "20 occurrences")
plot.responseCurves(data = df, niche = niche.range[12], n.occurrences = "100 occurrences")
plot.responseCurves(data = df, niche = niche.range[12], n.occurrences = "300 occurrences")
plot.responseCurves(data = df, niche = niche.range[12], n.occurrences = "1000 occurrences")

plot.responseCurves(data = df, niche = niche.range[13], n.occurrences = "20 occurrences")
plot.responseCurves(data = df, niche = niche.range[13], n.occurrences = "100 occurrences")
plot.responseCurves(data = df, niche = niche.range[13], n.occurrences = "300 occurrences")
plot.responseCurves(data = df, niche = niche.range[13], n.occurrences = "1000 occurrences")

plot.responseCurves(data = df, niche = niche.range[14], n.occurrences = "20 occurrences")
plot.responseCurves(data = df, niche = niche.range[14], n.occurrences = "100 occurrences")
plot.responseCurves(data = df, niche = niche.range[14], n.occurrences = "300 occurrences")
plot.responseCurves(data = df, niche = niche.range[14], n.occurrences = "1000 occurrences")

plot.responseCurves(data = df, niche = niche.range[15], n.occurrences = "20 occurrences")
plot.responseCurves(data = df, niche = niche.range[15], n.occurrences = "100 occurrences")
plot.responseCurves(data = df, niche = niche.range[15], n.occurrences = "300 occurrences")
plot.responseCurves(data = df, niche = niche.range[15], n.occurrences = "1000 occurrences")

plot.responseCurves(data = df, niche = niche.range[16], n.occurrences = "20 occurrences")
plot.responseCurves(data = df, niche = niche.range[16], n.occurrences = "100 occurrences")
plot.responseCurves(data = df, niche = niche.range[16], n.occurrences = "300 occurrences")
plot.responseCurves(data = df, niche = niche.range[16], n.occurrences = "1000 occurrences")

plot.responseCurves(data = df, niche = niche.range[17], n.occurrences = "20 occurrences")
plot.responseCurves(data = df, niche = niche.range[17], n.occurrences = "100 occurrences")
plot.responseCurves(data = df, niche = niche.range[17], n.occurrences = "300 occurrences")
plot.responseCurves(data = df, niche = niche.range[17], n.occurrences = "1000 occurrences")

plot.responseCurves(data = df, niche = niche.range[18], n.occurrences = "20 occurrences")
plot.responseCurves(data = df, niche = niche.range[18], n.occurrences = "100 occurrences")
plot.responseCurves(data = df, niche = niche.range[18], n.occurrences = "300 occurrences")
plot.responseCurves(data = df, niche = niche.range[18], n.occurrences = "1000 occurrences")

plot.responseCurves(data = df, niche = niche.range[19], n.occurrences = "20 occurrences")
plot.responseCurves(data = df, niche = niche.range[19], n.occurrences = "100 occurrences")
plot.responseCurves(data = df, niche = niche.range[19], n.occurrences = "300 occurrences")
plot.responseCurves(data = df, niche = niche.range[19], n.occurrences = "1000 occurrences")

plot.responseCurves(data = df, niche = niche.range[20], n.occurrences = "20 occurrences")
plot.responseCurves(data = df, niche = niche.range[20], n.occurrences = "100 occurrences")
plot.responseCurves(data = df, niche = niche.range[20], n.occurrences = "300 occurrences")
plot.responseCurves(data = df, niche = niche.range[20], n.occurrences = "1000 occurrences")

plot.responseCurves(data = df, niche = niche.range[21], n.occurrences = "20 occurrences")
plot.responseCurves(data = df, niche = niche.range[21], n.occurrences = "100 occurrences")
plot.responseCurves(data = df, niche = niche.range[21], n.occurrences = "300 occurrences")
plot.responseCurves(data = df, niche = niche.range[21], n.occurrences = "1000 occurrences")

plot.responseCurves(data = df, niche = niche.range[22], n.occurrences = "20 occurrences")
plot.responseCurves(data = df, niche = niche.range[22], n.occurrences = "100 occurrences")
plot.responseCurves(data = df, niche = niche.range[22], n.occurrences = "300 occurrences")
plot.responseCurves(data = df, niche = niche.range[22], n.occurrences = "1000 occurrences")

plot.responseCurves(data = df, niche = niche.range[23], n.occurrences = "20 occurrences")
plot.responseCurves(data = df, niche = niche.range[23], n.occurrences = "100 occurrences")
plot.responseCurves(data = df, niche = niche.range[23], n.occurrences = "300 occurrences")
plot.responseCurves(data = df, niche = niche.range[23], n.occurrences = "1000 occurrences")

plot.responseCurves(data = df, niche = niche.range[24], n.occurrences = "20 occurrences")
plot.responseCurves(data = df, niche = niche.range[24], n.occurrences = "100 occurrences")
plot.responseCurves(data = df, niche = niche.range[24], n.occurrences = "300 occurrences")
plot.responseCurves(data = df, niche = niche.range[24], n.occurrences = "1000 occurrences")

plot.responseCurves(data = df, niche = niche.range[25], n.occurrences = "20 occurrences")
plot.responseCurves(data = df, niche = niche.range[25], n.occurrences = "100 occurrences")
plot.responseCurves(data = df, niche = niche.range[25], n.occurrences = "300 occurrences")
plot.responseCurves(data = df, niche = niche.range[25], n.occurrences = "1000 occurrences")