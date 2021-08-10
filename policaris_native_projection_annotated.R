# Jesse Tabor
# April 15, 2020
# Hawaii non-native bee study using BIOMOD2
# Native range script

# This script creates species distribution models (SDMs) calibrated in the native range.
# We create the SDMs using species occurrence data (longitude, latitude) from gbif.org
# and current bioclimatic variables downloaded from Worldclim.org. We project these models 
# onto current climate scenarios in Hawaii to predict suitable habitat for non-native bees. 
# We also project models onto future climate scenarios in Hawaii to predict future suitable habitat.

# Set where the program saves and looks for data.
setwd("/Users/JTAdmin/Desktop/policaris")
list.files()

# Load libraries into the R environment.
library(rgbif)
library(biomod2)
library(ggplot2)
library(gridExtra)
library(ade4)
library(rgdal)
library(raster)
library(maptools)
library(rasterVis)
library(latticeExtra)
library(lattice)
library(sp)
library(dplyr)
library(devtools)
print("libraries loaded")

# Create a folder in your working directory and write your plots into that folder.
dir.create("mpoli_plots")

# Load species occurence data into R
bee <- read.table("mpoli2.csv",
                  header=TRUE, 
                  sep=",",
                  row.names=NULL)
print("species occurences loaded")

# Load your native range bioclimatic variables into R.
bio_1 <-raster("native_range_environment/mpoli_env/bioclim/bio_1.grd")
bio_2 <-raster("native_range_environment/mpoli_env/bioclim/bio_2.grd")
bio_3 <-raster("native_range_environment/mpoli_env/bioclim/bio_3.grd")
bio_4 <-raster("native_range_environment/mpoli_env/bioclim/bio_4.grd")
bio_5 <-raster("native_range_environment/mpoli_env/bioclim/bio_5.grd")
bio_6 <-raster("native_range_environment/mpoli_env/bioclim/bio_6.grd")
bio_7 <-raster("native_range_environment/mpoli_env/bioclim/bio_7.grd")
bio_8 <-raster("native_range_environment/mpoli_env/bioclim/bio_8.grd")
bio_9 <-raster("native_range_environment/mpoli_env/bioclim/bio_9.grd")
bio_10 <-raster("native_range_environment/mpoli_env/bioclim/bio_10.grd")
bio_11 <-raster("native_range_environment/mpoli_env/bioclim/bio_11.grd")
bio_12 <-raster("native_range_environment/mpoli_env/bioclim/bio_12.grd")
bio_13 <-raster("native_range_environment/mpoli_env/bioclim/bio_13.grd")
bio_14 <-raster("native_range_environment/mpoli_env/bioclim/bio_14.grd")
bio_15 <-raster("native_range_environment/mpoli_env/bioclim/bio_15.grd")
bio_16 <-raster("native_range_environment/mpoli_env/bioclim/bio_16.grd")
bio_17 <-raster("native_range_environment/mpoli_env/bioclim/bio_17.grd")
bio_18 <-raster("native_range_environment/mpoli_env/bioclim/bio_18.grd")
print("native range raster loaded")

# Stack native range rasters in raster stack.
bioclim_NA <- stack(c(bio_1, bio_10, bio_11,
                      bio_12, bio_13, bio_14,
                      bio_15, bio_16, bio_17,
                      bio_18, bio_2,
                      bio_3, bio_4, bio_5,
                      bio_6, bio_7, bio_8,
                      bio_9))
print("native range raster stacked")

# --------------------------------------------------------------------------------Principal component analysis (PCA)---------------------------------------------------------------------------------------------

# To reduce multicollinearity within the environmental variables, a principal component analysis is conducted to highlight the relationship between the target species 
# occurrences and the specific environmental combinations. Use the PCA plot to select a set of variables that are not too colinear (two variables pointing in
# orthogonal directions are independent, two variables pointing in the same or opposite directions are highly dependent, positively or negatively, respectively),
# and significantly contribute to the overall environmental variation (the longer the arrow, the more important the variable).

# Obtain identifiers of the cells where species occurs (lat/long).
points_bee <- data.frame(bee)[1:691, c("Longitude","Latitude")]
print("identifiers ready")

# Extract the value from each cell and put it in one central location.
bee_cell_id <- cellFromXY(subset(bioclim_NA,1),points_bee)
print("extracted cell id")

# Convert raster object into a data frame and remove non-defined area from the dataset, which gives values to every point.
bioclim_NA_df <- na.omit(as.data.frame(bioclim_NA))
print("removed non-defined area complete")

# Perform a principal component analysis (PCA) on environment to avoid model overfitting by decreasing the amount of variables used.
pca_NA <- dudi.pca(bioclim_NA_df, scannf = F, nf = 2)
print("pca complete")

# Write PCA plot to plots folder and determine which variables to use in SDMs
jpeg("/Users/JTAdmin/Desktop/policaris/mpoli_plots/mpoli_pca.jpeg", type="cairo")
s.corcircle(pca_NA$co,clabel = 1)
dev.off()

# ------------------------------------------------------------------BIOMOD_Modeling()------------------------------------------------------------------------------------------------------

# Create our initial native habitat/native range/native niche models that all other models will be built from.

# Stack your subset native range rasters into a new rasterstack.
bioclim_NA <- stack(c(bio_2,bio_4,bio_11,bio_12))
print("subset stack native rasters complete")

# Rearrange input data (species occurence data and native range raster) to make sure they can be used within biomod2. 
# The function allows to select pseudo-absences or background data in the case that true absences data are not available,
# or to add pseudo-absence data to an existing set of absence.
bee_data <- BIOMOD_FormatingData(resp.var = rep(1,nrow(bee)),
                                 expl.var = bioclim_NA,
                                 resp.xy = bee[,c('Longitude','Latitude')],
                                 resp.name = "mpoli.nat",
                                 PA.nb.rep = 3,PA.nb.absences = 10000,
                                 PA.strategy = 'random')
print("BIOMOD_FormattingData complete")

# Check to see if BIOMOD_FormatingData() function worked.
bee_data

# Plot bee_data and save to the plots folder
jpeg("/Users/JTAdmin/Desktop/policaris/mpoli_plots/bee_data.jpeg", type="cairo")
plot(bee_data)
dev.off()

# Parametrize and/or tune biomod’s single models options with BIOMOD_ModelingOptions() function.
bee_opt <- BIOMOD_ModelingOptions(GLM = list(type = 'quadratic',
                                             interaction.level = 1),
                                  GAM = list(algo ='GAM_mgcv'),
                                  GBM = list(n.trees = 1000),
                                  MARS = list(type = 'quadratic'))
print("BIOMOD_ModelingOptions complete")

# Next we run 9 different models on bee_data paired indavidually with 3 different pseudo-absence datasets
# (created with BIOMOD_FormatingData() function). We do this 4 times. This creates 108 total models.
# This function allows to calibrate and evaluate a range of species distribution models techniques
# run over a given species. Calibrations are made on the whole sample or a random subpart. The
# predictive power of the different models is estimated using a range of evaluation metrics.
bee_models <- BIOMOD_Modeling(data = bee_data, models = c("GLM","GAM","GBM","MARS","RF","ANN","FDA","CTA","MAXENT.Phillips"),
                              models.options = bee_opt, 
                              NbRunEval = 4, DataSplit = 80,
                              VarImport = 3, do.full.models = F,
                              modeling.id = "ex.2")
print("BIOMOD_Modeling complete")

# Get model evaluation scores 
bee_models_scores <- get_evaluations(bee_models)

# dim() function: Get or set the number of rows, columns, and layers of a Raster* object.
dim(bee_models_scores)

# Retrieve or set the dimnames of an object.
dimnames(bee_models_scores)

# model_scores_graph() function: This function is a graphic tool to represent evaluation scores of models
# produced with biomod2 according to 2 different evaluation methods. Models can be grouped in several ways
# (by algo, by CV run, ...) to highlight potential differences in models quality due to chosen models,
# or cross validation sampling bias. Points represent mean of evaluation score for a given condition and
# the lines represent the associated standard deviations.

# Model scores. Save the plot in the plots folder.
jpeg("/Users/JTAdmin/Desktop/policaris/mpoli_plots/native_model_scores.jpeg", type="cairo")
models_scores_graph(bee_models, by = "models",metrics = c("ROC","TSS"),
                    xlim = c(0.1,1), ylim = c(0.1,1))
dev.off()

# Get variable importance in the models
(bee_models_var_import <- get_variables_importance(bee_models))

# Cross-validation RUN scores. Save the plot in the plots folder.
jpeg("/Users/JTAdmin/Desktop/policaris/mpoli_plots/native_run_scores.jpeg", type="cairo")
models_scores_graph(bee_models, by = "cv_run",metrics = c("ROC","TSS"),
                    xlim = c(0.1,1), ylim = c(0.1,1))
dev.off()

# Pseudo-absence scores. Save the plot in the plots folder.
jpeg("/Users/JTAdmin/Desktop/policaris/mpoli_plots/native_pa_scores.jpeg", type="cairo")
models_scores_graph(bee_models, by = "data_set",metrics = c("ROC","TSS"),
                    xlim = c(0.1,1), ylim = c(0.1,1))
dev.off()

# Calculate the MEAN of variable importance by algorighm. This function tells us the importance of the individual
# variables in the model.
apply(bee_models_var_import, c(1,2), mean)

# Analyze how each environmental variable influence the species probability of prescence by creating response curve plots.
# These are graphical visualizations of the response curve of each varaible. How does each enviromental variable influence
# probability of prescence? Each line cooresponds to a diffrent model. To do this we first have to load the produced models.
bee_glm <- BIOMOD_LoadModels(bee_models, models = 'GLM')
bee_gbm <- BIOMOD_LoadModels(bee_models, models = 'GBM')
bee_rf <- BIOMOD_LoadModels(bee_models, models = 'RF')
bee_cta <- BIOMOD_LoadModels(bee_models, models = 'CTA')
bee_ann <- BIOMOD_LoadModels(bee_models, models = 'ANN')
bee_fda <- BIOMOD_LoadModels(bee_models, models = 'FDA')
bee_mars <- BIOMOD_LoadModels(bee_models, models = 'MARS')
bee_gam <- BIOMOD_LoadModels(bee_models, models = 'GAM')
bee_phi <- BIOMOD_LoadModels(bee_models, models = 'MAXENT.Phillips')

# biomod2::response.plot2() function creates the plot. Save the plot in the plots folder.
jpeg("/Users/JTAdmin/Desktop/policaris/mpoli_plots/glm.jpeg", type="cairo")
glm_eval_strip <- biomod2::response.plot2(models = bee_glm, Data = get_formal_data(bee_models,'expl.var'),
                                          show.variables = get_formal_data(bee_models,'expl.var.names'),
                                          do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE,
                                          display_title = FALSE, data_species = get_formal_data(bee_models, 'resp.var'))
dev.off()

# biomod2::response.plot2() function creates the plot. Save the plot in the plots folder.
jpeg("/Users/JTAdmin/Desktop/policaris/mpoli_plots/gbm.jpeg", type="cairo")
gbm_eval_strip <- biomod2::response.plot2(models = bee_gbm, Data = get_formal_data(bee_models, 'expl.var'),
                                          show.variables = get_formal_data(bee_models, 'expl.var.names'),
                                          do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE,
                                          display_title = FALSE, data_species = get_formal_data(bee_models, 'resp.var'))
dev.off()

# biomod2::response.plot2() function creates the plot. Save the plot in the plots folder.
jpeg("/Users/JTAdmin/Desktop/policaris/mpoli_plots/rf.jpeg", type="cairo")
rf_eval_strip <- biomod2::response.plot2(models = bee_rf, Data = get_formal_data(bee_models, 'expl.var'),
                                         show.variables = get_formal_data(bee_models, 'expl.var.names'),
                                         do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE,
                                         display_title = FALSE, data_species = get_formal_data(bee_models, 'resp.var'))
dev.off()

# biomod2::response.plot2() function creates the plot. Save the plot in the plots folder.
jpeg("/Users/JTAdmin/Desktop/policaris/mpoli_plots/gam.jpeg", type="cairo")
gam_eval_strip <- biomod2::response.plot2(models = bee_gam, Data = get_formal_data(bee_models, 'expl.var'),
                                          show.variables = get_formal_data(bee_models, 'expl.var.names'),
                                          do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE,
                                          display_title = FALSE, data_species = get_formal_data(bee_models, 'resp.var'))
dev.off()

# biomod2::response.plot2() function creates the plot. Save the plot in the plots folder.
jpeg("/Users/JTAdmin/Desktop/policaris/mpoli_plots/cta.jpeg", type="cairo")
cta_eval_strip <- biomod2::response.plot2(models = bee_cta, Data = get_formal_data(bee_models, 'expl.var'),
                                          show.variables = get_formal_data(bee_models, 'expl.var.names'),
                                          do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE,
                                          display_title = FALSE, data_species = get_formal_data(bee_models, 'resp.var'))
dev.off()

# biomod2::response.plot2() function creates the plot. Save the plot in the plots folder.
jpeg("/Users/JTAdmin/Desktop/policaris/mpoli_plots/ann.jpeg", type="cairo")
ann_eval_strip <- biomod2::response.plot2(models = bee_ann, Data = get_formal_data(bee_models, 'expl.var'),
                                          show.variables = get_formal_data(bee_models, 'expl.var.names'),
                                          do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE,
                                          display_title = FALSE, data_species = get_formal_data(bee_models, 'resp.var'))
dev.off()

# biomod2::response.plot2() function creates the plot. Save the plot in the plots folder.
jpeg("/Users/JTAdmin/Desktop/policaris/mpoli_plots/fda.jpeg", type="cairo")
fda_eval_strip <- biomod2::response.plot2(models = bee_fda, Data = get_formal_data(bee_models, 'expl.var'),
                                          show.variables = get_formal_data(bee_models, 'expl.var.names'),
                                          do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE,
                                          display_title = FALSE, data_species = get_formal_data(bee_models, 'resp.var'))
dev.off()

# biomod2::response.plot2() function creates the plot. Save the plot in the plots folder.
jpeg("/Users/JTAdmin/Desktop/policaris/mpoli_plots/mars.jpeg", type="cairo")
mars_eval_strip <- biomod2::response.plot2(models = bee_mars, Data = get_formal_data(bee_models, 'expl.var'),
                                           show.variables = get_formal_data(bee_models, 'expl.var.names'),
                                           do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE,
                                           display_title = FALSE, data_species = get_formal_data(bee_models, 'resp.var'))
dev.off()

# biomod2::response.plot2() function creates the plot. Save the plot in the plots folder.
jpeg("/Users/JTAdmin/Desktop/policaris/mpoli_plots/maxent.jpeg", type="cairo")
phi_eval_strip <- biomod2::response.plot2(models = bee_phi, Data = get_formal_data(bee_models, 'expl.var'),
                                          show.variables = get_formal_data(bee_models, 'expl.var.names'),
                                          do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE,
                                          display_title = FALSE, data_species = get_formal_data(bee_models, 'resp.var'))
dev.off()
print("response plots complete")

# Congratulations you have built 108 SDMs calibrated from the native range and
# analyzed the importance of algorithms, pseudo-absence, runs, and variables  in model performance. 
# Now we we will create an ensemble by averaging  all 108 SDMs into 1 final ensemble model. 

# ----------------------------------------------------------------------Ensemble modeling--------------------------------------------------------------------------------------------------

#BIOMOD_EnsembleModeling combines models and make ensemble predictions built with BIOMOD_Modeling.
#The ensemble predictions can also be evaluated against the original data given to BIOMOD_Modeling.
#Biomod2 proposes a range of options to build ensemble models and predictions and to assess the
#modeling uncertainty. The created ensemble models can then be used to project distributions over
#space and time as classical biomod2 models.

# In this case to reduce the number of outputs we only concider two "ensembling" options:
# Committee averaging and weighted mean. We also produce coefficent of variation that tell us the extent the models agree or disagree. 
# In this case we made a decision to mix  all models (i.e. all algorithms, all pseudo-absence sampling, all cross-validation runs) 
# to produce our ensemble models. TSS is used as the evaluation reference for committee building and defining weights. 
# In this case only models with a TSS greater than or equal to 0.8 are kept to build the final ensemble model. 
bee_ensemble_models <- BIOMOD_EnsembleModeling(modeling.output = bee_models,
                                               em.by = 'all',
                                               eval.metric = 'TSS',
                                               eval.metric.quality.threshold = 0.5,
                                               models.eval.meth = c('KAPPA','TSS','ROC'),
                                               prob.mean = FALSE,
                                               prob.cv = TRUE, 
                                               committee.averaging = TRUE,
                                               prob.mean.weight = TRUE,
                                               VarImport = 0)

# Now check the scores for the ensemble models
(bee_ensemble_models_scores <- get_evaluations(bee_ensemble_models))
print("native ensemble modeling complete")

# ----------------------------------------------------------------------------Current native range projection------------------------------------------------------------------------------

# For all the models currently implemented, BIOMOD_Projection() function is able to project potential distributions of
# species in other areas, other resolutions or other time scales. We now project all 108 current individual native range SDMs.
bee_models_proj_current <- BIOMOD_Projection(modeling.output = bee_models,
                                             new.env = bioclim_NA,
                                             proj.name = "current", 
                                             binary.meth = 'TSS', 
                                             do.stack = FALSE)
print("native BIOMOD_Projection complete")

# This function use projections of ‘individual models’ and ensemble models from BIOMOD_EnsembleModeling
# to build an ensemble of species’ projections over space and time. We now project current native range ensemble SDMs.
bee_ensemble_models_proj_current <- BIOMOD_EnsembleForecasting(EM.output = bee_ensemble_models,
                                                               projection.output = bee_models_proj_current,
                                                               binary.meth = "TSS",
                                                               output.format = ".img",
                                                               do.stack = FALSE)
print("native BIOMOD_EnsembleForcasting complete")

# Stack the native range ensemble SDMs.
stk_bee_ensemble_models_proj_current <- get_predictions(bee_ensemble_models_proj_current)

# Keep committee averaging and weighted mean ensamble models only.
stk_bee_ensemble_models_proj_current <- subset(stk_bee_ensemble_models_proj_current, grep("EMca|EMwmean", names(stk_bee_ensemble_models_proj_current)))

# Simplify the layer names for plotting conveniences.
names(stk_bee_ensemble_models_proj_current) <-sapply(strsplit(names(stk_bee_ensemble_models_proj_current), "_"), getElement, 2)

# Plot the committee averaging and weighted mean current native range ensemble models and write to plots folder.
jpeg("/Users/JTAdmin/Desktop/policaris/mpoli_plots/native_ensemble_projection.jpeg", type="cairo")
levelplot(stk_bee_ensemble_models_proj_current, main = "A. mellifera current ensemble native projection",
          col.regions = colorRampPalette(c("grey90","yellow4","green4"))(100))
dev.off()
print("current native range projection complete")

# -----------------------------------------------------Project current native range onto current Hawaii environment------------------------------------------------------------------------

# We project the "native niche" that the species fills in its "native range" onto the "Hawaiian environment" to predict suitable habitat in "Hawaii" for that species.

# Load our Hawaii range bioclimatic variables into R.
bio_1 <-raster("haw_raster_new/bio_1.grd")
bio_2 <-raster("haw_raster_new/bio_2.grd")
bio_3 <-raster("haw_raster_new/bio_3.grd")
bio_4 <-raster("haw_raster_new/bio_4.grd")
bio_5 <-raster("haw_raster_new/bio_5.grd")
bio_6 <-raster("haw_raster_new/bio_6.grd")
bio_7 <-raster("haw_raster_new/bio_7.grd")
bio_8 <-raster("haw_raster_new/bio_8.grd")
bio_9 <-raster("haw_raster_new/bio_9.grd")
bio_10 <-raster("haw_raster_new/bio_10.grd")
bio_11 <-raster("haw_raster_new/bio_11.grd")
bio_12 <-raster("haw_raster_new/bio_12.grd")
bio_13 <-raster("haw_raster_new/bio_13.grd")
bio_14 <-raster("haw_raster_new/bio_14.grd")
bio_15 <-raster("haw_raster_new/bio_15.grd")
bio_16 <-raster("haw_raster_new/bio_16.grd")
bio_17 <-raster("haw_raster_new/bio_17.grd")
bio_18 <-raster("haw_raster_new/bio_18.grd")
bio_19 <-raster("haw_raster_new/bio_19.grd")

# Stack Hawaii range rasters in raster stack.
bioclim_NA_current <- stack(c(bio_2,bio_4,bio_11,bio_12))
print("Hawaii raster stacked")

# For all the models currently implemented, BIOMOD_Projection() function is able to project potential distributions of
# species in other areas, other resolutions or other time scales. We now project all 108 current individual Hawaii range SDMs.
mpoli_predict <- BIOMOD_Projection(modeling.output = bee_models,
                                   new.env = bioclim_NA_current,
                                   proj.name = "mpoli_predict", 
                                   binary.meth = 'TSS', 
                                   do.stack = FALSE)
print("Hawaii mpoli BIOMOD_Projecion complete")

# This function use projections of ‘individual models’ and ensemble models from BIOMOD_EnsembleModeling
# to build an ensemble of species’ projections over space and time. We now project current Hawaii range ensemble SDMs.
bee_ensemble_models_predict <- BIOMOD_EnsembleForecasting(EM.output = bee_ensemble_models,
                                                          projection.output = mpoli_predict,
                                                          binary.meth = "TSS",
                                                          output.format = ".img",
                                                          do.stack = FALSE)
print("Hawaii mpoli BIOMOD_EnsembleForcasting complete")

# Stack the ensemble models.
stk_bee_ensemble_models_predict <- get_predictions(bee_ensemble_models_predict)

# Keep committee averaging and weighted mean ensemble models only.
stk_bee_ensemble_models_predict <- subset(stk_bee_ensemble_models_predict, grep("EMca|EMwmean", names(stk_bee_ensemble_models_predict)))

# Simplify the layer names for plotting conveniences.
names(stk_bee_ensemble_models_predict) <-sapply(strsplit(names(stk_bee_ensemble_models_predict), "_"), getElement, 2)

# Plot the committee averaging and weighted mean current Hawaii range ensemble models and write to plots folder.
jpeg("/Users/JTAdmin/Desktop/policaris/mpoli_plots/native_hawaii_projection.jpeg", type="cairo")
levelplot(stk_bee_ensemble_models_predict, main = "A. mellifera current Hawaii ensemble projection calibrated witih native range",
          col.regions = colorRampPalette(c("grey90","yellow4","green4"))(100))
dev.off()
print("Hawaii range hawaii projection complete")