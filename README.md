# Comparative-GLM-Model-Explorer
A Shiny App to inspect and compare GLM models for the freMTPL2freq dataset, implementing strong, modular visualization and model diagnostic checks, together with a short R script to generate example GLMs.

To use the model, simply run the app.R (through your preferred means) and put your models in a models/ folder in the same directory

Assumptions and Notes:
 - Models should be poisson, quasipoisson, negative binomial (such as from MASS), or inherit from zeroinfl (from library pscl)
 - Place .rds model files in a "models/" folder next to this app.R (use saveRDS(myModel, "filename.rds") to save the model)
 - Models should be saved with their model/data embedded. if myModel <- glm(...), then defaults should work
 - If myModel <- glm.nb(...), or myModel <- zeroinfl(...) you may need to manually attach data with myModel$data <- my_data
 - Using saveRDS(myModel, "filename.rds") should work for saving the objects properly
 - Make sure to include all columns and not just the ones you are fitting on, so that you can inspect residuals and fits across those predictors
 - Expected columns in the model data (case sensitive, in any order):
     IDpol, ClaimNb, Exposure, Area, VehPower, VehAge, DrivAge, BonusMalus, VehBrand, VehGas, Density, Region
 - Categorical predictors: Area, VehBrand, VehGas, Region
 - Continuous predictors : VehPower, VehAge, DrivAge, BonusMalus, Density
 - You should be able to easily generalize this project to any set of predictors by changing the list of expected predictors in the R script

