#model trainer
freq <- read_csv("your_file_path.csv")
library(MASS) #for negative binomial models
library(pscl) #for zero inflated models
library(splines)
library(dplyr)
library(forcats) #for aggregating categorical variables


#freq = freq[sample(1:nrow(freq), 100000), ] #include to make models train faster and be smaller
tiny_freq = freq[sample(1:nrow(freq), 5000), ]
a_toy_model = glm(ClaimNb ~ DrivAge + BonusMalus, family = poisson, data = tiny_freq, offset = log(Exposure))
saveRDS(a_toy_model, "a_toy_model.rds")


#simple, full models
# simple poisson model with all predictors
pois_full <- glm(ClaimNb ~ Area + VehPower + VehAge + DrivAge + BonusMalus + VehBrand + VehGas + Density + Region, family = poisson, data = freq, offset = log(Exposure))

# simple negbin model with all predictors
negbin_full <- glm.nb(ClaimNb ~ Area + VehPower + VehAge + DrivAge + BonusMalus + VehBrand + VehGas + Density + Region + offset(log(Exposure)), data = freq)
negbin_full$data <- freq

# simple quasipoisson model with all predictors
quasi_full <- glm(ClaimNb ~ Area + VehPower + VehAge + DrivAge + BonusMalus + VehBrand + VehGas + Density + Region, family = quasipoisson, data = freq, offset = log(Exposure))

saveRDS(pois_full, "pois_full.rds")
saveRDS(negbin_full, "negbin_full.rds")
saveRDS(quasi_full, "quasi_full.rds")

print("Built and saved full models")



#some of the categorical predictors can be aggregated as they are not statistically significantly different
#Most VehBrands are not statistically distinguishable, and Area has 3 visible groups
freq <- freq %>%
  mutate(Area = case_when(
    Area %in% c("A", "B", "C") ~ "ABC",
    Area %in% c("D", "E") ~ "DE",
    TRUE ~ Area #leaves area F separate
  )) %>%
  mutate(VehBrand = case_when(
    VehBrand %in% c("B4",  "B2",  "B12", "B3",  "B1",  "B6",  "B13","B11") ~ "B0",
    TRUE ~ VehBrand #leaves brand 5, 10, 14 separate
  ))


# Quadratic splines for DrivAge with one knot capture age well and yield significant coefficients
# Density has a large VIF and is well captured, and region is also generally well captured, and too random to justify adding intercepts 

# medium sized poisson model with DrivAge quadratic spline
pois_med_quad <- glm(ClaimNb ~ Area + VehPower + VehAge + bs(DrivAge, knots = 27.5, degree = 2) + BonusMalus + VehBrand + VehGas, family = poisson, data = freq, offset = log(Exposure))

# medium sized negbin model with DrivAge quadratic spline
negbin_med_quad <- glm.nb(ClaimNb ~ Area + VehPower + VehAge + bs(DrivAge, knots = 27.5, degree = 2) + BonusMalus + VehBrand + VehGas + offset(log(Exposure)), data = freq)
negbin_med_quad$data <- freq

#medium sized quasipoisson model with DrivAge quadratic spline
quasi_med_quad <- glm(ClaimNb ~ Area + VehPower + VehAge + bs(DrivAge, knots = 27.5, degree = 2) + BonusMalus + VehBrand + VehGas, family = quasipoisson, data = freq, offset = log(Exposure))

saveRDS(pois_med_quad, "pois_med_quad.rds")
saveRDS(negbin_med_quad, "negbin_med_quad.rds")
saveRDS(quasi_med_quad, "quasi_med_quad.rds")

print("Built and saved med_quad models")

#zero inflated models can work well with the data, but are harder to work with
#they can also struggle to train, at least on my pc
#use freq <- freq[sample(1:nrow(freq), 100000), ] to reduce the scale of the model
zero_model <- zeroinfl(ClaimNb ~ Area + VehPower + VehAge + bs(DrivAge, knots = 27.5, degree = 2) + BonusMalus + VehBrand + VehGas + offset(log(Exposure)), dist = "poisson", data = freq)
zero_model$data <- freq
saveRDS(zero_model, "zero_model.rds")
print("Built and saved zero-inflated model")



##medium sized models with a linear spline for driver age
# # medium sized poisson model with DrivAge linear spline
# pois_med_linear <- glm(ClaimNb ~ Area + VehPower + VehAge + bs(DrivAge, knots = 27.5, degree = 1) + BonusMalus + VehBrand + VehGas, family = poisson, data = freq, offset = log(Exposure))
# 
# # medium sized negbin model with DrivAge linear spline
# negbin_med_linear = glm.nb(ClaimNb ~ Area + VehPower + VehAge + bs(DrivAge, knots = 27.5, degree = 1) + BonusMalus + VehBrand + VehGas + offset(log(Exposure)), data = freq)
# negbin_med_linear$data <- freq
# 
# #medium sized poisson model with DrivAge linear spline
# quasi_med_linear = glm(ClaimNb ~ Area + VehPower + VehAge + bs(DrivAge, knots = 27.5, degree = 1) + BonusMalus + VehBrand + VehGas, family = quasipoisson, data = freq, offset = log(Exposure))
# 
# saveRDS(pois_med_linear, "pois_med_linear.rds")
# saveRDS(negbin_med_linear, "negbin_med_linear.rds")
# saveRDS(quasi_med_linear, "quasi_med_linear.rds")
# 
# print("Built and saved med_linear models")

# #medium sized models with cubic splines
# # medium sized poisson model with DrivAge cubic spline
# pois_med_cub <- glm(ClaimNb ~ Area + VehPower + VehAge + bs(DrivAge, knots = 27.5, degree = 3) + BonusMalus + VehBrand + VehGas, family = poisson, data = freq, offset = log(Exposure))
# 
# # medium sized negbin model with DrivAge cubic spline
# negbin_med_cub <- glm.nb(ClaimNb ~ Area + VehPower + VehAge + bs(DrivAge, knots = 27.5, degree = 3) + BonusMalus + VehBrand + VehGas + offset(log(Exposure)), data = freq)
# negbin_med_cub$data <- freq
# 
# #medium sized quasipoisson model with DrivAge cubic spline
# quasi_med_cub <- glm(ClaimNb ~ Area + VehPower + VehAge + bs(DrivAge, knots = 27.5, degree = 3) + BonusMalus + VehBrand + VehGas, family = quasipoisson, data = freq, offset = log(Exposure))
# 
# saveRDS(pois_med_cub, "pois_med_cub.rds")
# saveRDS(negbin_med_cub, "negbin_med_cub.rds")
# saveRDS(quasi_med_cub, "quasi_med_cub.rds")
# 
# print("Built and saved med_cub models")




# 
# #small models with minimal predictors
# #small poisson model
# pois_small <- glm(ClaimNb ~ Area + VehAge + bs(DrivAge, knots = 27.5, degree = 1) + BonusMalus, family = poisson, data = freq, offset = log(Exposure))
# 
# #small negbin model
# negbin_small <- glm.nb(ClaimNb ~ Area + VehAge + bs(DrivAge, knots = 27.5, degree = 1) + BonusMalus + offset(log(Exposure)), data = freq)
# negbin_small$data <- freq
# 
# #small quasipoisson model
# quasi_small <- glm(ClaimNb ~ Area + VehAge + bs(DrivAge, knots = 27.5, degree = 1) + BonusMalus, family = quasipoisson, data = freq, offset = log(Exposure))
# 
# saveRDS(pois_small, "pois_small.rds")
# saveRDS(negbin_small, "negbin_small.rds")
# saveRDS(quasi_small, "quasi_small.rds")
# 
# print("Built and saved small models")

