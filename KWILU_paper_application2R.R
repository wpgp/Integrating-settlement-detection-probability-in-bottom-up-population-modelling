#  This R script is designed to help interested users to replicate 
#  the Bayesian hierarchical population modelling methods used for Kwilu Province, DRC
#  Author: Dr Chris Nnanatu, WOrldPop, University of Southampton, October 2024

rm(list = ls()) # clear workspace
rm(list=ls()) #----Clear the workspace

# Load key libraries
packages <- c("raster", "haven", "sf","sp", "tmap","tidyverse",
              "lattice", "gridExtra", "devtools", "rlang", "DClusterm",
              "viridis", "tmaptools", "spdep", "ggplot2", "ggpubr",
              "psych", "knitr", "car", "MASS", "terra")


# install.packages("psych", type="binary")
#if(length(setdiff(packages, rownames(installed.packages()))) > 0) { 
#install.packages(setdiff(packages, rownames(installed.packages(type="binary")))) }


#Install INLA
#if(length(setdiff("INLA", rownames(installed.packages()))) > 0){
#  install.packages("INLA", type="binary", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
#}
library(INLA)
lapply(packages, library, character.only = TRUE) ##--access the libraries

#Specify Drive Path
path <- "//worldpop.files.soton.ac.uk/Worldpop/Projects/WP517763_GRID3_Phase2/Working/DRC/CHRIS_N/province/Kwilu/paper/application"

data_path <- paste0(path, "/data/")#---paths for survey data: 
results_path <- paste0(path , "/output/")

#------------------------------------------------------------
# Load data
#--------------------------------------------------------------
dim(dat <- read.csv(paste0(data_path,"Kwilu_EAS.csv"))) #pop data
# 14,283 EAs with 97 variables
str(dat)
names(dat) # view variable names

# Rename the level names of the GHSL-based settlement classes
class(dat$set_typ <- factor(dat$GHSL_SMOD))
levels(dat$set_typ) <- 1:length(levels(dat$set_typ))
table(dat$set_typ); table(dat$GHSL_SMOD)

# create/rename some key variables
dat <- dat %>%
  mutate(pop = PNLP_median_imputed_total,
         occ = ifelse(CIESIN_bcount>0,1,0),# create binary occupancy variable
         bldg = CIESIN_bcount,
         dens = PNLP_median_imputed_total/CIESIN_bcount) # density variable

table(dat$occ)# 916  EAs with zero buildings

# view the observation points
dat %>% ggplot(aes(x=long, y=lat))+
  geom_point(color="magenta", alpha=0.2)+
  theme_classic()+
  ylab("Latitude") +
  xlab("Longitude") 

# view the distribution of the population count by settlement types
 p1 <-  dat %>% arrange(desc(pop)) %>%
  mutate(settlement.type = set_typ) %>%
  ggplot(aes(x=long, y=lat, size=pop, fill=settlement.type)) +
  geom_point(alpha=0.5, shape=21, color="black") +
  scale_size(name="Observed \n Population \n Counts") +
   scale_fill_viridis(discrete=TRUE, option="B") +
  theme_classic() +
  ylab("Latitude") +
  xlab("Longitude") 
p1

# view the distribution of the building count by settlement types
  p2 <- dat %>% arrange(desc(bldg)) %>%
    mutate(settlement.type = set_typ) %>%
    ggplot(aes(x=long, y=lat, size=bldg, fill=settlement.type)) +
    geom_point(alpha=0.5, shape=21, color="black") +
    scale_size(name="Observed \n Building \n Counts") +
    scale_fill_viridis(discrete=TRUE, option="D") +
    theme_classic() +
    ylab("Latitude") +
    xlab("Longitude") 
p2

# Box plot for the population count
p3 <- dat %>% 
  mutate(settlement.type = set_typ) %>%
  dplyr::select(pop, settlement.type) %>%
  gather("settlement.type", "pop") %>%
  ggplot(aes(x = settlement.type, y = pop, color = settlement.type)) +  # ggplot function
  geom_boxplot() +
  theme_classic() +
  ylab("Observed Population Count") +
  xlab("Settlement type") 
p3

# Box plot for the building count
p4 <- dat %>% 
  mutate(settlement.type = set_typ) %>%
  dplyr::select(bldg, settlement.type) %>%
  gather("settlement.type", "bldg") %>%
  ggplot(aes(x = settlement.type, y = bldg, color = settlement.type)) +  # ggplot function
  geom_boxplot() +
  theme_classic() +
  ylab("Observed Building Count") +
  xlab("Settlement type") 
p4

# check summary stats of total population and building counts
summary(dat$pop)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 1.0   144.0   398.0   578.2   727.0 29435.0    2490

summary(dat$bldg)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00   12.99   58.17   89.50  125.49 4686.98 

length(dat$bldg[dat$bldg==0])  # 916 EAs without buildings

#---------------------------------------------------
# Carry out GLM-based covariate selection 
#---------------------------------------------------
#Scaling function to scale covariates
stdize <- function(x)
{ stdz <- (x - mean(x, na.rm=T))/sd(x, na.rm=T)
return(stdz) 
}

# Covs selection
covs <- dat  %>% 
  dplyr:: select(starts_with("x")) # x48 looks strange

# check covariates classes
apply(covs, 2, class) 
# Calculate mean and standard deviation of covariates
#cov_stats <- data.frame(Covariate = colnames(covs),
#                        Mean = apply(covs, 2, mean, na.rm = TRUE),
#                        Std_Dev = apply(covs, 2, sd, na.rm = TRUE))


#apply scaling function
covs <- apply(covs, 2, stdize) %>%    #z-score
  as_tibble()

# Select response variable for covariates selection
dat$ldens <- log(dat$dens)
dat$ldens[is.infinite(dat$ldens)] = NA
covs$ldens <- dat$ldens
dim(covs_selection <- covs %>% drop_na(ldens))

covs_selection <- data.frame(covs_selection)
# Covariate Selection -----------------------------------------------------
fit.dens  <- glm(ldens~ ., data = covs_selection, family = gaussian)
summary(fit.dens)

# Stepwise Regression both forward and backward
stepwr <- stepAIC(fit.dens, scale = 0,
                     direction = c("both"),
                     trace = 1, keep = NULL, steps = 1000, use.start = FALSE,
                     k = 2)
stepwr$call

# Refit model for the selected covariates
fitt <- glm(stepwr$formula, family = stepwr$family,
            data = stepwr$data)
summary(fitt)

# Calculate the VIF
viff = vif(fitt)
vifl5 <- viff[viff < 5] # VIF less than 5
paste(names(vifl5),
      "+ ", collapse="")

# Refit the selected model for covariates with < 5 vif
fitt2 <- glm(ldens ~ x9 + x10 + x11 + x13 + x14 + x16 + x18 + x20 + x21 + x48 + 
               x52 + x54 + x63 + x64 + x69 + x71 + x76 + x80, family = stepwr$family,
            data = stepwr$data)
summary(fitt2)

# extract the p-values
pvalues <- coef(summary(fitt2))[, "Pr(>|t|)"]

# Select significant covariates only 
 sig_covs <- names(which(pvalues < 0.05)[-1]) # intercept not required
paste(sig_covs,
      "+ ", collapse="")

# Refit the selected model for significant covariates with < 5 vif
fitt3 <- glm(ldens ~ x9 + x13 + x14 + x16 + x18 + x20 + x21 + x48 + 
               x52 + x54 + x63 + x64 + x69 + x71 + x76 + x80, family = stepwr$family,
             data = stepwr$data)
summary(fitt3)

# final checks for vif
viff3 = vif(fitt3)
vifl5b <- viff3[viff3 < 5] # VIF less than 5
paste(names(vifl5b),
      "+ ", collapse="")

# final best covs
# x9 + x13 + x14 + x16 + x18 + x20 + x21 + x48 + x52 + x54 + x63 
# + x64 + x69 + x71 + x76 + x80 +
mod_covs <- names(vifl5b)

#-------------------------------------------------------
# Prepare data for Bayesian Spatial Modelling
#-------------------------------------------------------
# scale covariates for modelling and set infinite values to NA
dat$dens[is.infinite(dat$dens)] = NA
dat[,mod_covs] <- apply(dat[,mod_covs], 2, stdize)


# Build Mesh
plot(dat$long, dat$lat) # view the points ro reconfirm

#  extract the coordinates
coords <- cbind(dat$long, dat$lat)


# Create the hull boundary and the mesh
non_convex_bdry_kwi <- inla.nonconvex.hull(coords, -0.035, -0.05, resolution = c(100, 100))
mesh <- inla.mesh.2d(boundary = non_convex_bdry_kwi, max.edge=c(0.2,1), 
                     offset = c(0.3, 0.32),
                     cutoff = 0.1)

# visualise mesh with points
par(mfrow=c(1,1))
plot(mesh)
points(coords, col="dark red", cex=0.4)


#  Build projector matrix A
A <-inla.spde.make.A(mesh=mesh,loc=as.matrix(coords));dim(A)

#  Create the SPDE object
spde <- inla.spde2.matern(mesh, alpha=2)

# specify the observation indices for estimation 
iset <- inla.spde.make.index(name = "s", spde$n.spde)

# function for extracting random effects
str_ranef <- function(dat, strat, st)
{
  uniq <- unique(strat[strat!=0])
  ranef <- rep(0, length(uniq))
  for(i in  1:nrow(dat))
  {
    
    for(j in 1:length(uniq))
    {
      if(strat[i]==uniq[j]) ranef[i] = st[j]
    }
    
  }
  ranef
}


# add observation level ids
dat$eps <- 1:nrow(dat) 
table(dat$set_typ <- factor(dat$set))

# select variables for stacking
covars <- dat[,c(mod_covs,"set_typ","eps")]; dim(covars)

# steck for occupancy
stk_occ <- inla.stack(data=list(y=dat$occ), #the response
                      
                      A=list(A,1),  #the A matrix; the 1 is included to make the list(covariates)
                      
                      effects=list(c(list(Intercept=1), #the Intercept
                                     iset),  #the spatial index
                                   #the covariates
                                   list(covars)),
                      #this is a quick name so you can call upon easily
                      tag='est')


# stack for density
dat$dens[is.infinite(dat$dens)] = NA
stk_dens <- inla.stack(data=list(y=dat$dens), #the response
                       
                       A=list(A,1),  #the A matrix; the 1 is included to make the list(covariates)
                       
                       effects=list(c(list(Intercept=1), #the Intercept
                                      iset),  #the spatial index
                                    #the covariates
                                    list(covars)),
                       #this is a quick name so you can call upon easily
                       tag='est')


# -----------------------------------------------------------------
# Models for occupancy
#------------------------------------------------------------------
# Model 1
f1occ <-  y ~ -1 + Intercept + x9 + x13 + x14 + x16 + x18 + x20 + x21 + x48 + 
  x52 + x54 + x63 + x64 + x69 + x71 + x76 + x80 +
  f(eps, model="iid") +  f(s, model=spde) + f(set_typ, model="iid") 

mod1occ <- inla(f1occ,
                family="binomial",
                control.compute = list(dic=T, waic=T, cpo=T, config=T),
                data = inla.stack.data(stk_occ),
                control.predictor = list(A = inla.stack.A(stk_occ), compute=T)) 

summary(mod1occ)

# Model 2
f2occ <-  y ~ -1 + Intercept + x9 + x13 + x14 + x16 + x18 + x20 + x21 + x48 + 
  x52 + x54 + x63 + x64 + x69 + x71 + x76 + x80 +
  f(eps, model="iid") +  f(set_typ, model="iid") 

mod2occ <- inla(f2occ,
                family="binomial",
                control.compute = list(dic=T, waic=T, cpo=T, config=T),
                data = inla.stack.data(stk_occ),
                control.predictor = list(A = inla.stack.A(stk_occ), compute=T)) 

summary(mod2occ)


# Model 3
f3occ <-  y ~ -1 + Intercept + x9 + x13 + x14 + x16 + x18 + x20 + x21 + x48 + 
  x52 + x54 + x63 + x64 + x69 + x71 + x76 + x80 +
  f(eps, model="iid") + f(s, model=spde)

mod3occ <- inla(f3occ,
                family="binomial",
                control.compute = list(dic=T, waic=T, cpo=T, config=T),
                data = inla.stack.data(stk_occ),
                control.predictor = list(A = inla.stack.A(stk_occ), compute=T)) 

summary(mod3occ)


# -----------------------------------------------------------------
# Models for Density
#------------------------------------------------------------------
# Model 1
f1dens <-  y ~ -1 + Intercept +  x9 + x13 + x14 + x16 + x18 + x20 + x21 + x48 + 
  x52 + x54 + x63 + x64 + x69 + x71 + x76 + x80 +
  f(eps, model="iid") +  f(s, model=spde) + f(set_typ, model="iid") 

mod1dens <- inla(f1dens,
                 family="gamma",
                 control.compute = list(dic=T, waic=T, cpo=T, config=T),
                 data = inla.stack.data(stk_dens),
                 control.predictor = list(A = inla.stack.A(stk_dens), compute=T)) 

summary(mod1dens)

# Model 2
#mod2b
f2dens <-  y ~ -1 + Intercept +  x9 + x13 + x14 + x16 + x18 + x20 + x21 + x48 + 
  x52 + x54 + x63 + x64 + x69 + x71 + x76 + x80 +
  f(eps, model="iid") + f(set_typ, model="iid") 

mod2dens <- inla(f2dens,
                 family="gamma",
                 control.compute = list(dic=T, waic=T, cpo=T, config=T),
                 data = inla.stack.data(stk_dens),
                 control.predictor = list(A = inla.stack.A(stk_dens), compute=T)) 

summary(mod2dens)


# Model 3
f3dens <-  y ~ -1 + Intercept +  x9 + x13 + x14 + x16 + x18 + x20 + x21 + x48 + 
  x52 + x54 + x63 + x64 + x69 + x71 + x76 + x80 +
  f(eps, model="iid") +  f(s, model=spde)

mod3dens <- inla(f3dens,
                 family="gamma",
                 control.compute = list(dic=T, waic=T, cpo=T, config=T),
                 data = inla.stack.data(stk_dens),
                 control.predictor = list(A = inla.stack.A(stk_dens), compute=T)) 

summary(mod3dens)
#index3dens <-inla.stack.index(stk_dens, "est")$data #---extract the data location indices
#pred3dens <- exp(mod3dens$summary.linear.predictor[index3dens,"mean"]) #--predicted mean count

#----------------------------------
# Model Selection based on DIC
#-----------------------------------
# select the best fit model of occupancy models
t(c(mod1=mod1occ$dic$dic,mod2=mod2occ$dic$dic,mod3=mod3occ$dic$dic)) 
#  mod1     mod2     mod3
# 4304.482 5012.666 4601.975 (Model 1 - best)

# select the best fit model of density models
t(c(mod1=mod1dens$dic$dic,mod2=mod2dens$dic$dic,mod3=mod3dens$dic$dic)) 
#    mod1     mod2     mod3
# 55685.92 57941.14 56206.84 (Model 1 - best)

#----------------------------------------------------------
## extract posterior estimates based on the best fit models
#----------------------------------------------------------

# occupancy
index1occ <-inla.stack.index(stk_occ, "est")$data #---extract the location indices
pred1occ <- plogis(mod1occ$summary.linear.predictor[index1occ,"mean"]) #--predicted mean count
summary(pred1occ) # check the support
dat$occ_prd <- pred1occ

# population density
index1dens <-inla.stack.index(stk_dens, "est")$data #---extract the data location indices
pred1dens <- exp(mod1dens$summary.linear.predictor[index1dens,"mean"]) #--predicted mean count
sum(dat$pop_prd <- pred1dens*dat$bldg) # check the predicted total 
dat$dens_prd <- pred1dens

# quick view of posterior distributions
hist(dat$occ_prd); hist(dat$dens_prd)

# Add adjusted totals by occupancy and modelled building counts
sum(dat$pop_prd2 <- pred1dens*dat$bldg*dat$occ_prd) # occupancy probability adjusted


# Check: Compare health area observed versus predicted total 
(dat_ha <- data.frame(dat %>% group_by(health_area_name) %>%
                            summarise(obs = sum(pop, na.rm=T),
                                      pred = round(sum(pop_prd, na.rm=T)),
                                      pred2 = round(sum(pop_prd2, na.rm=T)),
                                      bldg = round(sum(bldg, na.rm=T)))))

dat_ha <- dat_ha[order(dat_ha$health_area_name),]
apply(dat_ha[,2:5], 2, sum, na.rm=T)

# function for fit metrics
model_metrics <- function(obs,pred)
{
  residual = pred - obs
  MAE = mean(abs(residual), na.rm=T)
  MSE = mean(residual^2, na.rm=T)
  RMSE = sqrt(MSE)
  BIAS = mean(residual, na.rm=T)
  corr = cor(obs[!is.na(obs)],pred[!is.na(obs)])
  
  output <- list(MAE  = MAE,
                 RMSE = RMSE,
                 #BIAS = abs(BIAS),
                 corr = corr)
  return(output)
}
##-- Calculate Model fit metrics


met_dat <- data.frame(dat %>% 
                        dplyr::select(pop,pop_prd, 
                                      pop_prd2) %>%
                        drop_na())

# EA level
met1a <- unlist(model_metrics(met_dat$pop, 
                              met_dat$pop_prd))

met2a <- unlist(model_metrics(met_dat$pop, 
                              met_dat$pop_prd2)) # 

rbind(met1a, met2a)
#         MAE     RMSE      corr
# met1a 103.6604 224.1521 0.9686154
# met2a 105.9473 226.8958 0.9679107


# Health Area level
met1b <- unlist(model_metrics(dat_ha$obs, 
                              dat_ha$pred))

met2b <- unlist(model_metrics(dat_ha$obs, 
                              dat_ha$pred2))
rbind(met1b, met2b)
#      MAE    RMSE      corr
# met1b 989.7451 1525.16 0.9668249
# met2b 995.6305 1521.30 0.9675091


# Make scatter plots
# EA level
par(mfrow=c(2,2), mar=c(3,2,1,2))
plot(dat$pop, dat$pop_prd)
abline(a=0, b=1,col="red", lwd=2)
plot(dat$pop, dat$pop_prd2)
abline(a=0, b=1,col="red", lwd=2)


# Health area
#par(mfrow=c(1,2), mar=c(3,2,1,2))
plot(dat_ha$obs, dat_ha$pred)
abline(a=0, b=1,col="red", lwd=2)
plot(dat_ha$obs, dat_ha$pred2)
abline(a=0, b=1,col="red", lwd=2)


#-----------------------------------------------------------------------------
# prepare grid covariates for posterior simulation and predictions
#-----------------------------------------------------------------------------
dim(pred_covs_kwi <- readRDS(paste0(data_path,"Kwilu_covs_stack.rds")))# Load grid covariates
# 210,761 grid cells
pred_covs_kwi <- data.frame(pred_covs_kwi) # Ensure it's a data frame

# explore the grid data
dim(pred_covs_kwi)
head(pred_covs_kwi)
names(pred_covs_kwi)
covs_prd <- intersect(mod_covs, names(pred_covs_kwi))

# rename the GPS coordinates
pred_covs_kwi$lon <- pred_covs_kwi$x
pred_covs_kwi$lat <- pred_covs_kwi$y

# Check the grid cell values of the covariates for NAs
sapply(pred_covs_kwi[,covs_prd], function(x) length(x[is.na(x)]))#
#   x9   x14   x16   x18   x21   x48   x52   x63   x64   x69   x71   x76   x80 
#   0     0     0 28440 77366     0     0     0     0     0     0     0     0 

#  Build projector matrix Apred for prediction
coords_prd <- cbind(pred_covs_kwi$lon,pred_covs_kwi$lat)
Apred <-inla.spde.make.A(mesh=mesh,loc=as.matrix(coords_prd));dim(Apred)

# -------------------------------
# Scale the prediction covariates
pred_covs_kwi[,covs_prd] <- as.data.frame(apply(pred_covs_kwi[,covs_prd],2, stdize))

# check and drop covariates with unrealistically high values after scaling
apply(pred_covs_kwi[,covs_prd], 2, summary)
#------------------------------------------------
# Refit the best models with the final covariates
#-----------------------------------------------
# NB: Due to the high number of missing values, x18 and x21 were also dropped
# Only x16, x63, x64 and x80 were retained - cutt of point: max=4.6
#     The model was then updated. 
# Occupancy

f1occl <-  y ~ -1 + Intercept + x16 + x63 + x64 +
  f(eps, model="iid") +  f(s, model=spde) + f(set_typ, model="iid") 

mod1occl <- inla(f1occl,
                family="binomial",
                control.compute = list(dic=T, waic=T, cpo=T, config=T),
                data = inla.stack.data(stk_occ),
                control.predictor = list(A = inla.stack.A(stk_occ), compute=T)) 

summary(mod1occl)

# Population density
#f1densl <-  y ~ -1 + Intercept + x16  + x63 + x64 + x80 +
f1densl <-  y ~ -1 + Intercept + x16  + x63 + x64  +
  f(eps, model="iid") +  f(s, model=spde) + f(set_typ, model="iid") 

mod1densl <- inla(f1densl,
                 family="gamma",
                 control.compute = list(dic=T, waic=T, cpo=T, config=T),
                 data = inla.stack.data(stk_dens),
                 control.predictor = list(A = inla.stack.A(stk_dens), compute=T)) 

summary(mod1densl)
index1dens <-inla.stack.index(stk_dens, "est")$data #---extract the data location indices
pred1densl <- exp(mod1densl$summary.linear.predictor[index1dens,"mean"]) #--predicted mean count
sum(dat$pop_prdl <- pred1densl*dat$bldg) # check the predicted total 
dat$dens_prdl <- pred1densl

plot(dat$pop_prdl, dat$pop)


#-------------------------------------------------------------------
#  Run Posterior Simulation and Grid prediction simultaneously
#----------------------------------------------------------------
# Occupancy
simOcc <- function(model, dat, Aprediction, run)
{
  fixedeff  <- occ_hat  <- matrix(0, nrow=nrow(dat), ncol = run)
  inla.seed = as.integer(runif(1)*.Machine$integer.max)
  inla.seed =  106089# sampling seed for reproducible results
  set.seed(inla.seed)
  print(inla.seed)
  
  # draw posterior samples
  m1.samp <- inla.posterior.sample(run, 
                                   model, seed = inla.seed ,
                                   selection=list(
                                     x16=1,
                                     x63=1,
                                     x64=1
                                   ),
                                   num.threads="1:1")
  
  # extract spatial random effects
  sfield_nodes_mean <- model$summary.random$s['mean']
  field_mean <- (Aprediction%*% as.data.frame(sfield_nodes_mean)[, 1])
  
  # Run grid predictions with each simulated set of parameter values
  for(i in 1:run)
  {
    
    fixedeff[,i] <-  
      model$summary.fixed['Intercept', 'mean'] +
      m1.samp[[i]]$latent[1,] * dat[,'x16'] +
      m1.samp[[i]]$latent[2,] * dat[,'x63'] +
      m1.samp[[i]]$latent[3,] * dat[,'x64'] +
      
      
      dat$set_ranef + 
      rnorm(nrow(dat), 0, 1/model$summary.hyperpar$mean[1]) + 
      
      field_mean[,1]
    
    occ_hat[,i]<- exp(fixedeff[,i])/(1+exp(fixedeff[,i]))
  }
  
  # carry out inference
  dat$mean_occ_hat <- apply(occ_hat, 1, mean, na.rm=T) # mean density
  dat$mean_occ_hat <- apply(occ_hat, 1, mean, na.rm=T) # mean pop count
  dat$lower_occ_hat <- apply(occ_hat, 1, quantile, probs=c(0.025), na.rm=T) # lower
  dat$upper_occ_hat <- apply(occ_hat, 1, quantile, probs=c(0.975), na.rm=T) # upper
  dat$sd_occ_hat <- apply(occ_hat, 1, sd, na.rm=T) # standard deviation
  dat$cv_occ_hat <- dat$sd_occ_hat/dat$mean_occ_hat# coefficient of variation
  
  # store data
  output <- list(occ_hat = occ_hat,
                 est_data = dat)
}

run = 100
pred_covs_kwi$set_ranef  <- str_ranef(pred_covs_kwi, 
                                      pred_covs_kwi$DRC_GHSL_SMOD, 
                                      mod1occl$summary.random$set_typ$mean)
system.time(str(sim.occ <-  simOcc(mod1occl,pred_covs_kwi, Apred, run))) 
pred_covs_kwi$occ_hat <- sim.occ$est_data$mean_occ_hat
min(pred_covs_kwi$occ_hat); max(pred_covs_kwi$occ_hat)


# Density
simDens <- function(model, dat, Aprediction, run)
{
  fixedeff  <- dens_hat <- pop_hat <- pop_hat_sc1<- pop_hat_prob <- pop_hat_occ <- matrix(0, nrow=nrow(dat), ncol = run)
  inla.seed = as.integer(runif(1)*.Machine$integer.max)
  inla.seed =  1060897359 # sampling seed for reproducible results
  set.seed(inla.seed)
  print(inla.seed)
  
  # draw posterior samples
  m1.samp <- inla.posterior.sample(run, 
                                   model, seed = inla.seed ,
                                   selection=list(
                                     x16=1,
                                     x63=1,
                                     x64=1
                                   ),
                                   num.threads="1:1")
  
  # extract spatial random effects
  sfield_nodes_mean <- model$summary.random$s['mean']
  field_mean <- (Aprediction%*% as.data.frame(sfield_nodes_mean)[, 1])
  
  # Run grid predictions with each simulated set of parameter values
  for(i in 1:run)
  {
    
    fixedeff[,i] <-  
      model$summary.fixed['Intercept', 'mean'] +
      m1.samp[[i]]$latent[1,] * dat[,'x16'] +
      m1.samp[[i]]$latent[2,] * dat[,'x63'] +
      m1.samp[[i]]$latent[3,] * dat[,'x64'] +

      
      dat$set_ranef + 
      rnorm(nrow(dat), 0, 1/model$summary.hyperpar$mean[2]) + 
      
      field_mean[,1]
    
    dens_hat[,i]<- exp(fixedeff[,i])
    pop_hat[,i] <- dens_hat[,i]*dat$bldp
    
  }
  
# carry out inference
  dat$mean_dens_hat <- apply(dens_hat, 1, mean, na.rm=T) # mean density
  dat$mean_pop_hat <- apply(pop_hat, 1, mean, na.rm=T) # mean pop count
  dat$lower_pop_hat <- apply(pop_hat, 1, quantile, probs=c(0.025), na.rm=T) # lower
  dat$upper_pop_hat <- apply(pop_hat, 1, quantile, probs=c(0.975), na.rm=T) # upper
  dat$sd_pop_hat <- apply(pop_hat, 1, sd, na.rm=T) # standard deviation
  dat$cv_pop_hat <- dat$sd_pop_hat/dat$mean_pop_hat# coefficient of variation
  
# store data
  output <- list(pop_hat = pop_hat,
                 est_data = dat)
}

run = 100

# Scenario 1 - usual model
pred_covs_kwi$bldp <- round(pred_covs_kwi$COD_bld_count_20240617_unpublished)
pred_covs_kwi$set_ranef  <- str_ranef(pred_covs_kwi, 
                                       pred_covs_kwi$DRC_GHSL_SMOD, 
                                      mod1densl$summary.random$set_typ$mean)
system.time(str(sim.dens1 <-  simDens(mod1densl,pred_covs_kwi, Apred, run))) 
pred_covs_kwi$pop_hat1 <- sim.dens1$est_data$mean_pop_hat
sum(pred_covs_kwi$pop_hat1)

hist(pred_covs_kwi$pop_hat1)
max(pred_covs_kwi$pop_hat1)
min(pred_covs_kwi$pop_hat1)
max(pred_covs_kwi$pop_hat1/round(pred_covs_kwi$COD_bld_count_20240617_unpublished))


# Scenario 2 - occupancy probability adjsusted
pred_covs_kwi$bldp <- pred_covs_kwi$occ_hat*pred_covs_kwi$COD_bld_count_20240617_unpublished
pred_covs_kwi$set_ranef  <- str_ranef(pred_covs_kwi, 
                                      pred_covs_kwi$DRC_GHSL_SMOD, 
                                      mod1densl$summary.random$set_typ$mean)
system.time(str(sim.dens2 <-  simDens(mod1densl,pred_covs_kwi, Apred, run))) 
pred_covs_kwi$pop_hat2 <- sim.dens2$est_data$mean_pop_hat
sum(pred_covs_kwi$pop_hat2)

hist(pred_covs_kwi$pop_hat2)
max(pred_covs_kwi$pop_hat2)
min(pred_covs_kwi$pop_hat2)
max(pred_covs_kwi$pop_hat2/round(pred_covs_kwi$COD_bld_count_20240617_unpublished))

##  Health area and province totals
# Extract subsects of the posterior estimates for further processing
dim(data.count1 <- data.frame(cbind(sim.dens2$est_data[,c("lon","lat","health_area",
                                                   "health_area_name")], sim.dens1$pop_hat)))

dim(data.count2 <- data.frame(cbind(sim.dens1$est_data[,c("lon","lat","health_area",
                                                    "health_area_name")], sim.dens2$pop_hat)))


#-----Calculate Provincial Totalswith uncertainties
prov_est <- function(dat, run)
{
  p_hat <- dat[,5:(run+4)]
  tots <- apply(p_hat,2, sum, na.rm=T) #Col sums
  
  tot_sd  <- sd(tots, na.rm=T)
  
  tot_mean  <- mean(tots, na.rm=T)
  
  tot_lower <- quantile(tots, probs=c(0.025))
  tot_median <- quantile(tots, probs=c(0.5))
  tot_upper <- quantile(tots, probs=c(0.975))
  
  return(estimates <- data.frame(estimates=unlist(list(total=tot_mean, 
                                                       lower=tot_lower, median=tot_median, upper=tot_upper))))
}

(prov_count1 <- t(prov_est(data.count1, run))) 
(prov_count2 <- t(prov_est(data.count2, run))) 
prv_est <- data.frame(rbind(prov_count1,prov_count2))
rownames(prv_est) <- c("unadjusted", "adjusted")
prv_est

# Helath area population estimates
harea_est <- function(datr, run)
{
  names <-  unique(datr$health_area_name)
  uniR <-as.numeric(as.factor(unique(datr$health_area_name)))
  outR <- matrix(0, nrow=length(uniR), ncol=5)
  for(j in 1:length(uniR))
  {
    reg <- datr[datr$health_area_name==names[j],]
    rtots <- apply(reg[,5:(4+run)], 2, sum, na.rm=T)
    
    rtot_mean  <- mean(rtots, na.rm=T)
    rtot_sd <- sd(rtots, na.rm=T)
    
    rtot_lower <- quantile(rtots, probs=c(0.025))
    rtot_median <- quantile(rtots, probs=c(0.5))
    rtot_upper <- quantile(rtots, probs=c(0.975))
    rtot_uncert <- (rtot_upper - rtot_lower)/rtot_mean
    
    restimates <- round(c(rtot_mean, rtot_lower, rtot_median,rtot_upper, rtot_uncert),2)
    outR[j,] <- restimates
  }
  outR <- data.frame(outR)
  return(reg_est <- data.frame(id = uniR,
                               names = names,
                               total = outR[,1],
                               lower = outR[,2],
                               median = outR[,3],
                               upper = outR[,4],
                               uncertainty = outR[,5]))
}

dim(harea_count1 <- harea_est(data.count1[!is.na(data.count1$health_area),], run))
dim(harea_count2 <- harea_est(data.count2[!is.na(data.count2$health_area),], run))


####-----------------------------------------------------------------------------

names(dat)
dat_ha <- data.frame(dat %>% group_by(health_area_name) %>%
                       summarise(total1 = sum(pop, na.rm=T),
                                 build = sum(CIESIN_bcount)))

# rename mispelt health area name to the correct name 
dat_ha$health_area_name[1] = "All?gresse"

grd_bldg <- data.frame(pred_covs_kwi %>% group_by(health_area_name) %>%
  summarise(cbldg = sum(COD_bld_count_20240617_unpublished, na.rm=T)))

unique(pred_covs_kwi$health_area_name)

dat_ha$health_area_name

dt_count <- data.frame(health_area_name = harea_count1$names,
                       pred1 = harea_count1$total,
                       pred2 = harea_count2$total
)
sort(dt_count$health_area_name)
dim(jnd_dat <- full_join(dat_ha,dt_count,
                         by="health_area_name"))

jnd_dat <- merge(jnd_dat, grd_bldg, by = "health_area_name")
jnd_dat <- jnd_dat %>% mutate(obs = total1) 

# compare building counts
par(mfrow=c(1,1))
boxplot(jnd_dat$build, jnd_dat$cbldg)

# save the health area comparison data
write.csv(jnd_dat, paste0(results_path, "/health_area_data_validation_new.csv"),
          row.names=F)

par(mfrow=c(2,2))
plot(jnd_dat$obs, jnd_dat$pred1)
abline(0,1)
plot(jnd_dat$obs, jnd_dat$pred2)
abline(0,1)


model_metrics <- function(obs,pred)
{
  residual = pred - obs
  MAE = mean(abs(residual), na.rm=T)
  MSE = mean(residual^2, na.rm=T)
  RMSE = sqrt(MSE)
  #BIAS = mean(residual, na.rm=T)
  corr = cor(obs[!is.na(obs)],pred[!is.na(obs)])
  
  output <- list(MAE  = MAE ,
                 RMSE = RMSE,
                 #BIAS = abs(BIAS),
                 corr = corr)
  return(output)
}

# METRICS
met1 <- unlist(model_metrics(jnd_dat$obs, 
                              jnd_dat$pred1))
met2 <- unlist(model_metrics(jnd_dat$obs, 
                              jnd_dat$pred2))

rbind(met1, met2)
#       MAE     RMSE      corr
#met1 3092.093 4947.384 0.5941576 - base
#met2 2719.723 4543.286 0.6489489 - adjusted 


#-----------------------------------------------------------------
# Save the estimates of the best fit model 
# province Totals
write.csv(prov_count1, paste0(results_path, "/Province_total_counts_base.csv"), # unscaled 
          row.names=F)

write.csv(prov_count2, paste0(results_path, "/Province_total_counts_adjusted.csv"), # unscaled 
          row.names=F)

# Health Area Totals
write.csv(harea_count1, paste0(results_path, "/Health_area_total_counts_base.csv"), # unscaled 
          row.names=F)

write.csv(harea_count2, paste0(results_path, "/adjusted/Health_area_total_counts_adjusted.csv"), # unscaled 
          row.names=F)


#-------------------------------------------------------------------
#  Write and save raster files of the grid predictions
#-------------------------------------------------------------------
ref_coords <- cbind(sim.dens1$est_data$lon, sim.dens1$est_data$lat)
x <- as.matrix(ref_coords)


library(raster)
#-------------------------------------------
# base (not adjusted for occupancy probability)
#---------------------------------------------

# mean
z1b <- as.matrix(sim.dens1$est_data$mean_pop_hat)
#sum(sim.dens1$est_data$mean_pop_hat, na.rm=T)
h1b = rasterFromXYZ(cbind(x, z1b))

# lower
z1bL <- as.matrix(sim.dens1$est_data$lower_pop_hat)
#sum(sim.dens1$est_data$mean_pop_hat, na.rm=T)
h1bL = rasterFromXYZ(cbind(x, z1bL))

# upper
z1bU <- as.matrix(sim.dens1$est_data$upper_pop_hat)
#sum(sim.dens1$est_data$mean_pop_hat, na.rm=T)
h1bU = rasterFromXYZ(cbind(x, z1bU))

# write and save
writeRaster(h1b, filename=paste0(results_path, "/Kwilu_mean_base.tif"),
            overwrite=TRUE, options = c('COMPRESS' = 'LZW'))

writeRaster(h1bL, filename=paste0(results_path, "/Kwilu_lower_base.tif"),
            overwrite=TRUE, options = c('COMPRESS' = 'LZW'))

writeRaster(h1bU, filename=paste0(results_path, "/Kwilu_upper_base.tif"),
            overwrite=TRUE, options = c('COMPRESS' = 'LZW'))


#-----------------------------------------------
#                   adjusted 
#-----------------------------------------------
# mean
z2b <- as.matrix(sim.dens2$est_data$mean_pop_hat)
#sum(sim.dens1$est_data$mean_pop_hat, na.rm=T)
h2b = rasterFromXYZ(cbind(x, z2b))

# lower
z2bL <- as.matrix(sim.dens2$est_data$lower_pop_hat)
#sum(sim.dens1$est_data$mean_pop_hat, na.rm=T)
h2bL = rasterFromXYZ(cbind(x, z2bL))

# upper
z2bU <- as.matrix(sim.dens2$est_data$upper_pop_hat)
#sum(sim.dens1$est_data$mean_pop_hat, na.rm=T)
h2bU = rasterFromXYZ(cbind(x, z2bU))

# write and save
writeRaster(h2b, filename=paste0(results_path, "/Kwilu_mean_adjusted.tif"),
            overwrite=TRUE, options = c('COMPRESS' = 'LZW'))

writeRaster(h2bL, filename=paste0(results_path, "/Kwilu_lower_adjusted.tif"),
            overwrite=TRUE, options = c('COMPRESS' = 'LZW'))

writeRaster(h2bU, filename=paste0(results_path, "/Kwilu_upper_adjusted.tif"),
            overwrite=TRUE, options = c('COMPRESS' = 'LZW'))



##------------------------------------------------------------------
#    Model Cross Validation (K-fold)
#-------------------------------------------------------------------
# Function: Model fit metrics 
mod_metrics2 <- function(obs, pred)
{
  residual = pred - obs
  MAE = mean(abs(residual), na.rm=T)# Mean Absolute Error
  MSE = mean(residual^2, na.rm=T) # Mean Square Error
  RMSE = sqrt(MSE) # Root Mean Square Error
  BIAS = mean(residual, na.rm=T) # Bias
  CORR = cor(obs[!is.na(obs)], pred[!is.na(obs)]) # Correlation Coefficient
  output <- list(MAE  = MAE,
                 RMSE = RMSE,
                # BIAS = abs(BIAS),
                 CC=CORR)
  return(output)
}

#----Extract settement type effects
set_t <- function(dat, st)
{
  uniq <- unique(dat$set_typ)
  uniq[1]
  for(i in  1:nrow(dat))
  {
    
    for(j in 1:3)
    {
      if(dat$set_typ[i]==uniq[j]) dat$set_typ2[i] = st[j]
    }
    
  }
  dat$set_typ2
}


cross_validate <- function(dat, n.folds, mod1, mod2, formula1, formula2, A, seed)
{
  #--------------------------------------------
  # dat: the input survey data containing the all the variables
  # n.folds: number of test (k) folds to use
  # mod: the best model of the full or reference data
  # A: the projection  matrix used in training the full data model
  # seed: a random sample seed to make results reproducible 
  #--------------------------------------------
  seed = 1325
  set.seed(seed) 
  #dat <- dat2 # survey data
  N <- nrow(dat)
  dat$bld <- dat$bldg
  dat$obs <- dat$pop
  #n.folds <- 9
  #  k_fold
  ######
  
  table(ind_train <- factor(sample(x = rep(1:n.folds, each = floor(N / n.folds)),  # Sample IDs for training data
                                   size = N)))
  
  table(as.numeric(ind_train)) 
  dat$k_fold <- as.numeric(ind_train)
  coords = cbind(dat$lon, dat$lat)
  
  
  k_uniq <-sort(unique(dat$k_fold))
  
  
  #mod1 <- mod1densl
  #mod2 <- mod1occl
  #---------------------------------------------------------------
  #                   in-sample
  #---------------------------------------------------------------
  
  met_list_in <- list()
  pred_list_in <- list()
  for(i in 1:length(k_uniq))
  {
    #i =2
    
    print(paste0("in-sample cross-validation using fold ", i, sep=""))
    test_ind <- which(dat$k_fold==k_uniq[i])
    dim(test <- dat[test_ind, ]) #---test set for fold i
    
    
    train_coords <- coords
    test_coords <- coords[test_ind,]
    
    
   
    # spatial random effects based on the full data best model
    sfield_nodes_mean1 <- mod1$summary.random$s['mean']
    field_mean1 <- (A%*% as.data.frame(sfield_nodes_mean1)[, 1])
    
    ##--------
    fixed1 <-  
      mod1$summary.fixed['Intercept', 'mean'] +
      mod1$summary.fixed['x16', 'mean'] * test[,'x16'] +
      mod1$summary.fixed['x63', 'mean'] * test[,'x63'] +
      mod1$summary.fixed['x64', 'mean'] * test[,'x64'] +
      
      #rnorm(nrow(test), 0, 1/mod1$summary.hyperpar$mean[2]) + #---settlement type and region nested effects
      rnorm(nrow(test), 0, 1/mod1$summary.hyperpar$mean[5]) + #---settlement type random effect
      
      mod1$summary.random$eps['mean'][test_ind,1] + #--uncorrelated spatial random effects
      field_mean1[test_ind,1]
      dens_ht1 <- exp(fixed1)
      sum(pop_ht1 <- dens_ht1*test$bld)
    
    
      
      ##--------
      sfield_nodes_mean2 <- mod2$summary.random$s['mean']
      field_mean2 <- (A%*% as.data.frame(sfield_nodes_mean2)[, 1])
      
      fixed2 <-  
        mod2$summary.fixed['Intercept', 'mean'] +
        mod2$summary.fixed['x16', 'mean'] * test[,'x16'] +
        mod2$summary.fixed['x63', 'mean'] * test[,'x63'] +
        mod2$summary.fixed['x64', 'mean'] * test[,'x64'] +
        
        #rnorm(nrow(test), 0, 1/mod2$summary.hyperpar$mean[1]) + #---settlement type and region nested effects
        rnorm(nrow(test), 0, 1/mod2$summary.hyperpar$mean[4]) + #---settlement type random effect
        
        mod2$summary.random$eps['mean'][test_ind,1] + #--uncorrelated spatial random effects
        field_mean2[test_ind,1]
        test$occ_hat <- exp(fixed2)/(1+exp(fixed2))
        sum(test$pop_ht <- pop_ht1*test$occ_hat)
        
    # visualise samples
    #par(mfrow =c(1,2))
    #plot(shp)
    #points(train_coords, col="blue", pch =15, cex=0.6)
    #points(test_coords, col="orange", pch=15, cex=0.6)
    
    
    par(mfrow =c(1,1))
    plot(test$obs, test$pop_ht, xlab = "Observed", 
         ylab = "Predicted", col=c('blue','orange'),
         pch=c(16,16), cex.axis=1.5)
    abline(0,1)
    legend("topleft", c("Observed", "Predicted"), col=c("blue", "orange"), pch=c(16,16),
           bty="n", cex=1.5) 
    
    
    
    # calculate fit metrics
    met_in <- mod_metrics2(test$obs,  
                           test$pop_ht)
    
    met_list_in[[i]]<- unlist(met_in)
    pred_list_in[[i]] <- data.frame(obs = round(test$obs), pred = round(test$pop_ht),
                                    fold = rep(i, length(test$obs)),
                                    data = rep("insample", length(test$obs)))
  }
  met_list_in_dat <- do.call(rbind,met_list_in)
  metrics_in <- apply(met_list_in_dat, 2, mean)
  pred_list_in_dat <- do.call(rbind,pred_list_in)
  
  #-----------------------------------------------------------
  #               out - of -sample
  #-----------------------------------------------------------
  met_list_out <- list()
  pred_list_out <- list()
  for(i in 1:length(k_uniq))
  {
    #i=2
    print(paste0("out-of-sample cross-validation using fold ", i, sep=""))
    train_ind <- which(dat$k_fold!=k_uniq[i])
    test_ind <- which(dat$k_fold==k_uniq[i])
    dim(train <- dat[train_ind, ])#---train set for fold i
    dim(test <- dat[test_ind, ]) #---test set for fold i
    
    
    train_coords <- coords[train_ind,]
    test_coords <- coords[test_ind,]
    
    
    ###---Create projection matrices for training and testing datasets
    Ae<-inla.spde.make.A(mesh=mesh,loc=as.matrix(train_coords));dim(Ae) #training
    
    
    #####------------------------
    # population density
    covars_train <- train[,c("x16", "x63", "x64", "set_typ", "eps")]; dim(covars_train)
    stk_train_dens <- inla.stack(data=list(y=train$dens), #the response
                            
                            A=list(Ae,1),  #the A matrix; the 1 is included to make the list(covariates)
                            
                            effects=list(c(list(Intercept=1), #the Intercept
                                           iset),  #the spatial index
                                         #the covariates
                                         list(covars_train)
                            ), 
                            tag='train')
    
    
    
    # Occupancy
  stk_train_occ <- inla.stack(data=list(y=train$occ), #the response
                            
                            A=list(Ae,1),  #the A matrix; the 1 is included to make the list(covariates)
                            
                            effects=list(c(list(Intercept=1), #the Intercept
                                           iset),  #the spatial index
                                         #the covariates
                                         list(covars_train)
                            ), 
                            tag='train')
  
    ###---Rerun INLA for model test prediction
  # population density
  #mod1 = mod1densl
  #mod2 = mod1occl
  #formula1 = f1densl
    model1 <-inla(formula1, #the formula
                 data=inla.stack.data(stk_train_dens,spde=spde),  #the data stack
                 family= 'gamma',   #which family the data comes from
                 control.predictor=list(A=inla.stack.A(stk_train_dens),compute=TRUE),  #compute gives you the marginals of the linear predictor
                 control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                 verbose = FALSE) #can include verbose=TRUE to see the log of the model runs
    summary(model1)
    
    
  # Occupancy 
    #formula2 = f1occl
    model2 <-inla(formula2, #the formula
                 data=inla.stack.data(stk_train_occ,spde=spde),  #the data stack
                 family= 'binomial',   #which family the data comes from
                 control.predictor=list(A=inla.stack.A(stk_train_occ),compute=TRUE),  #compute gives you the marginals of the linear predictor
                 control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                 verbose = FALSE) #can include verwbose=TRUE to see the log of the model runs
    summary(model2)
    
    # Extract Spatial random effects from the full data best model
    sfield_nodes_mean1 <- mod1$summary.random$s['mean']
    field_mean1 <- (A%*% as.data.frame(sfield_nodes_mean1)[, 1])
    
    
    sfield_nodes_mean2 <- mod2$summary.random$s['mean']
    field_mean2 <- (A%*% as.data.frame(sfield_nodes_mean2)[, 1])
    
    
    ##--------
    fixed1 <-  
      model1$summary.fixed['Intercept', 'mean'] +
      model1$summary.fixed['x16', 'mean'] * test[,'x16'] +
      model1$summary.fixed['x63', 'mean'] * test[,'x63'] +
      model1$summary.fixed['x64', 'mean'] * test[,'x64'] +
      
      #rnorm(nrow(test), 0, 1/model$summary.hyperpar$mean[2]) + #---settlement type and region nested effects
      rnorm(nrow(test), 0, 1/model1$summary.hyperpar$mean[5]) + #---settlement type random effect
 
      mod1$summary.random$eps['mean'][test_ind,1] + #--uncorrelated spatial random effects
      field_mean1[test_ind,1]
      dens_ht1 <- exp(fixed1)
      sum(pop_ht1 <- dens_ht1*test$bld)
    
    
    # visualise samples
    # par(mfrow =c(1,2))
    #plot(shp)
    #points(train_coords, col="blue", pch =15, cex=0.6)
    #points(test_coords, col="orange", pch=15, cex=0.6)
      
      
      fixed2 <-  
        model2$summary.fixed['Intercept', 'mean'] +
        model2$summary.fixed['x16', 'mean'] * test[,'x16'] +
        model2$summary.fixed['x63', 'mean'] * test[,'x63'] +
        model2$summary.fixed['x64', 'mean'] * test[,'x64'] +
        
        #rnorm(nrow(test), 0, 1/model$summary.hyperpar$mean[2]) + #---settlement type and region nested effects
        rnorm(nrow(test), 0, 1/model2$summary.hyperpar$mean[4]) + #---settlement type random effect
        
        mod2$summary.random$eps['mean'][test_ind,1] + #--uncorrelated spatial random effects
        field_mean2[test_ind,1]
        test$occ_hat <- exp(fixed2)/(1+exp(fixed2))
        sum(pop_ht <- pop_ht1*test$occ_hat)
      
      
    
    par(mfrow =c(1,1))
    plot(test$obs, pop_ht, xlab = "Observed", 
         ylab = "Predicted", col=c('blue','orange'),
         pch=c(16,16), cex.axis=1.5)
    abline(0,1)
    legend("topleft", c("Observed", "Predicted"), col=c("blue", "orange"), pch=c(16,16),
           bty="n", cex=1.5) 
    
    
    
    # calculate fit metrics
    met_out <- mod_metrics2(test$obs,  
                            pop_ht)
    
    met_list_out[[i]]<- unlist(met_out)
    pred_list_out[[i]] <- data.frame(obs = test$obs, pred = pop_ht,
                                     fold = rep(i, length(test$obs)),
                                     data = rep("outsample", length(test$obs)))
  }
  met_list_out_dat <- do.call(rbind,met_list_out)
  metrics_out <- apply(met_list_out_dat, 2, mean) # fit metrics
  pred_list_out_dat <- do.call(rbind,pred_list_out)# predictions
  
  cv_mets <- rbind(metrics_in, metrics_out)
  output <- list( met_list_in_dat = met_list_in_dat,
                  met_list_out_dat = met_list_out_dat,
                  pred_dat = rbind(pred_list_in_dat, pred_list_out_dat),
                  cv_metrics = rbind(metrics_in, metrics_out))
}

dat <- dat
#dat$obs <- dat$pop # observed household size

(cross_val <- cross_validate(dat, n.folds = 9, 
                             mod1 = mod1densl, 
                             mod2 = mod1occl,
                             formula1 = f1densl,
                             formula2 = f1occl,
                             A, 
                             seed = 13235))

cross_val$met_list_in_dat  # in-sample metrics per fold 
cross_val$met_list_out_dat  # out-of-sample metrics per fold
cross_val$cv_metrics    # combined averaged metrics
#               MAE     RMSE        CC
#metrics_in  269.6381 519.9225 0.8008034
#metrics_out 263.3326 520.2749 0.8012217

cross_val$pred_dat  # combined prediction data

#-------------------------------------------------------------
# Scatter plots of model cross-validation results
#------------------------------------------------------------
pred.data <- cross_val$pred_dat
pred.data$Fold <- factor(pred.data$fold)
pred.data$data <- factor(pred.data$data,
                         levels = c("insample", "outsample"),
                         labels = c("In-Sample", "Out-of-Sample"))

table(pred.data$data)
names(pred.data)

library(ggpubr)
levels(pred.data$Fold) <- paste("Fold", 1:9)
plot_cval <- ggscatter(pred.data, x = "obs", y = "pred",
                       add = "reg.line",                         # Add regression line
                       facet.by = "data",
                       conf.int = TRUE,                          # Add confidence interval
                       color = "Fold", 
                       palette = "lancet",           # Color by groups "cyl"
                       shape = "Fold"                             # Change point shape by groups "cyl"
) 

rcval <-  ggpar(plot_cval, xlab="Observed counts", ylab="Predicted Counts",
                legend = "top", 
                legend.title = "Fold (k{=9}-fold)",size=20,
                font.legend=c(12),
                font.label = list(size = 12, face = "bold", color ="red"),
                font.x = c(12),
                font.y = c(12),
                font.main=c(14),
                font.xtickslab =c(12),
                font.ytickslab =c(12),
                # orientation = "reverse",
                xtickslab.rt = 45, ytickslab.rt = 45)
rcval


# Individually
plot_cval2 <- ggscatter(pred.data, x = "obs", y = "pred",
                       add = "reg.line",                         # Add regression line
                       facet.by = "Fold",
                       #conf.int = TRUE,                          # Add confidence interval
                       color = "data", 
                       palette = "lancet"#,           # Color by groups "cyl"
                       #shape = "Fold"                             # Change point shape by groups "cyl"
) 
rcval2 <-  ggpar(plot_cval2, xlab="Observed counts", ylab="Predicted Counts",
                 legend = "right", 
                 #legend.title = "Dataset",size=20,
                 legend.title = list(color = "Dataset", 
                                     linetype = "data", shape = "Fold"),
                 font.legend=c(12),
                 font.label = list(size = 12, face = "bold", color ="red"),
                 font.x = c(12),
                 font.y = c(12),
                 font.main=c(14),
                 font.xtickslab =c(12),
                 font.ytickslab =c(12),
                 # orientation = "reverse",
                 xtickslab.rt = 45, ytickslab.rt = 45)
rcval2


### Altogether
plot_cval3 <- ggscatter(pred.data, x = "obs", y = "pred",
          add = "reg.line",                         # Add regression line
          #facet.by = "Fold",
          #conf.int = TRUE,                          # Add confidence interval
          color = "data", 
          palette = "lancet",           # Color by groups "cyl"
          shape = "Fold"                             # Change point shape by groups "cyl"
) 

rcval3 <-  ggpar(plot_cval3, xlab="Observed counts", ylab="Predicted Counts",
                legend = "right", 
                #legend.title = "Dataset",size=20,
                legend.title = list(color = "Dataset", 
                                    linetype = "data", shape = "Fold"),
                font.legend=c(12),
                font.label = list(size = 12, face = "bold", color ="red"),
                font.x = c(12),
                font.y = c(12),
                font.main=c(14),
                font.xtickslab =c(12),
                font.ytickslab =c(12),
                # orientation = "reverse",
                xtickslab.rt = 45, ytickslab.rt = 45)
rcval3


# Visualise the gridded data
library(terra)
# base
kwi_grd.base <- rast(paste0(results_path, "/Kwilu_mean_base.tif"))
plot(kwi_grd.base, main="Gridded population counts \n for Kwilu",
     col=viridis(100))

# adjusted
kwi_grd.adj <- rast(paste0(results_path, "/Kwilu_mean_adjusted.tif"))
plot(kwi_grd.adj, main="Gridded population counts \n for Kwilu",
     col=viridis(100))


library(raster)
# load kwilu shapefile
dim(kwi_shp <- st_read(paste0(data_path,"Kwilu/Kwilu_4km.shp")))
kwi_sp <- as(st_geometry(kwi_shp), "Spatial")
plot(kwi_sp, col="black")
plot(kwi_grd.adj,  legend=T, col = viridis(100), 
     asp = NA, cex.axis=1.4, add=T)


## Zoom in
new.extent <- extent(c(18.82, 18.84, -5.06, -5.04))
plot(kwi_grd.adj, ext = new.extent, legend=T, 
     col = viridis(100),  asp = NA, cex.axis=1.4)


# plot the histogram
kwiHist<-hist(kwi_grd.adj,
              breaks=6,
              main="Histogram of Gridded Population \n Estimates for Kwilu Province",
              col="lightblue",  # changes bin color
              xlab= "Population Count")  # label the x-axis

kwi_grd.adj_df <- data.frame(kwi_grd.adj)
hist(log(kwi_grd.adj_df$Kwilu_mean_adjusted),
     breaks=100, prob = T,
     main="Histogram of Gridded Population \n Estimates for Kwilu Province",
     col="lightblue",  # changes bin color
     xlab= "Natural Log of Population Count")  # label the x-axis
lines(density(log(kwi_grd.adj_df$Kwilu_mean_adjusted)), lwd=3, col="magenta")



##

kwiHist$breaks
## [1]   0  50 100 150 200 250 300

kwiHist$counts

## [1] 9649 1978  401   75   15    5
plot(kwi_sp, col="black")
plot(kwi_grd.adj, 
     breaks = kwiHist$breaks, 
     col = terrain.colors(6),
     #col = viridis(6),
     #palette= viridis(100),
     main="Digital Surface Model (DSM) - HARV",
     add=T)
## Zoom in
new.extent <- extent(c(18.82, 18.84, -5.06, -5.04))
plot(kwi_grd.adj, ext = new.extent, legend=T, 
     col = terrain.colors(6),  asp = NA, cex.axis=1.4)

##
rm(dat_ha,r2p_merg,r2p, r2p_merg2,r2pFort, z1b,
   z1bL,z1bU, z2b, z2bL, z2bU)


## Overlay over basemaps
library(basemaps)
library(ggplot2)
library(tidyterra)
library(raster)
library(ggmap)
library(mapview)
#set_defaults(kwi_grd.adj, map_service = "esri", map_type = "world_imagery")
x <- basemap_raster(kwi_grd.adj, map_service = "esri", map_type = "world_imagery")
x_terr <- rast 

ggplot() +
  geom_spatraster_rgb(data = x_terr) +
  geom_spatraster(data = kwi_grd.adj) +
  geom_sf(data=kwi_shp, colour = "yellow", fill = NA)+
 scale_fill_gradient(
    limits = c(0,170), low = "red", high = "blue",
    na.value = NA, breaks = seq(0,170, by=30)
  ) +
  ggtitle("") +
  theme(
    panel.spacing.x = unit(0, "lines"), axis.title.x = element_blank(),
    axis.text = element_blank(), axis.ticks = element_blank(),
    axis.title.y = element_blank(), plot.title = element_text(
      hjust = 0.5,
      size = 18, face = "bold"
    ), strip.background = element_blank(),
    strip.text = element_text(size = 14, color = "black", face = "bold"),
    legend.text = element_text(size = 14, face = "bold"),
    legend.title = element_text(hjust = 0.1, size = 16, face = "bold"), 
    panel.background = element_blank()
  ) +
  guides(fill = guide_colourbar(
    title = "Population Count"#, title.hjust = 0.5,
   # barheight = 15, barwidth = 2
  ))

out_path <- "//worldpop.files.soton.ac.uk/Worldpop/Projects/WP517763_GRID3_Phase2/Working/DRC/CHRIS_N/new_data/Kwilu/output"


#load("//worldpop.files.soton.ac.uk/Worldpop/Projects/WP517763_GRID3_Phase2/Working/DRC/CHRIS_N/new_data/Kwilu/output/Kwilu.Rdata")
save.image("//worldpop.files.soton.ac.uk/Worldpop/Projects/WP517763_GRID3_Phase2/Working/DRC/CHRIS_N/new_data/Kwilu/output/Kwilu.Rdata")
