#==============================================================================#
#                                                                              #
#                                 SMR analysis                                 #
#                       Feral cats - Otway Ranges - 2017                       #
#   updated script - new secr version / fixed tn error / grids as sessions     #
#                                 Matthew Rees                                 #
#                                   23/5/2020                                  #
#                                                                              #
#==============================================================================#


## BACKGROUND
# this code builds MLE spatial mark-resight models (sighting-only) for feral cats in the Otway Ranges, VIC, AUS - from 1 survey in 2017. 
# we used camera-traps and natural coat markings to identify cats
# this is not the exact code used for https://doi.org/10.1016/j.biocon.2019.108287, but an updated script with the following changes:
# 1) blurry tabby cats correctly placed as tn - mark status unknown (thanks to Joanne Potts for picking this error up)
# 2) camera-trap grids seperated by being different sessions, rather than using a habitat mask with covariate (a more efficient approach)
# 3) "fastproximity = FALSE" arguement in models with time covariates as well as models which are being compared to time covariates (necessary due to change in default settings in the new secr package) 
# NOTE: these changes have no effect on the results - all results presented in Rees et al. 2019 are correct. 



#==============================================================================#
## LOAD PACKAGES / DATA

# secr version used: 4.2.2 / R version 3.6.1
library(secr)

#save.image("workspace.RData")
load("workspace.RData")

# note 1: to reduce temporal autocorrelation: all individuals / mark types detections have been reduced to a binary 0-1 per camera-trap per [24-hour] occasion. 
# note 2: in the previous version of this script tn cats were incorrectly put as tm - this did NOT effect the results because pID was fixed to 1 (excluding these cats from the models - as they should be).
# note 3: mark status unknown detections are usually included in the model as unmarked detections. Given the context of this study--underestimating cat density is the conservative route--we have kept these detections seperate - and left secr to discard these sightings. 
trapfile <- "data/trap.txt"                                     # trap coordinates and usage
captfile <- "data/capt.txt"                                     # identified, marked cats. 
tu <- as.matrix(read.csv("data/unmarked_dethist.csv")[, -1])    # unmarked cats
tn <- as.matrix(read.csv("data/unknown_dethist.csv")[, -1])     # marked status unknown (blurry tabby cats) - these cats do not get used in the density estimate, just reported in the model summary
#==============================================================================#




#==============================================================================#
## MAKE CAPTHIST FILE 
# note 1: we treat camera-traps as proximity detectors instead of count detectors to (a) reduce temporal autocorrelation (see above note 1), (b) because we had variable usage and (c) are modelling time-based covariates. 
# note 2: this data is currently formatted in one combined session - with 72 occasions overall. This is sighting-only mark-resight so all occasions need to be specified as sighting occasions (0). 
CH <- read.capthist(captfile = captfile, trapfile = trapfile, detector = "proximity", markocc = rep(0, 72))

# add un unmarked and unknown sightings 
# note - again, tn detections are often combined with tu detections, however, in our case it is conservative to underestimate cat density (see Rees et al. 2019)
MRCH <- addSightings(CH, unmarked = tu, nonID = NULL, uncertain = tn)

# split camera-trap into two sessions to distinguish them 
# note: in the previous version of this script - I seperated grids with a habitat mask covariate - but this is quicker & simpler
MRCH2 <- split(MRCH, grepl('T',rownames(traps(MRCH))), bytrap = TRUE)
# rename sessions
session(MRCH2) <- c('north',"south")
# drop unused occasions 
MRCH2[[1]] <- subset(MRCH2[[1]], occasions=13:72)
usage(traps(MRCH2[[1]])) <- usage(traps(MRCH2[[1]]))[,1:60]
MRCH2[[2]] <- subset(MRCH2[[2]], occasions=1:69)
usage(traps(MRCH2[[2]])) <- usage(traps(MRCH2[[2]]))[,1:69]

# check it's all g. 
verify(MRCH2[[1]])
verify(MRCH2[[2]])
summary(MRCH2, moves = TRUE)
par(mfrow = c(1,2))
plot(MRCH2, tracks = TRUE)
#==============================================================================#




#==============================================================================#
## MAKE MASK 
mask <- make.mask(traps(MRCH2), buffer = 3500, type = "trapbuffer", spacing = 346*0.6)
#==============================================================================#




#==============================================================================#
## FIT NULL MARK-RESIGHT MODELS WITH DIFFERENT DETECTOR FUNCTIONS
fit_MR_HN <- secr.fit(MRCH2, mask = mask, trace = TRUE, detectfn = 0, details = list(knownmarks = FALSE), fixed = list(pID = 1), ncores = 3)
fit_MR_HR <- secr.fit(MRCH2, mask = mask, trace = TRUE, detectfn = 1, details = list(knownmarks = FALSE), fixed = list(pID = 1), ncores = 3, start = fit_MR_HN)  # fails without specifying starting values
fit_MR_EX <- secr.fit(MRCH2, mask = mask, trace = TRUE, detectfn = 2, details = list(knownmarks = FALSE), fixed = list(pID = 1), ncores = 3)

# choose best detector function use for proceeding models based on lowest AICc model (if HR and EX have lower AICc, but by less than 2 of HN, use HN)
df_fits <- secrlist(fit_MR_HN, fit_MR_HR, fit_MR_EX)
AIC(df_fits)
## --> use HR detector function 
#==============================================================================#




#==============================================================================#
# ADJUST MODEL WITH BEST DETECTOR FUNCTION FOR OVERDISPERSION IN THE UNMARKED SIGHTINGS 
fit_MR_HR_adj <- secr.fit(MRCH2, mask = mask, trace = TRUE, detectfn = 1, details = list(knownmarks = FALSE, nsim = 10000, fastproximity = FALSE), fixed = list(pID = 1), ncores = 3, start = fit_MR_HR)
fit_MR_HR_adj$details$chat[1:2]
#==============================================================================#




#==============================================================================#
## MODEL EXPLANATORY COVARIATES

# linear time trend (T) 
fit_MR_HR_T <- secr.fit(MRCH2, mask = mask, trace = TRUE, detectfn = 1, fixed = list(pID = 1), details = list(knownmarks = F, chat = fit_MR_HR_adj$details$chat, fastproximity = FALSE), start = fit_MR_HR_adj, model = list(g0 ~ T))

# difference in density between grids (grid)
fit_MR_HR_grid  <- secr.fit(MRCH2, mask = mask, trace = FALSE, detectfn = 1, fixed = list(pID = 1), details = list(knownmarks = F, chat = fit_MR_HR_adj$details$chat, fastproximity = FALSE), start = fit_MR_HR_adj, model = list(D ~ session))

# difference in density between grids (grid)
fit_MR_HR_grid_T  <- secr.fit(MRCH2, mask = mask, trace = FALSE, detectfn = 1, fixed = list(pID = 1), details = list(knownmarks = F, chat = fit_MR_HR_adj$details$chat, fastproximity = FALSE), start = fit_MR_HR_adj, model = list(g0 ~ T, D ~ session))

# compare models using AICc
model_fits <- secrlist(fit_MR_HR_adj, fit_MR_HR_T, fit_MR_HR_grid, fit_MR_HR_grid_T)
AIC(model_fits)
# --> fit_MR_HR_T wins
#==============================================================================#



#==============================================================================#
# PLOT LINEAR TREND
all_predicted <- predict(model_fits$fit_MR_HR_T, newdata = data.frame(T = 1:72), realnames = "g0")
predicted_values <- unlist(sapply(all_predicted, function(x) x[2]))
lower_bound <- unlist(sapply(all_predicted, function(x) x[4]))
upper_bound <- unlist(sapply(all_predicted, function(x) x[5]))
time_sequence <- 1:72

par(mar=c(5.1,6,4.1,2.1), mfrow = c(1,1))
plot(predicted_values ~ time_sequence, type = "l", lwd = 3, col = "black",  ylim = c(0, 0.15), las = 1,  xlab = NA, ylab = NA, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.1)
lines(lower_bound ~ time_sequence, lty = 2, col = "gray50", lwd = 3)  # lty = 2 means dashed lines, lty = 1 is solid, lty = 3 is dotted
lines(upper_bound ~ time_sequence, lty = 2, col = "gray50", lwd = 3)
mtext(expression("Survey duration (days)"), side = 1, las = 1, line = 3, cex = 1.2)
mtext(expression("g0"), side = 2, las = 1, line = 4, cex = 1.2)
#==============================================================================#



# END 