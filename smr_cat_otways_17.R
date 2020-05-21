#==============================================================================#
#                                                                              #
#                                 SMR analysis                                 #
#                       Feral cats - Otway Ranges - 2017                       #
#            Updated script - new secr version and fixed tn error              #
#                                 Matthew Rees                                 #
#                                   0/5/2020                                   #
#                                                                              #
#==============================================================================#



## LOAD R PACKAGE (secr version used: 4.2.2 / R version 3.6.1)
library(secr)

## LOAD / SAVE DATA
#save.image("workspace.RData")
#load("workspace.RData")


# LOAD CAPTURE HISTORIES AND TRAPFILE 
trapfile <- "data/trap.txt"                                     # trap coordinates and usage
captfile <- "data/capt.txt"                                     # identified, marked cats
tu <- as.matrix(read.csv("data/unmarked_dethist.csv")[, -1])    # unmarked cats
tn <- as.matrix(read.csv("data/unknown_dethist.csv")[, -1])     # marked status unknown (blurry tabby cats)


## MAKE CAPTHIST FILE 
CH <- read.capthist(captfile = captfile, trapfile = trapfile, detector = "proximity", markocc = rep(0, 72))
# add un unmarked and unknown sightings 
MRCH <- addSightings(CH, unmarked = tu, nonID = NULL, uncertain = tn)
# remove repeat detections of an individual within 1 occasion - in order to reduce temporal autocorrelation
MRCH <- reduce(MRCH, dropunused = FALSE)
# split camera-trap into two sessions to distinguish them (I previously did this with a habitat mask - but this is simpler)
MRCH2 <- split(MRCH, grepl('T',rownames(traps(MRCH))), bytrap = TRUE)
# rename sessions
session(MRCH2) <- c('north',"south")

# drop unused occasions 
MRCH2[[1]] <- subset(MRCH2[[1]], occasions=13:72)
usage(traps(MRCH2[[1]])) <- usage(traps(MRCH2[[1]]))[,1:60]
MRCH2[[2]] <- subset(MRCH2[[2]], occasions=1:69)
usage(traps(MRCH2[[2]])) <- usage(traps(MRCH2[[2]]))[,1:69]
summary(MRCH2[[1]])
summary(MRCH2[[2]])

# check it looks right:
summary(MRCH2, moves = TRUE)
par(mfrow = c(1,2))
plot(MRCH2, tracks = TRUE)


## MAKE MASK 
mask <- make.mask(traps(MRCH2), buffer = 3500, type = "trapbuffer", spacing = 346*0.6)


## FIT NULL MARK-RESIGHT MODELS WITH DIFFERENT DETECTOR FUNCTIONS
fit_MR_HN <- secr.fit(MRCH2, mask = mask, trace = TRUE, detectfn = 0, details = list(knownmarks = FALSE), fixed = list(pID = 1), ncores = 3)
fit_MR_HR <- secr.fit(MRCH2, mask = mask, trace = TRUE, detectfn = 1, details = list(knownmarks = FALSE), fixed = list(pID = 1), ncores = 3, start = fit_MR_HN)  # fails without specifying starting values
fit_MR_EX <- secr.fit(MRCH2, mask = mask, trace = TRUE, detectfn = 2, details = list(knownmarks = FALSE), fixed = list(pID = 1), ncores = 3)


# CHOOSE BEST DETECTOR FUNCTION TO PROCEED WITH BASED ON LOWEST AICc VALUE:
df_fits <- secrlist(fit_MR_HN, fit_MR_HR, fit_MR_EX)
AIC(df_fits)
## --> use HR


# ADJUST MODEL WITH BEST DETECTOR FUNCTION FOR OVERDISPERSION IN THE UNMARKED SIGHTINGS 
fit_MR_HR_adj <- secr.fit(MRCH2, mask = mask, trace = TRUE, detectfn = 1, details = list(knownmarks = FALSE, nsim = 10000, fastproximity = FALSE), fixed = list(pID = 1), ncores = 3, start = fit_MR_HR)
fit_MR_HR_adj$details$chat[1:2]


## MODEL EXPLANATORY COVARIATES
# linear time trend (T) 
fit_MR_HR_T <- secr.fit(MRCH2, mask = mask, trace = TRUE, detectfn = 1, fixed = list(pID = 1), details = list(knownmarks = F, chat = fit_MR_HR_adj$details$chat, fastproximity = FALSE), start = fit_MR_HR_adj, model = list(g0 ~ T))
# difference in density between grids (grid)
fit_MR_HR_grid  <- secr.fit(MRCH2, mask = mask, trace = FALSE, detectfn = 1, fixed = list(pID = 1), details = list(knownmarks = F, chat = fit_MR_HR_adj$details$chat, fastproximity = FALSE), start = fit_MR_HR_adj, model = list(D ~ session))
#predict(model_fits$fit_MR_HR_grid, newdata = data.frame(name = 0:1), realnames = "D")  # 0 = southern grid, 1 = northern grid
# difference in density between grids (grid)
fit_MR_HR_grid_T  <- secr.fit(MRCH2, mask = mask, trace = FALSE, detectfn = 1, fixed = list(pID = 1), details = list(knownmarks = F, chat = fit_MR_HR_adj$details$chat, fastproximity = FALSE), start = fit_MR_HR_adj, model = list(g0 ~ T, D ~ session))

# compare models
model_fits <- secrlist(fit_MR_HR_adj, fit_MR_HR_T, fit_MR_HR_grid, fit_MR_HR_grid_T)
AIC(model_fits)





# T PLOT 

all_predicted <- predict(model_fits$fit_MR_HR_T, newdata = data.frame(T = 1:72), realnames = "g0")
predicted_values <- unlist(sapply(all_predicted, function(x) x[2]))
lower_bound <- unlist(sapply(all_predicted, function(x) x[4]))
upper_bound <- unlist(sapply(all_predicted, function(x) x[5]))
time_sequence <- 1:72

par(mar=c(5.1,6,4.1,2.1))

plot(predicted_values ~ time_sequence, type = "l", lwd = 3, col = "black",  ylim = c(0, 0.15), las = 1,  xlab = NA, ylab = NA, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.1)
lines(lower_bound ~ time_sequence, lty = 2, col = "gray50", lwd = 3)  # lty = 2 means dashed lines, lty = 1 is solid, lty = 3 is dotted
lines(upper_bound ~ time_sequence, lty = 2, col = "gray50", lwd = 3)

mtext(expression("Survey duration (days)"), side = 1, las = 1, line = 3, cex = 1.2)
mtext(expression("g0"), side = 2, las = 1, line = 4, cex = 1.2)


# END 


