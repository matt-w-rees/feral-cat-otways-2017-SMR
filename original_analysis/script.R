#==============================================================================#
#                                                                              #
#                                 SMR analysis                                 #
#                       Feral cats - Otway Ranges - 2017                       #
#                                 Matthew Rees                                 #
#                               1/7/2019 8:44:53                               #
#                                                                              #
#==============================================================================#



# LOAD R PACKAGE (secr version used: 3.2.1)
#install.packages("secr")
library(secr)


# LOAD CAPTURE HISTORIES AND TRAPFILE 
trapfile <- "sdh_identifiabletrap.txt"                          # trap coordinates and usage
captfile <- "sdh_identifiablecapt.txt"                          # identified, marked cats
tm <- as.matrix(read.csv("unidentified_dethist.csv")[, -1])     # unidentified, marked cats 
tu <- as.matrix(read.csv("unmarked_dethist.csv")[, -1])         # unmarked cats

## SET MODEL PARAMETERS
noccasions = 72                       
buffer = 3500                       
spacing <- 346*0.6      # 346 = sigma estimate 
# to run test script, you can increase this spacing for models to run faster (e.g. spacing <- 600)



## MAKE CAPTHIST FILE 
MRCH <- read.capthist(captfile = captfile, trapfile = trapfile, detector = "proximity", markocc = rep(0, noccasions))
MRCH <- addSightings(MRCH, tu, tm)
MRCH <- reduce(MRCH, dropunused = FALSE)        
# check it looks right:
summary(MRCH, moves = TRUE)
plot(MRCH, tracks = TRUE)

## MAKE HABITAT MASK (WITH THE DIFFERENT GRIDS AS A SPATIAL COVARIATE)
grid <- rgdal::readOGR (dsn = "./mask", layer = "treatment_mask") 
baitmask <- make.mask(traps(MRCH), buffer = buffer, type = "trapbuffer", spacing = spacing, poly = grid, keep.poly = FALSE)
baitmask <- addCovariates(baitmask, grid[1:2,])
# plot mask:
plot(baitmask, covariate = 'name', dots = FALSE, border = 100, title = 'treatment')
plotMaskEdge(baitmask, add = TRUE)
plot(traps(MRCH), add = TRUE)


## FIT NULL MARK-RESIGHT MODELS WITH DIFFERENT DETECTOR FUNCTIONS
fit_MR_HN <- secr.fit(MRCH, mask = baitmask, trace = TRUE, detectfn = 0, details = list(knownmarks = F), fixed = list(pID = 1))
fit_MR_HR <- secr.fit(MRCH, mask = baitmask, trace = TRUE, detectfn = 1, details = list(knownmarks = F), fixed = list(pID = 1), start = fit_MR_HN)  # fails without specifying starting values
fit_MR_EX <- secr.fit(MRCH, mask = baitmask, trace = TRUE, detectfn = 2, details = list(knownmarks = F), fixed = list(pID = 1))
# CHOOSE BEST DETECTOR FUNCTION TO PROCEED WITH BASED ON LOWEST AICc VALUE:
df_fits <- secrlist(fit_MR_HN, fit_MR_HR, fit_MR_EX)
save("df_fits.RData", df_fits)  
load("df_fits.RData")   
AIC(df_fits)


# ADJUST MODEL WITH BEST DETECTOR FUNCTION FOR OVERDISPERSION IN THE UNMARKED SIGHTINGS 
fit_MR_HR_adj <- secr.fit(MRCH, mask = baitmask, trace = FALSE, detectfn = 1, details = list(knownmarks = F, nsim = 10000), fixed = list(pID = 1), start = fit_MR_HR)
fit_MR_HR_adj$details$chat[1:2]


## MODEL EXPLANATORY COVARIATES

# linear time trend (T) 
fit_MR_HR_T <- secr.fit(MRCH, mask = baitmask, trace = FALSE, detectfn = 1, fixed = list(pID = 1), details = list(knownmarks = F, chat = fit_MR_HR_adj$details$chat), start = fit_MR_HR_adj, model = list(g0 ~ T))

# difference in density between grids (grid)
fit_MR_HR_grid  <- secr.fit(MRCH, mask = baitmask, trace = FALSE, detectfn = 1, fixed = list(pID = 1), details = list(knownmarks = F, chat = fit_MR_HR_adj$details$chat), start = fit_MR_HR_adj, model = list(D ~ name))
#predict(model_fits$fit_MR_HR_grid, newdata = data.frame(name = 0:1), realnames = "D")  # 0 = southern grid, 1 = northern grid


# DIFFRENT DENSITY BETWEEN GRIDS + LINEAR TIME TREND
fit_MR_HR_grid_T  <- secr.fit(MRCH, mask = baitmask, trace = FALSE, detectfn = 1, fixed = list(pID = 1), details = list(knownmarks = F, chat = fit_MR_HR_adj$details$chat), start = fit_MR_HR_adj, model = list(g0 ~ T, D ~ name))



# SAVE AND COMPARE MODELS USING AIC
model_fits <- secrlist(fit_MR_HR_adj, fit_MR_HR_T, fit_MR_HR_grid, fit_MR_HR_grid_T)
save("model_fits.RData", model_fits)  
load("model_fits.RData")   # note - use model_fits$fit_MR_HR_T etc. after models are loaded
AIC(model_fits)


PLOTS

#det fn plot
plot(model_fits$fit_MR_HR_T, xval = 0:700,  ylim = c(0,0.15), limits = T)


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


