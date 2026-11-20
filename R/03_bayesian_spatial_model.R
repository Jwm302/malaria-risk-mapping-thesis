# Transcribed from the appendix of my MCom thesis (Rhodes University, 2014). Original code written by me in R;
# not re-run since 2014. See README.md for what changed in transcription.
# Bayesian spatial GLM (spBayes spGLM, exponential covariance), MCMC diagnostics, prediction on a 20 km grid

#
# Spatial Analysis
#

# using full dataset keep only Stage 5 variables
dfFull = s.p[ which ( colnames (s.p) %in% c ( c ( "Pos" , "Examined" , "Lon" , "Lat" , "winterTemp" , "altitude" , "bio5" , "bio7" , "DstTClW" , "bio18" , "totRain" )))]

# using derivation dataset
dfFull = s.p[ derivIndex , ] [ which ( colnames (s.p[derivIndex ,]) %in% c ( c ( "Pos" , "Examined" , "Lon" , "Lat" , "winterTemp" , "altitude" , "bio5" , "bio7" , "DstTClW" , "bio18" , "totRain" )))]

# the following code demonstrates the spatial model for the full dataset the same workings apply when using only the derivation data

# transform DF to a SPDF
coordinates ( dfFull ) <- cbind ( "Lon" , "Lat" )

# standardise data
dfFull@data [ ,3:9] <- scale (dfFull@data[ ,3:9])

# run a non-spatial GLM to obtain starting values for the MH step
dfFull = s.p[ derivIndex , ] [ which ( colnames (s.p[derivIndex ,]) %in% c ( c ( "Pos" , "Examined" , "Lon" , "Lat" , "winterTemp" , "altitude" , "bio5" , "bio7" , "DstTClW" , "bio18" , "totRain" )))]

#starting values from non-spatial analysis for beta
fit <- glm (Pos/ Examined ~ winterTemp+altitude+bio5+bio7+DstTClW+bio18+totRain , weights =Examined, data =dfFull , family=binomial ( "logit" ))
beta.starting <- coefficients (fit)

# use the variance covariance matrix as the proposal ( tuning ) distribution for the MH step
beta.tuning <- t ( chol (vcov( fit )))

# get maximum Euclidean distance between sites
d.max <- max (iDist(dfFull@coords))

# this defines the number of simulations to be run in batches each of certain length as well the burn in period
n.batch <- 3500
batch.length <- 100
n.samples <- n.batch* batch.length
burn.in <- 0.8*n.samples

# 3 spatial GLMMs are run as follows
# notes : effective range of the spatial weights is controlled by phi and is roughly 3 / phi
# posterior inference is based on three MCMC chains each of length 350 000


m1.Full <- spGLM ( Pos ~ winterTemp+altitude+bio5+bio7+DstTClW+bio18+totRain ,
data =dfFull,
weights=dfFull@data$Examined , family = "binomial" , coords = dfFull@coords ,
starting = list ( "beta" = beta.starting , "phi" =3/ (0.5*d.max ), "sigma.sq" =1, "w" =0),
tuning = list ( "beta" = beta.tuning , "phi" =0.06, "sigma.sq" = 0.5 , "w" =0.5),
priors = list (phi.Unif = c (3/d.max ,3/ (0.1*d.max )) , "sigma.sq.IG" = c (2 ,1) ) ,
amcmc = list (n.batch = n.batch , batch.length =batch.length ,accept.rate = 0.43),
cov.model = "exponential" ,verbose=T,n.report =1)

m2.Full <- spGLM ( Pos ~ winterTemp+altitude+bio5+bio7+DstTClW+bio18+totRain ,
data =dfFull,
weights=dfFull@data$Examined , family = "binomial" , coords = dfFull@coords ,
starting = list ( "beta" = beta.starting , "phi" =3/ (0.5*d.max ), "sigma.sq" =1, "w" =0),
tuning = list ( "beta" = beta.tuning , "phi" =0.06, "sigma.sq" = 0.5 , "w" =0.5),
priors = list (phi.Unif = c (3/d.max ,3/ (0.1*d.max )) , "sigma.sq.IG" = c (2 ,1) ) ,
amcmc = list (n.batch = n.batch , batch.length =batch.length ,accept.rate = 0.43),
cov.model = "exponential" ,verbose=T,n.report =1)

m3.Full <- spGLM ( Pos ~ winterTemp+altitude+bio5+bio7+DstTClW+bio18+totRain ,
data =dfFull,
weights=dfFull@data$Examined , family = "binomial" , coords = dfFull@coords ,
starting = list ( "beta" = beta.starting , "phi" =3/ (0.5*d.max ), "sigma.sq" =1, "w" =0),
tuning = list ( "beta" = beta.tuning , "phi" =0.06, "sigma.sq" = 0.5 , "w" =0.5),
priors = list (phi.Unif = c (3/d.max ,3/ (0.1*d.max )) , "sigma.sq.IG" = c (2 ,1) ) ,
amcmc = list (n.batch = n.batch , batch.length =batch.length ,accept.rate = 0.43),
cov.model = "exponential" ,verbose=T,n.report =1)

# consolodate posteriors to perform convergence diagnostics on the fit of the spGLM MCMC chains
posteriors <- as.mcmc.list ( list (m1.Full $p.beta.theta.samples,m2.Full$p.beta.theta.samples , m3.Full$p.beta.theta.samples))

# compute Gelman diagnostics to assess convergence.These compare within-chain to
# between-chain variability.Values near 1 suggest full convergence .
print (gelman.diag (posteriors))


# define samples taken after burn in
sub.samps <- burn.in : n.samples

# credibility intervals for each simulated regression and variance parameter estimated
quantile (m1.Full $p.beta.theta.samples[sub.samps,1] , prob=c (0.025 , 0.975) )
quantile (m1.Full $p.beta.theta.samples[sub.samps,2] , prob=c (0.025 , 0.975) )
quantile (m1.Full $p.beta.theta.samples[sub.samps,3] , prob=c (0.025 , 0.975) )
quantile (m1.Full $p.beta.theta.samples[sub.samps,4] , prob=c (0.025 , 0.975) )
quantile (m1.Full $p.beta.theta.samples[sub.samps,5] , prob=c (0.025 , 0.975) )
quantile (m1.Full $p.beta.theta.samples[sub.samps,6] , prob=c (0.025 , 0.975) )
quantile (m1.Full $p.beta.theta.samples[sub.samps,7] , prob=c (0.025 , 0.975) )
quantile (m1.Full $p.beta.theta.samples[sub.samps,8] , prob=c (0.025 , 0.975) )
quantile (m1.Full $p.beta.theta.samples[sub.samps,9] , prob=c (0.025 , 0.975) )
quantile (m1.Full $p.beta.theta.samples[sub.samps,10] , prob=c (0.025 , 0.975) )

# odds ratio calculated for each simulated paramater ' s mean and credibility intervals
# calculated on odds ratio scale via exponentiation of paramater ' s coefficent

# Intercept odds ratio
intercept <- c ( exp ( mean (m1.Full $p.beta.theta.samples[sub.samps,1]) ) , exp ( quantile (m1.Full$p.beta.theta.samples[sub.samps,1] , prob=c (0.025 , 0.975) )))

# winterTemp odds ratio
winterTemp <- c ( exp ( mean (m1.Full $p.beta.theta.samples[sub.samps,2]) ) , exp ( quantile (m1.Full$p.beta.theta.samples[sub.samps,2] , prob=c (0.025 , 0.975) )))

# altitude odds ratio
altitude <- c ( exp ( mean (m1.Full $p.beta.theta.samples[sub.samps,3]) ) , exp ( quantile (m1.Full$p.beta.theta.samples[sub.samps,3] , prob=c (0.025 , 0.975) )))

# bio5 odds ratio
bio5 <- c ( exp ( mean (m1.Full $p.beta.theta.samples[sub.samps,4]) ) , exp ( quantile (m1.Full $ p.beta.theta.samples[sub.samps,4] , prob=c (0.025 , 0.975) )))

# bio7 odds ratio
bio7 <- c ( exp ( mean (m1.Full $p.beta.theta.samples[sub.samps,5]) ) , exp ( quantile (m1.Full $ p.beta.theta.samples[sub.samps,5] , prob=c (0.025 , 0.975) )))

# DstTClW odds ratio
DstTClW <- c ( exp ( mean (m1.Full $p.beta.theta.samples[sub.samps,6]) ) , exp ( quantile (m1.Full$p.beta.theta.samples[sub.samps,6] , prob=c (0.025 , 0.975) )))

# bio18 odds ratio
bio18 <- c ( exp ( mean (m1.Full $p.beta.theta.samples[sub.samps,7]) ) , exp ( quantile (m1.Full $p.beta.theta.samples[sub.samps,7] , prob=c (0.025 , 0.975) )))

# totRain odds ratio
totRain <- c ( exp ( mean (m1.Full $p.beta.theta.samples[sub.samps,8]) ) , exp ( quantile (m1.Full$p.beta.theta.samples[sub.samps,8] , prob=c (0.025 , 0.975) )))

# sigma . sq odds ratio
sigma.sq <- c ( exp ( mean (m1.Full $p.beta.theta.samples[sub.samps,9]) ) , exp ( quantile (m1.Full$p.beta.theta.samples[sub.samps,9] , prob=c (0.025 , 0.975) )))

# phi odds ratio
phi <- c ( exp ( mean (m1.Full $p.beta.theta.samples[sub.samps,10])) , exp ( quantile (m1.Full $ p.beta.theta.samples[sub.samps,10] , prob=c (0.025 , 0.975) )))


# spatial prediction #

# prepare grid for prediction of gridded sites across Botswana
# 20 km resolution

# read in prediction grid
grid20 <- read.table ( "Volumes/JUSTJUBBA/Spatial/gridPredBots20km.txt" )

# keep only relevant explanatory variables
pred.grid20 <- grid20 [ which ( colnames (grid20) %in% c ( "Lon" , "Lat" , "winterTemp" , "altitude" , "bio5" , "bio7" , "DstTClW" , "bio18" , "totRain" ))]

# convert DF to gridded SPDF and set correct projection
coordinates (pred.grid20 ) <- cbind ( "Lon" , "Lat" )
gridded (pred.grid20 ) <- TRUE
proj4string (pred.grid20) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 + towgs84=0,0,0"

pred.coords20 <- pred.grid20@coords
pred.covars20 <- scale (pred.grid20@data)

# run prediction command spPredict in order to get a prediction of malaria risk at each cell for all MCMC samples using one spGLM chain
m1.pred.Full.20 <- spPredict (m1.Full , pred.coords =pred.coords20, pred.covars=as.matrix ( cbind (1 , pred.covars20 )) , start=burn.in )

# mean and standard deviation of prediction probability at each cell
y.pred.grid.prob.mu <- apply (m1.pred.Full.20 $p.y.predictive.samples ,1 ,mean )
y.pred.grid.prob.sd <- apply (m1.pred.Full.20 $p.y.predictive.samples ,1 , sd )

# plot predicted mean probability of malaria risk at 20km resolution
# requires mba . surf to yield interpolation surfaces at 100x100 resolution on the x and y axis
res <- 100
surf <- mba.surf ( cbind (pred.grid20@coords , y.pred.grid.prob.mu), res , res , extend= TRUE, sp=TRUE) $ xyz.est
plot (botswana)
image.plot ( as.image.SpatialGridDataFrame(surf), asp=1.25, add =T)
plot (dfFull , add =T, pch = 20, cex = 0.3, col = ' black ' )
plot (botswana , add=T)
title (main="Predicted mean probability of malaria risk- 20km grid" )

# plot predicted standard deviation of malaria risk at 20km resolution
res <- 100
surf <- mba.surf ( cbind (pred.grid20@coords , y.pred.grid.prob.sd ), res, res, extend= TRUE, sp=TRUE) $ xyz.est
plot (botswana)
image.plot ( as.image.SpatialGridDataFrame(surf), asp=1.25, add =T)
plot (dfFull , add =T, pch = 20, cex = 0.3, col = ' black ' )
plot (botswana , add=T)
title (main="Predicted standard deviation of malaria risk- 20km grid" )
