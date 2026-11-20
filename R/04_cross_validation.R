# Transcribed from the appendix of my MCom thesis (Rhodes University, 2014). Original code written by me in R;
# not re-run since 2014. See README.md for what changed in transcription.
# Validation: mean error and mean absolute error, non-spatial vs spatial model

# Cross-validation calculations #

# Predicted vs observed prevelance at validation sites
# Notes : For this section of code the derivation dataset is used , i . e . dfDeriv.The same spatial code as above applies but inseatd of using the full dataset
# dfDeriv was used yielding 3 spatial chains given by : m1.Deriv , m2.Deriv , m3.Deriv .
# The code showing the spatial modelling for the derivation data is ommitted .

# non-spatial prediction accuracy

# derivation and validation subsets of the data
dfDeriv = s.p[ derivIndex , ] [ which ( colnames (s.p[derivIndex ,]) %in% c ( c ( "Pos" , "Examined" , "Lon" , "Lat" , "winterTemp" , "altitude" , "bio5" , "bio7" , "DstTClW" , "bio18" , "totRain" )))]
dfValid = s.p[- derivIndex , ] [ which ( colnames (s.p[derivIndex ,]) %in% c ( c ( "Pos" , "Examined" , "Lon" , "Lat" , "winterTemp" , "altitude" , "bio5" , "bio7" , "DstTClW" , "bio18" , "totRain" )))]

# non-spatial glm fit
fit.Deriv <- glm (Pos/ Examined ~ winterTemp+altitude+bio5+bio7+DstTClW+bio18+totRain , weights =Examined, data =dfDeriv, family = "binomial" )

# test prediction ability of non-spatial model by applying fitted model ( deriv data ) to validation sample ( predicted probabilities )
predValid <- predict (fit.Deriv, dfValid, type="response" , se =T)

# calculate mean absolute error (MAE ) and mean error (ME ) between observed and predicted using non-spatial model at validation sites on probability scale
obs <- dfValid$ Pos/ dfValid$ Examined
pred <- exp (predValid$ fit )/(1+exp (predValid$ fit ))
ME.nonSpatial <- mean (obs-pred )
MAE.nonSpatial <- mean ( abs (obs-pred ) )

# spatial prediction accuracy

# nx2 matrix of n prediction location coordinates
pred.coords <- dfValid@coords

# An nxq design matrix or data frame containing the covariates associated with pred . coords
pred.covars <- dfValid@data [ 3 : 9 ]

# holds the posterior predictive samples given model output m1.Deriv from spGLM function
m1.pred.valid <- spPredict (m1.Deriv , pred.coords =pred.coords ,pred.covars=as.matrix ( cbind (1 , pred.covars )) , start=burn.in )

# mean predicted probability of malaria at each site in validation subset
y.pred.valid.prob.mu <- apply (m1.pred.valid $p.y.predictive.samples ,1 ,mean )

# calculate ME and MAE between observed and predicted using spatial model at validation sites
obs <- dfValid$ Pos/ dfValid$ Examined
ME.spatial <- mean (obs-y.pred.valid.prob.mu)
MAE.spatial <- mean ( abs (obs-y.pred.valid.prob.mu))
