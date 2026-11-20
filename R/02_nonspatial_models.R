# Transcribed from the appendix of my MCom thesis (Rhodes University, 2014). Original code written by me in R;
# not re-run since 2014. See README.md for what changed in transcription.
# Non-spatial binomial GLMs: derivation / validation split, variable themes, bootstrap stepwise selection by stage

#
# Model Building : Non-spatial analysis #
#

# create a df with only the variables of interest
spatialVars = s.p[ which ( colnames (s.p) %in% c ( c ( "Pos" , "Examined" , "DstTClW" , "NDVI" , "altitude" , "SDTemp" , "SDNDVI" , "SDRain" , "q" , "rCIndex" , "summerTemp" , "winterTemp" , "summerNDVI" , "winterNDVI" , "summerRain" , "winterRain" , "totRain" ), rep (paste0("bio" ,1:19))))]

# randomly partition data keeping 85% for derivation data set
derivIndex <- createDataPartition ( spatialVars$ Pos , p = .85 , list =F, times = 1)

# create validation and derivation
spatialVarsDeriv <- spatialVars [ derivIndex ,]
spatialVarsTest <- spatialVars[- derivIndex , ]

# standardise variables
spatialVars [ ,3:36] <- scale (spatialVars[ ,3:36])

# 34 explanatory variables at start

# Univariate Logistic Regression - Model Building

AICList <-lapply ( c ( c ( "DstTClW" , "NDVI" , "altitude" , "SDTemp" , "SDNDVI" , "SDRain" , "q" , "rCIndex" , "summerTemp" , "winterTemp" , "summerNDVI" , "winterNDVI" , "summerRain" , "winterRain" , "totRain" ), rep (paste0("bio" ,1:19))),
function ( var ){
formula <- as.formula ( paste ( "Pos/Examined ~" , var ))
nonSpatialUniVar <- glm ( formula , data =spatialVarsDeriv, weights =Examined, family = binomial )
cbind ( summary (nonSpatialUniVar)$ aic , exp ( summary (nonSpatialUniVar)$ coef [, "Estimate" ]) , summary (nonSpatialUniVar)$ coef [, "Pr(>|z |)" ])
})
# http : / / rstudio-pubs-static.s3 . amazonaws . com/ 2989_ceae90d128554c728d5388439adf0661 . html access : 28 Feb


# put list into matrix
m <- matrix ( unlist (AICList) , ncol=6, byrow=TRUE)

# delete 1 st , 3 rd , 5 th column ( intercept details not required )
m <- m[ , - c (1 , 3 , 5) ]

# make col and row names nameable
dimnames (m) <- list ( rownames (m, do.NULL = FALSE, prefix = "row" ), colnames (m, do.NULL =FALSE, prefix ="col" ))

# name columns
colnames (m) <- c ( "AIC" , "OR" , "Pr(>|z |)" )

# name rows
rownames (m) <- c ( c ( "DstTClW" , "NDVI" , "altitude" , "SDTemp" , "SDNDVI" , "SDRain" , "q" , "rCIndex" , "summerTemp" , "winterTemp" , "summerNDVI" , "winterNDVI" , "summerRain" , "winterRain" , "totRain" ), rep (paste0("bio" ,1:19)))

# rank matrix from lowest AIC to highest
AICMatrix <- m[ order (m[ , 1 ] ) , ]

# preserve the order of univariate AIC rankings in columns of df
colTestCols <- paste0 (rownames (AICMatrix) [1:nrow (AICMatrix) ])
spatialVarsDerivColTest <- as.data.frame (spatialVarsDeriv[colTestCols])

# get column names for temp and rain related vars
tempThemeCols <- paste0 (c ( "bio1" , "bio2" , "bio3" , "bio4" , "bio5" , "bio6" , "bio7" , "bio8" , "bio9" , "bio10" , "bio11" , "summerTemp" , "winterTemp" , "SDTemp" ))
rainThemeCols <- paste0 (c ( "bio12" , "bio13" , "bio14" , "bio15" , "bio16" , "bio17" , "bio18" , "bio19" , "totRain" , "summerRain" , "winterRain" , "SDRain" , "q" , "rCIndex" ) )
NDVIThemCols <- paste0 (c ( "NDVI" , "summerNDVI" , "winterNDVI" , "SDNDVI" ))   # closing bracket missing in the printed appendix

# to preserve the order of columns by AIC rank remove all non-theme related vars from original df leaving only AIC ordered themed vars
tempTheme <- spatialVarsDerivColTest[- which ( colnames (spatialVarsDerivColTest) %in%c ( c (rainThemeCols) , c ( "NDVI" , "summerNDVI" , "winterNDVI" , "SDNDVI" , "DstTClW" , "altitude" )))]
rainTheme <- spatialVarsDerivColTest[- which ( colnames (spatialVarsDerivColTest) %in% c ( c (tempThemeCols) , c ( "NDVI" , "summerNDVI" , "winterNDVI" , "SDNDVI" , "DstTClW" , "altitude" )))]
NDVITheme <- spatialVarsDerivColTest[- which ( colnames (spatialVarsDerivColTest) %in% c ( c (tempThemeCols) , c (rainThemeCols) , c ( "DstTClW" , "altitude" )))]

### STAGE 2 ###

# function checks for multicolinearity in each theme
# set exact = F for appropriate asmyptotic methods to handle presence of ties
corrFunctionTheme <- function (varX , varDf){

# create matrix to store multicolinearity test results per variable
tst <- matrix ( data =NA, nrow = ncol (varDf) , ncol =4)

# name columns
dimnames (tst) <- list ( rownames (tst , do.NULL = FALSE, prefix = "row" ), colnames (tst , do.NULL = FALSE, prefix = "col" ))
colnames (tst) <- c ( "Upper p" , "Lower p" , "Two.sided p" , "rho" )

for (i in 1:ncol (varDf)){
tst [ i ,1] = cor.test(varX, varDf[,i], method = "spearm" ,alternative="g" ,exact =F)$p.value
tst [ i ,2] = cor.test(varX, varDf[,i], method = "spearm" ,alternative="l" ,exact =F)$p.value
tst [ i ,3] = cor.test(varX, varDf[,i], method = "spearm" ,alternative="two.sided" ,exact=F)$p.value
tst [ i ,4] = cor.test(varX, varDf[,i], method = "spearm" ,alternative="two.sided" ,exact=F)$ estimate
rownames (tst)[i] = names (varDf)[ i ]
}
return (tst)
}

# test multicollinearity among all variables and rank from lowest to highest according to Spearman ' s r and check if r > 0.85
# criteria for excluding a variable : keep variable with lowest AIC from univariate analysis in the presence of collinearity
# starting with lowest ranked AIC variable in each theme
# so that everything correlated to it will have a higher AIC and can be removed
# at each round the variable tested , with lowest AIC , is put into a list

# initialize lists
tempThemeKept <- list ()
tempThemeRemoved <- list ()
rainThemeKept <- list ()
rainThemeRemoved <- list ()
NDVIThemeKept <- list ()
NDVIThemeRemoved <- list ()

# runs until one var remains in themed df

# Temperature theme

while ( ncol (tempTheme) > 1){

# perform correlation test
corrTemp <- corrFunctionTheme (varX = tempTheme [ , 1 ] , varDf = tempTheme )
corrTempRanked <- corrTemp [ order ( abs (corrTemp[ ,4]) ) ,]
corrTempVars <- corrTempRanked [ which ( abs (corrTempRanked [ ,4]) > 0.85) , arr.ind =TRUE]

# whether or not correlated tested var always kept - lowest AIC
tempThemeKept [ length (tempThemeKept)+1] <- list ( colnames (tempTheme) [1])

# add lowest ranking AIC var to list

# else if statement needed becasue row names disappear when only 1 row remains in matrix

if ( length ( nrow (corrTempVars)) < 1) {
tempThemeRemoved [ length (tempThemeRemoved)+1] <- list (NA)
tempTheme <- tempTheme[- which ( colnames ( tempTheme ) %i n% colnames (tempTheme) [1]) ]
} else{
if ( nrow (corrTempVars) == 2) {
ind = which ( rownames (corrTempVars) != colnames (tempTheme) [1])
tempThemeRemoved [ length (tempThemeRemoved)+1] <- rownames (corrTempVars) [ ind ]
tempTheme <- tempTheme[- which ( colnames (tempTheme) %in% rownames ( corrTempVars ) ) ]
} else{
tempThemeRemoved [ length (tempThemeRemoved)+1] <- list ( rownames ( corrTempVars[- which ( rownames (corrTempVars) %in% colnames (tempTheme) [1]) ,]))
tempTheme <- tempTheme[- which ( colnames (tempTheme) %in% rownames ( corrTempVars ) ) ]
}
}

}


# Rain theme

while ( ncol (rainTheme) > 1){

# perform correlation test
corrRain <- corrFunctionTheme (varX = rainTheme [ , 1 ] , varDf = rainTheme )
corrRainRanked <- corrRain [ order ( abs (corrRain[ ,4])) ,]
corrRainVars <- corrRainRanked [ which ( abs (corrRainRanked [ ,4]) > 0.85) , arr.ind =TRUE]

# whether or not correlated tested var always kept - lowest AIC
rainThemeKept [ length (rainThemeKept)+1] <- list ( colnames (rainTheme) [1])

# add lowest ranking AIC var to list

# else if statement needed because row names disappear when only 1 row remains in matrix

if ( length ( nrow (corrRainVars)) < 1) {
rainThemeRemoved [ length (rainThemeRemoved)+1] <- list (NA)
rainTheme <- rainTheme[- which ( colnames ( rainTheme ) %in% colnames (rainTheme) [1]) ]
} else{
if ( nrow (corrRainVars) == 2) {
ind = which ( rownames (corrRainVars) != colnames (rainTheme) [1])
rainThemeRemoved [ length (rainThemeRemoved)+1] <- rownames (corrRainVars)[ ind ]
rainTheme <- rainTheme[- which ( colnames (rainTheme) %in% rownames ( corrRainVars ) ) ]
} else{
rainThemeRemoved [ length (rainThemeRemoved)+1] <- list ( rownames ( corrRainVars[- which ( rownames (corrRainVars) %in% colnames (rainTheme) [1]) ,]))
rainTheme <- rainTheme[- which ( colnames (rainTheme) %in% rownames ( corrRainVars ) ) ]
}
}

}


# NDVI theme

while ( ncol (NDVITheme) > 1){

# perform correlation test
corrNDVI <- corrFunctionTheme (varX = NDVITheme[ , 1 ] , varDf = NDVITheme )
corrNDVIRanked <- corrNDVI [ order ( abs (corrNDVI[ ,4]) ) ,]
corrNDVIVars <- corrNDVIRanked [ which ( abs (corrNDVIRanked [ ,4]) > 0.85) , arr.ind =TRUE]

# whether or not correlated tested var always kept - lowest AIC
NDVIThemeKept [ length (NDVIThemeKept)+1] <- list ( colnames (NDVITheme) [ 1 ] )

# add lowest ranking AIC var to list

# else if statement needed because row names disappear when only 1 row remains in matrix

if ( length ( nrow (corrNDVIVars)) < 1) {
NDVIThemeRemoved [ length (NDVIThemeRemoved)+1] <- list (NA)
NDVITheme <- NDVITheme[- which ( colnames ( NDVITheme ) %i n% colnames (NDVITheme) [ 1 ] ) ]
} else{
if ( nrow (corrNDVIVars) == 2) {
ind = which ( rownames (corrNDVIVars) != colnames (NDVITheme) [ 1 ] )
NDVIThemeRemoved [ length (NDVIThemeRemoved)+1] <- rownames (corrNDVIVars) [ ind ]
NDVITheme <- NDVITheme[- which ( colnames (NDVITheme) %in% rownames ( corrNDVIVars ) ) ]
} else{
NDVIThemeRemoved [ length (NDVIThemeRemoved)+1] <- list ( rownames ( corrNDVIVars[- which ( rownames (corrNDVIVars) %in% colnames (NDVITheme) [1]) ,]))
NDVITheme <- NDVITheme[- which ( colnames (NDVITheme) %in% rownames ( corrNDVIVars ) ) ]
}
}

}

# left with tempThemeKept , rainThemeKept , NDVIThemeKept , DstTClW , altitude


### STAGE 3 ###

spatialVarsDerivDf = spatialVarsDeriv [ which ( colnames (spatialVarsDeriv) %in% c ( c ( tempThemeKept , rainThemeKept , NDVIThemeKept) , c ( "DstTClW" , "altitude" , "Pos" , "Examined" )))]

fit.Stage3 <- glm (Pos/ Examined ~ bio9+SDTemp+bio5+summerTemp+bio3+bio7+rCIndex+bio13+ bio19+bio14+totRain+bio18+summerNDVI+DstTClW+altitude , weights =Examined, data = spatialVarsDerivDf , family="binomial" )

# this will yield a basic candidate model ( Candidate List : )
bootGLM.Stage3 <- boot.stepAIC( fit.Stage3 , spatialVarsDerivDf , direction = "backward" ,alpha=0.05,B=1000)

fitTest.1 <- glm (Pos/ Examined ~ bio9 , weights =Examined, data = spatialVarsDerivDf , family="binomial" )
summary (fitTest.1)

fitTest.2 <- glm (Pos/ Examined ~ bio9 + altitude , weights =Examined, data =spatialVarsDerivDf, family="binomial" )
summary (fitTest.2)


fitTest.3 <- glm (Pos/ Examined ~ bio9 + altitude + bio5 , weights = Examined , data =spatialVarsDerivDf, family="binomial" )
summary (fitTest.3)

fitTest.4 <- glm (Pos/ Examined ~ bio9+altitude+bio5+bio7 , weights =Examined, data =spatialVarsDerivDf, family="binomial" )
summary (fitTest.4)

fitTest.5 <- glm (Pos/ Examined ~ bio9+altitude+bio5+bio7+summerTemp, weights = Examined , data =spatialVarsDerivDf, family="binomial" )
summary (fitTest.5) # exclude summerTemp


fitTest.6 <- glm (Pos/ Examined ~ bio9+altitude+bio5+bio7+summerNDVI, weights = Examined , data =spatialVarsDerivDf, family="binomial" )
summary (fitTest.6) # exclude summerNDVI


fitTest.7 <- glm (Pos/ Examined ~ bio9+altitude+bio5+bio7+SDTemp , weights = Examined , data =spatialVarsDerivDf, family="binomial" )
summary (fitTest.7) # exclude SDTemp


fitTest.8 <- glm (Pos/ Examined ~ bio9+altitude+bio5+bio7+DstTClW, weights = Examined , data =spatialVarsDerivDf, family="binomial" )
summary (fitTest.8)

fitTest.9 <- glm (Pos/ Examined ~ bio9+altitude+bio5+bio7+DstTClW+bio18 , weights = Examined , data =spatialVarsDerivDf, family="binomial" )
summary (fitTest.10)

fitTest.10 <- glm (Pos/ Examined ~ bio9+altitude+bio5+bio7+DstTClW+bio18+rCIndex , weights =Examined, data =spatialVarsDerivDf, family="binomial" )
summary (fitTest.11) # exclude rCIndex

fitTest.11 <- glm (Pos/ Examined ~ bio9+altitude+bio5+bio7+DstTClW+bio18+bio3 , weights =Examined, data =spatialVarsDerivDf, family="binomial" )
summary (fitTest.12) # exclude bio3

fitTest.12 <- glm (Pos/ Examined ~ bio9+altitude+bio5+bio7+DstTClW+bio18+totRain , weights =Examined, data =spatialVarsDerivDf, family="binomial" )
summary (fitTest.13)

fitTest.13 <- glm (Pos/ Examined ~ bio9+altitude+bio5+bio7+DstTClW+bio18+totRain+ bio13 , weights =Examined, data =spatialVarsDerivDf, family="binomial" )
summary (fitTest.14) # exclude bio13


fitTest.14 <- glm (Pos/ Examined ~ bio9+altitude+bio5+bio7+DstTClW+bio18+totRain+ bio19 , weights =Examined, data =spatialVarsDerivDf, family="binomial" )
summary (fitTest.15) # exclude bio19


fitTest.15 <- glm (Pos/ Examined ~ bio9+altitude+bio5+bio7+DstTClW+bio18+totRain+ bio14 , weights =Examined, data =spatialVarsDerivDf, family="binomial" )
summary (fitTest.16) # exclude bio14

# Stage 4: bio9+altitude+bio5+bio7+DstTClW+bio18+totRain


### STAGE 5 ###

# further tests based on enviro themes of previously excluded variables
# criteria: bring variables that were exlcuded ( high correlation + AIC ranking ) into candidate model ( consider bootGLM + manual stepwise ) and re-assess based on
# frequency of selection in bootsrapped samples ( and if variable was not excluded by stepAIC () algorithm )

# kept bio9 ahead of bio6+bio1+winterTemp+bio11

fit.1 <- glm (Pos/ Examined ~ bio9+altitude+bio5+bio7+DstTClW+bio18+totRain+bio6+bio1+ winterTemp+bio11 , weights =Examined, data =spatialVarsDeriv, family = "binomial" )

bootGLM.1 <- boot.stepAIC( fit.1 , spatialVarsDeriv , direction = "backward" ,alpha= 0.05 , B = 1000) # keep winterTemp instead of bio9

# kept bio5 ahead of bio10
fit.2 <- glm (Pos/ Examined ~winterTemp+altitude+bio5+bio7+DstTClW+bio18+totRain+bio10 , weights =Examined, data =spatialVarsDeriv, family = "binomial" )

bootGLM.2 <- boot.stepAIC( fit.2 , spatialVarsDeriv , direction = "backward" ,alpha= 0.05 , B = 1000) # keep bio5

# kept bio7 ahead of bio2
fit.3 <- glm (Pos/ Examined ~winterTemp+altitude+bio5+bio7+DstTClW+bio18+totRain+bio2 , weights =Examined, data =spatialVarsDeriv, family = "binomial" )

bootGLM.3 <- boot.stepAIC( fit.3 , spatialVarsDeriv , direction = "backward" ,alpha= 0.05 , B = 1000) # keep bio7

# kept totRain ahead of bio12
fit.4 <- glm (Pos/ Examined ~winterTemp+altitude+bio5+bio7+DstTClW+bio18+totRain+bio12 , weights =Examined, data =spatialVarsDeriv, family = "binomial" )

bootGLM.4 <- boot.stepAIC( fit.4 , spatialVarsDeriv , direction = "backward" ,alpha= 0.05 , B = 1000) # keep totRain

fit.Stage5 <- glm (Pos/ Examined ~winterTemp+altitude+bio5+bio7+DstTClW+bio18+totRain , weights =Examined, data =spatialVarsDeriv, family = "binomial" )

bootGLM.Stage5 <- boot.stepAIC( fit.Stage5 , spatialVarsDeriv , direction = "backward" , alpha = 0.05 , B = 1000)

# End Stage 5 model : winterTemp+altitude+bio5+bio7+DstTClW+bio18+totRain

## end non-spatial model building
