# Transcribed from the appendix of my MCom thesis (Rhodes University, 2014). Original code written by me in R;
# not re-run since 2014. See README.md for what changed in transcription.
# Builds the spatial database: survey sites, distance to surface water, WorldClim climate layers and MODIS NDVI

#
# Compliling of Spatial Database #
#

# set workspace
setwd ( "/Volumes/JUSTJUBBA/Spatial" )

# check if package is installed , install if not load otherwise
packageInstallLoad <- function (x){
for (i in x){
# require returns TRUE invisibly if it was able to load package
if ( ! require (i ,character.only = TRUE ) ){
# If package was not able to be loaded then re- install
install.packages (i ,dependencies =TRUE )
# Load package after installing
require (i ,character.only = TRUE )
}
}
}

# packages needed for spatial db comilation and non-spatial model building
packageInstallLoad (c ( "raster" , "RMySQL" , "sp" , "rgdal" , "gdalUtils" , "maptools" , "maps" , "plyr" , "stats" , "glm2" , "bootStepAIC" , "texreg" , "xtable" , "tables" , "data.table" , "diagram" , "caret" , "bootStepAIC" ))

# packages needed for spatial analysis
packageInstallLoad (c ( "spBayes" , "MBA" , "fields" , "raster" , "coda" , "fields" ))

# connect to MySQL
con <- dbConnect (MySQL() ,
user = Sys.getenv("MYSQL_USER") ,   # credentials removed for publication
password = Sys.getenv("MYSQL_PASSWORD") ,
dbname = "MySql" ,
host = "localhost" )

# delete table if exists
dbSendQuery ( con , 'DROP TABLE IF EXISTS botsTable ; ' ) #decimal (9 ,6)

# create table
# MARA Botswana dataset
dbSendQuery ( con , 'CREATE TABLE botsTable ( Lat float (10 ,8) ,
Lon float (10 ,8) ,
Start_Mnth int ,
Start_Yr int ,
End_Mnth int ,
End_Yr int ,
AgeGroup_Lower int ,
AgeGroup_Upper int ,
Numb_Positive int ,
Numb_Examined int ) ; ' )

# import data from csv file into newly created table - save csv where root user has default access
dbSendQuery ( con , 'LOAD DATA LOCAL INFILE "/usr/local/BotsMARA.csv"
INTO TABLE botsTable
FIELDS TERMINATED BY ","
LINES TERMINATED BY "\n"
IGNORE 1 LINES ( Lat ,
Lon ,
Start_Mnth ,
Start_Yr ,
End_Mnth ,
End_Yr ,
AgeGroup_Lower ,
AgeGroup_Upper ,
Numb_Positive ,
Numb_Examined ) ; ' )

# select all childern between 1 and 15
rs <- dbSendQuery ( con , 'SELECT *
FROM botsTable
WHERE AgeGroup_Upper <= 15 ' )


# retrieve data from MySQL
sitesMultiple <- fetch (rs , n = 6000)

# clear previous MySQL transaction from memory
dbClearResult ( rs )

# close connection to MySQL
dbDisconnect (con)


# clean up data #

# average over multiple sites to get unique sets
sitesDf <- ddply ( sitesMultiple , . ( Lon , Lat ) , summarise , Month = round ( mean (Start_Mnth )) , Pos = round ( mean (Numb_Positive )) , Examined = round ( mean (Numb_Examined ) ) )


sitesDf <- ddply ( sitesMultiple , . ( Lon , Lat ) , summarise , Month = round ( mean (Start_Mnth )) , Pos = round ( mean (Numb_Positive )) , Examined = round ( mean (Numb_Examined ) ) )

# extract coords
coordsSitesDf <- sitesDf [ ,c ( "Lon" , "Lat" )]

# remove any duplicates
if ( any ( duplicated (coordsSitesDf)))
sitesDf <- sitesDf [ ! duplicated (coordsSitesDf) ,]

# remove obs where 0 ppl examined
sitesDf <- subset (sitesDf , sitesDf$ Examined > 0)

# extract coords from clean data
coordsSitesDf <- as.matrix (sitesDf[,c ( "Lon" , "Lat" )])

# Perform calculations to find distance to closest surface water body and find where these are #

# read shapefile into SpatialPointsDataFrame
surfWaterSpdf <- readShapePoints ("Surface_water/gns_swb/gns_swb" )

# convert SPDF to DF
surfWaterDf <- as.data.frame (surfWaterSpdf)

# 129 sites
locs <- SpatialPoints ( sitesDf [ ,1:2] , proj4string=CRS("+proj=longlat +datum=WGS84" ))

# 46591 water sources
src <- SpatialPoints (surfWaterDf [ ,1:2] , proj4string=CRS("+proj=longlat +datum=WGS84" ))

# distances of 46576 water sources per site
distances <- lapply (1: length (locs), function (i) spDistsN1(src , locs[i], longlat=TRUE) )

# get min distance per site in km
sitesDf$DstTClW <- sapply (distances , min )

# get index of min distance sites to find their coords
minPos <- sapply (distances , which.min )

sitesDf$LonWater <- surfWaterDf [ minPos , 1 ]
sitesDf$ LatWater <- surfWaterDf [ minPos , 2 ]

# get map boundary of Botswana explore the sample site data #

# get boundary shape for Botswana
data (wrld_simpl )
botswana <- wrld_simpl [ wrld_simpl$NAME == "Botswana" ,]

# extract coords
coordsSitesDf <- sitesDf [ ,c ( "Lon" , "Lat" )]

# convert DF to SPDF
sitesSpdf <- SpatialPointsDataFrame ( coordsSitesDf , sitesDf )

# assign projection to SPDF
projection ( sitesSpdf ) <- projection (botswana)

# only keep sample sites that are in Botswana
sitesSpdf <- sitesSpdf [botswana ,] # 122 obs

# plot sample sites representing proportion of malaria cases out of no . examined
# over the maximum of this ratio by the size of the circle
plot (botswana)
plot (sitesSpdf , add =T, asp = 1, cex = 4 * sitesSpdf$ Pos/ sitesSpdf$ Examined/max ( sitesSpdf$ Pos/ sitesSpdf$ Examined ) , pch = 1)

# Notes : There seems to be some association btw big circles and other big ones
# that are close together.Informal test showing that circles close together exhibit similiar malaria intensitie although .
# Patterns of attribute seem not to be random .


# WorldClim climate layers & MODIS NDVI raster processing #

# get and process modis data . hdf files
# download all available monthly images for the years 2000 - 2013. MODIS PRODUCT : MOD13A3 , Terra , Vegetation Indices , Tile , 1000m, Monthly
ModisDownload (x="MOD13A3" ,h=c (21 ,22) ,v=c (8 ,9) , dates=c ( ' 2000.01.01 ' , ' 2013.12.31 ' ), mosaic=F, proj=F )

load ( "ModisLP.RData" )
source ( "ModisDownload.R" )

# set wd to where . hdf files live
setwd ( "/Volumes/JUSTJUBBA/Spatial/NDVI_MODIS_Botswana" )

# put each . hdf file in list
out.files <- list.files(getwd () , pattern="hdf$" ,full.names =F)

# get list of subdatasets from . hdf files , choose subdataset 1 - mean monthly NDVI
sdsList <- sapply (X= out.files , FUN = function (out.files ){get_subdatasets (out.files ) [1]})

# see which raster data has more rows/ columns between MODIS & WorldClim- must resample to grid with smallest no . rows/ cols
gdalwarp ( srcfile= sdsList [ 1 ] , t_srs="+proj=longlat +datum=WGS84 +no_defs" ,dstfile= "/Volumes/JUSTJUBBA/Spatial/NDVIBotsTest.tiff" ,te=c (bbox(botswana) [1] ,bbox( botswana ) [ 2 ] , bbox( botswana ) [ 3 ] , bbox( botswana ) [ 4 ] ) )

NDVIRast <- raster ("NDVIBotsTest.tiff" )
worldClimRast <- crop ( raster ("WorldClim_botswana/tmean1_36.tif" ), extent(botswana))

# compare resolution
ncell (NDVIRast) ; ncell (worldClimRast)

# Notes : NDVI has less cells so re-align WorldClim rasters to that of NDVI rasters
# Use NDVI as the model raster to get correct dimensions for WorldClim and NDVI processing

# get sub dataset for each month
# then remove the null entries to get length of list

janList <- llply (sdsList , function (x){ if ( substr (x, 32, 34) %in% c ( "001" , "002" )){ return (x)}})
janList<-janList [ ! sapply (janList , is.null )]

febList <- llply (sdsList , function (x){ if ( substr (x, 32, 34) %in% c ( "032" , "033" )){ return (x)}})
febList<-febList [ ! sapply (febList , is.null )]

marList <- llply (sdsList , function (x){ if ( substr (x, 32, 34) %in% c ( "061" , "062" )){ return (x)}})
marList<-marList [ ! sapply (marList , is.null )]

aprList <- llply (sdsList , function (x){ if ( substr (x, 32, 34) %in% c ( "091" , "092" )){ return (x)}})
aprList<-aprList [ ! sapply (aprList , is.null )]

mayList <- llply (sdsList , function (x){ if ( substr (x, 32, 34) %in% c ( "122" , "123" )){ return (x)}})
mayList<-mayList [ ! sapply (mayList , is.null )]

junList <- llply (sdsList , function (x){ if ( substr (x, 32, 34) %in% c ( "152" , "153" )){ return (x)}})
junList<-junList [ ! sapply (junList , is.null )]

julList <- llply (sdsList , function (x){ if ( substr (x, 32, 34) %in% c ( "182" , "183" )){ return (x)}})
julList <- julList [ ! sapply (julList , is.null )]

augList <- llply (sdsList , function (x){ if ( substr (x, 32, 34) %in% c ( "213" , "214" )){ return (x)}})
augList <-augList [ ! sapply (augList , is.null )]

sepList <- llply (sdsList , function (x){ if ( substr (x, 32, 34) %in% c ( "244" , "245" )){ return (x)}})
sepList <-sepList [ ! sapply (sepList , is.null )]

octList <- llply (sdsList , function (x){ if ( substr (x, 32, 34) %in% c ( "274" , "275" )){ return (x)}})
octList <-octList [ ! sapply (octList , is.null )]

novList <- llply (sdsList , function (x){ if ( substr (x, 32, 34) %in% c ( "305" , "306" )){ return (x)}})
novList <-novList [ ! sapply (novList , is.null )]

decList <- llply (sdsList , function (x){ if ( substr (x, 32, 34) %in% c ( "335" , "336" )){ return (x)}})
decList <-decList [ ! sapply (decList , is.null )]
# put above lists into another list for processing
NDVIMonths <- list (janList , febList , marList , aprList , mayList, junList , julList , augList , sepList , octList , novList , decList )

# loop through each month , looping again through all the years (2000 - 2013) of data reprojecting , merging tiles and converting to tiff
# format for each available month and year in one step

for (j in 1 : 12){
i <- 0
while (i < 14){
if (i < 10 & length (NDVIMonths [[ j ]] [ grep (paste0("A200" , "" ,i),NDVIMonths[[j]])]) >0){
gdalwarp ( srcfile= NDVIMonths [[ j ]] [ grep (paste0("A200" , "" ,i),NDVIMonths[[j]])], t_srs="+proj=longlat +datum=WGS84 +no_defs" ,dstfile=paste0("NDVI" ,j,"200" ,i,".tiff" ), te = c (bbox(NDVIRast) [1] ,bbox(NDVIRast) [2] ,bbox(NDVIRast) [3] , bbox (NDVIRast) [ 4 ] ) , tr = c (xres(NDVIRast), yres(NDVIRast)))
i <- i+1
}
else if (i < 10){
print (paste0("no data available for month" ,j, "and year 200" ,i))
i <- i+1
}
else if (i > 9 & length (NDVIMonths [[ j ]] [ grep (paste0("A20" , "" ,i),NDVIMonths[[j ]]) ]) > 0){
gdalwarp ( srcfile= NDVIMonths [[ j ]] [ grep (paste0("A20" , "" ,i),NDVIMonths[[j]])],t_srs="+proj=longlat +datum=WGS84 +no_defs" ,dstfile=paste0("NDVI" ,j,"20" , i, ".tiff" ),te = c (bbox(NDVIRast) [1] ,bbox(NDVIRast) [2] ,bbox(NDVIRast) [3] , bbox (NDVIRast) [ 4 ] ) , tr = c (xres(NDVIRast), yres(NDVIRast)))
i <- i+1
}
else if (i > 9 & i<14){
print (paste0("no data available for month" ,j, "and year 20" ,i))
i <- i+1
}
}
}

# (1) intialise empty stack for each climate and environmental factor
# (2) for each climate stack except NDVI merge monthly tiles , do appropriate raster calculations and crop to extent of model NDVI stack then add each monthly raster layer to the stack in a loop ( no
# loop needed for altitude )
# (3) for the NDVI montly layers apply NDVIRasterFunction to each monthly tiff image for each year , function involves calculating the mean NDVI value for each month across all years and cropping to
# same extent as botswana using NDVI as model raster

# Notes : temperature layers must be divided by 10 and NDVI layers must be divided by 10000

# put all WorldClim images in a list for processing
climFilelist <- list.files("/Volumes/JUSTJUBBA/Spatial/WorldClim_Botswana" ,pattern="tif$" ,full.names=T )

# initialize all climatic and environmental stacks
meanTempStack <- stack ()
maxTempStack <- stack ()
minTempStack <- stack ()
rainStack <- stack ()
bioStack <- stack ()
NDVIStack <- stack ()

NDVIRasterFunction<- function (y){
x <- stack ( list.files("/Volumes/JUSTJUBBA/Spatial/NDVI_MODIS_Botswana/" ,pattern= y, full.names=T ) )
x <- calc (x, function (z) z* 0.0001)
x <- mean (x, na.rm =T)
x <- crop (x , extent (NDVIRast) )
return (x)
}

for (i in 1:12){
NDVIStack <- stack (NDVIStack , NDVIRasterFunction( paste0 ("NDVI" ,i,20)))
}

for (i in 1 : 12){
meanTempStack <- stack (meanTempStack , crop( raster ( climFilelist [ grep (paste0("tmean" , i, "_36" ), climFilelist)])/ 10 , extent (NDVIRast) ) )

minTempStack <- stack (minTempStack , crop( raster ( climFilelist [ grep (paste0("tmin" ,i, "_36" ), climFilelist)])/ 10 , extent (NDVIRast) ) )

maxTempStack <- stack (maxTempStack , crop( raster ( climFilelist [ grep (paste0("tmax" ,i, "_36" ), climFilelist)])/ 10 , extent (NDVIRast) ) )

rainStack <- stack ( rainStack , crop( raster ( climFilelist [ grep (paste0("prec" ,i,"_36" ), climFilelist)]), extent(NDVIRast)))
}

for (i in 1 : 19){
if (i < 12){
bioStack <- stack (bioStack , crop( raster ( climFilelist [ grep (paste0("bio" ,i,"_36" ), climFilelist ) ])/ 10 , extent (NDVIRast) ) )
}
else{
bioStack <- stack (bioStack , crop( raster ( climFilelist [ grep (paste0("bio" ,i,"_36" ), climFilelist ) ]) , extent(NDVIRast)))
}
}
altitudeLayer <- crop ( raster ("/Volumes/JUSTJUBBA/Spatial/WorldClim_Botswana/alt_36.tif" ), extent(NDVIRast))


# resample WolrdClim layers to that of model NDVI layers
meanTempStack <- resample (meanTempStack , NDVIStack)
minTempStack <- resample (minTempStack , NDVIStack)
maxTempStack <- resample (maxTempStack , NDVIStack)
rainStack <- resample ( rainStack , NDVIStack)
bioStack <- resample ( bioStack , NDVIStack)
altitudeLayer <- resample ( altitudeLayer , NDVIStack)

# caclulate min annual temp
minAnnualTempLayer <- (minTempStack [[1]] + minTempStack [[2]] + minTempStack [[3]] + minTempStack [[ 4 ]] + minTempStack [[ 5 ]] + minTempStack [[ 6 ]] + minTempStack [[ 7 ]] + minTempStack [[ 8 ]] + minTempStack [[ 9 ]] + minTempStack [[10]] + minTempStack [[11]] + minTempStack [[12]] ) / 12

# name each layer in stack indicating appropriate month
names (NDVIStack) <- rep (paste0("NDVI" ,1:12))
names (rainStack) <- rep (paste0("rain" ,1:12))
names (meanTempStack) <- rep (paste0("meanTemp" ,1:12))
names (maxTempStack) <- rep (paste0("maxTemp" ,1:12))
names (minTempStack) <- rep (paste0("minTemp" ,1:12))
names (bioStack) <- rep (paste0("bio" ,1:19))
names (minAnnualTempLayer) <- "minAnnualTemp"
names (altitudeLayer) <- "altitude"


# create multilayered spatial database - last step #

# write SpatialPointsDataFrame values to polygon that are inside Botswana
writeOGR ( sitesSpdf , dsn = "getwd ()" ,layer=' sitesBots ' ,driver=' ESRI Shapefile ' , overwrite = T)

# read polygon as a new SPDF object
sitesPolyPoints = readOGR("getwd ()" , "sitesBots" )

# add coords onto @data component of SPDF
sitesPolyPoints@data = cbind (sitesPolyPoints@data , sitesPolyPoints@coords)

# extract raster values at matching coords and add to @data component for each climate layer in stack
# then give layers column names and append them to the SPDF

rainLayers <- raster ::: extract (rainStack , as (sitesPolyPoints , "SpatialPoints" ))
colnames (rainLayers) <- rep (paste0("rain" ,1:12))
sitesPolyPoints@data <- cbind (sitesPolyPoints@data , rainLayers)

meanTempLayers <- raster ::: extract (meanTempStack, as (sitesPolyPoints , "SpatialPoints" ))
colnames (meanTempLayers) <- rep (paste0("meanTemp" ,1:12))
sitesPolyPoints@data <- cbind (sitesPolyPoints@data , meanTempLayers)

maxTempLayers <- raster ::: extract (maxTempStack, as (sitesPolyPoints , "SpatialPoints" ))
colnames (maxTempLayers) <- rep (paste0("maxTemp" ,1:12))
sitesPolyPoints@data <- cbind (sitesPolyPoints@data , maxTempLayers)

minTempLayers <- raster ::: extract (minTempStack , as (sitesPolyPoints , "SpatialPoints" ))
colnames (minTempLayers) <- rep (paste0("minTemp" ,1:12))
sitesPolyPoints@data <- cbind (sitesPolyPoints@data , minTempLayers)

NDVILayers <- raster ::: extract (NDVIStack, as (sitesPolyPoints , "SpatialPoints" ))
colnames (NDVILayers) <- rep (paste0("NDVI" ,1:12))
sitesPolyPoints@data <- cbind (sitesPolyPoints@data , NDVILayers)

bioLayers <- raster ::: extract (bioStack , as (sitesPolyPoints , "SpatialPoints" ))
colnames (bioLayers) <- rep (paste0("bio" ,1:19))
sitesPolyPoints@data <- cbind (sitesPolyPoints@data , bioLayers)

minAnnualTemp <- raster ::: extract (minAnnualTempLayer , as (sitesPolyPoints , "SpatialPoints" ))
sitesPolyPoints@data <- cbind (sitesPolyPoints@data , minAnnualTemp)

altitude <- raster ::: extract (altitudeLayer ,as (sitesPolyPoints , "SpatialPoints" ))
sitesPolyPoints@data <- cbind (sitesPolyPoints@data , altitude)


# assign SPDF to DF for caclulation convenience
s.p <- as.data.frame (sitesPolyPoints)

# subset meanTemp variables
s.pTemp <- subset (s.p, select = rep (paste0("meanTemp" ,1:12)))
s.pRain <- subset (s.p, select = rep (paste0("rain" ,1:12)))
s.pNDVI <- subset (s.p, select = rep (paste0("NDVI" ,1:12)))

# caclulate annual std deviation of monthly variables
s.p$SDTemp <- apply (s.pTemp, 1, sd )
s.p$SDRain <- apply (s.pRain, 1, sd )
s.p$SDNDVI <- apply (s.pNDVI, 1, sd )

s.p$ totTemp <- apply (s.pTemp, 1, sum )
s.p$ totRain <- apply (s.pRain, 1, sum )
s.p$totNDVI <- apply (s.pNDVI, 1, sum )

# summer ( months 12 , 1, 2, 3) and winter caclulations (4 ,5 ,6 ,7 ,8 ,9 ,10)

s.p$summerTemp <- apply (s.p[c ( "meanTemp12" , "meanTemp1" , "meanTemp2" , "meanTemp3" )], 1, mean )
s.p$ winterTemp <- apply (s.p[c ( "meanTemp4" , "meanTemp5" , "meanTemp6" , "meanTemp7" , "meanTemp8" , "meanTemp9" , "meanTemp10" )], 1, mean )
s.p$summerRain <- apply (s.p[c ( "rain12" , "rain1" , "rain2" , "rain3" )], 1, mean )
s.p$ winterRain <- apply (s.p[c ( "rain4" , "rain5" , "rain6" , "rain7" , "meanTemp8" , "rain9" , "rain10" )], 1, mean )
s.p$summerNDVI <- apply (s.p[c ( "NDVI12" , "NDVI1" , "NDVI2" , "NDVI3" )], 1, mean )
s.p$ winterNDVI <- apply (s.p[c ( "NDVI4" , "NDVI5" , "NDVI6" , "NDVI7" , "NDVI8" , "NDVI9" , "NDVI10" )], 1, mean )

meanTempCols <- rep (paste0("meanTemp" ,1:12))
rainCols <- rep (paste0("rain" ,1:12))
NDVICols <- rep (paste0("NDVI" ,1:12))

# rainfall concentration index and mean peak month around which rainfall is concentrated calculations #

angle <- rep (paste0("angle" ,1:12))
s.pRain[angle] <- NA

for (i in 1:12){
s.pRain[angle][ i] <- (i* 2 * pi )/ 12
}

r1 <- rep (paste0("r1" ,1:12))
r2 <- rep (paste0("r2" ,1:12))
s.pRain[r1] <- NA
s.pRain[r2] <- NA

for (i in 1:12){
s.pRain[r1][ i] <- s.pRain[rainCols][ i]* cos (s.pRain[angle ][ i ])
s.pRain[r2][ i] <- s.pRain[rainCols][ i]* sin (s.pRain[angle ][ i ])
}

s.pRain$ r <- sqrt (( apply (s.pRain[r1] , 1, sum ))^2 + (apply (s.pRain[r2] , 1, sum ))^2)

# concentration index
s.p$ rCIndex <- (100* s.pRain$ r)/ s.p$ totRain

# mean peak month around which rainfall is concentrated
s.p$q <- atan ( apply (s.pRain[r2] , 1, sum ) / apply (s.pRain[r1] , 1, sum ))

# save workspace and write DF containing spatial database to table for later use
save.image ( ' spatial_Bots.RData ' )
write.table (s.p, "/Volumes/JUSTJUBBA/Spatial/s.pBots" )

# create prediction grid covering Botswana

# Notes : The same code principles apply as with the compilation of the spatial database above except at a lower resolution .
# Instead of extracting values at sample points all of the raster data is used so
# that an environmental or climatic value is present for each grid cell.All raster images
# are stored in stacks and the same process as above is used to get a corresponding attribute for each cell .
# Therefore this code is omitted .

# a prediction grid prepared using similar raster calculations as above
# covering Botswana at a 20km resolution was written to a table called gridPredBots20km :

write.table (grid20 , "/Volumes/JUSTJUBBA/Spatial/gridPredBots20km" )
