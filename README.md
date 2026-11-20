# Malaria risk mapping with spatial GLMMs (MCom thesis code, 2014)

R code from my MCom Statistics thesis at Rhodes University (2012–2014): **spatial generalised linear mixed models for
malaria risk mapping in Botswana**. I wrote all of it myself in R, before AI coding assistants existed.

## What the code does

| Script | Step |
|---|---|
| [`R/01_spatial_database.R`](R/01_spatial_database.R) | Builds the spatial database: malaria survey sites loaded into MySQL, distance from each site to the closest surface water body, WorldClim climate layers (rainfall, temperature, bioclimatic variables) and MODIS NDVI, resampled to a common grid and extracted at each site |
| [`R/02_nonspatial_models.R`](R/02_nonspatial_models.R) | Non-spatial binomial GLMs of the proportion positive: an 85 / 15 derivation / validation split, candidate variables grouped into themes (temperature, rainfall, vegetation), correlated variables removed within each theme, and bootstrap stepwise selection by AIC in stages |
| [`R/03_bayesian_spatial_model.R`](R/03_bayesian_spatial_model.R) | Bayesian spatial GLM with a Gaussian process random effect (`spBayes::spGLM`, exponential covariance, adaptive MCMC), convergence checks over several chains, and prediction of malaria risk on a 20 km grid with `spPredict` |
| [`R/04_cross_validation.R`](R/04_cross_validation.R) | Validation on the held-out sites: mean error and mean absolute error of the non-spatial and the spatial model |

## About this copy

- **Transcribed from the thesis appendix (a PDF).** The original `.R` files are not available. The text was extracted from the
  typeset appendix and repaired where typesetting had damaged it: minus signs (`<-`), spaces inside names (`as . data . frame`),
  letters spread apart (`f i l e s`), split numbers (`[[ 1 0 ]]`). Every script parses in R. Comments are kept as written,
  typos included.
- **One change to the code itself:** a closing bracket missing from the printed appendix (marked in `02_nonspatial_models.R`).
- **Credentials removed:** the local MySQL user name and password now read from environment variables.
- **Not re-run since 2014.** Several packages it uses (`rgdal`, `maptools`) have since been retired from CRAN, so it is a
  record of the work, not a pipeline to run today.
- **No data included.** The survey data (`BotsMARA.csv`), climate, NDVI and surface-water layers are not part of this repository.
