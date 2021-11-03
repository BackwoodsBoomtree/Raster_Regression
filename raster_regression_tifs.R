library(terra)

x_dir    <- "//isilondata.rccc.ou.edu/eomfdata/TROPOMI/gridded/official-Philipp/4326/fixed/2018"
y_dir    <- "G:/MCD43C4/0.20/8-day_nirv"
out_dir  <- "C:/Russell/R_Scripts/TROPOMI_2/regression_rasters/8-day/0.20/MCD43"
out_name <- "MCD43_NIRvvsSIF_8-day"

rastlm <- function(x) {
  full <- length(x)
  half <- full / 2
  
  if (all(is.na(x[1:half])) || all(is.na(x[(half + 1):full]))){ 
    
    return(c(NA,NA,NA,NA,NA,NA))
    
  } else { 
    reg       <- lm(x[1:half] ~ x[(half +1):full])
    s         <- summary(reg)
    r2        <- s$r.squared
    pval      <- s$coefficients[8]
    slope     <- s$coefficients[2]
    intercept <- s$coefficients[1]
    rmse      <- sqrt(mean(s$residuals^2))
    n         <- nobs(reg)
    
    return(c(r2, pval, slope, intercept, rmse, n)) 
  }
}

rast_reg <- function(x_dir, y_dir, out_dir, out_name) {
  
  y <- stack(list.files(y_dir, full.names = TRUE, pattern = "*.tif"))
  x <- stack(list.files(x_dir, full.names = TRUE, pattern = "*.tif"))
  
  y <- mask(y, x)
  x <- mask(x, y)
  
  # Combine bricks into single brick. lm convention is y~x
  yx <- stack(y, x)
  
  beginCluster(12)
  lm.result <- clusterR(yx, calc, args = list(fun = rastlm))
  endCluster()
  
  writeRaster(lm.result[[1]], paste0(out_dir, "/", out_name, "_Rsquare.tif"), overwrite=TRUE)
  writeRaster(lm.result[[2]], paste0(out_dir, "/", out_name, "_Pval.tif"), overwrite=TRUE)
  writeRaster(lm.result[[3]], paste0(out_dir, "/", out_name, "_Slope.tif"), overwrite=TRUE)
  writeRaster(lm.result[[4]], paste0(out_dir, "/", out_name, "_Intercept.tif"), overwrite=TRUE)
  writeRaster(lm.result[[5]], paste0(out_dir, "/", out_name, "_RMSE.tif"), overwrite=TRUE)
  writeRaster(lm.result[[6]], paste0(out_dir, "/", out_name, "_Nobs.tif"), overwrite=TRUE)
  
  remove(x, y, yx) # get it out of memory
  
}

rast_reg(x_dir, y_dir, out_dir, out_name)