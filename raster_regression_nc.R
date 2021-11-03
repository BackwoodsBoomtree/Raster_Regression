library(terra)

in_file  <- "G:/TROPOMI/esa/gridded/TROPOMI.SIF.201805-202108.0.20deg.modisLike.esa.735nm.nc"
x_name   <- "NIRv"
y_name   <- "SIF_743"
out_dir  <- "G:/Russell/Projects/SLUE/raster_regressions"
out_name <- "TROPOMI_SIF_vs_NIRv_0.20_conus_full_record"

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

rast_reg <- function(in_file, x_name, y_name, out_dir, out_name) {
  
  data <- sds(in_file)
  y <- data[[y_name]]
  x <- data[[x_name]]
  
  y <- mask(y, x)
  x <- mask(x, y)
  
  # Combine bricks into single brick. lm convention is y~x
  yx <- c(y, x)
  
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

rast_reg(in_file, x_name, y_name, out_dir, out_name)