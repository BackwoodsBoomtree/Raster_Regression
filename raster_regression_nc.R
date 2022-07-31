
library(raster)

in_file  <- "G:/TROPOMI/esa/gridded/1deg/monthly/ebf/TROPOMI.ESA.SIF.2018-2021.EBF.monthly.1deg.clearsky.cold.nc"
x_name   <- "NIRv_RAD"
y_name   <- "SIF_743"
out_dir  <- "G:/SIF_comps/figs/raster_regressions/clearsky_cold"
out_name <- "TROPOMI_SIF743_vs_NIRv_RAD_monthly_1deg_clearsky_cold_n30_2018-2021"
f_name   <- "n" # Filter by value. Example, error, std, or n. If none use NA.
f_thresh <- 30 # Values >= will be kept

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
  
  y <- brick(in_file, varname = y_name)
  x <- brick(in_file, varname = x_name)
  
  y <- mask(y, x)
  x <- mask(x, y)
  
  # Filter as needed
  if (!is.na(f_name)) {
    print(paste0("Filtering data using: ", f_name))
    print(paste0("Filter threshold is: ", f_thresh))
    
    f <- brick(in_file, varname = f_name)
    f[f < f_thresh] <- NA
    
    y <- mask(y, f)
    x <- mask(x, f)

  }
  
  # Combine bricks into single brick. lm convention is y~x
  yx <- stack(y, x)
  
  beginCluster(12)
  lm.result <- clusterR(yx, calc, args = list(fun = rastlm))
  endCluster()
  
  writeRaster(lm.result[[1]], paste0(out_dir, "/", out_name, "_Rsquare.tif"), overwrite = TRUE)
  writeRaster(lm.result[[2]], paste0(out_dir, "/", out_name, "_Pval.tif"), overwrite = TRUE)
  writeRaster(lm.result[[3]], paste0(out_dir, "/", out_name, "_Slope.tif"), overwrite = TRUE)
  writeRaster(lm.result[[4]], paste0(out_dir, "/", out_name, "_Intercept.tif"), overwrite = TRUE)
  writeRaster(lm.result[[5]], paste0(out_dir, "/", out_name, "_RMSE.tif"), overwrite = TRUE)
  writeRaster(lm.result[[6]], paste0(out_dir, "/", out_name, "_Nobs.tif"), overwrite = TRUE)
  
  remove(x, y, yx) # get it out of memory
  
}

rast_reg(in_file, x_name, y_name, out_dir, out_name)