
library(raster)
library(ncdf4)

x_list  <- list.files("G:/TROPOMI/esa/gridded/1deg/sun_global/clearsky", pattern = "*.nc", full.names = TRUE)
y_list  <- list.files("G:/SIF_OCO2_005", pattern = "*.nc", full.names = TRUE)
x_name   <- "SIF_Corr_743"
y_name   <- "sif_ann"
out_dir  <- "G:/SIF_comps/figs/raster_regressions/clearsky/sun_global"
out_name <- "TROPOMI_SIF_Corr_vs_Sun_SIF_1deg_clearsky_n30_2019-2020"
f_name   <- NA # Filter by value. Example, error, std, or n. If none use NA.
f_thresh <- 30 # Values >= will be kept

# Build bricks from list
for (i in 1:length(x_list)) {
  x_t <- raster(x_list[1], varname = x_name)
  y_t <- raster(y_list[1], varname = y_name)
  
  if (i == 1) {
    x_stack <- x_t
    y_stack <- y_t
  } else {
    x_stack <- addLayer(x_stack, x_t)
    y_stack <- addLayer(y_stack, y_t)
  }
}

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

rast_reg <- function(y_stack, x_stack, out_dir, out_name) {
  
  y <- brick(y_stack)
  x <- brick(x_stack)

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

rast_reg(y_stack, x_stack, out_dir, out_name)