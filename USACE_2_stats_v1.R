setwd('~/Documents/psu/ESM 566/snowpack-project')


usace_2_sage <- read.csv('Data/USACE_2_allhourlydata.csv')

usace_2 <- subset(usace_2_sage, select = c(X.1,
                                           X,
                                           TIMESTAMP,
                                           RECORD,
                                           Temp_Avg, 
                                           RH,
                                           Flux_density_Avg,
                                           Flux_Density_2_Avg,
                                           WS_ms,
                                           WindDir,
                                           D_snow_temp_Avg,
                                           VW_Avg,
                                           VW_2_Avg,
                                           VW_3_Avg))

# calculate snow depth by subtracting the surface distance measure from 5 meters
usace_2$snow_depth <- (4.62 - usace_2$D_snow_temp_Avg)

# remove NAs

# periods where snow depth is decreasing for six hours
# defining one of these periods: 
# the endpoint is less than the starting point?
# monotonically decreases the whole time?

# for each hour, if snow depth in each of the following five hours is less than the one before,
# save that chunk of time as one of the periods into the new data frame

usace_2_sixhourdecreases <- data.frame()

#### My initial too strict implementation of six hour period as strictly monotonic
for (i in 1:nrow(usace_2)) {
  na <- 0
  for (j in 0:5) {
    if (is.na(usace_2[i+j,]$snow_depth)) {
      na <- 1
    }
  }
  if (na == 1) {
    next
  }
  if (usace_2[i,]$snow_depth > usace_2[i+1,]$snow_depth) {
    if (usace_2[i+1,]$snow_depth > usace_2[i+2,]$snow_depth) {
      if (usace_2[i+2,]$snow_depth > usace_2[i+3,]$snow_depth) {
        if (usace_2[i+3,]$snow_depth > usace_2[i+4,]$snow_depth) {
          if (usace_2[i+4,]$snow_depth > usace_2[i+5,]$snow_depth) {
            rbind(usace_2_sixhourdecreases, usace_2[i:i+5,])
          }
        }
      }
    }
  }
}

