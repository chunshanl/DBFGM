### Process raw sst data in .nc format
### Output data of shape n * T * p
library(ncdf4)
library(anytime)

setwd("./Helper_functions")
source('Helper_functions_sst.R')
setwd("..")

### Select locations
selected_points = c(-130,30,
                    -160, 30,
                    150, 30,
                    170, 30,
                    -120,10,
                    -150, 10,
                    160, 10,
                    179, 10,
                    -90, -10,
                    -110, -10,
                    -140, -10,
                    -170, -10,
                    -90, -30,
                    -110, -30,
                    -140, -30,
                    -170, -30)
locations = matrix(selected_points, ncol = 2, byrow = T)

# ### Load .nc data from different years, filter locations, save as csv
# Raw data is downloaded from the ERA5 reanalysis
setwd('./sst_data')
filenames <- list.files(pattern="*.nc", full.names=TRUE)
for (cur_filename in filenames){
  nc = nc_open(cur_filename)

  print(nc$var$sst$dim[[1]]$name)
  longitude = nc$var$sst$dim[[1]]$vals
  print(max(longitude))

  print(nc$var$sst$dim[[2]]$name)
  #print(nc$var$sst$dim[[1]]$vals)
  latitude = nc$var$sst$dim[[2]]$vals
  # plot(nc$var$sst$dim[[2]]$vals)

  hour_since_1900 = nc$var$sst$dim[[3]]$vals
  timestamp = as.POSIXct("1900-01-01 00:00:00.0", tz = 'UTC' ) + hour_since_1900 * 3600
  print(max(timestamp))
  print(min(timestamp))

  data = ncvar_get(nc, 'sst')   # longitude, latitude, time

  nc_close(nc)

  ## filter data
  data_filtered = c()
  temp1 = 1:length(longitude)
  temp2 = 1:length(latitude)
  for (i in 1:dim(locations)[1]){
    ind1 = temp1[locations[i,1] == longitude]
    ind2 = temp2[locations[i,2] == latitude]
    data_filtered = c(data_filtered, data[ind1, ind2,])
  }
  data_filtered = matrix(c(data_filtered, hour_since_1900), nrow = length(hour_since_1900))

  ## Save data
  cur_name = paste(format(timestamp[1], format = '%Y'), '.csv', sep = "")
  write.table(data_filtered, file=cur_name,
              row.names = F, col.names = F, sep = ",")

}

### Combine all files from different years
data_all = data.frame(matrix(ncol = 17, nrow = 0))
# setwd('./sst_data')
filenames <- list.files(pattern="*.csv", full.names=TRUE)
for (cur_name in filenames){
  cur_data = read.csv(cur_name, header = F)
  data_all = rbind(data_all, cur_data)
}
cur_name = paste('sst_update.csv')
write.table(data_all, file=cur_name,
            row.names = F, col.names = F, sep = ",")

### Read and reshape data --------
sst_input = read.csv('sst_data/sst_update.csv')
p = dim(sst_input)[2] -1
# Get time
data_time = sst_input[,p + 1]
timestamp = as.POSIXct("1900-01-01 00:00:00.0", tz = 'UTC' ) + data_time * 3600
# Get temperature
sst_input = sst_input[,1:p]

## Moving average
smooth_data = TRUE
if (smooth_data){
  sst_input_smoothed = sst_input
  for (j in 1:p){
    sst_input_smoothed[[j]] = moving_average(sst_input[[j]], window_length = 13)
  }
  #plot(sst_input_smoothed[1:366,1],type = 'l')
  #lines(sst_input[1:366,1],type = 'l', col = 2)
  sst_input = sst_input_smoothed
}

## Reshape data into n*T*p
data_time = anytime(timestamp)
n = 2022-1979 + 1
T_data = 366
cur_year = format(data_time, '%Y')  # char
col_num = as.numeric(cur_year) - 1978
cur_day = (data_time - anytime(cur_year))/24
row_num = round(as.numeric(cur_day)) + 1

Y = array(0, dim = c(n, T_data, p))
for (i in 1:nrow(sst_input)){
  Y[col_num[i], row_num[i],] = as.numeric(sst_input[i,])
}
Y = Y[,-(365:366),]

## Start from March
Y_chop = array(NA, c(n-1, 364, p))
ind_start = (31+29)
ind_end = 364*n - (364 - 31 - 28)
for(j in 1:p){
  y_vec = c(t(Y[,,j]))
  y_vec = y_vec[ind_start:ind_end]
  y_matrix = matrix(y_vec, nrow = 364, ncol = n-1)
  Y_chop[,,j] = t(y_matrix)
  
}
Y = Y_chop
j=2
matplot(t(Y[,,j]),xlab = '', ylab = '', type = 'l')

# Standardize data
Y_sd = Y
n = dim(Y)[1]
for (i in 1:n){
  for(j in 1:p){
    Y_sd[i,,j] = (Y[i,,j] - mean(Y[i,,j]))/sd(Y[i,,j])
  }
}
data = list()
data$Y = Y_sd

save(data, file= 'sst_processed_update.Rdata')

## Plot data -----------------------------------------
i = 1
j = 2
the_plot <- function()
{
  matplot(t(Y[,,j]), type = 'l', xlab = 'Date since March 1st', 
          ylab = 'SST (degrees K)')
}
png(
  "test.png",
  width     = 6,
  height    = 3.25,
  units     = "in",
  res       = 1200,
  pointsize = 4
)
par(
  mgp=c(6,2,0),
  mar      = c(10, 8, 3, 3),
  xaxs     = "i",
  yaxs     = "i",
  cex.axis = 2.5,
  cex.lab  = 2.7
)
the_plot()
dev.off()


