##### EXAMPLE OF TeWA PREDICTION ######

#### Reading functions and data####
library("WaveletComp")
library("TTR")
library("ggplot2")
library("forecast")
library("TSPred") 
library("WaveletArima")

setwd("C:/Users/rober/Dropbox/Rutinas/")

# Uploading the TeWA functions:
source("TeWA_functions.r")

# Uploading the raw predictors:
rawpredictors = readRDS("predictores.rds")
smtpredictors = readRDS("smtpredictores.rds")
already_smoothed = TRUE

# Uploading the observed time series of the target variable to forecast it: 
# It is necessary to use the same length than 'rawpredictors'
serieobs = readRDS("serieobs.rds")

# Setting time horizon ahead for the forecast
horizon = 90

# Optional: Setting scale of the smoothing of seasonal signal
# It depends on the horizon and on the (quasiperiodic) signal/noise ratio of the target variable
nsmooth = round(horizon/4) 

# Optional: Setting the preferred time windows to analyse ARIMA
# It depends on the horizon and on the (quasiperiodic) signal/noise ratio of the target variable
training = 4.5*horizon

# Optional: If exists, we can also upload the real observations of the period to validate
# It is necessary to use the same length than the value of 'horizon'
validation_obs = readRDS("validacion_obs.rds")



####  Before forecasting: It is recommended to smooth the atmospheric indices to reduce noise #### 
if(!already_smoothed)
{
  atmospherics = c("ULMO", "WeMO", "MO", "AJSL", "GJSL")
  spredictors= rawpredictors
  for(i in 1:dim(rawpredictors)[2])
    if(colnames(rawpredictors)[i] %in% atmospherics)
      smtpredictors[,i] = smooth_with_wavelet(rawpredictors[,i])
  
  saveRDS(smtpredictors, "smtpredictores.rds")
}





####  Now, we can apply the smoothed predictors #### 

out = teWA_example(serieobs, smtpredictors,  horizon=horizon, nsmooth=nsmooth, training=training, Show=T, validation_obs=NULL)



####  Finally, we represent the prediction #### 

plot(out$pred, axes=T, main="Final prediction", xlab="Time", ylab="Anomaly", xlim=c(start(out$pred$mean)[1]-training,start(out$pred$mean)[1]+horizon))
lines(out$fitted ,col="indianred",lwd=2)
lines(ts(validation_obs,start=start(out$pred$mean)))
abline(h=0, lty=2,col="gray")
box(col="darkblue")


