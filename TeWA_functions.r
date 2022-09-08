

teWA_example = function(serieobs, spredictors, horizon, nsmooth=NULL, despl=0, Show=T, Plot=T, 
                        validation_obs=NULL, training=NULL, ylim=NULL)
{
  require(WaveletComp)
  
  if(is.null(training))
    training = round(4.5*horizon)
  if(is.null(nsmooth))
    nsmooth = round(horizon/4) 
  if(is.null(ylim))
    ylim = c(min(serieobs,na.rm=T),max(serieobs,na.rm=T))
  
  teleconnection_only = forecast_indexes(serieobs, spredictors, nsmooth=nsmooth, horizon=horizon)
  
  if(is.null(validation_obs))
    validation_obs = ts(rep(NA, horizon))
  mytimeserie = ts(c(serieobs,validation_obs),start=start(serieobs))
  
  
  recorte = NULL
  for(i in 1:dim(teleconnection_only$model)[1])
    recorte = cbind(recorte, tail(spredictors[,as.character(teleconnection_only$model[i,"indice"])],training+teleconnection_only$model[i,"lag"])[1:(training+horizon)])
  
  st0 = lm(tail(serieobs,training) ~ head(recorte,training))
  summary(st0)
  prediction_tele = st0$coefficients[1] + apply(st0$coefficients[-1]*t(recorte),2,sum)
  
  
  serieobs_aux = tail(serieobs,training)-st0$fitted.values
  serieobs_aux[is.na(serieobs_aux )] = 0
  
  
  val = data.frame(val=serieobs_aux)
  my.w = analyze.wavelet(val,"val",loess.span = 0,
                         dt = 1, dj = 1/50,
                         lowerPeriod = 4, upperPeriod = 100,
                         make.pval = TRUE, n.sim = 10)
  
  
  wavelet_fitted = reconstruct(my.w, plot.waves = FALSE, lwd = c(1,2),plot=F,
                               legend.coords = "bottomleft", ylim = c(-1.8, 1.8))
  
  
  wavelet_reconstruct = wavelet_fitted$series$val.r
  
  
  model3 = auto.arima(wavelet_reconstruct)
  model2=auto.arima(serieobs_aux)
  
  future3=forecast(model3, h=horizon)
  future2=forecast(model2, h=horizon)
  
  
  prediction_tele2 = prediction_tele + c(future2$x, future2$mean)
  prediction_tele3 = prediction_tele + c(future3$x, future3$mean)
  
  st0 = lm(tail(serieobs,training) ~ head(recorte,training))
  pval = summary(st0)$coefficients[-1,4]
  lg.pval = pval < 0.05
  prediction_tele = st0$coefficients[1] + apply(st0$coefficients[-1][lg.pval]*t(recorte[,lg.pval]),2,sum,na.rm=T)
  
  if(Show & !is.na(validation_obs[1]) & length(validation_obs)==length(teleconnection_only$simulation))
  {
    lg = !is.na(serieobs*teleconnection_only$simulation)
    cor_training = cor(teleconnection_only$simulation[lg],serieobs[lg])
    cor_validation = cor(tail(teleconnection_only$prediction,horizon),validation_obs)
    cat("R_long_training_tele = ",round(cor_training,2),", ", "R_short_validation_tele = ", round(cor_validation,2),"\n")
     
  }
  
  
  
  
  
  my.wx = analyze.wavelet(val,
                          dt = 1, dj = 1/20,
                          lowerPeriod = 6, upperPeriod = horizon,
                          make.pval = TRUE, n.sim = 10)
  maximum.level = 1.001*max(my.wx$Power.avg)
  #wt0 = wt.avg(my.wx, maximum.level = maximum.level, label.avg.axis = F, show.legend =F, plot=F)
  my.wx$Power.avg
  my.wx$Power.avg.pval
  my.wx$Period
  
  
  if(despl>0)
  {
    despl = as.numeric(head(tail(future3$mean,horizon),1)) + as.numeric(head(tail(prediction_tele,horizon),1)) -  as.numeric(tail(serieobs, 1))
    despl = despl/2
  }
  
  Apred = Arima(tail(serieobs,2*horizon+10), order=c(min(c(22,round(horizon/4))),1,1))
  Ensemb = forecast(Apred, horizon)
  Ensemb$mean = ts(future3$mean + as.numeric(tail(prediction_tele,horizon)) - despl, end=end(Ensemb$mean))
  Ensemb$lower[,1] = future3$lower[,1] + as.numeric(tail(prediction_tele,horizon)- despl)
  Ensemb$lower[,2] = future3$lower[,2] + as.numeric(tail(prediction_tele,horizon)- despl)
  Ensemb$upper[,1] = future3$upper[,1]  + as.numeric(tail(prediction_tele,horizon)- despl)
  Ensemb$upper[,2] = future3$upper[,2]  + as.numeric(tail(prediction_tele,horizon)- despl)
  Ensemb$x = tail(serieobs,training)
  Ensemb$fitted = wavelet_reconstruct+st0$fitted.values
  
  if(Show)
  {
    plot(tail(serieobs,training+horizon),type="l", ylab="Predictors", ylim=c(1,16), xlim=c(start(tail(serieobs,training+horizon))[1]-20,end(tail(serieobs,training+horizon))[1]) ,col="white", axes=F)
    for(i in 1:dim(spredictors)[2])
    {
      lines(ts(head(spredictors[,i],training)*0.15+i+1,start=start(tail(serieobs,training))),type="l")
      text(end(serieobs)+5,i+1,colnames(spredictors)[i],cex=0.5, adj = 0)
    }
    box(col="gray10", lty=2)
    
    plot(tail(serieobs,training+horizon),type="l", ylab="Observed time series (anomaly)",ylim=c(-2,2), xlim=c(start(tail(serieobs,training+horizon))[1]-20,end(tail(serieobs,training+horizon))[1]) ,axes=F)
    abline(h=0, lty=2,col="gray")
    box(col="gray10", lty=2)
    
    
    
    plot(tail(mytimeserie,training+horizon),type="l", ylab = "Selected best-delayed predictors", ylim=c(1,16), col="white", axes=F)
    abline(h=0, lty=2,col="gray")
    escogidos = paste(teleconnection_only$model[,1])
    for(i in 1:dim(spredictors)[2])
    {
      lines(ts(head(spredictors[,i],training)*0.15+i+1,start=start(tail(serieobs,training))),type="l",col="gray")
      text(end(serieobs)+5,i+1,colnames(spredictors)[i],cex=0.5, adj = 0,col="gray")
      
      if(colnames(spredictors)[i] %in% escogidos)
      {
        j = which(escogidos == colnames(spredictors)[i])[1]
        M = training-teleconnection_only$model[j,"lag"]
        lines(ts(head(spredictors[,escogidos[j]],training)*0.15+i+1,start=start(tail(serieobs,training))),type="l",col="blue")
        lines(ts(head(spredictors[,escogidos[j]],M)*0.15+i+1,start=start(tail(serieobs,training))),type="l",col="indianred")
        points(end(serieobs)[1]-teleconnection_only$model[j,"lag"],spredictors[M,escogidos[j]]*0.15+i+1, pch=20,col="blue",cex=0.5)
        text(end(serieobs)+5,i+1,teleconnection_only$model[j,"indice"],cex=0.5, adj = 0)
      }
    }
    box(col="violet")
    
    par(xpd=T)
    plot(tail(mytimeserie,training+horizon),type="l", ylim=c(1,16),col="white", axes=F,ylab="Final predictors", main="Building the multi-regression model",xlab="")
    rect(start(tail(mytimeserie,horizon))[1]-training-20,0,start(tail(mytimeserie,horizon))[1],17,col="papayawhip",border="papayawhip")
    rect(start(tail(mytimeserie,horizon))[1],0,start(tail(mytimeserie,horizon))[1]+100,17,col="lightblue1",border="lightblue1")
    text(start(tail(mytimeserie,horizon))[1]-training/2,-1,"Training period",col="indianred",cex=0.6)
    text(start(tail(mytimeserie,horizon))[1]+48.1,-0.5,"Prediction",col="blue",cex=0.6,)
    text(start(tail(mytimeserie,horizon))[1]+48.1,-1.5,"period",col="blue",cex=0.6)
    escogidos = paste(teleconnection_only$model[,1])
    par(xpd=F)
    for(i in 1:dim(spredictors)[2])
    {
      if(colnames(spredictors)[i] %in% escogidos)
      {
        j = which(escogidos == colnames(spredictors)[i])[1]
        lines(ts(head(recorte[,j],training+horizon)*0.15+i+1,start=start(tail(serieobs,training))),type="l",col="blue")
        lines(ts(head(recorte[,j],training)*0.15+i+1,start=start(tail(serieobs,training))),type="l",col="indianred")
        points(start(tail(mytimeserie,horizon))[1],recorte[training,j]*0.15+i+1, pch=20,col="blue",cex=0.5)
        
      }
    }
    abline(v=start(tail(mytimeserie,horizon)),lty=2,col="darkblue",lwd=0.8)
    
    readline(prompt="Press [enter] to continue")
    
    plot(tail(mytimeserie ,training+horizon), ylab="Anomaly", type="l", ylim=ylim,axes=T,col="white")
    lines(tail(serieobs, training+horizon))
    lines(ts(tail(prediction_tele,horizon),start=start(tail(mytimeserie ,horizon))),col="blue", ylim=c(-1.9,2.4),axes=F,lwd=2)
    lines(ts(st0$fitted.values,start=start(tail(serieobs,training))),col="indianred",lwd=2)
    abline(v=start(tail(mytimeserie ,horizon)),lty=2,col="darkblue")
    abline(h=0, lty=2,col="gray")
    #text(start(tail(serieobs,N))+N/3,1.5,"Fitted model",cex=0.8,col="indianred")
    box(col="violet")
    legend("topleft",legend=c("Fitted model"),lwd=2,col="indianred",text.col="indianred",cex=0.7, bty="n")
    
    readline(prompt="Press [enter] to continue")
    
    plot(tail(mytimeserie,training+horizon),col="white", axes=T,ylab="Residual anomaly", ylim = ylim)
    lines(ts(serieobs_aux,start=start(tail(serieobs,training))),col="darkgreen")
    abline(h=0, lty=2,col="gray")
    box(col="yellow3")
    
    legend("topleft",legend=c("Residual signal"),lwd=2,col="darkgreen",text.col="darkgreen",cex=0.7, bty="n")
    
    readline(prompt="Press [enter] to continue")
    
    
    wt.image(my.w, color.key = "quantile", n.levels = 250, label.time.axis=T,label.period.axis=T,
             periodlab = "", timelab = "",
             legend.params = list(lab = "wavelet power levels", mar = 5.1,cex.axis=0.5,cex=0.5))
    
    readline(prompt="Press [enter] to continue")
    
    plot(my.wx$Power.avg,my.wx$Period,typ="l", axes=T,log="y", ylab="Period", xlab="Power")
    lines(my.wx$Power.avg[my.wx$Power.avg.pval<0.05],my.wx$Period[my.wx$Power.avg.pval<0.05],typ="l", col="red",lwd=2)
    box(col="yellow3")
    legend("bottom",legend=c("Amplitude with p-vale < 0.05"),col="red", lty=1, lwd=2,cex=0.6, bty="n")
    
    readline(prompt="Press [enter] to continue")
    
    plot(future3, axes=T, main="", xlab="Time", ylab="Residual anomaly", ylim=ylim)
    lines(ts(wavelet_fitted$series$val,start=start(future3)))
    lines(ts(wavelet_reconstruct,start=start(future3)),col="indianred",lwd=2)
    abline(h=0, lty=2,col="gray")
    box(col="yellow3")
    
    readline(prompt="Press [enter] to continue")
    fitted = ts(wavelet_reconstruct+st0$fitted.values,end=start(Ensemb$mean))
    
    plot(Ensemb, axes=T, main="", xlab="Time",  ylab="Total anomaly", ylim=ylim, xlim=c(start(Ensemb$upper)[1]-training,start(Ensemb$upper)[1]+horizon))
    lines(fitted,col="indianred",lwd=2)
    
    validation = ts(validation_obs,start=start(Ensemb$mean))
    lines(ts(validation_obs,start=start(Ensemb$mean)))
    
    abline(h=0, lty=2,col="gray")
    box(col="darkblue")
    title("Final prediction")
    
  }
  
  return(list(pred=Ensemb, fitted = ts(wavelet_reconstruct+st0$fitted.values,end=start(Ensemb$mean))))
}


#-------------------------------------------------------------------------#
#-------------------------------------------------------------------------#


smooth_with_wavelet = function(val, lowerPeriod =3, upperPeriod = 100)
{
  my = analyze.wavelet(data.frame(val=val),"val",loess.span = 0,
                       dt = 1, dj = 1/50,
                       lowerPeriod = lowerPeriod, upperPeriod = upperPeriod,
                       make.pval = TRUE, n.sim = 10)
  my.rec = reconstruct(my, plot.waves = FALSE, lwd = c(1,2),plot=F)
  
  
  return(my.rec$series$val.r)
}

#-------------------------------------------------------------------------#
#-------------------------------------------------------------------------#

forecast_indexes = function(serieobs, spredictors, nsmooth, horizon, plot=F, print=F, order=0)
{
  model_lags = data.frame(indice=NULL,lag=NULL,cor=NULL,coef=NULL)
  model_inds = NULL
  model_futs = NULL
  if(plot) 
    par(mfrow=c(3,5))
  for(ipred in 1:dim(spredictors)[2])
  {
    #Vamos a calcular el suavizado tipo ma(ULMO[lg_fech_ULMO,2][lg.ajuste],nsmooth)
    #indice = ma(get(indices_emp[ipred])[get(paste0("lg_fech_",indices_emp[ipred])),2][lg.ajuste],nsmooth)
    
    indice = as.ts(SMA(spredictors[,ipred],nsmooth)[,1])
    out = search_lag_corr(serieobs, indice, horizon=horizon, main=colnames(spredictors)[ipred], plot=plot, control=T)
    if(print)
    {
      print(colnames(spredictors)[ipred])
      print(out)
    }
    
    n = length(indice)
    if(length(out$lag)>0)
      for(ilag in 1:length(out$lag))
        if(out$lag[ilag] > horizon)
        {
          lag = out$lag[ilag]
          indice_lag = c(rep(NA,lag),indice[-c((n-lag+1):n)])
          indice_fut = c(rep(NA,lag),indice)[1:(n+horizon)]
          model_lags = rbind(model_lags, data.frame(indice = colnames(spredictors)[ipred], lag=out$lag[ilag], cor=out$cor[ilag]))
          model_inds = cbind(model_inds,indice_lag)
          model_futs = cbind(model_futs,indice_fut)
        }
    
  }
  colnames(model_inds) = paste(model_lags$indice,model_lags$lag,sep="-")
  colnames(model_futs) = paste(model_lags$indice,model_lags$lag,sep="-")
  st0 = lm(serieobs ~ 0+model_inds)
  simulation = apply(st0$coefficients*t(model_inds),2,sum)
  prediction = apply(st0$coefficients*t(model_futs),2,sum)
  simulation = ts(simulation, start=start(serieobs))
  prediction = tail(ts(prediction, start=start(serieobs)),horizon)
  
  
  Apred = forecast_arima(tail(serieobs,2*horizon+10), horizon, order=c(order,1,1))
  inercial = Apred$mean
  #inercial = as.numeric(tail(serieobs,1))
  
  ###if(class(Wfore)!="try-error")
  ###  inercial = as.numeric(wavelet_pred)
  
  peso = 1/seq(1,horizon,1)^(1-1/seq(1,horizon,1)^0.15)
  prediction = inercial*peso + prediction*(1-peso)
  
  #peso = rev(1/seq(1,length(simulation),1)^(1-1/seq(1,length(simulation),1)^0.15))
  #simulation = as.numeric(tail(serieobs,1))*peso + simulation*(1-peso)
  
  prediction = ts(prediction, start=end(serieobs))
  simulation = ts(simulation, start=start(serieobs))
  
  return(list(model=model_lags, simulation=simulation, prediction=prediction, arima=Apred, factor=factor))
}


search_lag_corr = function(serieobs, indice, horizon, plot=F, signif.lev=0.999999, control=F, main="serieobs vs Indice")
{
  umbral = qnorm((1 + signif.lev)/2)/sqrt(sum(!is.na(serieobs)))
  
  ccf0 = ccf(as.numeric(serieobs), as.numeric(indice), lag.max=2*horizon, plot=plot, main=main, na.action = na.pass)
  
  if(plot)
    abline(h=c(-umbral,umbral),col="red",lty=2)
  d0 = diff(ccf0$acf); 
  selecciona_lags = d0>=0 & c(d0[-1],1) <= 0 & ccf0$acf[-1] > 0 |  d0<=0 & c(d0[-1],-1) >= 0 & ccf0$acf[-1] <  0
  out = data.frame(lag=ccf0$lag[-1], cor= ccf0$acf[-1])[ ccf0$lag[-1] >= -2 & selecciona_lags & abs(ccf0$acf[-1]) > umbral, ]
  
  if(control & dim(out)[1]==0)
    out = data.frame(lag=0, cor = ccf0$acf[ccf0$lag==0])
  
  if(control)
    out[out$lag<0,] = data.frame(lag=0, cor = ccf0$acf[ccf0$lag==0])
  
  if(plot)
    abline(v=out$lag,col="green",lty=2)
  
  
  return(out)
}
forecast_arima = function(x, horizon, order=0)
{
  fit =  try(Arima(x, order=order),silent=T)
  if(class(fit)[1]=="try-error")
    fit =  try(Arima(x, order=c(0,1,1)),silent=T)
  
  Arima.pred = ts(forecast(fit, horizon),start=end(x))
  
  return(Arima.pred)
  
}
