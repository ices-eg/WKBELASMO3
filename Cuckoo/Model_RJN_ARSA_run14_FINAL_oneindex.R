
rm(list=ls(all=TRUE))


library(spict)


dir.create("Model_results_ARSA_run14_oneindex")
dir.create("Model_results_ARSA_run14_oneindex/Plots")


############# SAVE DATA ###########

xtab<-function(x,caption='Table X.', file=stdout(), width='"100%"', cornername='', dec=rep(1,ncol(x))){
  nc<-ncol(x)
  lin<-paste('<table width=',width,'>', sep='')
  lin<-c(lin,sub('$','</td></tr>',sub('\\. |\\.$','.</b> ',
                                      sub('^', paste('<tr><td colspan=',nc+1,'><b>',sep=''), caption))))
  hr<-paste('<tr><td colspan=',nc+1,'><hr noshade></td></tr>', sep='')
  lin<-c(lin,hr)
  cnames<-colnames(x)
  cnames<-paste(sub('$','</b></td>',sub('^','<td align=right><b>',cnames)), collapse='\t')
  lin<-c(lin,paste('<tr>',paste('<td align=left><b>',cornername,'</b></td>',sep=''),cnames,'</tr>'))
  lin<-c(lin,hr)
  rnames<-sub('$','</b></td>',sub('^','<tr> <td align=left><b>',rownames(x)))
  #x<-sapply(1:ncol(x),function(i)sub('NA','  ',format(round(x[,i],dec[i]))))
  x<-sapply(1:ncol(x),function(i)sub('NA','  ',formatC(round(x[,i],dec[i]),digits=dec[i], format='f')))
  for(i in 1:nrow(x)){
    thisline<-paste(rnames[i],paste(sub('$','</td>',sub('^','<td align=right>',x[i,])), collapse='\t'),'</tr>', sep='')
    lin<-c(lin,thisline)
  }
  lin<-c(lin,hr)
  lin<-c(lin,'</table><br>\n')
  writeLines(lin,con=file)
}

######################################


## Load data - catch data

rjn.landings <- read.csv("rjn.27.9a_landings.csv", header=TRUE, sep=";")
names(rjn.landings)

rjn.index <- read.csv("ARSA_indices.csv", header=TRUE, sep=";", dec=",")
names(rjn.index)


## remove 1998-2000

rjn.index <- rjn.index[rjn.index$Year>2000,]

rjn.landings <- rjn.landings[rjn.landings$Year>2000,]


## create a list with these objects


## ARSA

rjn_data <- list(timeC = rjn.landings$Year, obsC = rjn.landings$Total_landings,
                 timeI = c(rjn.index$Year+3/12, rjn.index$Year+11/12),
                 obsI = c(rjn.index$Q1_B_kg, rjn.index$Q4_B_kg))



## plot series

plotspict.data(rjn_data)


## Plot data with basic model fitting

plotspict.ci(rjn_data)



## Prior fixing, if needed


### fix parameter to Scahaeffer production curve

rjn_data$phases$logn <- -1


# process error
rjn_data$priors$logsdb <- c(log(0.15), 0.5, 1)



## disable noise ratios logalpha and logbeta

rjn_data$priors$logalpha <- c(1, 1, 0)
rjn_data$priors$logbeta <- c(1, 1, 0)


## UNCERTAINTY - ajustar para este caso (este exemplo ? do lagostim)

rjn_data$stdevfacC <- c(rep(3,8), rep(1,6), rep(3,8)) # Capturas: incerteza nos 5 primeiros anos (2003-2007)e em 2010


## r prior


rjn_data$priors$logr <- c(log(0.41),0.5,1)



## Fitting the model

rjn_model <- fit.spict(rjn_data)
rjn_model


## plot model results

plot(rjn_model)


plotspict.priors(rjn_model)



### b/k estimate
exp(rjn_model$value[38])



#save Model parameter estimates w 95% CI
xtab(sumspict.parest(rjn_model),caption="ARSA_run14_Parameter estimates",cornername="Parameter",file="Model_results_ARSA_run14_oneindex/ARSA_run14_Parameter_estimates.xls",dec=rep(4,ncol(sumspict.parest(rjn_model))))
xtab(sumspict.states(rjn_model),caption="ARSA_run14_states of the model",cornername="",file="Model_results_ARSA_run14_oneindex/ARSA_run14_states of the model.xls",dec=rep(4,ncol(sumspict.states(rjn_model))))
xtab(sumspict.predictions(rjn_model),caption="ARSA_run14_predictions",cornername="",file="Model_results_ARSA_run14_oneindex/ARSA_run14_predictions.xls",dec=rep(4,ncol(sumspict.predictions(rjn_model))))


##check model adequacy (SPiCT guidelines)

rjn_model$opt$convergence #should equal 0

all(is.finite(rjn_model$sd)) #should be TRUE

calc.bmsyk(rjn_model) # should be between 0.1 and 0.9

calc.om(rjn_model)  # should not span more than 1 order of magnitude
xtab(calc.om(rjn_model),caption="ARSA_run14_magnitude",cornername="",file="Model_results_ARSA_run14_oneindex/ARSA_run14_magnitude.xls",dec=rep(4,ncol(calc.om(rjn_model))))


rjn_model_ini <- check.ini(rjn_model, ntrials=30)
rjn_model_ini$check.ini$resmat #the estimates should be the same for all initial values
xtab(rjn_model_ini$check.ini$resmat,caption="ARSA_run14_initial values",cornername="",file="Model_results_ARSA_run14_oneindex/ARSA_run14_initial values.xls",dec=rep(4,ncol(rjn_model_ini$check.ini$resmat)))



## Check residuals

rjn_model <- calc.osa.resid(rjn_model)
plotspict.diagnostic(rjn_model)


#covariance between the model parameters (fixed effects)
#rjn_model$cov.fixed
#correlation between the model parameters (fixed effects)
cov2cor(rjn_model$cov.fixed)
xtab(cov2cor(rjn_model$cov.fixed),caption="ARSA_run14_correlation between parameters",cornername="",file="Model_results_ARSA_run14_oneindex/ARSA_run14_correlation between parameters.xls",dec=rep(4,ncol(cov2cor(rjn_model$cov.fixed))))


## Retrospective analysis

retro_rjn_model <- retro(rjn_model, nretroyear=5)

plotspict.retro(retro_rjn_model)


##Hindcast
hc.rjn_model = hindcast(rjn_model,npeels =5)

plotspict.hindcast(hc.rjn_model,legend.pos = "topleft")

tiff(filename = "Model_results_ARSA_run14_oneindex/Plots/hindcast.tiff")
plotspict.hindcast(hc.rjn_model,legend.pos = "topleft")
dev.off()


## plot results

tiff(filename = "Model_results_ARSA_run14_oneindex/Plots/Fit_results_1.tiff")
plot(rjn_model)
dev.off()


## other plots


tiff(filename = "Model_results_ARSA_run14_oneindex/Plots/plot_data.tiff")
plotspict.data(rjn_data)
dev.off()

tiff(filename = "Model_results_ARSA_run14_oneindex/Plots/Priors.tiff")
plotspict.priors(rjn_model)
dev.off()


tiff(filename = "Model_results_ARSA_run14_oneindex/Plots/Biomass.tiff")
plotspict.biomass(rjn_model)
dev.off()

tiff(filename = "Model_results_ARSA_run14_oneindex/Plots/FishingM.tiff")
plotspict.f(rjn_model)
dev.off()

tiff(filename = "Model_results_ARSA_run14_oneindex/Plots/BBmsy.tiff")
plotspict.bbmsy(rjn_model)
dev.off()

tiff(filename = "Model_results_ARSA_run14_oneindex/Plots/FFmsy.tiff")
plotspict.ffmsy(rjn_model)
dev.off()


tiff(filename = "Model_results_ARSA_run14_oneindex/Plots/Diagnostics.tiff")
plotspict.diagnostic(rjn_model)
dev.off()

tiff(filename = "Model_results_ARSA_run14_oneindex/Plots/Retrospective.tiff")
plotspict.retro(retro_rjn_model)
dev.off()



## Tables


tab1 <- sumspict.parest(rjn_model);
xtab(tab1,caption="Parameter estimates",cornername="Parameter",file="Model_results_ARSA_run14_oneindex/Parameter_estimates.xls",dec=rep(4,ncol(tab1)))

tab2 <- sumspict.srefpoints(rjn_model);
xtab(tab2,caption="Stochastic reference points",cornername="Reference points",file="Model_results_ARSA_run14_oneindex/Reference_points.html",dec=rep(4,ncol(tab2)))

tab3 <- sumspict.states(rjn_model);
xtab(tab3,caption="Estimated states",cornername="",file="Model_results_ARSA_run14_oneindex/Estimated_states.html",dec=rep(4,ncol(tab3)))

tab4 <- sumspict.predictions(rjn_model);
xtab(tab4,caption="Forecast",cornername="",file="Model_results_ARSA_run14_oneindex/Forecast.html",dec=rep(4,ncol(tab4)))

tab5 <- get.par("logBBmsy",rjn_model,exp=TRUE)
tab5_<- tab5[grep(".",rownames(tab5), fixed=TRUE, invert=TRUE),] 
xtab(tab5_,caption="B/Bmsy",cornername="",file="Model_results_ARSA_run14_oneindex/BBmsy.html",dec=rep(4,ncol(tab5_)))
xtab(tab5,caption="B/Bmsy",cornername="",file="Model_results_ARSA_run14_oneindex/BBmsy_all.html",dec=rep(4,ncol(tab5)))
#write.csv(tab5_, "BBmsy.csv", dec=",")


tab6 <- get.par("logFFmsy",rjn_model,exp=TRUE)
tab6_<- tab6[grep(".9375",rownames(tab6), fixed=TRUE, invert=FALSE),] 
xtab(tab6,caption="F/Fmsy",cornername="",file="Model_results_ARSA_run14_oneindex/FFmsy.html",dec=rep(4,ncol(tab6)))
xtab(tab6_,caption="F/Fmsy",cornername="",file="Model_results_ARSA_run14_oneindex/FFmsy_end year.html",dec=rep(4,ncol(tab6_)))

tab6b <- get.par("logF",rjn_model,exp=TRUE)
xtab(tab6b,caption="F",cornername="",file="Model_results_ARSA_run14_oneindex/F_results.html",dec=rep(4,ncol(tab6b)))


tab6c <- get.par("logCpred",rjn_model,exp=TRUE)
xtab(tab6c,caption="Catch",cornername="",file="Model_results_ARSA_run14_oneindex/Catch.html",dec=rep(4,ncol(tab6c)))





### FORECAST


## The management period should be set before running check.inp or fitting the model
## Here the assumption is that the managment starts with one intermediate year
## The default is that the management period starts immediatelly after the last data point

## Time point to evaluate the model states - usually the end of the management period
rjn_data$maninterval <- c(2024, 2025)
rjn_data$maneval <- 2025


fit <- fit.spict(rjn_data)


fit <- add.man.scenario(fit, "F=0", ffac = 0)
fit <- add.man.scenario(fit, "F=Fsq", ffac = 1)
fit <- add.man.scenario(fit, "F=Fmsy", breakpointB = c(1/3, 1/2))
fit <- add.man.scenario(fit, "F=Fmsy_C_fractile", fractiles = list(catch = 0.35), breakpointB = c(1/3, 1/2))
fit <- add.man.scenario(fit, "F=Fmsy_C_fractile_20", fractiles = list(catch = 0.20), breakpointB = c(1/3, 1/2))
fit <- add.man.scenario(fit, "F=Fmsy_C_fractile_10", fractiles = list(catch = 0.10), breakpointB = c(1/3, 1/2))
fit <- add.man.scenario(fit, "F=Fmsy_C_fractile_15", fractiles = list(catch = 0.15), breakpointB = c(1/3, 1/2))

res <- sumspict.manage(fit, include.unc=TRUE,include.abs = TRUE)
res

write.csv(res, "Model_results_ARSA_run14_oneindex/Forecast_selected_scenarios.csv")


plotspict.hcr(fit)

plotspict.catch(fit)

fit$man


