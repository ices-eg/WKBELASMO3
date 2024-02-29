
rm(list=ls(all=TRUE))


#library(icesAdvice)
#library(icesTAF)
library(TMB)
library(spict)


dir.create("Model_results_Scenario50")
dir.create("Model_results_Scenario50/Plots")

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
## Load data ####


rjh <- read.csv("rjh.27.9a_input.csv", header=TRUE, sep=";")
names(rjh)
colnames(rjh)[1] <- "Year"

names(rjh)

## nao considerar desembarques iniciais
#rjh <- rjh[rjh$Year>2007,] 


## BIOMASS INDICES TO CONSIDER

PT_LPUE <- rjh$LPUE
PT_LPUE <- PT_LPUE[-(1:8)]


## TIME

YEAR <- rjh$Year

YEAR_LPUE <- seq(2008, 2022, 1)


## create a list with these objects, biomass index set at the middle of the year (+0.5; remove to set at the beggining of the year)

## CENARIO 1 (PT_LPUE)

rjh_data <- list(timeC = YEAR, obsC = rjh$Landings,
                 timeI = list(YEAR_LPUE+0.5),
                 obsI = list(PT_LPUE))



## plot series

plotspict.data(rjh_data)


## Plot data with basic model fitting (linear regression)
# two top plots come from plotspict.data, with the dashed horizontal line representing a guess of MSY.
# guess comes from a linear regression between the index and the catch divided by the index (middle row,
# left). This regression is expected to have a negative slope. A similar plot can be made showing catch versus
# catch/index (middle row, right) to approximately find the optimal effort (or effort proxy). The proportional
# increase in the index as a function of catch (bottom row, right) should show primarily positive increases
# in index at low catches and vice versa. Positive increases in index at large catches could indicate model
# violations. In the current plot these are not seen.

plotspict.ci(rjh_data)



## UNCERTAINTY TO INPUT DATA (if needed)

#inp$stdevfacC <- c(rep(2,9), rep(1,27)) # Capturas: incerteza nos 9 primeiros anos 
#rjh_data$stdevfacC <- c(rep(6,5), rep(1,14)) # Capturas: incerteza nos 5 primeiros anos (2003-2007)e em 2010
rjh_data$stdevfacC <- c(rep(2,8), rep(1,15)) #capturas: incerteza para o periodo 2000-2007
# rjh_data$stdevfacI <- c(rep(1,11), 2,1,1) # indice: incerteza em 2017 e 2018
# rjh_data$stdevfacI <- c(rep(2,14) # indice: incerteza em 2017 e 2018






## Prior fixing ####

#if needed...


#fix parameter to Schaeffer production curve (initial parameter)

rjh_data$phases$logn <- -1 # for Schaeffer

#rjh_data$priors$logn<-c(log(2),0.2,1)  # for Tighter Schaefer


# add prior to the initial depletion (rjh_data$priors$logbkfrac = c(log(0.5),0.5,1) and c(log(0.3),0.5,1))


#rjh_data$priors$logbkfrac = c(log(0.5),0.5,1) #initial deploection of 0.2, cv


# prior no r, valor de acordo com leslie matrix

rjh_data$priors$logr <- c(log(0.22),0.5, 1) 


# alpha e Beta 
rjh_data$priors$logalpha  <- c(1, 1, 0)
rjh_data$priors$logbeta   <- c(1, 1, 0)

# priors for logsdb (to correct overfitting of the production curve; when sdb estimate is very low)

rjh_data$priors$logsdb <- c(log(0.15), 0.5, 1) # prior value provided by Tobias (pers. comm) in WKLIFE XII



## Fitting the model ####

rjh_model <- fit.spict(rjh_data)

## summary results of the model

rjh_model
summary(rjh_model)


write.csv(capture.output(summary(rjh_model)), "Model_results_Scenario50/run_1fit.csv")

#Extrair o valor de bkfrac:

exp(rjh_model$value[34])#0.2722013

# retrieve the score for "objective function at optimum"
rjh_model$opt$objective 
# 9.273344



## plot model results

plot(rjh_model)

plotspict.priors(rjh_model)

#save Model parameter estimates w 95% CI
xtab(sumspict.parest(rjh_model),caption="Scenario50_Parameter estimates",cornername="Parameter",file="Model_results_Scenario50/Scenario50_Parameter_estimates.xls",dec=rep(4,ncol(sumspict.parest(rjh_model))))
xtab(sumspict.states(rjh_model),caption="Scenario50_states of the model",cornername="",file="Model_results_Scenario50/Scenario50_states of the model.xls",dec=rep(4,ncol(sumspict.states(rjh_model))))
xtab(sumspict.predictions(rjh_model),caption="Scenario50_predictions",cornername="",file="Model_results_Scenario50/Scenario50_predictions.xls",dec=rep(4,ncol(sumspict.predictions(rjh_model))))



# Estimates of deterministic reference points with 95% CIs. These are the reference points
# one would derive if stochasticity were ignored. Can be extracted with sumspict.drefpoints(res).
# Estimates of stochastic reference points with 95% CIs. These are the reference points of
# the stochastic model. The column ?rel.diff.Drp? shows the relative difference when compared to the
# deterministic reference points. The information can be extracted with sumspict.srefpoints(res).



##check model adequacy ####

# Assessment convergence if the model converged should equal 0 
rjh_model$opt$convergence # =>0

##If all variance parameters of the model parameters are finite should be TRUE 
all(is.finite(rjh_model$sd)) #=> TRUE

# Check for realistic production curve
#The shape of the production curve should not be too skewed, should be between 0.1 and 0.9 
calc.bmsyk(rjh_model) #=> 0.5

#Check for assessment uncertainty
## High assessment uncertainty can indicate a lack of contrast in the input data or violation of the ecological model assumptions. The main variance parameters (logsdb, logsdc, logsdi, logsdf) should not be unrealistically high.should not span more than 1 order of magnitude 
calc.om(rjh_model) # => 1,1
xtab(calc.om(rjh_model),caption="Scenario50_magnitude",cornername="",file="Model_results_Scenario50/Scenario50_magnitude.xls",dec=rep(4,ncol(calc.om(rjh_model))))

#Check if Initial values do not influence the parameter estimates:
#the estimates should be the same for all initial values =>0o

rjh_model_ini <- check.ini(rjh_model, ntrials=30)
rjh_model_ini$check.ini$resmat 

xtab(rjh_model_ini$check.ini$resmat,caption="Scenario50_initial values",cornername="",file="Model_results_Scenario50/Scenario50_initial values.xls",dec=rep(4,ncol(rjh_model_ini$check.ini$resmat)))


## Check residuals ####

#Test for no violation of model assumptions
rjh_model <- calc.osa.resid(rjh_model)
sumspict.diagnostics(rjh_model) # as to be higher than 0.05
plotspict.diagnostic(rjh_model)



#covariance between the model parameters (fixed effects)
#rjh_model$cov.fixed

#correlation between the model parameters (fixed effects)
cov2cor(rjh_model$cov.fixed)
xtab(cov2cor(rjh_model$cov.fixed),caption="Scenario50_correlation between parameters",cornername="",file="Model_results_Scenario50/Scenario50_correlation between parameters.xls",dec=rep(4,ncol(cov2cor(rjh_model$cov.fixed))))


## Retrospective analysis ####

# mohns_rho: if -0.2 < mohns_rho < 0.2 => OK
mohns_rho(retro(rjh_model,nretroyear=3), what=c("BBmsy","FFmsy")) # => OK
# BBmsy       FFmsy 
#0.6899073 -0.1342973

retro_rjh_model <- retro(rjh_model, nretroyear=3)

plotspict.retro(retro_rjh_model)

# Retrospective analysis of fixed effects (with 95% CI)
plotspict.retro.fixed(retro_rjh_model)  

# Time to B_MSY under different scenarios about F
plotspict.tc(retro_rjh_model)



## plot results

tiff(filename = "Model_results_Scenario50/Plots/Fit_results_1.tiff")
plot(rjh_model)
dev.off()


##some plots

tiff(filename = "Model_results_Scenario50/Plots/Biomass.tiff")
plotspict.biomass(rjh_model)
dev.off()

tiff(filename = "Model_results_Scenario50/Plots/FishingM.tiff")
plotspict.f(rjh_model)
dev.off()

tiff(filename = "Model_results_Scenario50/Plots/BBmsy.tiff")
plotspict.bbmsy(rjh_model)
dev.off()

tiff(filename = "Model_results_Scenario50/Plots/FFmsy.tiff")
plotspict.ffmsy(rjh_model)
dev.off()


tiff(filename = "Model_results_Scenario50/Plots/Diagnostics.tiff")
plotspict.diagnostic(rjh_model)
dev.off()

tiff(filename = "Model_results_Scenario50/Plots/Retrospective.tiff")
plotspict.retro(retro_rjh_model)
dev.off()


## EXTRACAO DE TABELAS



tab1 <- sumspict.parest(rjh_model);
xtab(tab1,caption="Parameter estimates",cornername="Parameter",file="Model_results_Scenario50/Parameter_estimates.xls",dec=rep(4,ncol(tab1)))

tab2 <- sumspict.srefpoints(rjh_model);
xtab(tab2,caption="Stochastic reference points",cornername="Reference points",file="Model_results_Scenario50/Reference_points.html",dec=rep(4,ncol(tab2)))

tab3 <- sumspict.states(rjh_model);
xtab(tab3,caption="Estimated states",cornername="",file="Model_results_Scenario50/Estimated_states.html",dec=rep(4,ncol(tab3)))

tab4 <- sumspict.predictions(rjh_model);
xtab(tab4,caption="Forecast",cornername="",file="Model_results_Scenario50/Forecast.html",dec=rep(4,ncol(tab4)))


tab5 <- get.par("logBBmsy",rjh_model,exp=TRUE)
tab5_<- tab5[grep(".",rownames(tab5), fixed=TRUE, invert=TRUE),] 
xtab(tab5_,caption="B/Bmsy",cornername="",file="Model_results_Scenario50/BBmsy.html",dec=rep(4,ncol(tab5_)))
xtab(tab5,caption="B/Bmsy",cornername="",file="Model_results_Scenario50/BBmsy_all.html",dec=rep(4,ncol(tab5)))
#write.csv(tab5_, "BBmsy.csv", dec=",")


tab6 <- get.par("logFFmsy",rjh_model,exp=TRUE)
tab6_<- tab6[grep(".9375",rownames(tab6), fixed=TRUE, invert=FALSE),] 
xtab(tab6,caption="F/Fmsy",cornername="",file="Model_results_Scenario50/FFmsy.html",dec=rep(4,ncol(tab6)))
xtab(tab6_,caption="F/Fmsy",cornername="",file="Model_results_Scenario50/FFmsy_end year.html",dec=rep(4,ncol(tab6_)))

tab6b <- get.par("logF",rjh_model,exp=TRUE)
xtab(tab6b,caption="F",cornername="",file="Model_results_Scenario50/F_results.html",dec=rep(4,ncol(tab6b)))


tab6c <- get.par("logCpred",rjh_model,exp=TRUE)
xtab(tab6c,caption="Catch",cornername="",file="Model_results_Scenario50/Catch.html",dec=rep(4,ncol(tab6c)))



### FORECAST ####


## The management period should be set before running check.inp or fitting the model
## Here the assumption is that the managment starts with one intermediate year
## The default is that the management period starts immediatelly after the last data point

## Time point to evaluate the model states - usually the end of the management period
rjh_data$maninterval <- c(2024, 2025)
rjh_data$maneval <- 2025


fit <- fit.spict(rjh_data)


fit <- add.man.scenario(fit, "F=0", ffac = 0)
fit <- add.man.scenario(fit, "F=Fsq", ffac = 1)
fit <- add.man.scenario(fit, "F=Fmsy", breakpointB = c(1/3, 1/2))
fit <- add.man.scenario(fit, "F=Fmsy_C_fractile", fractiles = list(catch = 0.35), breakpointB = c(1/3, 1/2))
fit <- add.man.scenario(fit, "F=Fmsy_C_fractile_20", fractiles = list(catch = 0.20), breakpointB = c(1/3, 1/2))
fit <- add.man.scenario(fit, "F=Fmsy_C_fractile_10", fractiles = list(catch = 0.10), breakpointB = c(1/3, 1/2))
fit <- add.man.scenario(fit, "F=Fmsy_C_fractile_15", fractiles = list(catch = 0.15), breakpointB = c(1/3, 1/2))

res <- sumspict.manage(fit, include.unc=TRUE,include.abs = TRUE)
res

write.csv(res, "Model_results_Scenario50/Forecast_selected_scenarios.csv")

plotspict.hcr(fit)

plotspict.catch(fit)

fit$man



# hindcast ####


hindcasthc.rjh_model = hindcast(rjh_model,npeels =3)

plotspict.hindcast(hindcasthc.rjh_model,legend.pos = "topleft")

tiff(filename = "Model_results_Scenario50/Plots/hindcast.tiff")
plotspict.hindcast(hindcasthc.rjh_model,legend.pos = "topleft")
dev.off()


#Caspers code from https://community.ices.dk/ExpertGroups/benchmarks/2024/WKBELASMO3/2022%20Meeting%20Documents/07.%20Software/Casper_mrci_function.R

source("Casper_mrci_function.R")

