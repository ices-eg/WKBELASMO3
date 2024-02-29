library(spict)
# rjc <- read.csv("rjh.27.9a_input.csv", header=T, sep=",", dec=".")
# inp <- list()

# inp$obsC  <- rjc$Landings
# inp$timeC <- rjc$YEAR        
# inp$obsI  <- rjc$PT_LPUE
# inp$timeI <- rjc$YEAR+ 0.5

# rjc_data$stdevfacC <- c(rep(2,8), rep(1,15)) #capturas: incerteza para o periodo 2000-2007


# inp$stdevfacC <- c(rep(2,20), rep(1,13))
# inp$stdevfacI <- c(rep(2,6), rep(1,3), 2, rep(1,3), rep(2,2), rep(1,18))
inp<-rjh_data #Put you input data here

inp <- check.inp(inp)
plotspict.data(inp, qlegend = TRUE)

#inp$ini$logn <- log(2); inp$phases$logn <- -1  # Fixing n to resemble the Schaefer production model

#inp$priors$logr      <- c(log(0.27),0.5,1)   # intrinsic biomass growth (obtained from jbleslie function (R package JABBA))

# inp$priors$logsdb     <- c(log(0.15), 0.5, 1)  # prior value provided by Tobias (pers. comm) in WKLIFE XII
# 
# inp$priors$logalpha <- c(1, 1, 0)         #disabled
# inp$priors$logbeta  <- c(1, 1, 0)         #disabled

# fit9t <- fit.spict(inp)  
# summary(fit9t)


plot(fit)

retr <- retro(fit)

plotspict.retro(retr)

## Function to calculate CI for mohn's rho
mrci<-function(spictfit,n, npeels=3,mc.cores=1,seed=1,probs=c(0.025,0.5,0.975)){

    do.one<-function(i){
        print(i)
        set.seed(seed+i)
        sim <- sim.spict(spictfit,nobs=1) ## obs, nobs is not used
        cat("refitting\n")
        tmpfit <- fit.spict(sim)
        if(tmpfit$opt$convergence==0){
            cat("running retro\n")
            try(retr <- retro(tmpfit,npeels,mc.cores=mc.cores))
            mr <- c(NA,NA)
            try(mr <- mohns_rho(retr))
            
            return(mr)
        }
        else return( c(NA,NA))
    }
    res = lapply(1:n,do.one)
    resmat = do.call(rbind,res)
    
    apply(resmat,2,FUN=quantile,probs=probs,na.rm=TRUE)

}

test = mrci(fit,n=400,mc.cores=1)

test
