library(rstan)
library(data.table)
library(lubridate)
library(gdata)
library(EnvStats)
library(tictoc)
tic()
countries <- c( "Germany_2nd") #Introduce the name of the country to study. 
#The dataset must have the form indicated in line 22 (i.e. COVID-19-country.csv)
N2 = 169 # Increase this for a further forecast. Number of days of the simulation.


args = commandArgs(trailingOnly=TRUE)
if(length(args) == 0) {
  args = 'base'
} 
StanModel = args[1]

print(sprintf("Running %s",StanModel))

## Reading all data
if((length(countries) == 1) && (file.exists(paste0('data/COVID-19-',countries,'.csv'))) ){
   d=read.csv(paste0('data/COVID-19-',countries,'.csv'))
}else{
  d=read.csv('data/COVID-19-up-to-date.csv')   
}

## get CFR
cfr.by.country = read.csv("data/weighted_fatality.csv")
cfr.by.country$country = as.character(cfr.by.country[,2])
cfr.by.country$country[cfr.by.country$country == "United Kingdom"] = "United_Kingdom"   

serial.interval = read.csv("data/serial_interval.csv")
SI=serial.interval$fit
n_SI = length(SI)
if (n_SI < N2){
  SI[(n_SI+1):N2] = rep(0.,(N2-n_SI))
}

covariates = read.csv('data/int_ger2.csv', stringsAsFactors = FALSE) #name of the interventions file


covariates <- covariates[1:length(covariates$Country), c(1,2,3,4,5,6, 7, 8)]  
inf.to.death = read.csv("data/infection_to_death.csv")  
ITD=inf.to.death$fit
n_ITD = length(ITD)
if (n_ITD < N2){
  ITD[(n_ITD+1):N2] = rep(0.,(N2-n_ITD))
}


p <- ncol(covariates) - 1
forecast = 0


RMSD_deaths <- c()
RMSD_cumdeaths <- c()


for(Country in countries) {
  #for ( q in 1:2){
  for (q in 1:length(covariates$Country)){
    
    dates = list()
    reported_cases = list()
    stan_data = list(M=length(countries),N=NULL,p=p,x1=poly(1:N2,2)[,1],x2=poly(1:N2,2)[,2],
                     y=NULL,covariate1=NULL,covariate2=NULL,covariate3=NULL,covariate4=NULL,covariate5=NULL,covariate6=NULL,covariate7=NULL,deaths=NULL,f=NULL,
                     N0=6,cases=NULL,LENGTHSCALE=7,SI=SI,
                     EpidemicStart = NULL) 
    deaths_by_country = list()
    
    CFR=cfr.by.country$weighted_fatality[cfr.by.country$Region..subregion..country.or.area.. == Country]
    
    d1=d[d$Countries.and.territories==Country,]
    d1$date = as.Date(d1$DateRep,format='%d/%m/%Y')
    d1$t = decimal_date(d1$date) 
    d1=d1[order(d1$t),]
    index = which(d1$Cases>0)[1]
    n_deaths_0 = 10  # defines when to start the analysis; default 10, but 1 seems better for small-size countries of regions
    #index1 = which(cumsum(d1$Deaths)>=n_deaths_0)[1] 
    #index2 = index1-30   
    index1 = 10
    index2 = 0
    
    print(sprintf("First non-zero cases is on day %d, and 30 days before %d deaths is day %d",index,n_deaths_0, index2))   
    d1=d1[index2:nrow(d1),]
    if(length(countries)==1) {
      stan_data$EpidemicStart = as.array(index1+1-index2)
    } else {
      stan_data$EpidemicStart = c(stan_data$EpidemicStart,index1+1-index2)
    }
    
    dates[[Country]] = d1$date
    # hazard estimation
    N = length(d1$Cases)
    print(sprintf("%s has %d days of data",Country,N))
    forecast = N2 - N
    if(forecast < 0) {
      print(sprintf("%s: %d", Country, N))
      print("ERROR!!!! increasing N2")
      N2 = N
      forecast = N2 - N
    }
    
    
    times <- d1$date   
    times_forecast <- times[N] + 1:forecast  
    times[(N+1):N2] <- times_forecast[1:forecast]
    
    
    f=CFR * ITD   #  infection to death probability read from a file
    
    
    
    
    
    y=c(as.vector(as.numeric(d1$Cases)),rep(-1,forecast))
    reported_cases[[Country]] = as.vector(as.numeric(d1$Cases))
    deaths=c(as.vector(as.numeric(d1$Deaths)),rep(-1,forecast))
    cases=c(as.vector(as.numeric(d1$Cases)),rep(-1,forecast))
    deaths_by_country[[Country]] = as.vector(as.numeric(d1$Deaths))
    
    cum_cases = cumsum(cases)
    cum_cases[(N+1):N2] = rep(cum_cases[N], N2-N)
    cum_deaths = cumsum(deaths)
    cum_deaths[(N+1):N2] = rep(cum_deaths[N], N2-N)
    
    covariates1 <- covariates[q, 2:8]
    #Names of the columns for each intervention (start and finish)
    covariate1 = (times >= as.Date(covariates1[1,"other_interv1"]))*1
    covariate2 = (times >= as.Date(covariates1[1,"other_interv1_off"]))*(-1)
    covariate1 = covariate1 + covariate2
    
    covariate3 = (times >= as.Date(covariates1[1,"major_intervention"]))*1
    covariate6 = (times >= as.Date(covariates1[1,"major_intervention_off"]))*(-1)
    covariate3 = covariate3 + covariate6
    
    covariate4 = (times >= as.Date(covariates1[1,"other_interv2"]))*1
    covariate5 = (times >= as.Date(covariates1[1,"other_interv2_off"]))*(-1)
    covariate4 = covariate4 + covariate5
    
    covariate7 = (times >= as.Date(covariates1[1,"interv_dummy"]))*1
    
    
    d1["other_interv1"] <- covariate1[1:N]
    d1["other_interv1_off"]   <- covariate2[1:N]
    d1["major_intervention"]  <- covariate3[1:N]
    d1["other_interv2"]  <- covariate4[1:N]
    d1["other_interv2_off"]  <- covariate5[1:N]
    d1["major_intervention_off"]   <- covariate6[1:N]
    d1["interv_dummy"]   <- covariate7[1:N]
    
    
    ## append data
    stan_data$N = c(stan_data$N,N)
    stan_data$y = c(stan_data$y,y[1]) # just the index case!
    stan_data$covariate1 = cbind(stan_data$covariate1,covariate1)
    stan_data$covariate2 = cbind(stan_data$covariate2,covariate2)
    stan_data$covariate3 = cbind(stan_data$covariate3,covariate3)
    stan_data$covariate4 = cbind(stan_data$covariate4,covariate4)
    stan_data$covariate5 = cbind(stan_data$covariate5,covariate5)
    stan_data$covariate6 = cbind(stan_data$covariate6,covariate6)
    stan_data$covariate7 = cbind(stan_data$covariate7,covariate7) 
    stan_data$f = cbind(stan_data$f,f)
    stan_data$deaths = cbind(stan_data$deaths,deaths)
    stan_data$cases = cbind(stan_data$cases,cases)
    stan_data$inf_to_death = cbind(stan_data$inf_to_death,ITD)   
    stan_data$CFR = CFR  
    
    
    stan_data$cum_cases = cbind(stan_data$cum_cases,cum_cases)
    stan_data$cum_deaths = cbind(stan_data$cum_deaths,cum_deaths)
    
    stan_data$N2=N2
    stan_data$x=1:N2
    if(length(stan_data$N) == 1) {
      stan_data$N = as.array(stan_data$N)
    }
    
    stan_data$covariate7 = 0 # models should only take 6 covariates
    
    
    
    stan_data$y = t(stan_data$y)
    options(mc.cores = parallel::detectCores())
    rstan_options(auto_write = TRUE)
    m = stan_model(paste0('stan-models/',StanModel,'.stan'))
    
    
    
    #fit = sampling(m,data=stan_data,iter=4000,warmup=2000,chains=8,thin=4,control = list(adapt_delta = 0.90, max_treedepth = 10))
    fit = sampling(m,data=stan_data,iter=200,warmup=100,chains=4,thin=4,control = list(adapt_delta = 0.90, max_treedepth = 10))   # for speed
    
    
    
    out = rstan::extract(fit)
    prediction = out$prediction
    estimated.deaths = out$E_deaths
    estimated.deaths.cf = out$E_deaths0
    
    RMSD_deaths[q] <- sqrt(1/N2*(sum(((apply(estimated.deaths, FUN = median, MARGIN = 2)-deaths)[1:nrow(d1)])^2))) #mejor que mediana para no coger valores extremos
    RMSD_cumdeaths[q] <- sqrt(1/N2*(sum(cumsum(((apply(estimated.deaths, FUN = median, MARGIN = 2)))-cum_deaths)[1:nrow(d1)]^2)))
    
    JOBID = Sys.getenv("PBS_JOBID")
    if(JOBID == "")
      JOBID = as.character(abs(round(rnorm(1) * 1000000)))
    print(sprintf("Jobid = %s",JOBID))
    
    filename <- paste0('results/',countries,'cov-',q) #We save the results file wuth the name of the country and the index of the covariates
    save.image(paste0(filename, ".Rdata"))
    
    save(fit,prediction,dates,reported_cases,deaths_by_country,countries,estimated.deaths,estimated.deaths.cf,out,covariates, RMSD_deaths, RMSD_cumdeaths, file=paste0(filename,'-stanfit.Rdata'))
    
    # to visualize results
    #system(paste0("Rscript plot-3-panel.r ", filename,'.Rdata'))
    #source("pred.r")
    toc()
  }
}
save(RMSD_deaths, RMSD_cumdeaths, file=paste0('results/',countries,'/RMSD_models.Rdata'))

