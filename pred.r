library(rstan)
library(data.table)
library(lubridate)
library(gdata)
library(EnvStats)
library(matrixStats)

Country0 = "Germany_2nd"
if(length(Country0) > 1) {
  print("ERROR!!!! this script is valid only for 1 country at a time")
}
for (pp in 1:1 ){
  s_resultsFile = "none"
  s_resultsFile = paste0("results/Germany_2nd/best_pred501.Rdata") 
  pr = 1# Spain, conf 14M, no contention no worklockdown,5May20 CFR=0.009262, 4000 iter
  #    s_resultsFile = "results/base-722907.Rdata"   # Rioja
  #    s_resultsFile = "results/base-879695.Rdata"   # RiojaNORES, excluding data from residences with cases
  #    s_resultsFile = "results/base-288586.Rdata"   # Iceland
  #    s_resultsFile = "results/base-358177.Rdata"   # Germany
  #    s_resultsFile = "results/base-1945735.Rdata"  # France
  #    s_resultsFile = "results/base-384766.Rdata"   # Italy
  #    s_resultsFile = "results/base-1442193.Rdata"  .# United_Kingdomº
  
  
  if (file.exists(paste0(s_resultsFile))){   # first will try to load data from previous results file
    load(s_resultsFile)

   # if (any(countries == Country0)){
    #  m = which(countries==Country0)
    #}else{  
    #  print("ERROR!!! Input country not found in results file")
    #}
    m <- 1
    set.seed(as.numeric(Sys.time()))
    I=length( prediction[,1,m] )
    predictionI = matrix(0, I, N2)
    E_deathsI = matrix(0, I, N2)
    RtI = matrix(0,I, N2)
    
    predictionI[,1:N2] = out$prediction[,1:N2,m]
    E_deathsI[,1:N2] = out$E_deaths[,1:N2,m]
    RtI[,1:N2] = out$Rt[,1:N2,m]
    
    
    active = list()
    if (any(names(d1) == "total.active")){
      active[[Country0]] = as.vector(as.numeric(d1$total.active))    # only if this data exists
    }else{
      active[[Country0]] = rep(0.,N)     # empty vector, just to plot predicted active
    }
    
    
    
    
    
    
  }else{     # if not results file, predicted data will be generated
    
    
    N2 = 300   # number of days to analyze
    
    ##### set of parameters: number of cases first 5 days, R0, Rt reduction after intervention0, idem. intervention1, idem. intervention2, n deaths to define starting time, CFR
    params = c(18.3501510,2.5188597,1., 0.1635222,1.,1,0.009262 )  # RiojaNORES, excl residences, from base-879695.Rdata, median E_deaths set
    params = c(11.42053, 1.833035,1.,0.129883, 1.,1, 0.00556 )  # Iceland, fitting deaths and accumulated deaths, from base-288586.Rdata, median E_deaths set
    
    Country = Country0
    ## Reading all data
    # if COVID-19-"Country".csv exits, read it... if not, read the general file
    if (file.exists(paste0('data/COVID-19-',Country,'.csv'))){
      d=read.csv(paste0('data/COVID-19-',Country,'.csv'))
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
    covariates = read.csv('data/int_ger2.csv', stringsAsFactors = FALSE)
    covariates <- covariates[1:length(covariates$Country), c(1,2,3,4,5,6, 7, 8)]   
    inf.to.death = read.csv("data/infection_to_death.csv")  
    ITD=inf.to.death$fit
    n_ITD = length(ITD)
    if (n_ITD < N2){
      ITD[(n_ITD+1):N2] = rep(0.,(N2-n_ITD))
    }  
    
    
    p <- ncol(covariates) - 1
    
    dates = list()
    reported_cases = list()
    active = list()
    stan_data = list(M=length(Country),N=NULL,p=p,x1=poly(1:N2,2)[,1],x2=poly(1:N2,2)[,2],
                     y=NULL,covariate1=NULL,covariate2=NULL,covariate3=NULL,covariate4=NULL,covariate5=NULL,covariate6=NULL,covariate7=NULL,deaths=NULL,f=NULL,
                     N0=6,cases=NULL,LENGTHSCALE=7,
                     EpidemicStart = NULL) 
    deaths_by_country = list()
    
    
    CFR=params[7]   
    population = 1000.*cfr.by.country$population[cfr.by.country$Region..subregion..country.or.area.. == Country0]
    
    covariates1 <- covariates[pr, 2:8]
    
    d1=d[d$Countries.and.territories==Country,]
    d1$date = as.Date(d1$DateRep,format='%d/%m/%Y')
    d1$t = decimal_date(d1$date) 
    d1=d1[order(d1$t),]
    index = which(d1$Cases>0)[1]
    n_deaths_0 = params[6]
    index1 = which(cumsum(d1$Deaths)>=n_deaths_0)[1] 
    index2 = index1-30
    
    print(sprintf("First non-zero cases is on day %d, and 30 days before %d deaths is day %d",index,n_deaths_0, index2))   
    d1=d1[index2:nrow(d1),]
    stan_data$EpidemicStart = c(stan_data$EpidemicStart,index1+1-index2)
    
    
    
    
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
    
    
    
    
    
    f=CFR * ITD[1:N2]   # infection to death distribution read from a file
    
    
    
    
    y=c(as.vector(as.numeric(d1$Cases)),rep(-1,forecast))
    
    reported_cases[[Country]] = as.vector(as.numeric(d1$Cases))
    deaths=c(as.vector(as.numeric(d1$Deaths)),rep(-1,forecast))
    cases=c(as.vector(as.numeric(d1$Cases)),rep(-1,forecast))
    deaths_by_country[[Country]] = as.vector(as.numeric(d1$Deaths))
    
    if (any(names(d1) == "total.active")){
      active[[Country]] = as.vector(as.numeric(d1$total.active))    # only if this data exists
    }else{
      active[[Country]] = rep(0.,length(d1$Cases))     # empty vector, just to plot predicted active
    }
    
    
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
    stan_data$deaths = cbind(stan_data$deaths,deaths)
    stan_data$cases = cbind(stan_data$cases,cases)
    
    
    stan_data$N2=N2
    stan_data$x=1:N2
    if(length(stan_data$N) == 1) {
      stan_data$N = as.array(stan_data$N)
    }
    
    
    stan_data$covariate7 = 0 # models should only take 6 covariates
    
    
    
    stan_data$y = t(stan_data$y)
    
    
    
    # predictive model 
    
    I=100 # number of iterations   
    predictionI = matrix(0, I, stan_data$N2)
    E_deathsI = matrix(0, I, stan_data$N2)
    RtI = matrix(0,I, stan_data$N2)
    
    for(i_iter in (1:I)){
      
      prediction = matrix(0,stan_data$N2,stan_data$M)
      E_deaths  = matrix(0,stan_data$N2,stan_data$M)
      Rt = matrix(0,stan_data$N2,stan_data$M)
      alpha = rep(0,6)
      predictionN0=params[1]   # number of cases in the first N0 days  
      Rt1= params[2]  # initial Rt   
      alpha[1] = log( params[3] )  # effect of intervention 1   
      alpha[3] = log( params[4] )  # effect of major intervention   
      alpha[4] = log( params[5] )  # effect of intervention 2   
      n_discrete = 1000000 # when number of cases < n_discrete, a random discrete distribution will be generated from the probability curve
      if(I!=1 && n_discrete ==0){    # add noise when using several iterations
        r_noise = 0.01 
        predictionN0 = rnorm(1, params[1], sd=abs(r_noise*params[1]))
        if (predictionN0 < 1){predictionN0=1}
        Rt1= rnorm(1,params[2],sd=abs(r_noise*params[2]))     
        if (Rt1 < 0.){Rt1=0.}
        alpha[1] = rnorm( 1, log(params[3]), sd=abs(r_noise*log(params[3])) )  
        alpha[3] = rnorm( 1, log(params[4]), sd=abs(r_noise*log(params[4])) )    
        alpha[4] = rnorm( 1, log(params[5]), sd=abs(r_noise*log(params[5])) )    
      }   
      for (m in 1:stan_data$M){
        prediction[1:stan_data$N0,m] = rep(predictionN0,stan_data$N0)  
        Rt[,m] = Rt1 * exp(stan_data$covariate1[,m] * alpha[1] +    
                             stan_data$covariate3[,m] * alpha[3] + stan_data$covariate4[,m] * alpha[4] )
        random_infected = matrix(0,N2,N2)  
        for (i in (1:stan_data$N0)){    
          if(prediction[i,m] < n_discrete){		  
            random_infected[i,] = table(cut(sample(N2,prediction[i,m],prob=SI[1:N2],replace=T), breaks=c(seq(0,N2,1)) ))   # discrete infected precomputed 
          }
        }                  
        maxcumpred = sum(prediction[1:stan_data$N0,m])	  
        for (i in (stan_data$N0+1):stan_data$N2) {
          convolution=0
          for(j in 1:(i-1)) {
            if(prediction[j,m] < n_discrete){
              convolution = convolution + random_infected[j,i-j]   # discrete infected precomputed 
            }else{		  
              convolution = convolution + prediction[j, m]*SI[i-j]   # continuous predicted infected distribution 
            }
          }
          prediction[i, m] = Rt[i,m] * convolution
          if ((maxcumpred + prediction[i,m] ) > population){
            prediction[i,m] = population - maxcumpred 
          }
          maxcumpred = maxcumpred + prediction[i,m]
          if(prediction[i,m] < n_discrete){
            random_infected[i,] = table(cut(sample(N2,prediction[i,m],prob=SI[1:N2],replace=T), breaks=c(seq(0,N2,1)) ))   # discrete infected precomputed 
          }
        }
        
        E_deaths[1, m]= 1e-9
        random_deaths = matrix(0,N2,N2)  # discrete death precomputed 
        if(CFR*prediction[1,m] < n_discrete){     
          random_deaths[1,] = table(cut(sample(N2,round(CFR*prediction[1,m]),prob=ITD[1:N2],replace=T), breaks=c(seq(0,N2,1)) ))   # discrete death precomputed
        }
        for (i in 2:stan_data$N2){
          E_deaths[i,m]= 0
          for(j in 1:(i-1)){
            if(CFR*prediction[j,m] < n_discrete){
              E_deaths[i,m] = E_deaths[i,m] + random_deaths[j,i-j]   # discrete death precomputed 
            }else{	  
              E_deaths[i,m] = E_deaths[i,m] + CFR*prediction[j,m]*ITD[i-j]   # continuous estimated death distribution 
            }
          }
          if(CFR*prediction[i,m] < n_discrete){
            random_deaths[i,] = table(cut(sample(N2,round(CFR*prediction[i,m]),prob=ITD[1:N2],replace=T), breaks=c(seq(0,N2,1)) ))   # discrete death precomputed
          }
        }
      }
      
      
      
      predictionI[i_iter,] = prediction[,m]
      E_deathsI[i_iter,] = E_deaths[,m]
      RtI[i_iter,] = Rt[,m]
      
      
    }
    
    
    
  } 
  
  
  
  
  # infection-to-recovery curve
  h_ITR=rep(0.,100)
  
  mean1_ITR = 5.1; cv1_ITR = 0.86; # infection to onset
  mean2_ITR = 24.7  ; cv2_ITR = 0.35  # onset to recovery
  x1_ITR = rgammaAlt(5e6,mean1_ITR,cv1_ITR) # infection-to-onset 
  x2_ITR = rgammaAlt(5e6,mean2_ITR,cv2_ITR) # onset-to-recovery
  fc_ITR = ecdf(x1_ITR+x2_ITR)
  convolution_ITR = function(u) (fc_ITR(u))
  
  h_ITR[1] = (convolution_ITR(1.5) - convolution_ITR(0))
  for(i_ITR in 2:length(h_ITR)) {
    h_ITR[i_ITR] = (convolution_ITR(i_ITR+.5) - convolution_ITR(i_ITR-.5)) / (1-convolution_ITR(i_ITR-.5))
  }
  s_ITR = rep(0,100)
  s_ITR[1] = 1
  for(i_ITR in 2:100) {
    s_ITR[i_ITR] = s_ITR[i_ITR-1]*(1-h_ITR[i_ITR-1])
  }
  ITR =  s_ITR * h_ITR
  
  
  if (N2 <= 100){
    ITR = ITR[1:N2]
  }else{
    ITR[101:N2] = rep(0.,length(101:N2))
  }
  
  
  
  
  predactiveI = matrix(0, I, N2)
  predrecoveredI = matrix(0, I, N2)
  cumpredrecoveredI = matrix(0, I, N2)
  n_discrete = 0 # when number of cases < n_discrete, a random discrete distribution will be generated from the probability curve
  for(i_iter in (1:I)){
    
    predrecoveredI[i_iter,1] = 1e-9
    random_recov = matrix(0,N2,N2)  # discrete recovered precomputed 
    if((1-CFR)*predictionI[i_iter,1] < n_discrete){
      random_recov[1,] = table(cut(sample(N2,round((1-CFR)*predictionI[i_iter,1]),prob=ITR[1:N2],replace=T), breaks=c(seq(0,N2,1)) ))   # discrete recovered precomputed 
    }
    for (i_recov in 2:N2){
      predrecoveredI[i_iter,i_recov] = 0
      for (j_recov in 1:(i_recov-1)){
        if((1-CFR)*predictionI[i_iter,j_recov] < n_discrete){
          predrecoveredI[i_iter,i_recov] = predrecoveredI[i_iter,i_recov] + random_recov[j_recov,i_recov-j_recov]   # discrete recovered precomputed 
        }else{
          predrecoveredI[i_iter,i_recov] = predrecoveredI[i_iter,i_recov] + (1-CFR)*predictionI[i_iter,j_recov]*ITR[i_recov-j_recov]   # continuous estimated recovered distribution
        }
      }
      if((1-CFR)*predictionI[i_iter,i_recov] < n_discrete){
        random_recov[i_recov,] = table(cut(sample(N2,round((1-CFR)*predictionI[i_iter,i_recov]),prob=ITR[1:N2],replace=T), breaks=c(seq(0,N2,1)) ))   # discrete death precomputed 
      }
    }
    
    
    cumpredrecoveredI[i_iter,] = cumsum(predrecoveredI[i_iter,])
  }
  
  
  
  
  
  cumpreddeathsI = matrix(0, I, N2)
  for (i_acum in (1:length(E_deathsI[,1]))){
    cumpreddeathsI[i_acum,] = cumsum(E_deathsI[i_acum,1:N2])
  }
  
  
  
  
  
  cumpredictionI = matrix(0, I, N2)
  for (i_acum in (1:length(predictionI[,1]))){
    cumpredictionI[i_acum,] = cumsum(predictionI[i_acum,1:N2])
    
    predactiveI[i_acum,] = cumpredictionI[i_acum,] - cumpredrecoveredI[i_acum,] - cumpreddeathsI[i_acum,]
  }
  
  
  
  
  
  
  if(I==1){
    
    predicted_cases <- predictionI[,1:N2]
    predicted_cases_li <- predictionI[,1:N2]
    predicted_cases_ui <- predictionI[,1:N2]
    predicted_cases_li2 <-  predictionI[,1:N2]
    predicted_cases_ui2 <- predictionI[,1:N2]
    
    estimated_deaths  <- E_deathsI[,1:N2]
    estimated_deaths_li  <- E_deathsI[,1:N2]
    estimated_deaths_ui   <- E_deathsI[,1:N2]
    estimated_deaths_li2  <- E_deathsI[,1:N2]
    estimated_deaths_ui2  <- E_deathsI[,1:N2]
    
    rt <- RtI[,1:N2]
    rt_li <- RtI[,1:N2] 
    rt_ui <-  RtI[,1:N2]
    rt_li2 <-  RtI[,1:N2]
    rt_ui2 <-  RtI[,1:N2]
    
    predactive <- predactiveI[,1:N2]
    predactive_li <- predactiveI[,1:N2]
    predactive_ui <- predactiveI[,1:N2]
    predactive_li2 <-  predactiveI[,1:N2]
    predactive_ui2 <- predactiveI[,1:N2]
    
    cumpredicted_cases <- cumpredictionI[,1:N2]
    cumpredicted_cases_li <- cumpredictionI[,1:N2]
    cumpredicted_cases_ui <- cumpredictionI[,1:N2]
    cumpredicted_cases_li2 <-  cumpredictionI[,1:N2]
    cumpredicted_cases_ui2 <- cumpredictionI[,1:N2]
    
    cumpredicted_deaths <- cumpreddeathsI[,1:N2]
    cumpredicted_deaths_li <- cumpreddeathsI[,1:N2]
    cumpredicted_deaths_ui <- cumpreddeathsI[,1:N2]
    cumpredicted_deaths_li2 <-  cumpreddeathsI[,1:N2]
    cumpredicted_deaths_ui2 <- cumpreddeathsI[,1:N2]
    
    
    
  }else{	
    
    singleI = FALSE # to choose a single set of parameters from the fitted model, use FALSE for taking all
    seleMODE = "MEDIAN"  # "RANDOM"  to choose a random set of parameters from the fitted model; "MEDIAN" to choose the set that provdes the median E_deaths
    if (file.exists(paste0(s_resultsFile)) && singleI==TRUE){
      if (seleMODE == "RANDOM") { i_sele = round(runif(1, min=1, max=N2), digits=0) }
      if (seleMODE == "MEDIAN") {
        i_sele=1
        i_max = 9999999
        median_deaths <-  colQuantiles(E_deathsI[,1:N2], probs=.5)
        for (ii_sele in 1:length(E_deathsI[,1]) ) {
          if(   sum( abs( E_deathsI[ii_sele,1:N2] - median_deaths )  )  < i_max){
            i_max = sum( abs( E_deathsI[ii_sele,1:N2] - median_deaths ))
            i_sele = ii_sele
          }
        }
      }
      print("selected set of parameters")
      print(i_sele)
      predicted_cases <- predictionI[i_sele,1:N2]
      predicted_cases_li <- predictionI[i_sele,1:N2]
      predicted_cases_ui <- predictionI[i_sele,1:N2]
      predicted_cases_li2 <-  predictionI[i_sele,1:N2]
      predicted_cases_ui2 <- predictionI[i_sele,1:N2]
      
      estimated_deaths  <- E_deathsI[i_sele,1:N2]
      estimated_deaths_li  <- E_deathsI[i_sele,1:N2]
      estimated_deaths_ui   <- E_deathsI[i_sele,1:N2]
      estimated_deaths_li2  <- E_deathsI[i_sele,1:N2]
      estimated_deaths_ui2  <- E_deathsI[i_sele,1:N2]
      
      rt <- RtI[i_sele,1:N2]
      rt_li <- RtI[i_sele,1:N2]
      rt_ui <-  RtI[i_sele,1:N2]
      rt_li2 <-  RtI[i_sele,1:N2]
      rt_ui2 <-  RtI[i_sele,1:N2]
      
      predactive <- predactiveI[i_sele,1:N2]
      predactive_li <- predactiveI[i_sele,1:N2]
      predactive_ui <- predactiveI[i_sele,1:N2]
      predactive_li2 <-  predactiveI[i_sele,1:N2]
      predactive_ui2 <- predactiveI[i_sele,1:N2]
      
      cumpredicted_cases <- cumpredictionI[i_sele,1:N2]
      cumpredicted_cases_li <- cumpredictionI[i_sele,1:N2]
      cumpredicted_cases_ui <- cumpredictionI[i_sele,1:N2]
      cumpredicted_cases_li2 <-  cumpredictionI[i_sele,1:N2]
      cumpredicted_cases_ui2 <- cumpredictionI[i_sele,1:N2]
      
      cumpredicted_deaths <- cumpreddeathsI[i_sele,1:N2]
      cumpredicted_deaths_li <- cumpreddeathsI[i_sele,1:N2]
      cumpredicted_deaths_ui <- cumpreddeathsI[i_sele,1:N2]
      cumpredicted_deaths_li2 <-  cumpreddeathsI[i_sele,1:N2]
      cumpredicted_deaths_ui2 <- cumpreddeathsI[i_sele,1:N2]
      
      
      
    }else{
      
      
      # prediction for n iterations
      predicted_cases <- colQuantiles(predictionI[,1:N2], probs=.5)             
      predicted_cases_li <- colQuantiles(predictionI[,1:N2], probs=.025)
      predicted_cases_ui <- colQuantiles(predictionI[,1:N2], probs=.975)
      predicted_cases_li2 <- colQuantiles(predictionI[,1:N2], probs=.25)
      predicted_cases_ui2 <- colQuantiles(predictionI[,1:N2], probs=.75)
      
      estimated_deaths <-  colQuantiles(E_deathsI[,1:N2], probs=.5)       
      estimated_deaths_li <- colQuantiles(E_deathsI[,1:N2], probs=.025)
      estimated_deaths_ui <- colQuantiles(E_deathsI[,1:N2], probs=.975)
      estimated_deaths_li2 <- colQuantiles(E_deathsI[,1:N2], probs=.25)
      estimated_deaths_ui2 <- colQuantiles(E_deathsI[,1:N2], probs=.75)
      
      rt <- colQuantiles(RtI[,1:N2],probs=.5)    
      rt_li <- colQuantiles(RtI[,1:N2],probs=.025)
      rt_ui <- colQuantiles(RtI[,1:N2],probs=.975)
      rt_li2 <- colQuantiles(RtI[,1:N2],probs=.25)
      rt_ui2 <- colQuantiles(RtI[,1:N2],probs=.75)
      
      predactive <- colQuantiles(predactiveI[,1:N2], probs=.5)   
      predactive_li <- colQuantiles(predactiveI[,1:N2], probs=.025)
      predactive_ui <- colQuantiles(predactiveI[,1:N2], probs=.975)
      predactive_li2 <- colQuantiles(predactiveI[,1:N2], probs=.25)
      predactive_ui2 <- colQuantiles(predactiveI[,1:N2], probs=.75)
      
      cumpredicted_cases <- colQuantiles(cumpredictionI[,1:N2], probs=.5)   
      cumpredicted_cases_li <- colQuantiles(cumpredictionI[,1:N2], probs=.025)
      cumpredicted_cases_ui <- colQuantiles(cumpredictionI[,1:N2], probs=.975)
      cumpredicted_cases_li2 <- colQuantiles(cumpredictionI[,1:N2], probs=.25)
      cumpredicted_cases_ui2 <- colQuantiles(cumpredictionI[,1:N2], probs=.75)
      
      cumpredicted_deaths <- colQuantiles(cumpreddeathsI[,1:N2], probs=.5)
      cumpredicted_deaths_li <- colQuantiles(cumpreddeathsI[,1:N2], probs=.025)
      cumpredicted_deaths_ui <- colQuantiles(cumpreddeathsI[,1:N2], probs=.975)
      cumpredicted_deaths_li2 <-  colQuantiles(cumpreddeathsI[,1:N2], probs=.25)
      cumpredicted_deaths_ui2 <- colQuantiles(cumpreddeathsI[,1:N2], probs=.75)
      
      
    }
  }
  
  
  
  # to visualize results
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(rstan)
  library(data.table)
  library(lubridate)
  library(gdata)
  library(bayesplot)
  library(EnvStats)
  library(matrixStats)
  library(scales)
  library(gridExtra)
  library(ggpubr)
  library(cowplot)
  source("utils/geom-stepribbon.r")
  
  
  
  i=1  # this script is valid only for 1 country at a time
  
  # to read result files with old intervention names
  
  if(any(names(covariates) == "confinement")){ 
    
    covariates$major_intervention <- covariates$confinement
    covariates$major_intervention_off <- covariates$confinement_off
    covariates$other_interv1 <- covariates$contention
    covariates$other_interv1_off <- covariates$contention_off
    covariates$other_interv2 <- covariates$work_lockdown
    covariates$other_interv2_off <- covariates$work_lockdown_off
    covariates$interv_dummy <- covariates$school_universities
    
    covariates <- covariates[5,c(1,9:15)]
    
  }
  
  covariates_country <- covariates[, 2:8]
  
  # Remove covariates not used
  covariates_country$interv_dummy = NULL
  covariates_country_long <- gather(covariates_country[], key = "key",
                                    value = "value")
  covariates_country_long$x <- rep(NULL, length(covariates_country_long$key))
  un_dates <- unique(covariates_country_long$value)
  
  for (k in 1:length(un_dates)){
    idxs <- which(covariates_country_long$value == un_dates[k])
    max_val <- round(max(rt_ui)) + 0.3      
    for (j in idxs){
      covariates_country_long$x[j] <- max_val
      max_val <- max_val - 0.3
    }
  }
  
  
  covariates_country_long$value <- as_date(covariates_country_long$value)
  covariates_country_long$country <- rep(Country,
                                         length(covariates_country_long$value))
  
  
  
  
  
  
  
  
  reported_cases[[i]][(N+1):N2] <- rep(0,forecast)   
  deaths_by_country[[i]][(N+1):N2] <- rep(0,forecast)   
  active[[i]][(N+1):N2] <- rep(0,forecast)    # if this exists
  data_country <- data.frame("time" = times,
                             "country" = rep(Country, N2),  
                             "reported_cases" = reported_cases[[i]],
                             "reported_cases_c" = cumsum(reported_cases[[i]]),
                             "predicted_cases_c" = cumpredicted_cases,
                             "predicted_min_c" = cumpredicted_cases_li,   
                             "predicted_max_c" = cumpredicted_cases_ui,
                             "predicted_cases" = predicted_cases,
                             "predicted_min" = predicted_cases_li,
                             "predicted_max" = predicted_cases_ui,
                             "predicted_min2" = predicted_cases_li2,
                             "predicted_max2" = predicted_cases_ui2,
                             "predicted_min2_c" = cumpredicted_cases_li2,
                             "predicted_max2_c" = cumpredicted_cases_ui2,
                             "deaths" = deaths_by_country[[i]],
                             "deaths_c" = cumsum(deaths_by_country[[i]]),
                             "estimated_deaths_c" =  cumpredicted_deaths,
                             "death_min_c" = cumpredicted_deaths_li,
                             "death_max_c"= cumpredicted_deaths_ui,
                             "death_min2_c" = cumpredicted_deaths_li2,
                             "death_max2_c" = cumpredicted_deaths_ui2,
                             "estimated_deaths" = estimated_deaths,
                             "death_min" = estimated_deaths_li,
                             "death_max"= estimated_deaths_ui,
                             "death_min2" = estimated_deaths_li2,
                             "death_max2"= estimated_deaths_ui2,
                             "rt" = rt,
                             "rt_min" = rt_li,
                             "rt_max" = rt_ui,
                             "rt_min2" = rt_li2,
                             "rt_max2" = rt_ui2,
                             "active" = active[[i]],
                             "predactive" = predactive,
                             "predactive_min" = predactive_li,
                             "predactive_max" = predactive_ui,
                             "predactive_min2" = predactive_li2,
                             "predactive_max2" = predactive_ui2
  )
  
  
  
  data_country$reported_cases_c[(N+1):N2] <- rep(0,forecast) 
  data_country$deaths_c[(N+1):N2] <- rep(0,forecast)
  
  
  
  
  JOBID = Sys.getenv("PBS_JOBID")
  if(JOBID == "")
    JOBID = as.character(abs(round(rnorm(1) * 1000000)))
  print(sprintf("Jobid = %s",JOBID))
  
  save.image(paste0('results/pred/','pred2-',JOBID,'.Rdata'))
  
  save(prediction,dates,reported_cases,deaths_by_country,Country,estimated_deaths,covariates,stan_data,data_country,file=paste0('results/pred/','pred-',JOBID,'-stanfit.Rdata'))
  
  
  
  
  
  # p1 plot with daily detected cases (PCR)
  
  data_cases_95 <- data.frame(data_country$time, data_country$predicted_min,
                              data_country$predicted_max)
  names(data_cases_95) <- c("time", "cases_min", "cases_max")
  data_cases_95$key <- rep("nintyfive", length(data_cases_95$time))
  data_cases_50 <- data.frame(data_country$time, data_country$predicted_min2,
                              data_country$predicted_max2)
  names(data_cases_50) <- c("time", "cases_min", "cases_max")
  data_cases_50$key <- rep("fifty", length(data_cases_50$time))
  data_cases <- rbind(data_cases_95, data_cases_50)
  levels(data_cases$key) <- c("ninetyfive", "fifty")
  
  p1 <- ggplot(data_country) +
    geom_bar(data = data_country, aes(x = time, y = reported_cases),
             fill = "coral4", stat='identity', alpha=0.5) +
    geom_line(data = data_country, aes(x = time, y = predicted_cases),
              col = "black", alpha = 0.5) +
    geom_ribbon(data = data_cases,
                aes(x = time, ymin = cases_min, ymax = cases_max, fill = key)) +
    xlab("") +
    ylab("Daily number of infections") +
    scale_x_date(date_breaks = "months", labels = date_format("%e %b")) +
    scale_fill_manual(name = "", labels = c("50%", "95%"),
                      values = c(alpha("deepskyblue4", 0.55),
                                 alpha("deepskyblue4", 0.45))) +
    theme_pubr() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8),
          legend.position = "None") +
    guides(fill=guide_legend(ncol=1))
  
  
  
  # p2 plot with daily deaths
  
  data_deaths_95 <- data.frame(data_country$time, data_country$death_min,
                               data_country$death_max)
  names(data_deaths_95) <- c("time", "death_min", "death_max")
  data_deaths_95$key <- rep("nintyfive", length(data_deaths_95$time))
  data_deaths_50 <- data.frame(data_country$time, data_country$death_min2,
                               data_country$death_max2)
  names(data_deaths_50) <- c("time", "death_min", "death_max")
  data_deaths_50$key <- rep("fifty", length(data_deaths_50$time))
  data_deaths <- rbind(data_deaths_95, data_deaths_50)
  levels(data_deaths$key) <- c("ninetyfive", "fifty")
  
  p2 <-   ggplot(data_country, aes(x = time)) +
    geom_bar(data = data_country, aes(y = deaths, fill = "reported"),
             fill = "coral4", stat='identity', alpha=0.5) +
    geom_line(data = data_country, aes(x = time, y = estimated_deaths),
              col = "black", alpha = 0.5) +
    geom_ribbon(
      data = data_deaths,
      aes(ymin = death_min, ymax = death_max, fill = key)) +
    scale_x_date(date_breaks = "months", labels = date_format("%e %b")) +
    scale_fill_manual(name = "", labels = c("50%", "95%"),
                      values = c(alpha("deepskyblue4", 0.55),
                                 alpha("deepskyblue4", 0.45))) +
    theme_pubr() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8),
          legend.position = "None") +
    guides(fill=guide_legend(ncol=1))
  
  
  
  # p3 plot with daily Rt values
  
  plot_labels <- c("Major intervention",
                   "Major intervention off",
                   "Other interv. 1",
                   "Other interv. 1 off",
                   "Other interv. 2",
                   "Other interv. 2 off")
  # labels seem to be in alphabetic order, probably because variables are automatically sorted in alphabetic order
  
  # Plotting interventions
  data_rt_95 <- data.frame(data_country$time,
                           data_country$rt_min, data_country$rt_max)
  names(data_rt_95) <- c("time", "rt_min", "rt_max")
  data_rt_95$key <- rep("nintyfive", length(data_rt_95$time))
  data_rt_50 <- data.frame(data_country$time, data_country$rt_min2,
                           data_country$rt_max2)
  names(data_rt_50) <- c("time", "rt_min", "rt_max")
  data_rt_50$key <- rep("fifty", length(data_rt_50$time))
  data_rt <- rbind(data_rt_95, data_rt_50)
  levels(data_rt$key) <- c("ninetyfive", "fifth")
  
  p3 <- ggplot(data_country) +
    geom_stepribbon(data = data_rt, aes(x = time, ymin = rt_min, ymax = rt_max,
                                        group = key,
                                        fill = key)) +
    geom_line(data = data_country, aes(x = time, y = rt),
              col = "black", alpha = 0.5) +
    geom_hline(yintercept = 1, color = 'black', size = 0.1) +
    geom_segment(data = covariates_country_long,
                 aes(x = value, y = 0, xend = value, yend = max(x)),
                 linetype = "dashed", colour = "grey", alpha = 0.75) +
    geom_point(data = covariates_country_long, aes(x = value,
                                                   y = x,
                                                   group = key,
                                                   shape = key,
                                                   col = key), size = 2) +
    xlab("") +
    ylab(expression(R[t])) +
    scale_fill_manual(name = "", labels = c("50%", "95%"),
                      values = c(alpha("seagreen", 0.75), alpha("seagreen", 0.5))) +
    scale_shape_manual(name = "Interventions", labels = plot_labels,
                       values = c(21, 22, 23, 24, 25, 12)) +
    scale_colour_discrete(name = "Interventions", labels = plot_labels) +
    scale_x_date(date_breaks = "months", labels = date_format("%e %b"),
                 limits = c(data_country$time[1],
                            data_country$time[length(data_country$time)])) +
    theme_pubr() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8)) +
    theme(legend.position="right")
  
  
  
  
  # p4 plot of acummulated detected cases
  
  data_cumcases_95 <- data.frame(data_country$time, data_country$predicted_min_c,
                                 data_country$predicted_max_c)
  names(data_cumcases_95) <- c("time", "cumcases_min", "cumcases_max")
  data_cumcases_95$key <- rep("nintyfive", length(data_cumcases_95$time))
  data_cumcases_50 <- data.frame(data_country$time, data_country$predicted_min2_c,
                                 data_country$predicted_max2_c)
  names(data_cumcases_50) <- c("time", "cumcases_min", "cumcases_max")
  data_cumcases_50$key <- rep("fifty", length(data_cumcases_50$time))
  data_cumcases <- rbind(data_cumcases_95, data_cumcases_50)
  levels(data_cumcases$key) <- c("ninetyfive", "fifty")
  
  p4 <- ggplot(data_country) +
    geom_bar(data = data_country, aes(x = time, y = reported_cases_c),
             fill = "coral4", stat='identity', alpha=0.5) +
    geom_line(data = data_country, aes(x = time, y = predicted_cases_c),
              col = "black", alpha = 0.5) +
    geom_ribbon(data = data_cumcases,
                aes(x = time, ymin = cumcases_min, ymax = cumcases_max, fill = key)) +
    xlab("") +
    ylab("Cumulative number of infections") +
    scale_x_date(date_breaks = "months", labels = date_format("%e %b")) +
    scale_fill_manual(name = "", labels = c("50%", "95%"),
                      values = c(alpha("deepskyblue4", 0.55),
                                 alpha("deepskyblue4", 0.45))) +
    theme_pubr() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8),
          legend.position = "None") +
    guides(fill=guide_legend(ncol=1))
  
  
  # p5 plot of acummulated deaths
  
  data_cumdeaths_95 <- data.frame(data_country$time, data_country$death_min_c,
                                  data_country$death_max_c)
  names(data_cumdeaths_95) <- c("time", "cumdeaths_min", "cumdeaths_max")
  data_cumdeaths_95$key <- rep("nintyfive", length(data_cumdeaths_95$time))
  data_cumdeaths_50 <- data.frame(data_country$time, data_country$death_min2_c,
                                  data_country$death_max2_c)
  names(data_cumdeaths_50) <- c("time", "cumdeaths_min", "cumdeaths_max")
  data_cumdeaths_50$key <- rep("fifty", length(data_cumdeaths_50$time))
  data_cumdeaths <- rbind(data_cumdeaths_95, data_cumdeaths_50)
  levels(data_cumdeaths$key) <- c("ninetyfive", "fifty")
  
  p5 <- ggplot(data_country) +
    geom_bar(data = data_country, aes(x = time, y = deaths_c),
             fill = "coral4", stat='identity', alpha=0.5) +
    geom_line(data = data_country, aes(x = time, y = estimated_deaths_c),
              col = "black", alpha = 0.5) +
    geom_ribbon(data = data_cumdeaths,
                aes(x = time, ymin = cumdeaths_min, ymax = cumdeaths_max, fill = key)) +
    xlab("") +
    ylab("Cumulative number of deaths") +
    scale_x_date(date_breaks = "months", labels = date_format("%e %b")) +
    scale_fill_manual(name = "", labels = c("50%", "95%"),
                      values = c(alpha("deepskyblue4", 0.55),
                                 alpha("deepskyblue4", 0.45))) +
    theme_pubr() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8),
          legend.position = "None") +
    guides(fill=guide_legend(ncol=1))
  
  
  
  # p6 plot of total active cases
  
  data_active_95 <- data.frame(data_country$time, data_country$predactive_min,
                               data_country$predactive_max)
  names(data_active_95) <- c("time", "active_min", "active_max")
  data_active_95$key <- rep("nintyfive", length(data_active_95$time))
  data_active_50 <- data.frame(data_country$time, data_country$predactive_min2,
                               data_country$predactive_max2)
  names(data_active_50) <- c("time", "active_min", "active_max")
  data_active_50$key <- rep("fifty", length(data_active_50$time))
  data_active <- rbind(data_active_95, data_active_50)
  levels(data_active$key) <- c("ninetyfive", "fifty")
  
  p6 <-   ggplot(data_country, aes(x = time)) +
    geom_bar(data = data_country, aes(y = active, fill = "reported"),
             fill = "coral4", stat='identity', alpha=0.5) +
    geom_line(data = data_country, aes(x = time, y = predactive),
              col = "black", alpha = 0.5) +
    geom_ribbon(
      data = data_active,
      aes(ymin = active_min, ymax = active_max, fill = key)) +
    scale_x_date(date_breaks = "months", labels = date_format("%e %b")) +
    scale_fill_manual(name = "", labels = c("50%", "95%"),
                      values = c(alpha("deepskyblue4", 0.55),
                                 alpha("deepskyblue4", 0.45))) +
    theme_pubr() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8),
          legend.position = "None") +
    guides(fill=guide_legend(ncol=1))
  
  
  
  
  
  
  
  
  
  
  
  p <- plot_grid(p1, p2, p3, p4, p6, p5, ncol = 3, rel_widths = c(1, 1, 2))
  
  
  
  filename_pred <- paste0('pred-',JOBID)
  
  save_plot(filename = paste0("figures/Germany_2nd/", "secondbest_pred2", ".pdf"),
            p, base_width = 14, base_height = 7)
  print(pr)
  
}


