library(pomp)
library(plyr)
library(reshape2)
library(magrittr)
library(ggplot2)
library(foreach)
library(doParallel)

registerDoParallel()
set.seed(998468235L,kind="L'Ecuyer")
stopifnot(packageVersion("pomp")>="1.4.8")


daturl <- "http://kingaa.github.io/pomp/vignettes/twentycities.rda"
datfile <- file.path("/Users/cliffordallotey/Downloads","twentycities.rda") #Change this to be downloaded online  !!!
download.file(daturl,destfile=datfile,mode="wb")
load(datfile)

demog$town = factor(demog$town)
measles$town = factor(measles$town)

"creating City datasets"
for (names in levels(demog$town)) {
  tmp<- subset(demog, town == names)
  tmp<-tmp[,-1]
  tmp %>% subset(year>=1944 & year<1964) %>% 
    summarize(
      time=seq(from=min(year),to=max(year),by=1/12),
      pop=predict(smooth.spline(x=year,y=pop),x=time)$y,
      birthrate=predict(smooth.spline(x=year+0.5,y=births),x=time-4)$y
    ) -> covar
  
  assign( paste0(names,"_covar"),covar)
}


"Cases"
for (names in levels(measles$town)) {
  tmp <- subset(measles, town == names)
  tmp %>% 
    dcast(date~"cases", fun.aggregate = sum) %>%
    mutate(year=as.integer(format(date,"%Y"))) %>%
    subset(year>=1944 & year<1965) %>%
    mutate(time=(julian(date,origin=as.Date("1944-01-01")))/365.25+1944) %>%
    subset(time>1944 & time<1965, select=c(time,cases)) -> tmp
  assign(paste0(names,"_cases"),tmp)
}


## ----mles,include=FALSE--------------------------------------------------Parameter estimates from Aaron's model.
read.csv(text="
         town,loglik,loglik.sd,mu,delay,sigma,gamma,rho,R0,amplitude,alpha,iota,cohort,psi,S_0,E_0,I_0,R_0,sigmaSE
Bedwellty,-1125.1,0.14,0.02,4,57.9,146,0.311,24.7,0.16,0.937,0.0396,0.351,0.951,0.0396,2.64e-05,2.45e-05,0.96,0.0611
Birmingham,-3239.3,1.55,0.02,4,45.6,32.9,0.544,43.4,0.428,1.01,0.343,0.331,0.178,0.0264,8.96e-05,0.000335,0.973,0.0611
Bradford,-2586.6,0.68,0.02,4,45.6,129,0.599,32.1,0.236,0.991,0.244,0.297,0.19,0.0365,7.41e-06,4.59e-06,0.964,0.0451
Bristol,-2681.6,0.5,0.02,4,64.3,82.6,0.626,26.8,0.203,1.01,0.441,0.344,0.201,0.0358,9.62e-06,5.37e-06,0.964,0.0392
Cardiff,-2364.9,0.73,0.02,4,39,143,0.602,34.4,0.223,0.996,0.141,0.267,0.27,0.0317,1.01e-05,9.21e-06,0.968,0.0539
Consett,-1362.9,0.73,0.02,4,42.6,172,0.65,35.9,0.2,1.01,0.0731,0.31,0.406,0.0322,1.83e-05,1.97e-05,0.968,0.0712
Dalton.in.Furness,-726.1,0.3,0.02,4,73.6,257,0.455,28.3,0.203,0.989,0.0386,0.421,0.818,0.0387,2.23e-05,2.36e-05,0.961,0.0779
Halesworth,-318.6,0.51,0.02,4,49.6,210,0.754,33.1,0.381,0.948,0.00912,0.547,0.641,0.0526,1.99e-05,2.82e-05,0.947,0.0748
Hastings,-1583.7,0.21,0.02,4,56.3,74.1,0.695,34.2,0.299,1,0.186,0.329,0.396,0.0233,5.61e-06,3.4e-06,0.977,0.0955
Hull,-2729.4,0.39,0.02,4,42.1,73.9,0.582,38.9,0.221,0.968,0.142,0.275,0.256,0.0371,1.2e-05,1.13e-05,0.963,0.0636
Leeds,-2918.6,0.23,0.02,4,40.7,35.1,0.666,47.8,0.267,1,1.25,0.592,0.167,0.0262,6.04e-05,3e-05,0.974,0.0778
Lees,-548.1,1.1,0.02,4,45.6,244,0.612,29.7,0.153,0.968,0.0311,0.648,0.681,0.0477,2.66e-05,2.08e-05,0.952,0.0802
Liverpool,-3403.1,0.34,0.02,4,49.4,39.3,0.494,48.1,0.305,0.978,0.263,0.191,0.136,0.0286,0.000184,0.00124,0.97,0.0533
London,-3804.9,0.16,0.02,4,28.9,30.4,0.488,56.8,0.554,0.976,2.9,0.557,0.116,0.0297,5.17e-05,5.14e-05,0.97,0.0878
Manchester,-3250.9,0.66,0.02,4,34.4,56.8,0.55,32.9,0.29,0.965,0.59,0.362,0.161,0.0489,2.41e-05,3.38e-05,0.951,0.0551
Mold,-296.5,0.25,0.02,4,67.4,301,0.131,21.4,0.271,1.04,0.0145,0.436,2.87,0.064,2.61e-05,2.27e-05,0.936,0.0544
Northwich,-1195.1,2.25,0.02,4,45.6,147,0.795,30.1,0.423,0.948,0.0602,0.236,0.402,0.0213,1.32e-05,1.58e-05,0.979,0.0857
Nottingham,-2703.5,0.53,0.02,4,70.2,115,0.609,22.6,0.157,0.982,0.17,0.34,0.258,0.05,1.36e-05,1.41e-05,0.95,0.038
Oswestry,-696.1,0.49,0.02,4,37.3,168,0.631,52.9,0.339,1.04,0.0298,0.263,0.476,0.0218,1.56e-05,1.61e-05,0.978,0.0699
Sheffield,-2810.7,0.21,0.02,4,54.3,62.2,0.649,33.1,0.313,1.02,0.853,0.225,0.175,0.0291,6.04e-05,8.86e-05,0.971,0.0428
         ",stringsAsFactors=F) -> mles
mles$phi<-0.5
## ----mle-----------------------------------------------------------------


#Creating parameter vectors
for (name in levels(measles$town)) {
  mles %>% subset(town==paste0(name)) -> mle
  paramnames <- c("R0","mu","sigma","gamma",
                  "rho","amplitude","phi",
                  "S_0","E_0","I_0","R_0")
  mle %>% extract(paramnames) %>% unlist() %>% eval() -> theta
  assign(paste0(name,"_pars"),theta)
}


# created the Datasets required for this work
###############################################################################################################'
###############################################################################################################'
########################       Creating pomp Body parts      ##################################################'
###############################################################################################################'
###############################################################################################################'



initz <- Csnippet("
                   double m = pop/(S_0+E_0+I_0+R_0);
                 S = m*S_0;
                 E = m*E_0;
                 I = m*I_0;
                 R = m*R_0;
                 H = 0;
                 ")

#' SEIR Model, the extra compartment is for reporting the number of infected
#' The model is scaled with total popultion 1 and includes seasonal forcing
seasonal.sir.ode <- Csnippet("
                              double beta0;
                             // transmission rate
                             beta0 = R0*(gamma+mu)*(sigma+mu)/sigma;
                             double Beta = beta0*(1+amplitude*cos(2*M_PI*(t-phi)))/pop;
                             

                             DS =  birthrate -  Beta*S*I  - mu*S;
                             DE =  Beta*S*I - (sigma+mu)*E;
                             DI =  sigma*E - (gamma+mu)*I;
                             DR =  gamma*I  - mu*R;
                             DH =  gamma*I; 
                             ") #Seasonally forced snippet

toEst <- Csnippet("
                  Tmu = log(mu);
                  Tsigma = log(sigma);
                  Tgamma = log(gamma);
                  TR0 = log(R0);
                  Trho = logit(rho);
                  Tamplitude = logit(amplitude);
                  Tphi = logit(phi);
                  ")

fromEst <- Csnippet("
                    Tmu = exp(mu);
                    Tsigma = exp(sigma);
                    Tgamma = exp(gamma);
                    TR0 = exp(R0);
                    Trho = expit(rho);
                    Tamplitude = expit(amplitude);
                    Tphi = expit(phi);
                    ")
#' Measurement process, calling on both rmeasure and dmeasure on H
#' The Data is  binomial distributed 
#' 

dmeas <- Csnippet("
  double m = rho*H;
                  double v = m*(1.0-rho);
                  lik = dnorm(cases,m,sqrt(v),give_log);
                  ")

rmeas <- Csnippet("
                  double m = rho*H;
                  double v = m*(1.0-rho);
                  cases = rnorm(m,sqrt(v));
                  ")
# dmeas <- Csnippet("lik = dpois(cases,rho*H,give_log);
#                    ")
# 
# rmeas <- Csnippet("cases = rpois(rho*H);")

##############################################################################################'
#################################################################################################'
##### Main ######################################################################################'
######################################## Object #################################################'
#################################################################################################'
########################### to ##################################################################'
#################################################################################################'
#################################################################################################'
#################################################################### run ########################'
#################################################################################################'
#################################################################################################'

stew(file="DeterministicPOMP_M.rda",{
  
  registerDoParallel()
  #####################################################################################################'
  
  #' Parameters to be estimated
  names<-levels(demog$town)
 
  
  for (name in names) {
    #Assigning to parest the right name associated parameters
    parest <- get(paste0(name,"_pars"))
    
    #######################################################################'
    #                            POMP BODY
    #######################################################################'
    get(paste0(name,"_cases")) %>% 
       pomp(
         times="time",
         t0=with(get(paste0(name,"_cases")),2*time[1]-time[2]),
         skeleton=vectorfield(seasonal.sir.ode),
         rmeasure=rmeas,
         covar=get(paste0(name,"_covar")),
         tcovar="time",
         toEstimationScale=toEst,
         fromEstimationScale=fromEst,
         dmeasure=dmeas,
         initializer=initz,
         statenames=c("S","E","I","R","H"),
         zeronames=c("H"),
         paramnames=c("R0","amplitude","gamma","mu","phi","sigma","S_0","E_0","R_0","I_0","rho")
    ) -> m1
    
    
    
    
    
    coef(m1)<- parest
    
    
    
    #######################################################################'
    ####################################################'
    
    #first  fitting iteration 
    firstFit <- traj.match(m1 , start=coef(m1), est = c("R0","amplitude","gamma","mu","phi","sigma"),
                           method ="subplex",transform = T)
    parest<-coef(firstFit)
    secondFit <- traj.match(m1 , start=parest, est = c("R0","amplitude","gamma","mu","phi","sigma"),
                            method ="subplex",transform = T)
    
    
    
    theta<-coef(m1) <- coef(secondFit)
    #Saving model
    assign(paste0(name,"_Model"),m1)
    
    
    # Point estimate is in secondFit.
    ###' Confidence intervals
    parnames <- names(coef(m1))
  
    
    
    
    
    #For loop for Confidence intervals
    for (parr in parnames) {
      estpars <- setdiff(names(theta),c(paste0(parr),"S_0","E_0","R_0","I_0","rho"))# parameter space to be profiled on
      
      
      
      theta.t <- partrans(m1,theta,"toEstimationScale")
      theta.t.hi <- theta.t.lo <- theta.t
      #parspace
      theta.t.lo[estpars] <- theta.t[estpars]-log(5) #Lower bound for parspace to be profiled on
      theta.t.hi[estpars] <- theta.t[estpars]+log(5) #Upper bound for parspace to be profiled on
      #estspace
      FROM <- theta.t[paste0(parr)]-log(5) #Lower bound for parspace to be profiled on
      TO <- theta.t[paste0(parr)]+log(5) #Upper bound for parspace to be profiled on
      
      
      
      profileDesign(
        assign(paste0(parr),seq(from=FROM ,to=TO ,length=50)),# 2 was 20
        lower=theta.t.lo,upper=theta.t.hi,nprof=40            # 4 was 40
      ) -> pd 
      names(pd)[1]<-paste0(parr)
      
      pd <- as.data.frame(t(partrans(m1,t(pd),"fromEstimationScale")))

      
      ######################################################################################################'
      ######################################################################################################'
      #########################  starters   created   for maximisation ######################################'
      ######################################################################################################'
      
      ###########################################################################################'
      ######################################################################################################'
      
      dtaf<-data.frame()# to store the info
      
      foreach (p=iter(pd,"row"),
               .combine=rbind,
               .errorhandling="remove",
               .inorder=FALSE,
               .options.mpi=list(chunkSize=1,seed=1598260027L,info=TRUE)
      ) %dopar% {
        p <- unlist(p)
        
        tic <- Sys.time()
        
        library(magrittr)
        library(plyr)
        library(reshape2)
        library(pomp)
        
        options(stringsAsFactors=FALSE)
        
        
        
        m1 %>% 
          traj.match( start=unlist(p), est = estpars,
                      method ="subplex",transform = T) -> mf
        
        ##################################################################################################'       
        
        
        
        
        ll <- logLik(mf)
        toc <- Sys.time()
        etime <- toc-tic
        units(etime) <- "hours"
        
        data.frame(as.list(coef(mf)),
                   loglik = ll,
                   etime = as.numeric(etime))
        
      }->dtat
      
      dtaf <- rbind(dtaf,dtat)
      # Table of info
      
      assign(paste0(name,"_",parr,"-profile1"),dtaf)
      
      
      # Confidence interval 95%
      maxloglik <- max(dtaf[["loglik"]])
      cutoff <- maxloglik-qchisq(p=0.95,df=1)/2
      subset(dtaf,loglik>cutoff)-> CI
      assign(paste0(name,parr,"_CI"),CI)
      
    }
    
  }
  
  
  #################################################################################################'
  #
  
})->Realdata
################################################################################################'
#################################################################################################'
#################################################################################################'
#################################################################################################











































