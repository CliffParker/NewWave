library(pomp)
library(plyr)
library(reshape2)
library(magrittr)
library(ggplot2)
library(foreach)
library(doParallel)

registerDoParallel()

set.seed(998468235L,kind="L'Ecuyer")
theme_set(theme_bw())
stopifnot(packageVersion("pomp")>="1.4.8")
set.seed(594709947L)

daturl <- "http://kingaa.github.io/pomp/vignettes/twentycities.rda"
datfile <- file.path("/Users/cliffordallotey/Downloads","twentycities.rda")
download.file(daturl,destfile=datfile,mode="wb")
load(datfile)

measles %>% 
  dcast(date~"cases", fun.aggregate = sum) %>%
  mutate(year=as.integer(format(date,"%Y"))) %>%
  subset(year>=1948 & year<1957) %>%
  mutate(time=(julian(date,origin=as.Date("1948-01-01")))/365.25+1948) %>%
  subset(time>1948 & time<1957, select=c(time,cases)) -> dat




demog %>% 
  dcast(year ~"pop", fun.aggregate = sum, value.var = c("pop")) -> popsum
demog %>% 
  dcast(year ~"births", fun.aggregate = sum, value.var = c("births")) -> birthsum
join(popsum,birthsum, by = "year") -> mydem

mydem %<>% subset(year>=1948 & year<1957) 




dat %>% ggplot(aes(x=time,y=cases))+geom_line()

mydem%>% melt(id="year") %>%
  ggplot(aes(x=year,y=value))+geom_point()+
  facet_wrap(~variable,ncol=1,scales="free_y")


#Now, we smooth the covariates. Note that we delay the entry of newborns into the susceptible pool.
mydem%>% 
  summarize(
    time=seq(from=min(year),to=max(year),by=1/12),
    pop=predict(smooth.spline(x=year,y=pop),x=time)$y,
    birthrate=predict(smooth.spline(x=year+0.5,y=births),x=time-4)$y
  ) -> covar

plot(pop~time,data=covar,type='l')
points(pop~year,data=mydem)

plot(birthrate~time,data=covar,type='l')
points(births~year,data=mydem)

plot(birthrate~I(time-4),data=covar,type='l')
points(births~I(year+0.5),data=mydem)

###############################################################################################################'
###############################################################################################################'
######################## Work begins here   ##################################################'
###############################################################################################################'
###############################################################################################################'

rproc <- Csnippet("
                  double beta, br, seas, foi, dw, births;
                  double rate[6], trans[6];
                  
                  // cohort effect
                  if (fabs(t-floor(t)-251.0/365.0) < 0.5*dt) 
                  br = cohort*birthrate/dt + (1-cohort)*birthrate;
                  else 
                  br = (1.0-cohort)*birthrate;
                  
                  // term-time seasonality
                  t = (t-floor(t))*365.25;
                  if ((t>=7&&t<=100) || (t>=115&&t<=199) || (t>=252&&t<=300) || (t>=308&&t<=356))
                  seas = 1.0+amplitude*0.2411/0.7589;
                  else
                  seas = 1.0-amplitude;
                  
                  // transmission rate
                  beta = R0*(gamma+mu)*seas;
                  // expected force of infection
                  foi = beta*pow(I+iota,alpha)/pop;
                  // white noise (extrademographic stochasticity)
                  dw = rgammawn(sigmaSE,dt);
                  
                  rate[0] = foi*dw/dt;  // stochastic force of infection
                  rate[1] = mu;             // natural S death
                  rate[2] = sigma;        // rate of ending of latent stage
                  rate[3] = mu;             // natural E death
                  rate[4] = gamma;        // recovery
                  rate[5] = mu;             // natural I death
                  
                  // Poisson births
                  births = rpois(br*dt);
                  
                  // transitions between classes
                  reulermultinom(2,S,&rate[0],dt,&trans[0]);
                  reulermultinom(2,E,&rate[2],dt,&trans[2]);
                  reulermultinom(2,I,&rate[4],dt,&trans[4]);
                  
                  S += births   - trans[0] - trans[1];
                  E += trans[0] - trans[2] - trans[3];
                  I += trans[2] - trans[4] - trans[5];
                  R = pop - S - E - I;
                  W += (dw - dt)/sigmaSE;  // standardized i.i.d. white noise
                  C += trans[4];           // true incidence
                  ")

initlz <- Csnippet("
                   double m = pop/(S_0+E_0+I_0+R_0);
                   S = nearbyint(m*S_0);
                   E = nearbyint(m*E_0);
                   I = nearbyint(m*I_0);
                   R = nearbyint(m*R_0);
                   W = 0;
                   C = 0;
                   ")
dmeas <- Csnippet("
                  double m = rho*C;
                  double v = m*(1.0-rho+psi*psi*m);
                  double tol = 1.0e-18;
                  if (cases > 0.0) {
                  lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)-pnorm(cases-0.5,m,sqrt(v)+tol,1,0)+tol;
                  } else {
                  lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)+tol;
                  }
                  ")
rmeas <- Csnippet("
                  double m = rho*C;
                  double v = m*(1.0-rho+psi*psi*m);
                  double tol = 1.0e-18;
                  cases = rnorm(m,sqrt(v)+tol);
                  if (cases > 0.0) {
                  cases = nearbyint(cases);
                  } else {
                  cases = 0.0;
                  }
                  ")
#
dat <- as.data.frame(London_cases.rda)

dat %>% 
  pomp(t0=with(dat,2*time[1]-time[2]),
       time="time",
       rprocess=euler.sim(rproc,delta.t=1/365.25),
       initializer=initlz,
       dmeasure=dmeas,
       rmeasure=rmeas,
       covar=covar,
       tcovar="time",
       zeronames=c("C","W"),
       statenames=c("S","E","I","R","C","W"),
       paramnames=c("R0","mu","sigma","gamma","alpha","iota",
                    "rho","sigmaSE","psi","cohort","amplitude",
                    "S_0","E_0","I_0","R_0")
  ) -> m1


toEst <- Csnippet("
  Tmu = log(mu);
                  Tsigma = log(sigma);
                  Tgamma = log(gamma);
                  Talpha = log(alpha);
                  Tiota = log(iota);
                  Trho = logit(rho);
                  Tcohort = logit(cohort);
                  Tamplitude = logit(amplitude);
                  TsigmaSE = log(sigmaSE);
                  Tpsi = log(psi);
                  TR0 = log(R0);
                  to_log_barycentric (&TS_0, &S_0, 4);
                  ")

fromEst <- Csnippet("
                    Tmu = exp(mu);
                    Tsigma = exp(sigma);
                    Tgamma = exp(gamma);
                    Talpha = exp(alpha);
                    Tiota = exp(iota);
                    Trho = expit(rho);
                    Tcohort = expit(cohort);
                    Tamplitude = expit(amplitude);
                    TsigmaSE = exp(sigmaSE);
                    Tpsi = exp(psi);
                    TR0 = exp(R0);
                    from_log_barycentric (&TS_0, &S_0, 4);
                    ")

pomp(m1,toEstimationScale=toEst,
     fromEstimationScale=fromEst,
     statenames=c("S","E","I","R","C","W"),
     paramnames=c("R0","mu","sigma","gamma","alpha","iota",
                  "rho","sigmaSE","psi","cohort","amplitude",
                  "S_0","E_0","I_0","R_0")) -> m1

# A plot of Cases population and birth rate
m1 %>% as.data.frame() %>% 
  melt(id="time") %>%
  ggplot(aes(x=time,y=value))+
  geom_line()+
  facet_grid(variable~.,scales="free_y")


theta= c(R0 =56.8 , mu = .02, sigma = 28.9, gamma= 30.4, alpha= 0.976 ,iota = 2.9 ,
         rho = .488, sigmaSE = .0878, psi =.116, cohort =.557, amplitude= .554,
         S_0 = .0297, E_0 = 5.17e-05, I_0=5.14e-5 ,R_0= .97)
coef(m1)<- theta


#A single simulation od the pomp object
plot(simulate(m1))
# A look at he simulation object
#simulate(m1,nsim=10,as.data.frame=TRUE)->x
simulate(m1,params=theta,nsim=100,as.data.frame=TRUE,include.data=TRUE)->x

#' A sense of the variation in simulations based on the chosen parameter
#' Simulation is done a 100 times and infomation is used to give the variation in simulation.
m1 %>% 
  simulate(params=theta,nsim=100,as.data.frame=TRUE,include.data=TRUE) %>%
  subset(select=c(time,sim,cases)) %>%
  mutate(data=sim=="data") %>%
  ddply(~time+data,summarize,
        p=c(0.05,0.5,0.95),q=quantile(cases,prob=p,names=FALSE)) %>%
  mutate(p=mapvalues(p,from=c(0.05,0.5,0.95),to=c("lo","med","hi")),
         data=mapvalues(data,from=c(TRUE,FALSE),to=c("data","simulation"))) %>%
  dcast(time+data~p,value.var='q') %>%
  ggplot(aes(x=time,y=med,color=data,fill=data,ymin=lo,ymax=hi))+
  geom_ribbon(alpha=0.2)


# A number of simulations displayed independently
m1 %>% 
  simulate(params=theta,nsim=9,as.data.frame=TRUE,include.data=TRUE) %>%
  ggplot(aes(x=time,y=cases,group=sim,color=(sim=="data")))+
  guides(color=FALSE)+
  geom_line()+facet_wrap(~sim,ncol=2)



# Likelihood estimation
library(foreach)
library(doParallel)

registerDoParallel()

set.seed(998468235L,kind="L'Ecuyer")
#' Four times estimation of likelihood
#' pfilter does the likelihood estimation
foreach(i=1:4,
        .packages="pomp",
        .options.multicore=list(set.seed=TRUE)
) %dopar% {
  pfilter(m1,Np=10000,params=theta)
} -> pfs
logmeanexp(sapply(pfs,logLik),se=TRUE)




#' Parameter estimation

fit<-mif2(m1, Nmif = 100, start = theta, Np = 1000, 
          rw.sd = rw.sd( R0 = ivp(.52) ,
                         mu = ivp(.5) ,
                         rho =  ivp(.5) ,
                         amplitude = ivp(.5)) , transform = T,
          cooling.type = "hyperbolic", cooling.fraction.50 = .05,
          tol = 1e-17, max.fail = Inf, verbose = getOption("verbose"))

fit<-continue(fit, Nmif = 100, 
          rw.sd = rw.sd( R0 = ivp(.52) ,
                         mu = ivp(.5) ,
                         rho =  ivp(.5) ,
                         amplitude = ivp(.5)))

#Log liklihood
logLik(fit)

#AIC 2(#pars - loglik)
AIC(fit)# write a function that calculates this
Aic <- function(fit) {
  l = logLik(fit)
  k = length(coef(fit))
      return(2*(k - l))  
}

Aic(fit)




































#####################################################################################################'
#####################################################################################################'
#####################################################################################################'
#####################################################################################################'
#####################################################################################################'
#####################################################################################################'
#####################################################################################################'
#####################################################################################################'
# Thu Mar  9 19:14:13 2017 ------------------------------
 
# names <- c("a","b","c")
#  a_1 =1
#  b_1 =2
#  c_1 =3
# 
# 
# 
# for (name in names) {
#   #tmp <-paste0(name,"_1")
#  assign("parest", get(paste0(name,"_1")))           #paste0(name,"_1")// gives me name
#   print(eval(parest))
# }
# 
#  
#  for (name in names) {
#    parest<-as.name(paste0(name,"_1"))
#    #assign("parest", get(paste0(name,"_1")))           #paste0(name,"_1")// gives me name
#    print(parest)
#  }
#----
#' which parameters are of interest? rewrite the rw.sd.
#' 



#' 
#' Parameters to be estimated
estpar <- c("R0")
names<-c("London")

for (name in names) {
  #Assigning to parest the right name associated parameters
  parest <- get(paste0(name,"_pars"))
  
  #first fit with larger sd's of rw
  firstFit <- mif2(m1, Nmif = 2, start = parest, Np = 2, 
                   rw.sd = rw.sd( R0 = ivp(1) ,
                                  mu = ivp(1) ,
                                  rho =  ivp(1) ,
                                  amplitude = ivp(1)) , transform = T,
                   cooling.type = "hyperbolic", cooling.fraction.50 = .05,
                   tol = 1e-17, max.fail = Inf, verbose = getOption("verbose"))
  
  # Second fit with smaller sd's of rw
  secondFit<-continue(firstFit, Nmif = 2, 
                      rw.sd = rw.sd( R0 = ivp(.05) ,
                                     mu = ivp(.05) ,
                                     rho =  ivp(.05) ,
                                     amplitude = ivp(.05)))
  
  theta<-coef(m1) <- coef(secondFit)
  # Point estimate is in secondFit.
  
  
  ###' Confidence intervals
  ###'  Dont use in parnames restrict it
  parnames <- names(coef(m1))
  
  #For loop for Confidence intervals
  for (parr in parnames) {
    #theta.t is a copy of theta
    theta.t <- partrans(m1,theta,"toEstimationScale") #parameters on toEstimation scale
    
    pd <- parmat(theta.t, nrep = 3)# A parameter replication is done to obtain a parameter matrix
    pd<-as.data.frame(t(pd))# the matix is transformed into a dataset
    #the column corresponding to whats to be estimated is replaced
    
    pd[[parr]] <- seq(from= (pd[[parr]][1]-log(2)) ,to=pd[[parr]][1]+log(2), length.out = 3) 
    
    # Transformation back to scale
    pd <- as.data.frame(t(partrans(m1,t(pd),"fromEstimationScale")))
    
    dtaf<-data.frame()
    ############# for each row it estimates the likelihood
    foreach (p=iter(pd,"row"),
             .combine=rbind,
             .errorhandling="remove",
             .inorder=FALSE,
             .options.mpi=list(chunkSize=1,seed=1598260027L,info=TRUE)
    ) %dopar% {
      p <- unlist(p)
      
      ## Runs 10 particle filters to assess Monte Carlo error in likelihood
      foreach(i=1:3,
              .packages="pomp",
              .options.multicore=list(set.seed=TRUE)
      ) %dopar% {
        pfilter(m1,Np=10,params=p)
      } -> pf
      logmeanexp(sapply(pf,logLik),se=TRUE)->ll

      nfail <- sapply(pf,getElement,"nfail")

      toc <- Sys.time()
      etime <- toc-tic
      units(etime) <- "hours"

      data.frame(as.list(p),
                 loglik = ll[1],
                 loglik.se = ll[2],
                 nfail.min = min(nfail),
                 nfail.max = max(nfail),
                 etime = as.numeric(etime))->dtat
      dtaf <- rbind(dtaf,dtat)
      
    }
    # Table of info
    assign(paste0(name,parr,"_table"),dtaf)
    # Confidence interval 95%
    maxloglik <- max(dtaf[["loglik"]])
    cutoff <- maxloglik-qchisq(p=0.95,df=1)/2
    subset(dtaf,loglik>cutoff)-> CI
    assign(paste0(name,parr,"_CI"),CI)
    
  }
  
}


#################################################################################################'
#################################################################################################'
#################################################################################################'
#################################################################################################'
#################################################################################################'
#################################################################################################'
#################################################################################################'
#################################################################################################'
#################################################################################################'
#################################################################################################'
#################################################################################################'
#################################################################################################'
stew(file="box_search_global.rda",{
  n_global <- getDoParWorkers()
  t_global <- system.time({
    mf1 <- mifs_local[[1]]
    guesses <- as.data.frame(apply(params_box,1,function(x)runif(300,x[1],x[2])))
    results_global <- foreach(guess=iter(guesses,"row"), 
                              .packages='pomp', 
                              .combine=rbind,
                              .options.multicore=list(set.seed=TRUE),
                              .export=c("mf1","fixed_params")
    ) %dopar% 
    {
      mf <- mif2(mf1,start=c(unlist(guess),fixed_params))
      mf <- mif2(mf,Nmif=100)
      ll <- replicate(10,logLik(pfilter(mf,Np=100000)))
      ll <- logmeanexp(ll,se=TRUE)
      c(coef(mf),loglik=ll[1],loglik=ll[2])
    }
  })
},seed=1270401374,kind="L'Ecuyer")
results_global <- as.data.frame(results_global)
results <- rbind(results,results_global[names(results)])
write.csv(results,file="bsflu_params.csv",row.names=FALSE)

library(plyr)
all <- ldply(list(guess=guesses,result=subset(results,loglik>max(loglik)-50)))
pairs(~Beta+mu_I+rho,data=all,col=ifelse(all$.id=="guess",grey(0.5),"red"),pch=16)
#################################################################################################'
#################################################################################################'
#################################################################################################'
#################################################################################################

bake("sigmaSE-profile1.rds",{
  
 
}) -> sigmaSE_prof







#################################################################################################'
#################################################################################################'
#################################################################################################'
#################################################################################################'
#################################################################################################'

