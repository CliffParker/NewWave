

#' This is my personall pomp model from studying the course
#' its a stochastic model and its different from the exact model used in the paper. 
#' But im using the euler method derived on the normal SEIR setup with seasonal forcing.
#' i changed the scale by a factor of 1e+7 but still not working
#' The initial parameters are the same as those used in the previous models

#############################'
library(ggplot2)
library(pomp)
library(plyr)
library(reshape2)
library(magrittr)


#Old Dataset
# Data = read.table(file.choose())
# Data = Data[1:469,2]
# time <- seq(0, 9, by = 9/((469)*7) )[1:((469)*7)]
# timesD<- cbind(t=0:3282,time) 
# ii <- which (timesD[,1] %in% seq(0, 3282, by = 7))
# time <-timesD[ii,2]
# DataT <- as.data.frame(cbind(time = time, data = Data)) # Appending time and Data


daturl <- "http://kingaa.github.io/pomp/vignettes/twentycities.rda"
datfile <- file.path("/Users/cliffordallotey/Downloads","twentycities.rda")
download.file(daturl,destfile=datfile,mode="wb")
load(datfile)

demog$town = factor(demog$town)
measles$town = factor(measles$town)

"creating City datasets"
for (names in levels(demog$town)) {
  tmp<- subset(demog, town == names)
  assign( paste0(names,"_demog.rda"),tmp)
}

for (names in levels(measles$town)) {
  tmp <- subset(measles, town == names)
  assign(paste0(names,"_cases.rda"),tmp)
}

London_cases.rda %>% 
  dcast(date~"cases", value.var = "cases") %>%
  mutate(year=as.integer(format(date,"%Y"))) %>%
  subset(year>=1944 & year<1965) %>%
  mutate(time=(julian(date,origin=as.Date("1944-01-01")))/365.25+1944) %>%
  subset(time>1944 & time<1965, select=c(time,cases)) -> dat


# Csnippet for model structure
rproc <- Csnippet("
                  double Beta = beta0*(1+beta1*cos(2*M_PI*(t-phi)))/pop;
                  double t8 = rbinom(R,1-exp(-mud*dt));  // natural R death
                  double rate[6], trans[6];
                  double births;

                  
                  rate[0] = Beta*I;  //  infection
                  rate[1] = mud;             // natural S death
                  rate[2] = lambda;        // rate of ending of latent stage
                  rate[3] = mud;             // natural E death
                  rate[4] = gamma;        // recovery
                  rate[5] = mud;             // natural I death
                  
                  //births
                  births = rbinom(S+E+I+R,1-exp(-mub*dt));
                  
                  // transitions between classes
                  reulermultinom(2,S,&rate[0],dt,&trans[0]);
                  reulermultinom(2,E,&rate[2],dt,&trans[2]);
                  reulermultinom(2,I,&rate[4],dt,&trans[4]);
                  
                  S += births - trans[0] - trans[1];
                  E += trans[0] - trans[2] - trans[3];
                  I += trans[2] - trans[4] - trans[5];
                  R += trans[4] - t8;
                 
                  H += trans[4];           // true incidence
                  "
                  )

# The above uses true incidence for reporting

init <- Csnippet("
                   double m = pop/(S_0+E_0+I_0+R_0);
                   S = nearbyint(m*S_0);
                   E = nearbyint(m*E_0);
                   I = nearbyint(m*I_0);
                   R = nearbyint(m*R_0);
                   H = 0;
                   ")


dmeas <- Csnippet("lik = dbinom(cases,H,rho,give_log);")

rmeas <- Csnippet("cases = rbinom(H,rho);")

###########################################################'
########## POMP OBJECT          ######################'
###########################################################'
###########################################################'
###########################################################'
pomp(data=dat,
     times="time",
     t0=with(dat,time[1]),
     rprocess = euler.sim(rproc,delta.t=1/365),
     rmeasure=rmeas,
     dmeasure=dmeas,
     initializer=init,
     statenames=c("S","E","I","R","H"),
     paramnames=c("beta0","beta1","gamma","mub","mud","phi","lambda","S_0","E_0","R_0","I_0","rho","pop")
) -> Stoch
###########################################################'
pomp(Stoch,zeronames="H") -> Stoch


toEst <- Csnippet("
                  Tmub = log(mub);
                  Tmud = log(mud);
                  Tlambda = log(lambda);
                  Tgamma = log(gamma);
                  Tbeta0 = log(beta0);
                  Tpop = log(pop);
                  Trho = logit(rho);
                  Tbeta1 = logit(beta1);
                  Tphi = logit(phi);
                  ")

fromEst <- Csnippet("
                    Tmub = exp(mub);
                    Tmud = exp(mud);
                    Tlambda = exp(lambda);
                    Tgamma = exp(gamma);
                    Tbeta0 = exp(beta0);
                    Tpop = exp(pop);
                    Trho = expit(rho);
                    Tbeta1 = expit(beta1);
                    Tphi = expit(phi);
                    ")

pomp(Stoch,toEstimationScale=toEst,
     fromEstimationScale=fromEst,
     statenames=c("S","E","I","R","H"),
     paramnames=c("beta0","beta1","gamma","mub","mud",
                  "phi","va","lambda","S_0","E_0","I_0","R_0","rho","pop")) -> Stoch






params1 <- c(mub=1/50, beta0= 403.1, mud =1/50, lambda = 365/8, gamma = 365/14,
             va = 0,  beta1 = .2, phi = .5,S_0 =1/23,I_0=1e-4, E_0=2e-4, R_0 =1 -1/23 -3e-4 , rho=0.5, pop = 1e+7)

coef(Stoch) <- params1
# Simulation 
simulate(Stoch,params = params1, nsim=1,states=TRUE,as.data.frame=TRUE) -> x

(x[,,1000:1100])
plot(simulate(Stoch))



#' A sense of the variation in simulations based on the chosen parameter
#' Simulation is done a 100 times and infomation is used to give the variation in simulation.
Stoch %>% 
  simulate(params=params1,nsim=100,as.data.frame=TRUE,include.data=TRUE) %>%
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
Stoch %>% 
  simulate(params=params1,nsim=9,as.data.frame=TRUE,include.data=TRUE) %>%
  ggplot(aes(x=time,y=cases,group=sim,color=(sim=="data")))+
  guides(color=FALSE)+
  geom_line()+facet_wrap(~sim,ncol=2)



# Likelihood estimation working with warning
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
  pfilter(Stoch,Np=10000,params=params1)
} -> pfs
logmeanexp(sapply(pfs,logLik),se=TRUE)







# Error here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Parameter Estimation
fit<-mif2(Stoch, Nmif = 100, start = params1, Np = 1000, 
          rw.sd = rw.sd( beta0 = ivp(.52) ,
                  beta1 = ivp(.5) ,
                  phi =  ivp(.5) ,
                  pop = ivp(.5)) , transform = T,
     cooling.type = "hyperbolic", cooling.fraction.50 = .05,
     tol = 1e-17, max.fail = Inf, verbose = getOption("verbose"))

save(list="fit", file= file.path("/Users/cliffordallotey/Data/fit.rda"))

#pars=c("beta0","beta1","gamma","mub","mud","phi","lambda","S_0","E_0","I_0","rho","pop")

save()

load("fit.rda")
x= coef(fit)

coef(Stoch)<-x
plot(simulate(Stoch))













































