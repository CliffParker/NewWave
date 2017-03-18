library(pomp)
library(ggplot2)

#time Dataset
Data = read.table(file.choose())
Data = Data[1:469,2]/5.1e+07
time <- seq(0, 9, by = 9/((469)*7) )[1:((469)*7)]
timesD<- cbind(t=0:3282,time) 
ii <- which (timesD[,1] %in% seq(0, 3282, by = 7))
time <-timesD[ii,2]
DataT <- as.data.frame(cbind(time = time, data = Data)) # Appending time and Data


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#' Setting Initial Condition
#' S+E+I+R = 1
init2 <- Csnippet("
                  S = S_0;
                  E = E_0;
                  I = I_0;
                  R = 1-S_0-E-0-I_0;
                  ")


#' SEIR Model, the extra compartment is for reporting the number of infected
#' The model is scaled with total popultion 1 and includes seasonal forcing
seasonal.sir.ode <- Csnippet("
   double Beta = beta0*(1+beta1*cos(2*M_PI*(t-phi)));
                             DS =  mub*(1-va) -  Beta*S*I  - mud*S;
                             DE =  Beta*S*I - (lambda+mud)*E;
                             DI =  lambda*E - (gamma+mud)*I;
                             DR =  gamma*I + (mub*va) - mud*R;
                             DH =  Beta*S*I; 
                             ") #Seasonally forced snippet
# $$$$$$$$
#' Measurement process, calling on both rmeasure and dmeasure on H
#' The Data is  binomial distributed 
#' 
dmeas <- Csnippet("lik = dbinom(cases,H,rho,give_log);")
rmeas <- Csnippet("cases = rbinom(H,rho);")

#' This this the POMP object the most important part 
pomp(data=data.frame(time=time,cases=Data),
     times="time",t0=-20,
     skeleton=vectorfield(seasonal.sir.ode),
     rmeasure=rmeas,
     dmeasure=dmeas,
     initializer=init2,
     statenames=c("S","E","I","R","H"),
     paramnames=c("beta0","beta1","gamma","mub","mud","phi","va","lambda","S_0","E_0","I_0","rho")
) -> seas.sir


pomp(seas.sir,zeronames="H") -> seas.sir

###################################################'
########## Rescalling of parameters ###############'
###################################################'
#######

toEst <- Csnippet("
                  Tbeta0 = log(beta0);
                  Tgamma = log(gamma);
                  Tmub = log(mub);
                  Tmud = log(mud);
                  Tlambda = log(lambda);
                  Tva = log(va);
                  Tbeta1 = logit(beta1);
                  Tphi = logit(phi);
                  Trho = logit(rho);
                  ")

fromEst <- Csnippet("
                    Tbeta0 = exp(beta0);
                    Tgamma = exp(gamma);
                    Tmub = exp(mub);
                    Tmud = exp(mud);
                    Tlambda = exp(lambda);
                    Tva = exp(va);
                    Tbeta1 = expit(beta1);
                    Tphi = expit(phi);
                    Trho = expit(rho);
                    ")

pomp(seas.sir,toEstimationScale=toEst,
     fromEstimationScale=fromEst, paramnames=c("beta0","beta1","gamma","mub","mud","phi","va","lambda","S_0","E_0","I_0","rho")
     ) -> seas.sir

###########
###################################################'
###################################################'

params1 <- c(mub=1/50, beta0= 403.1, mud =1/50, lambda = 365/8, gamma = 365/14,
             va = 0,  beta1 = .2, phi = .5,S_0=1/23,I_0=1e-4, E_0=0, rho=0.9)

coef(seas.sir) <- params1
trajectory(seas.sir,params=params1,as=TRUE) -> x

library(ggplot2)
ggplot(x,mapping=aes(x=time,y=I))+geom_path()
ggplot(x,mapping=aes(x=S,y=I))+geom_path()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#One Parameter Bifurication Diagram

p <- parmat(coef(seas.sir),1000)
dim(p); dimnames(p)

p["beta0",] <- seq(from=50,to=800,length=1000)
y <- trajectory(seas.sir,params=p,times=time)
matplot(p["beta0",],y["I",,],pch=".",col='black',xlab='beta0',ylab='I',log='x')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

##################################################' 
##################################################' 
############# Cost functions    ##################' 
##################################################' 
# No longer needed match.trajectory used below  ##' 
##################################################' 

# Sum of square error
sse <- function (params) {
  x <- trajectory(seas.sir,params=params)
  discrep <- x["H",,]-obs(seas.sir)
  sum(discrep^2)
}

# Focused on Beta0
f1 <- function (beta) {
  params <- c(mub=1/50, beta0= beta, mud =1/50, lambda = 365/8, gamma = 365/14,
              va = 0,  beta1 = .2, phi = .5,S_0=1/23,I_0=1e-4, E_0=0, rho=0.9)
  
  sse(params)
}
beta <- seq(from=390.8,to=900,by=.1)
SSE <- sapply(beta,f1)
beta.hat <- beta[which.min(SSE)]
plot(beta,SSE,type='l')
abline(v=beta.hat,lty=2)


grid <- expand.grid(beta0=seq(from=350,to=700,by = 1),
                    mub=1/50, mud =1/50, lambda = 365/8, gamma = 365/14,
                    va = 0,  beta1 = seq(from = 0, to = 1, by = .1), phi = .5,S_0=1/23,I_0=1e-4, E_0=0,rho=0.9)


x <- trajectory(seas.sir,params=t(grid),as.data.frame=TRUE)

library(plyr)


join(x,as.data.frame(seas.sir),by="time") -> x
ddply(x,~traj,summarize,sse=sum((cases-H)^2)) -> x
cbind(grid,x) -> grid

# lots of information here.
ggplot(data = grid, mapping = aes(x=beta0,y=beta1, z=sse, fill=sse )) + geom_tile()+geom_contour(binwidth=2)



#' Using optims optimizer
f2 <- function (par) {
  params <- c(mub=1/50, beta0= par[1], mud =1/50, lambda = 365/8, gamma = 365/14,
              va = 0,  beta1 = par[2], phi = par[3] ,S_0=1/23,I_0=1e-4, E_0=0, rho=0.9)
  
  sse(params)
}
optim(fn=f2,par=c(475,.71,.5)) -> fit2
fit2

#' Maximum likelihood
#' Modelling the data? which distribution will do?

loglik.normal <- function (params) {
  x <- trajectory(seas.sir,params=params)
  sum(dnorm(x=obs(seas.sir),mean=x["H",,],
            sd=params["sigma"],log=TRUE))
}


##############################################################'
###########Parameter estimation using traj.match##############'
##############################################################


traj.match(seas.sir, start=coef(seas.sir), est = c("beta0","beta1","phi"),
 method ="subplex",transform = FALSE)->fit

#'method = c("Nelder-Mead","subplex","SANN","BFGS",
#'"sannbox","nloptr")












































