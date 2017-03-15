
#Required Packages
library(deSolve)
require(pracma)
require(FME)
require(stats4)
require(bbmle)

#Data reading
Data = read.table(file.choose())
Data = Data[1:469,2]/5.1e+07
summary(Data)
head(Data)

#time Dataset


time <- seq(0, 37, by = 37/((469+1456)*7) )[1:((469+1456)*7)]
timesD<- cbind(t=0:13474,time) 
ii <- which (timesD[,1] %in% seq(0, 13474, by = 7))[-(1:1456)]


DataT <- cbind(time = timesD[ii,2], L= Data) # Appending time and Data

#Reporting Function
report=function(out){
  I= out[,4]
  C=c()
  for (i in 0:(468+1456) ){
    C[i+1]=0
    for (j in 1:7) {
      x=(7*i)+j
      C[i+1]=C[i+1] + (365/14)*I[x]*(1/(365.25*7))  
    }
  }
  return(C)
}


"Sensitivity Analysis"
########################

################################################################################################
#Initial Conditions
s = 1/23
e = 0
i = 1e-4
r = 1-1/23-1e-4
#SEIRR Gives the Model realisation weekly and appends the Dataframe with Observation (L)

times <- seq(0, 37, by = 37/((469+1456)*7) )[1:((469+1456)*7)]

parms <- c(mub=1.5/60, beta= 600, mud = 1/60, lambda = 365/10, gamma = 365/14,
           va = 0, s = 1/23, e = 0, i = 1e-4, r = 1-1/23-1e-4 , beta2 = .2, phi = .5)

state <- c(S = s , E = e , I = i, R = r)
  
sir_rhs=function(t,state,parms){
    
    
    with(as.list(c(state, parms)),{
      #rates of change
      dS <-  mub*(1-va) -  (S*I*beta*(1 + beta2*cos(2*pi*(t - phi))))  - mud*S
      dE <- (S*I*beta*(1 + beta2*cos(2*pi*(t - phi)))) - (lambda+mud)*E
      dI <-  lambda*E - (gamma+mud)*I
      dR <- gamma*I + (mub*va) - mud*R
 
      # return the rate of change
      return(list(c(dS, dE, dI, dR)))
    })
}
  

##############################################################################################



# Model Cost
SEIRRcost<- function(P = c(mub=1.5/60, beta= 600, mud = 1/60, lambda = 365/10, gamma = 365/14,
                           va = 0, s = 1/23, e = 0, i = 1e-4, r = 1-1/23-1e-4 , beta2 = .2, phi = .5), parset =names (P)){
  parms[parset] <- P
  
  out <- ode( state,times,  sir_rhs, parms )
  new=as.data.frame(out)
  
  
  time <- seq(0, 37, by = 37/((469+1456)*7) )[1:((469+1456)*7)]
  timesD<- cbind(t=0:13474,time) 
  ii <- which (timesD[,1] %in% seq(0, 13474, by = 7))
  
  Dat <- cbind(new[ii,], L=report(out))[-(1:1456),]

  
  
  cost <- modCost(Dat, DataT)
  return(cost)
}


#Parameter Estimation
Fit <- modFit(f = SEIRRcost, p = par1 )


SEIRR<-function(pars){
  #parms[names(P)] <- P
  #pars = as.vector(pars)
  #Initial Conditions
  s = as.vector(pars[7])
  e = as.vector(pars[8])
  i = as.vector(pars[9])
  r = as.vector(pars[10])
  n = s+e+i+r 
  
  state <- c(S = s , E = e , I = i, R = r, N = n)
  
  sir_rhs=function(t,state,pars){
    
    
    with(as.list(c(state, pars)),{
      #rates of change
      dS <-  pars[1]*(1-pars[6]) -  (S*I*pars[2]*(1 + pars[11]*cos(2*pi*(t - pars[12]))))  - pars[3]*S
      dE <- (S*I*pars[2]*(1 + pars[11]*cos(2*pi*(t - pars[12])))) - (pars[4]+pars[3])*E
      dI <-  pars[4]*E - (pars[5]+pars[3])*I
      dR <- pars[5]*I + (pars[1]*pars[6]) - pars[3]*R
      dN <- pars[1] - pars[3]*N  
      
      # return the rate of change
      return(list(c(dS, dE, dI, dR, dN)))
    })
  }
  
  
  times <- seq(0, 37, by = 37/((469+1456)*7) )[1:((469+1456)*7)]
  
  out <- ode(y = state, times = times, func = sir_rhs, parms = pars)
  
  new=as.data.frame(out)
  
  
  time <- seq(0, 37, by = 37/((469+1456)*7) )[1:((469+1456)*7)]
  timesD<- cbind(t=0:13474,time) 
  ii <- which (timesD[,1] %in% seq(0, 13474, by = 7))
  
  Dat <- cbind(new[ii,], L=report(out))[-(1:1456),]
  return(Dat)
}

Sfun<-sensFun(func = SEIRR, parms = parms, senspar = c(1:6,11,12))

summary(Sfun)
plot(Sfun, which = c("L","S","E","I","R","N"), xlab ="time", lwd = 2)

# Identifiability of parameters
ident <- collin(Sfun)
ident
plot(ident, log = "y")

## 20 = magical number above which there are identifiability problems
abline(h = 20, col = "red")
