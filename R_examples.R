# univariate volatility modelling

library(rugarch)
p = read.csv('index.csv')
## We multiply returns by 100 and de-mean them
y=diff(log(p$Index))*100
y=y-mean(y)
## GARCH(1,1)
spec1 = ugarchspec(variance.model = list( garchOrder = c(1, 1)),
                   mean.model = list( armaOrder = c(0,0),include.mean = FALSE))
res1 = ugarchfit(spec = spec1, data = y)
## ARCH(1)
spec2 = ugarchspec(variance.model = list( garchOrder = c(1, 0)),
                   mean.model = list( armaOrder = c(0,0),include.mean = FALSE))
res2 = ugarchfit(spec = spec2, data = y)
## tGARCH(1,1)
spec3 = ugarchspec(variance.model = list( garchOrder = c(1, 1)),
                   mean.model = list( armaOrder = c(0,0),include.mean = FALSE),
                   distribution.model = "std")
res3 = ugarchfit(spec = spec3, data = y)
## plot(res) shows various graphical analysis, works in command line


# Multivariate volatility and correlations modelling
## read data

p = read.csv('stocks.csv')
y=apply(log(p),2,diff)    # calculate returns
y = y[,1:2]          # consider first two stocks
y[,1] = y[,1]-mean(y[,1]) # subtract mean
y[,2] = y[,2]-mean(y[,2])
TT = dim(y)[1]
plot(p[,3],type = "l")

## EWMA model
## create a matrix to hold covariance matrix for each t
EWMA = matrix(nrow=TT,ncol=3)
lambda = 0.94
S = cov(y)           # initial (t=1) covar matrix
EWMA[1,] = c(S)[c(1,4,2)]      # extract var and covar
for (i in 2:dim(y)[1]){
  S = lambda*S+(1-lambda)*  y[i-1,] %*% t(y[i-1,])
  EWMA[i,] = c(S)[c(1,4,2)]
}
EWMArho = EWMA[,3]/sqrt(EWMA[,1]*EWMA[,2]) # calculate correlations
print(head(EWMArho))
print(tail(EWMArho))
plot(EWMArho,type="l")


## DCC-model
xspec = ugarchspec(mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
                   variance.model = list(garchOrder=c(1,1), model="sGARCH"), distribution.model="norm")
uspec = multispec(replicate(2, xspec))
spec = dccspec(uspec = uspec, dccOrder = c(1, 1), distribution = 'mvnorm')
res = dccfit(spec, data = y)
H=res@mfit$H
DCCrho=vector(length=dim(y)[1])
for(i in 1:dim(y)[1]){
  DCCrho[i] =  H[1,2,i]/sqrt(H[1,1,i]*H[2,2,i])
}
plot(DCCrho, type="l")

# compare EWMA and DCC correlations
matplot(cbind(EWMArho,DCCrho),type='l',las=1,lty=1,col=2:3,ylab="")
mtext("Correlations",side=2,line=0.3,at=1,las=1,cex=0.8)
legend("bottomright",c("EWMA","DCC"),lty=1,col=2:3,bty="n",cex=0.7)



# Relationship between EWMArho and rhoDCC
plot(DCCrho, EWMArho, main = "Relationship between EWMArho and rhoDCC")




## GO-GARCH
library(rmgarch)
spec = gogarchspec(mean.model = list(armaOrder = c(0, 0),
                                     include.mean =FALSE),
                   variance.model = list(model = "sGARCH",
                                         garchOrder = c(1,1)) ,
                   distribution.model =  "mvnorm"
)
fit = gogarchfit(spec = spec, data = y)
show(fit)
# the GARCH coefficients of the independent factors
coef(fit)
# a news-impact surface plot#
ni = nisurface(fit, type = "cov", pair = c(1, 2), factor = c(1,2), plot = TRUE)
# the time varying correlation array
mc = rcor(fit)
plot(mc[1,2,], type = "l")


## GO-GARCH example 2
library(gogarch)
## Not run:
library(vars)## Boswijk / van der Weide (2009)
data(BVDWSTOXX)

BVDWSTOXX <- zoo(x = BVDWSTOXX[, -1], order.by = BVDWSTOXX[, 1])
BVDWSTOXX <- window(BVDWSTOXX, end = as.POSIXct("2007-12-31"))
BVDWSTOXX <- diff(log(BVDWSTOXX))
sectors <- BVDWSTOXX[, c("AutoParts", "Banks", "OilGas")]
sectors <- apply(sectors, 2, scale, scale = FALSE)
gogmm <- gogarch(sectors, formula = ~garch(1,1), estby = "mm",lag.max = 100)



### RISK MEASURES
# Load the data
p = read.csv('stocks.csv')
y=apply(log(p),2,diff)  
y = y[,1] 



# Specify the probability p
p <- 0.05

# Assume we have a portfolio value of 1000 USD
portfolio <- 1000

# Sort the values in y using sort()
ys <- sort(y)

# Plot them
plot(ys, type = "l", main = "Sorted returns", las = 1)
# Number of observations
n <- length(y)

# Round up
quant <- ceiling(n*p)
ys[quant]



# Visually finding the 5% quantile
plot(ys, type = "l", main = "Sorted returns", las = 1)

# Adding the segments
segments(x0 = quant, y0 = -0.5, y1 = ys[quant], lty = "dashed", lwd = 2, col ="red")
segments(x0 = -1000, y0 = ys[quant], x1 = quant, lty = "dashed", lwd = 2, col ="red")

# Adding ticks
axis(1, at = quant, label = quant)
axis(2, at = ys[quant], label = round(ys[quant],3), las = 1)

-quantile(y,  probs = 0.05)*portfolio



# Use it to calculate VaR, scaling for the portfolio value
VaR <- -ys[quant] * portfolio

# Normal VaR

sigma = sd(y) # estimate volatility
VaR3 = -sigma * qnorm(p) * portfolio_value
print(VaR3)


# For ES, we get the mean of observations up to the 5th quantile
ES <- -mean(ys[1:quant]) * portfolio
ES


###### Multivariate HS

# Load the data
p <- read.csv("stocks.csv")
# Create a new data frame 

n <- nrow(p)

y <- ((p[2:n, ] - p[1:(n-1), ])/p[1:(n-1), ])
y <- as.matrix(y)
#y=apply(log(p),2,diff)  


w = matrix(c(0.3,0.2,0.5)) # vector of portfolio weights
## The number of columns of the left matrix must be the same as the number of rows of the right matrix
yp = y %*% w         # obtain portfolio returns
yps = sort(yp)

p=0.01
y1 = y[,1]           # select one asset
op = ceiling(length(y1)*p) # p percent smallest, rounded up
portfolio_value=1000
VaR2 = -yps[op]*portfolio_value
print(VaR2)

sigma = sqrt(t(w) %*% cov(y) %*% w)[1] # portfolio volatility
## Note: the trailing [1] is to convert a single element matrix to float
VaR4 = -sigma * qnorm(p)*portfolio_value
print(VaR4)

# Normal ES
sigma = sd(y1)
ES2 = sigma*dnorm(qnorm(p))/p * portfolio_value
print(ES2)

# Direct integration when computing ES
VaR = -qnorm(p)
integrand = function(q){q*dnorm(q)}
ES = -sigma*integrate(integrand,-Inf,-VaR)$value/p*portfolio_value
print(ES)

# rolling estimation window
# HS
WE=20
for (t in seq(length(y1)-5,length(y1))){
  t1=t-WE+1
  window= y1[t1:t] # estimation window
  sigma=sd(window)
  VaR6 = -sigma * qnorm(p) * portfolio_value
  print(VaR6)
}

# EWMA
lambda = 0.94
s11 = var(y1) # initial variance, using unconditional
for (t in 2:length(y1)){
  s11 = lambda * s11  + (1-lambda) * y1[t-1]^2
}
VaR7 = -qnorm(p) * sqrt(s11) * portfolio_value
print(VaR7)


# Portfolio EWMA
s = cov(y)           # initial covariance
for (t in 2:dim(y)[1]){
  s = lambda*s + (1-lambda)*y[t-1,] %*% t(y[t-1,])
}
sigma = sqrt(t(w) %*% s %*% w)[1] # portfolio vol
## Note: trailing[1] is to convert single element matrix to float
VaR8 = -sigma * qnorm(p) * portfolio_value
print(VaR8)


# Univariate with GARCH
library(rugarch)
spec = ugarchspec(variance.model = list( garchOrder = c(1, 1)),
                  mean.model = list( armaOrder = c(0,0),include.mean = FALSE))
res = ugarchfit(spec = spec, data = y1)
omega = res@fit$coef['omega']
alpha = res@fit$coef['alpha1']
beta = res@fit$coef['beta1']
sigma2 = omega + alpha * tail(y1,1)^2 + beta * tail(res@fit$var,1)
VaR9 = -sqrt(sigma2) * qnorm(p) * portfolio_value
names(VaR9)="VaR"
print(VaR9)




###### BACKTESTING and STRESSTESTING

p = read.csv('index.csv')
y = diff(log(p$Index)) # get returns


# set up backtest
WE = 1000  # estimation window length
p = 0.01   # probability
l1 = ceiling(WE*p)   # HS quantile
portfolio_value = 1  # portfolio value
VaR = matrix(nrow=length(y),ncol=4) # matrix for forecasts
## EWMA setup
lambda = 0.94
s11 = var(y)
for(t in 2:WE) {
  s11=lambda*s11+(1-lambda)*y[t-1]^2
}


# Estimate VaR using several methods


## this will take several minutes to run
library(rugarch)
spec = ugarchspec(variance.model = list( garchOrder = c(1, 1)),
                  mean.model = list( armaOrder = c(0,0),include.mean = FALSE))
start_time <- Sys.time()
for (t in (WE+1):length(y)){
  t1 = t-WE;         # start of the data window
  t2 = t-1;          # end of the data window
  window = y[t1:t2]  # data for estimation
  s11 = lambda*s11 + (1-lambda)*y[t-1]^2
  VaR[t,1] = -qnorm(p) * sqrt(s11) * portfolio_value    # EWMA
  VaR[t,2] = - sd(window) * qnorm(p)*portfolio_value    # MA
  ys = sort(window)
  VaR[t,3] = -ys[l1]*portfolio_value     # HS
  res = ugarchfit(spec = spec, data = window,solver="hybrid")
  omega = res@fit$coef['omega']
  alpha = res@fit$coef['alpha1']
  beta = res@fit$coef['beta1']
  sigma2 = omega + alpha * tail(window,1)^2 + beta * tail(res@fit$var,1)
  VaR[t,4] = -sqrt(sigma2) * qnorm(p) * portfolio_value # GARCH(1,1)
}
end_time <- Sys.time()
print(end_time - start_time)


# plot the analysis
for (i in 1:4){
  VR = sum(y[(WE+1):length(y)]< -VaR[(WE+1):length(y),i])/(p*(length(y)-WE))
  s = sd(VaR[(WE+1):length(y),i])
  cat(i,"VR",VR,"VaR vol",s,"\n")
}
matplot(cbind(y,VaR),type='l',col=1:5,las=1,ylab="",lty=1:5)
legend("topleft",legend=c("Returns","EWMA","MA","HS","GARCH"),lty=1:5,col=1:5,bty="n")

# this is how to code a simple version of the Bernoulli coverage test

bern_test=function(p,v){
  lv=length(v)
  sv=sum(v)
  al=log(p)*sv+log(1-p)*(lv-sv)
  bl=log(sv/lv)*sv +log(1-sv/lv)*(lv-sv)
  return(-2*(al-bl))
}


# and this is the independence test of Christoferson
ind_test=function(V){
  J=matrix(ncol=4,nrow=length(V))
  for (i in 2:length(V)){
    J[i,1]=V[i-1]==0 & V[i]==0
    J[i,2]=V[i-1]==0 & V[i]==1
    J[i,3]=V[i-1]==1 & V[i]==0
    J[i,4]=V[i-1]==1 & V[i]==1
  }
  V_00=sum(J[,1],na.rm=TRUE)
  V_01=sum(J[,2],na.rm=TRUE)
  V_10=sum(J[,3],na.rm=TRUE)
  V_11=sum(J[,4],na.rm=TRUE)
  p_00=V_00/(V_00+V_01)
  p_01=V_01/(V_00+V_01)
  p_10=V_10/(V_10+V_11)
  p_11=V_11/(V_10+V_11)
  hat_p=(V_01+V_11)/(V_00+V_01+V_10+V_11)
  al = log(1-hat_p)*(V_00+V_10) + log(hat_p)*(V_01+V_11)
  bl = log(p_00)*V_00 + log(p_01)*V_01 + log(p_10)*V_10 + log(p_11)*V_11
  return(-2*(al-bl))
}


# applying this to S&P500
W1=WE+1
ya=y[W1:length(y)]
VaRa=VaR[W1:length(y),]
m=c("EWMA","MA","HS","GARCH")
for (i in 1:4){
  q= y[W1:length(y)]< -VaR[W1:length(y),i]
  v=VaRa*0
  v[q,i]=1
  ber=bern_test(p,v[,i])
  ind=ind_test(v[,i])
  cat(i,m[i], "\n",
      "Bernoulli - ","Test statistic:",ber,"  p-value:",1-pchisq(ber,1),"\n",
      "Independence - ", "Test statistic:",ind,"  p-value:",1-pchisq(ind,1),"\n \n")
}

# let's turn to expected shortfall across models

VaR2 = matrix(nrow=length(y), ncol=2)    # VaR forecasts for 2 models
ES = matrix(nrow=length(y), ncol=2)      # ES forecasts for 2 models
for (t in (WE+1):length(y)){
  t1 = t-WE;
  t2 = t-1;
  window = y[t1:t2]
  s11 = lambda * s11  + (1-lambda) * y[t-1]^2
  VaR2[t,1] = -qnorm(p) * sqrt(s11) * portfolio_value # EWMA
  ES[t,1] = sqrt(s11) * dnorm(qnorm(p)) / p
  ys = sort(window)
  VaR2[t,2] = -ys[l1] * portfolio_value  # HS
  ES[t,2] = -mean(ys[1:l1]) * portfolio_value
}


# and see how was the ratio of predicted vs realized
ESa = ES[W1:length(y),]
VaRa = VaR2[W1:length(y),]
for (i in 1:2){
  q = ya <= -VaRa[,i]
  nES = mean(ya[q] / -ESa[q,i])
  cat(i,"nES",nES,"\n")
}


################### Options and Bonds
# Black-Scholes formula
bs = function(X, P, r, sigma, T){
  d1 = (log(P/X) + (r + 0.5*sigma^2)*(T))/(sigma*sqrt(T))
  d2 = d1 - sigma*sqrt(T)
  Call = P*pnorm(d1,mean=0,sd=1)-X*exp(-r*(T))*pnorm(d2,mean=0,sd=1)
  Put = X*exp(-r*(T))*pnorm(-d2,mean=0,sd=1)-P*pnorm(-d1,mean=0,sd=1)
  Delta.Call = pnorm(d1, mean = 0, sd = 1)
  Delta.Put = Delta.Call - 1
  Gamma = dnorm(d1, mean = 0, sd = 1)/(P*sigma*sqrt(T))
  return(list(Call=Call,Put=Put,Delta.Call=Delta.Call,Delta.Put=Delta.Put,Gamma=Gamma))
}


f = bs(X = 90, P = 100, r = 0.05, sigma = 0.2, T = 0.5)
print(f)


# Pseudo randon numbers
x = seq(-3, 3, by = 0.1)
plot(x, pnorm(x), type = "l")

# random numbers in R
set.seed(12) # set seed
S = 10
runif(S)
rnorm(S)
rt(S,4)

# Bond pricing
yield = c(5.00, 5.69, 6.09, 6.38, 6.61,
          6.79, 6.94, 7.07, 7.19, 7.30) # yield curve
T = length(yield)    # number of time periods
r = 0.07   # initial yield rate
Par = 10   # par value
coupon = r * Par     # coupon payments
cc = rep(coupon, T)  # vector of cash flows
cc[T] = cc[T] + Par  # add par to cash flows
P = sum(cc/((1+yield/100)^(1:T)))     # calculate price
print(P)



# Let's simulate yields
set.seed(12)         # set seed
sigma = 1.5          # daily yield volatiltiy
S = 8      # number of simulations
r = rnorm(S, 0, sigma) # generate random numbers
ysim = matrix(nrow=length(yield),ncol=S)
for (i in 1:S) ysim[,i]=yield+r[i]
matplot(ysim,type='l')


# Let's now simulate Bond prices
SP = rep(NA, length = S)
for (i in 1:S){      # S simulations
  SP[i] = sum(cc/((1+ysim[,i]/100)^(1:T)))
}
SP = SP-(mean(SP) - P) # correct for mean
par(mfrow=c(1,2), pty="s")
barplot(SP)
hist(SP,probability=TRUE)
x=seq(6,16,length=100)
lines(x, dnorm(x, mean = mean(SP), sd = sd(SP)))
S = 50000
r = rnorm(S, 0, sigma) # generate random numbers
ysim = matrix(nrow=length(yield),ncol=S)
for (i in 1:S) ysim[,i]=yield+r[i]
SP = rep(NA, length = S)
for (i in 1:S){      # S simulations
  SP[i] = sum(cc/((1+ysim[,i]/100)^(1:T)))
}
SP = SP-(mean(SP) - P) # correct for mean
par(mfrow=c(1,2), pty="s")
barplot(SP)
hist(SP,probability=TRUE)
x=seq(6,16,length=100)
lines(x, dnorm(x, mean = mean(SP), sd = sd(SP)))


# let's value  options using Black-Scholes

P0 = 50    # initial spot price
sigma = 0.2          # annual volatility
r = 0.05   # annual interest
TT = 0.5   # time to expiration
X = 40     # strike price
f = bs(X = X, P = P0, r = r, sigma = sigma, T = TT) # analytical call price
## This calculation uses the Black-Scholes pricing function (Listing 6.1/6.2)
print(f)


# And now let's value options using MC, we're simulating the final prices

set.seed(12)         # set seed
S = 1e6    # number of simulations
F = P0*exp(r*TT)     # futures price
ysim = rnorm(S,-0.5*sigma^2*TT,sigma*sqrt(TT)) # sim returns, lognorm corrected
F = F*exp(ysim)      # sim futures price
SP = F-X   # payoff
SP[SP<0] = 0         # set negative outcomes to zero
fsim = SP*exp(-r*TT)           # discount
call_sim = mean(fsim)          # simulated price
print(call_sim)

# plot results
par(mfrow=c(1,2), pty="s")
hist(F, probability=TRUE, ylim=c(0,0.06))
x = seq(min(F), max(F), length=100)
lines(x, dnorm(x, mean = mean(F), sd = sd(SP)))
hist(fsim, nclass=100, probability=TRUE)



# Simulate VaR
set.seed(1)          # set seed
S = 1e7    # number of simulations
s2 = 0.01^2          # daily variance
p = 0.01   # probability
r = 0.05   # annual riskfree rate
P = 100    # price today
ysim = rnorm(S,r/365-0.5*s2,sqrt(s2)) # sim returns
Psim = P*exp(ysim)   # sim future prices
q = sort(Psim-P)     # simulated P/L
VaR1 = -q[p*S]
print(VaR1)



# Simulate option VaR

TT = 0.25  # time to expiration
X = 100    # strike price
sigma = sqrt(s2*250)           # annual volatility
f = bs(X, P, r, sigma, TT)     # analytical call price
fsim = bs(X,Psim,r,sigma,TT-(1/365)) # sim option prices
q = sort(fsim$Call-f$Call)     # simulated P/L
VaR2 = -q[p*S]
print(VaR2)


# example
X1 = 100
X2 = 110
f1 = bs(X1, P, r, sigma, TT)
f2 = bs(X2, P, r, sigma, TT)
f2sim = bs(X2, Psim, r, sigma, TT-(1/365))
f1sim = bs(X1, Psim, r, sigma, TT-(1/365))
q = sort(f1sim$Call + f2sim$Put + Psim-f1$Call - f2$Put-P);
VaR3 = -q[p*S]
print(VaR3)

# Simulate two asset returns
library (MASS)
set.seed(12)         # set seed
mu = rep(r/365, 2)   # return mean
Sigma = matrix(c(0.01, 0.0005, 0.0005, 0.02),ncol=2) # covariance matrix
y = mvrnorm(S,mu,Sigma)        # simulated returns

# VaR of portfolio of two assets
K=2
P = c(100, 50)       # prices
x = rep(1, 2)        # number of assets
Port = P %*% x       # portfolio at t
Psim = matrix(t(matrix(P,K,S)),ncol=K)*exp(y) # simulated prices
PortSim = Psim %*% x           # simulated portfolio value
q = sort(PortSim-Port[1,1])    # simulated P/L
VaR4 = -q[S*p]
print(VaR4)

# two asset case with an option
f = bs(X = P[2], P = P[2], r = r, sigma = sigma, T = TT)
fsim = bs(X = P[2], P = Psim[,2], r = r, sigma = sigma, T = TT-(1/365))
q = sort(fsim$Call + Psim[,1] - f$Call - P[1]);
VaR5 = -q[p*S]
print(VaR5)
