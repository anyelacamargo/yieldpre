library('minpack.lm');
library("optimx")
library(ggplot2);


optima = function(start1, w, fc)
{
  ## Optim
  w.optx <- optimx(par=start1, fn=fc, mydata=w, 
                   control=list(all.methods=TRUE, save.failures=TRUE, maxit=100))
  return(w.optx)
}

rhsfit <- function(varray, mydata) {
  sum((mydata$y - (varray[1] + (varray[2] * mydata$x) + (varray[3] *(exp(-varray[4]*mydata$x)))))^2)
}

# Effect of Soil wate
# Richter et al. Agricultural and forest, 109 (2001)
soilWaterStress = function ()
{
  Oth = 0.02; #(Sinclair (1986))
  fdt = 6:12; #(Richter)
  Orel = DAW / AWC;
  fstress = (2 / (1 + exp(-fdt * (Orel - Oth)))) - 1;
  return(fstress);
  
}

# Quantitative effect of fstress on y (f on A. Qi)
effectofSWSfitYield = function(fstress)
{
  p = list()
  v0 = 0.005866;
  umin0 = -0.000169;
  p$umin = umin0*(2 - fstress);
  p$v = v0*(1 + 0.1*(1 - fstress));
  return(p);
}
  

apparentGrowthRate = function()
{
  
  dydt = By - dy

}


relativeGrowthRate = function()
{
  tmax = t0 - (1/k) * log(-umin / (u0 - umin));
  u = umin * (1 -exp(-k*(t-tmax)));
  plot(t, u, xlab= 'time')
}


relativeGrowthRate = function()
{
  
  ymax = y0 * exp((u0/k) - ((umin/k) * (log(-umin / (u0 - umin)))));
  y = ymax * exp((umin*(t-tmax)) - ((umin / k) * (1 - (exp(-k*(t-tmax))))));
  plot(t, y, xlab= 'time')
}

roCrelativeGrowthRate = function()
{
  
  dudt = -k(u - umin)
  
}

# Extract parameters from Model I 
# model = Fitted model
# Return umin and tmax
extractParam = function(model)
{
  p = list();
  param = data.frame(summary(model)$parameters);
  p$umin = param['B.B',1];
  p$tmax = (1/param['v.v',1]) * log((param['v.v',1] * param['C.C',1]) / param['B.B',1]);
  return(p);
}



# Foliage Cover, Model I
# ydat = yield
# tdat = accumulated temperature
# Wekker and Jaggard. Annals of Botany 79 (1997) 
fitYield = function(ydat, tdat)
{
  wdata = data.frame(y=ydat, x=tdat);
  start1 = c(A = 1, B = 1, C = 0, v=0);
  eunsc = ydat ~ A + B*tdat + C*(exp(-v*tdat));
  optx = optima(start1, wdata, rhsfit);
  o = nlsLM(eunsc, 
            start=list(A = start1[1], B = start1[2], C = start1[3], r = start1[4]), trace = TRUE);
  
  
  plot(tdat, ydat, xlab = 'T accummulated temperature', ylab = 'foliage cover', main='Model I');
  lines(tdat, predict(o), col='red')
  cor(ydat, predict(o))
  return(o)
}


y = c(0.05, 0.1, 0.38, 0.6, 0.90, 0.98, 0.97);
t = c(300, 490, 510, 680, 800, 1250, 1600);
fm= fitYield(y, t);




#y0=0.000015; t0=100; u0=0.0050; umax=0.0051; ymax=0.90; kt=1.0; umin=-0.0050; k=0.00346 
#t = 0:2000;