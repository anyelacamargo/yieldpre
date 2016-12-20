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
  p$umin = coef(model)[2];
  p$tmax = (1/coef(model)[4]) * log((coef(model)[4] * coef(model)[3]) / coef(model)[2]);
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
            start=list(A = start1[1], B = start1[2], C = start1[3], v = start1[4]), trace = TRUE);
  
  A = coef(o)[1]; B = coef(o)[2]; C = coef(o)[3]; v = coef(o)[4] 
  plot(tdat, ydat, xlab = 'T accummulated temperature', 
       ylab = 'foliage cover', main='Model I', ylim=range(-0.4,1), xlim=range(0,2500), pch=21, bg='red');
  lines(tdat, predict(o), col='green')
  cor(ydat, predict(o))
  
  x = 0:2000;
  yeq = ydat ~ A + B*x + C*(exp(-v*x));
  y = eval(parse(text=yeq))
  lines(x,y);
  return(o)
}


y = c(0.745882353, 0.731678201, 0.689653979, 0.565986159, 0.553529412, 0.444429066, 0.380415225, 0.261262976
);
t = c(1478.125, 912.5,821.875,734.375,734.375,528.125,465.625,443.75);
fm= fitYield(y, t);




#y0=0.000015; t0=100; u0=0.0050; umax=0.0051; ymax=0.90; kt=1.0; umin=-0.0050; k=0.00346 
#t = 0:2000;