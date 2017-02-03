library('minpack.lm');
library("optimx")
library('ggplot2');


optima = function(start1, w, fc)
{
  ## Optim
  w.optx <- optimx(par=start1, fn=fc, mydata=w, 
                   control=list(all.methods=TRUE, save.failures=TRUE, maxit=1000))
  return(w.optx)
}

rhsfit <- function(varray, mydata) 
{
  y.pred = varray[1] + (varray[2] * mydata$x) + (varray[3] *(exp(-varray[4]*mydata$x)));
  RSS = sum((mydata$y -  y.pred)^2);
}

# Effect of Soil water stress
# Richter et al. Agricultural and forest, 109 (2001)
# DAW = Daily available water content Qi et al
# AWC = Available water content
soilWaterStress = function (DAW, AWC)
{
  Oth = 0.02; #(Sinclair (1986))
  fdt = 6:12; #(Richter)
  Orel = DAW / AWC;
  fstress = (2 / (1 + exp(-fdt * (Orel - Oth)))) - 1;
  return(fstress);
  
}

# Rooting depth (D)
# Dsowing = Sowing depth fixed
# T = daily air temperature above 3C
# l0 = sum of the lenght of epycol
# B0 = initial relative growth
# S Rate of change with B0 to 0 = 0.002715
# T0 = 120
rootingDepth = function()
{
  Dsowing = 0.02;
  l0 = 0.0491;
  B0 = 0.00935;
  S = 0.002715;
  T0 = 120;

  D = Dsowing + l0*exp((B0/S)*(1-exp(-S*(T-T0))));
  
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
  eunsc = ydat ~ A + B*tdat + C*(exp(-v*tdat)); # Equation
  
  # Transform to linear
  # log(y -A)  = log(Bt) + log(C) - vt;
  
  A.0 = min(wdata$y) * 0.5;
  model.0 = lm(log(y - A.0) ~ x, data=wdata);
  start1 = c(A = A.0, B = 1, C = exp(coef(model.0)[[1]]), v=coef(model.0)[[2]]); 
  
  # v is K in Werker
  optx = optima(start1, wdata, rhsfit);
  o = nlsLM(eunsc, 
            start=list(A = start1[1], B = start1[2], C = start1[3], v = start1[4]), trace = TRUE);
  
  A = coef(o)[1]; B = coef(o)[2]; C = coef(o)[3]; v = coef(o)[4] 
  plot(tdat, ydat, xlab = 'T accummulated temperature', 
       ylab = 'foliage cover', main='Model I', ylim=range(-0.4,1), xlim=range(0,2500), pch=21, bg='red');
  lines(tdat, predict(o), col='green')
  
  x = 0:2000;
  yeq = ydat ~ A + B*x + C*(exp(-v*x));
  y = eval(parse(text=yeq))
  lines(x,y);
  return(o)
}


ydat = c(0.75, 0.73, 0.69, 0.57, 0.55, 0.44, 0.38, 0.26);
tdat = c(1478.125, 912.5,821.875,734.375,734.375,528.125,465.625,443.75);
fm= fitYield(ydat, tdat);