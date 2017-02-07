library('minpack.lm');
library("optimx")
library('ggplot2');
source('C:/Anyela/M/senescence_disease/generic.R');

# create cumulative sum
# data = data to use
# variable: variable to cumulate
cumSum = function(data, variable)
{
  
  return(sumsum(data[[variable]]));
  
}




# Intercepted radiation use efficiency (Eq 7). A.Qi et al
interceptedRadiationUE = function(measurements_table)
{
  
  e0 = measurements_table[measurements_table$Variables=='e0', 'Value'];
  l = measurements_table[measurements_table$Variables=='l', 'Value'];
  #e = e0*(Ea/Ep)  exp(-l*W);
  
  
}


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
PredYield = function(ydat, tdat)
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

# Eq 17
# Soil water potential in the rooting zone
calQfcRooting = function(measurements_table)
{
  b = measurements_table[measurements_table$Variables=='b', 'Value'];
  Qfcroot = -5*(Q /calQfc(measurements_table))^b;
  return(Qfcroot);
}

# Eq 18
# Q
# SMD = Soil moisture deficit
calQ = function()
{
  
  Q = calQfc(measurements_table) - (SMD / D);
  return(Q);  
}

# Eq 19
# SSE = 
calSMD = function()
{
  
  R = cumsum(avgdata$Precipitation);
  
  SMD = SMD + SSE + Ea - R;
  
}



#Eq 16
calQfc = function(measurements_table)
{
  a1 = measurements_table[measurements_table$Variables=='a1', 'Value'];
  a2 = measurements_table[measurements_table$Variables=='a2', 'Value'];
  b = measurements_table[measurements_table$Variables=='b', 'Value'];
  Qfc = a2 * (a1/5) ^ (1/b)
  return(Qfc);
  
}

#Qrel = AWC / Qfc;
calQrel = function(measurements_table)
{
  
  Qrel = AWC / calQfc(measurements_table);
  
  
}

#Eq 13
# DAily crop actual transpiration
calEa = function()
{
  Ea = min(Ep, Emax);
  return(Ea);
  
}

# Eq 12
# Daily Crop potential transpiration
# f = 
# ETgrass = Daily potential evapotranspiration
calEp = function(measurements_table)
{
  
  ET_grass = measurements_table[measurements_table$Variables=='ET_grass', 'Value'];
  Ep = (1.25*cal) * ETgrass;
}


# Eq 11
# SSE
calSSE = function(measurements_table, avgdata)
{
  
  SSE = 1.5*(1 - calf(measurements_table, avgdata))
  return(SSE);
}


#Eq 10
# Net crop yield
calY = function(measurements_table, avgdata)
{
  k = measurements_table[measurements_table$Variables=='k', 'Value'];
  
  Y = W - (1/k)*log((k*W) + 1);
  return(Y);
  
}


#Eq 7
# Intercepted radiation use efficiency
cale = function(measurements_table, avgdata)
{
  e0 = measurements_table[measurements_table$Variables=='e0', 'Value'];
  l = measurements_table[measurements_table$Variables=='l', 'Value'];
  W = measurements_table[measurements_table$Variables=='W', 'Value'];
  Ea = calEa(measurements_table, avgdata);
  Ep = calEp(measurements_table, avgdata);
  W = calW(measurements_table, avgdata);
  
  e = e0*(Ea/Ep)*exp(-s*W);
}


#Eq 6
# Crop dry matter
calW = function(measurements_table, avgdata)
{
  dW = e*calf(measurements_table, avgdata) * S;
  return(W)
}


#Eq 5
# 
# D = rooting depth
calD = function(measurements_table, avgdata)
{
  
  Dsowing = measurements_table[measurements_table$Variables=='Dsowing', 'Value'];
  l0 = measurements_table[measurements_table$Variables=='l0', 'Value'];
  B0 = measurements_table[measurements_table$Variables=='B0', 'Value'];
  l0 = measurements_table[measurements_table$Variables=='l0', 'Value'];
  s = measurements_table[measurements_table$Variables=='s', 'Value'];
  T0 = measurements_table[measurements_table$Variables=='T0', 'Value'];
  T = cumsum(avgdata$HCAirTemperature);
  
  D = Dsowing + l0*exp((B0/s))*(1-exp(-s*(T - T0)));
  plot(T, D, xlab = 'T accummulated temperature', 
       ylab = 'Root depth', main='Model I', pch=21, bg='red');
  lines(tdat, predict(o), col='green'); 
  return(D);
  
}


# Eq 4
calv = function()
{
  v0 = measurements_table[measurements_table$Variables=='v0', 'Value'];
  v =v0(1+0.1(1-calfstress()));
  
}

# eq 3
calumin = function()
{
  µmin0 = measurements_table[measurements_table$Variables=='µmin0', 'Value'];
  umin = µmin0(2-calfstress());
  return(umin);
  
}

# Effect of Soil water stress
# Richter et al. Agricultural and forest, 109 (2001)
# Qfc = Daily available water content Qi et al
# AWC = Available water content at field capacity
# Eq 2 fstress
calfstress = function (DAW, AWC)
{
  Qfc = measurements_table[measurements_table$Variables=='Qfc', 'Value'];
  Qrel = measurements_table[measurements_table$Variables=='Qrel', 'Value'];
  Oth = 0.02; #(Sinclair (1986))
  fdt = 6:12; #(Richter)
  Qrel = AWC / Qfc;
  fstress = (2 / (1 + exp(-fdt * (Qrel - Oth)))) - 1;
  return(fstress);
  
}


# Calculate crop canopy cover
# Eq 1
calf = function(measurements_table, avgdata)
{
  f0 = measurements_table[measurements_table$Variables=='f0', 'Value'];
  umin = measurements_table[measurements_table$Variables=='µmin', 'Value'];
  T0 = measurements_table[measurements_table$Variables=='T0', 'Value'];
  u0 = measurements_table[measurements_table$Variables=='µ0', 'Value'];
  v = measurements_table[measurements_table$Variables=='v', 'Value'];
  T = cumsum(avgdata$HCAirTemperature);
  
  for(i in seq(-0.00001,0.00012, by=0.00001))
  {
    umin =i
    f = f0*exp(umin*(T-T0) + (((u0 - umin)/v) * (1-exp(-v*(T-T0)))));
    
    plot(T, f, xlab = 'T accummulated temperature', 
         ylab = 'Crop foliage cover', main='Model I', pch=21, bg='red');
    if(readline(i) =='q') { break()}
  }
  return(f);
  
}

#ydat = c(0.75, 0.73, 0.69, 0.57, 0.55, 0.44, 0.38, 0.26);
#tdat = c(1478.125, 912.5,821.875,734.375,734.375,528.125,465.625,443.75);

measurements_table = read.table('BroomsBarnParam.csv', header=TRUE, sep=',');

rawdata = read.table('BroomsBarnData.csv', sep=',', header=TRUE)
avgdata = averageValues(rawdata[, c('SolarRadiation', 'Precipitation','HCAirTemperatureOp1')], 
                        list(date=rawdata$date), rawdata);
colnames(avgdata)[2] = 'HCAirTemperature';

fc = calf(measurements_table, avgdata)
fm= fitYield(ydat, tdat);
