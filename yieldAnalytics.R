rm(list=ls()); # Delete files
cat("\014");
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
                   control=list(all.methods=TRUE, save.failures=TRUE, maxit=10000))
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
PredYieldDummy = function(ydat, tdat)
{
  wdata = data.frame(y=ydat, x=tdat);
  eunsc = ydat ~ A + B*tdat + C*(exp(-v*tdat)); # Equation
  
  # Transform to linear
  # log(y -A)  = log(Bt) + log(C) - vt;
  
  A.0 = min(wdata$y) * 0.5;
  model.0 = lm(log(y - A.0) ~ x, data=wdata);
  start1 = c(A = A.0, B = 0, C = exp(coef(model.0)[[1]]), v=coef(model.0)[[2]]); 
  
  # v is K in Werker
  optx = optima(start1, wdata, rhsfit);
  o = nlsLM(eunsc, 
            start=list(A = start1[1], B = start1[2], C = start1[3], v = start1[4]), trace = TRUE);
  
  A = coef(o)[1]; B = coef(o)[2]; C = coef(o)[3]; v = coef(o)[4] 
  
  
  plot(tdat, ydat, xlab = 'T accummulated temperature', 
       ylab = 'foliage cover', main='Model I', pch=21, bg='red', xlim=c(0, 2500), ylim=c(0, 1), col='red');
  points(tdat, predict(o), col='green')
  
  x = 0:2500;
  yeq = ydat ~ A + B*x + C*(exp(-v*x));
  y = eval(parse(text=yeq))
  lines(x,y);
  
  #get from linear
  umin = B;
  tmax = (1/v)*log((v*C)/B);
  log(ymax) = A + umin*((1/v) + tmax);
  ymax = log(A) + log(umin*((1/v) + tmax));
  #
  
  y0 = 0.000015;
  t0 =  100;
  umin = -0.0050;
  k = 0.003;
  u0 = 0.005;
  ymax = 0.90;
  tmax = t0 - (1/(k*log(-umin/(u0-umin))));
  
  y = y0*exp(umin*(tdat-t0) + ( ((u0-umin)/k) * (1-exp(-k*(tdat-t0))) ) );
  y1 = ymax*exp(umin*(tdat-tmax) - ((umin/k)*(1-exp(-k*(tdat-tmax)))));
  
  plot(tdat, y1)
  return(o)
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
  start1 = c(A = A.0, B = 0, C = exp(coef(model.0)[[1]]), v=coef(model.0)[[2]]); 
  
  # v is K in Werker
  optx = optima(start1, wdata, rhsfit);
  o = nlsLM(eunsc, 
            start=list(A = start1[1], B = start1[2], C = start1[3], v = start1[4]), trace = TRUE);
  
  A = coef(o)[1]; B = coef(o)[2]; C = coef(o)[3]; v = coef(o)[4] 
  
  umin = B;
  tmax = 1/v*log((v*C)/D);
  
  
  plot(tdat, ydat, xlab = 'T accummulated temperature', 
       ylab = 'foliage cover', main='Model I', pch=21, bg='red', xlim=c(0, 2500), ylim=c(0, 1), col='red');
  points(tdat, predict(o), col='green')
  
  x = 0:2500;
  yeq = ydat ~ A + B*x + C*(exp(-v*x));
  y = eval(parse(text=yeq))
  lines(x,y);
  
  y0 = 0.000015;
  t0 =  100;
  umin = -0.0050;
  k = 0.003;
  u0 = 0.005;
  ymax = 0.90;
  tmax = t0 - (1/(k*log(-umin/(u0-umin))));
  
  y = y0*exp(umin*(tdat-t0) + ( ((u0-umin)/k) * (1-exp(-k*(tdat-t0))) ) );
  y1 = ymax*exp(umin*(tdat-tmax) - ((umin/k)*(1-exp(-k*(tdat-tmax)))));
  
  plot(tdat, y1)
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
calcQ = function(Qfc, SMD, D)
{
  
  Q = Qfc - (SMD / D);
  return(Q);  
}

# Eq 19
# SSE = 
calcSMD = function(SMD, SSE, Ea, R)
{
    print(SSE)
    SMD = SMD + SSE + Ea - R;
    return(SMD);
}



# Eq 17
# psi_soil
calcpsi_soil = function(measurements_table, Q, Qfc)
{
  b = measurements_table[measurements_table$Variables=='b', 'Value'];
  psi_soil = -5*((Q/Qfc)^(-b));
  return(psi_soil);
  
}

#Eq 16
calcQfc = function(measurements_table)
{
  a1 = measurements_table[measurements_table$Variables=='a1', 'Value'];
  a2 = measurements_table[measurements_table$Variables=='a2', 'Value'];
  b = measurements_table[measurements_table$Variables=='b', 'Value'];
  Qfc = a2 * (a1/5) ^ (1/b)
  return(Qfc);
  
}

calcQpwp = function(measurements_table)
{
  a1 = measurements_table[measurements_table$Variables=='a1', 'Value'];
  a2 = measurements_table[measurements_table$Variables=='a2', 'Value'];
  b = measurements_table[measurements_table$Variables=='b', 'Value'];
  Qpwp = a2 * (a1/1500) ^ (1/b);
  print(Qpwp)
  return(Qpwp);
  
}

#Qrel = AWC / Qfc;
calcQrel = function(Q, Qfc, Qpwp)
{
  
  Qrel = (Q - Qpwp) / (Qfc - Qpwp); 
  #Qrel =  Q / Qfc;
  return(Qrel);
  
}


#Eq14_15
calcEmax = function(measurements_table, psi_soil, Qfc, D, Q)
{
  psi_crop = measurements_table[measurements_table$Variables=='psi_crop', 'Value'];
  c1 = measurements_table[measurements_table$Variables=='c1', 'Value'];
  c2 = measurements_table[measurements_table$Variables=='c2', 'Value'];
  b = measurements_table[measurements_table$Variables=='b', 'Value'];
  
  Emax = (psi_soil - psi_crop) / (c1+c2/D*( ( (Q/Qfc)^(-1*(2*b+3)) ) -1 ) );
  return(Emax)
  
}


#Eq 13
# DAily crop actual transpiration
calcEa = function(Ep, Emax)
{
  Ea = min(Ep, Emax);
  return(Ea);
  
}

# Eq 12
# Daily Crop potential transpiration
# f = 
# ETgrass = Daily potential evapotranspiration
calcEp = function(measurements_table, avgdata)
{
  
  ET_grass = measurements_table[measurements_table$Variables=='ET_grass', 'Value'];
  Ep = (1.25*avgdata) * ET_grass;
  return(Ep);
}


# Eq 11
# SSE
calcSSE = function(measurements_table, avgdata)
{
  
  SSE = 1.5*(1 - avgdata);
  return(SSE);
}


#Eq 10
# Net crop yield
calcY = function(measurements_table, avgdata)
{
  k = measurements_table[measurements_table$Variables=='k', 'Value'];
  
  Y = W - (1/k)*log((k*W) + 1);
  return(Y);
  
}


#Eq 7
# Intercepted radiation use efficiency
calce = function(measurements_table, avgdata)
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
calcW = function(measurements_table, avgdata)
{
  dW = e*calf(measurements_table, avgdata) * S;
  return(W)
}


#Eq 5
# 
# D = rooting depth

calcD = function(measurements_table, T)
{
  Dsowing = measurements_table[measurements_table$Variables=='Dsowing', 'Value'];
  l0 = measurements_table[measurements_table$Variables=='l0', 'Value'];
  B0 = measurements_table[measurements_table$Variables=='B0', 'Value'];
  s = measurements_table[measurements_table$Variables=='s', 'Value'];
  T0 = measurements_table[measurements_table$Variables=='T0', 'Value'];
  T = T;
  
  D = Dsowing + l0*exp((B0/s))*(1-exp(-s*(T - T0)));
  # plot(T, D, xlab = 'T accummulated temperature', 
  #      ylab = 'Root depth', main='Model I', pch=21, bg='red');
  # lines(tdat, predict(o), col='green'); 
  return(D);
  
}


# Eq 4
calcv = function(measurements_table, fstress)
{
  v0 = measurements_table[measurements_table$Variables=='v0', 'Value'];
  v = v0*(1+0.1*(1 - fstress));
  return(v);
  
}

# eq 3
calcumin = function(measurements_table, fstress)
{
  umin0 = measurements_table[measurements_table$Variables=='umin0', 'Value'];
  umin = umin0*(2 - fstress);
  return(umin);
  
}

# Effect of Soil water stress
# Richter et al. Agricultural and forest, 109 (2001)
# Qfc = Daily available water content Qi et al
# AWC = Available water content at field capacity
# Eq 2 fstress
calcfstress = function (fdt, Qrel)
{
  Oth = 0.02; #(Sinclair (1986))
  fstress = (2 / (1 + exp(-fdt * (Qrel - Oth)))) - 1;
  
  return(fstress);
  
}


# Calculate crop canopy cover
# Eq 1
calcf = function(measurements_table, T, umin, v)
{
  
  f0 = measurements_table[measurements_table$Variables=='f0', 'Value'];
  #umin = measurements_table[measurements_table$Variables=='µmin', 'Value'];
  T0 = measurements_table[measurements_table$Variables=='T0', 'Value'];
  u0 = measurements_table[measurements_table$Variables=='µ0', 'Value'];
  #v = measurements_table[measurements_table$Variables=='v', 'Value'];
  
  f = f0*exp(umin*(T-T0) + (((u0 - umin)/v) * (1-exp(-v*(T-T0)))));
  return(f);
  
}

# calculate Fdt
calcfdt = function(T)
{
  
  if(T <= 1700) { return(6)}
  else{
    return(12);
  }

}


#Determination of solar radiation from temperature data
# Rs = Krs*sqrt(Tmax-Tmin)*Ra
# Krs for coastal locations is 0.19
calcRs = function()
{
  Rs =  0.75*min((0.19*sqrt(Tmax-Tmin)*Ra),Rso);
  
}


calcRa = function()
{
  
  Gsc = 0.0820;
  n = 246;# number of days of a target day.
  d = 20;# degrees;
  
  delta =  0.409*sin(2*pi*(n)/365 - 1.39);
  dr = dr = 1 + 0.033*cos(2*pi*(n)/365);
  phi = (pi/180)*(-d);
  Ra = ((12*60)/pi)*Gsc*dr*((omega2-omega1)*sin(phi)*sin(delta) + 
                             cos(phi)*cos(delta) * (sin(omega2)-sin(omega1)));
  
}

# Calculate all equations
calAllEq = function(measurements_table, avgdata)
{

  v = list();
  avgdataCum = data.frame(date = avgdata$date, SolarRadiation = cumsum(avgdata$SolarRadiation),
                          Precipitation = cumsum(avgdata$Precipitation), 
                          HCAirTemperatureOp1 = cumsum(avgdata$HCAirTemperatureOp1));
  avgdataCum = avgdataCum[5:nrow(avgdataCum),];
  avgdataCum = avgdataCum[order(c(as.Date(avgdataCum$date, format='%d/%m/%Y'))),];
  avgdataCum$date = as.Date(avgdataCum$date, format='%d/%m/%Y');
  
  T0 = measurements_table[measurements_table$Variables=='T0', 'Value'];
  
  avgdataCum$HCAirTemperatureOp1[1] = 0;
  avgdataCum$Precipitation[1] = 0;
  
  v$date = avgdataCum$date;
  
  for(i in 1:nrow(avgdataCum))
  {
    #SMD = -i
    # Don't depend on anything else
    v$thermaltime[i] = avgdataCum$HCAirTemperatureOp1[i] + T0;
    v$rainfall[i] = avgdataCum$Precipitation[i] + 3.00;
    v$D[i] = calcD(measurements_table, v$thermaltime[i]); # rooting depth; Eq 5
    v$Qfc[i] = calcQfc(measurements_table); # Qfc Eq 16
    v$Qpwp[i] = calcQpwp(measurements_table);
    v$fdt[i] = calcfdt(v$thermaltime[i]);
    
    if(i == 1) { v$Q[i] = v$Qfc[i]}
    else{
        v$Q[i] = calcQ(v$Qfc[i], v$SMD[i], v$D[i]); # Eq 18 Q = Qfc
    }
    
    v$Qrel[i] = calcQrel(v$Q[i], v$Qfc[i], v$Qpwp[i]); #Qrel
    
      
    # Depend on other values
    
    if(i == 1) { v$SMD[i] = 0}
    else{
      
      v$SMD[i] = calcSMD(v$SMD[i-1], v$SSE[i], v$Ea[i], v$rainfall[i]); # Eq 12
    }
    
    # v$Q[i] = calcQ(v$Qfc[i], v$SMD[i], v$D[i]); # Eq 18
    # v$Qrel[i] = calcQrel(v$Q[i], v$Qfc[i]); #Qrel
    # v$psi_soil[i] = calcpsi_soil(measurements_table, v$Q[i], v$Qfc[i]); # Eq 17
    v$Emax[i] = calcEmax(measurements_table, v$psi_soil[i], v$Qfc[i], v$D[i], v$Q[i]); # Eq 14 and 15
    v$Ea[i] = calcEa(v$Ep[i], v$Emax[i]); # Eq 13
    
    # Qrel fro previous circle
    v$fstress[i] = calcfstress(v$fdt[i], v$Qrel[i] ); # Eq 2
    
    v$v[i] = calcv(measurements_table, v$fstress[i]); # Eq 4
    v$umin[i] =  calcumin(measurements_table, v$fstress[i]); # Eq 3
    
    if( i == 1 ) { v$f[i] = 0;}
    else{
    v$f[i] = calcf(measurements_table, v$thermaltime[i], v$umin[i], v$v[i]); # Eq 1
    }
    
    v$SSE[i] = calcSSE(measurements_table, v$f[i]); # Eq 11
    v$Ep[i] = calcEp(measurements_table, v$f[i]); # Eq 12
  }
  
  
  #plot(T, f, xlab = 'T accummulated temperature', 
  #ylab = 'Crop foliage cover', main='Model I', pch=21, bg='red');
  
  return(data.frame(v));
}

#ydat = c(0.75, 0.73, 0.69, 0.57, 0.55, 0.44, 0.38, 0.26);
#tdat = c(1478.125, 912.5,821.875,734.375,734.375,528.125,465.625,443.75);
# 1985
ydat = c(0.058064516,	0.090322581,	0.341935484,	0.6,	0.858064516,	0.929032258,	0.929032258)
tdat = c(318.8405797,	478.2608696,	579.7101449,	695.6521739,	927.536,	1275.362,	1579.710);


measurements_table = read.table('BroomsBarnParam.csv', header=TRUE, sep=',');

rawdata = read.table('BroomsBarnData.csv', sep=',', header=TRUE)
avgdata1 = averageValues(rawdata[, c('HCAirTemperatureOp1')], 
                        list(date=rawdata$date), funtype='mean');

avgdata2 = aggregate(rawdata[, c('SolarRadiation', 'Precipitation')], 
                         list(date=rawdata$date), 'sum');

avgdata = merge(avgdata1, avgdata2, by='date');
rm(avgdata1, avgdata2);
colnames(avgdata)[2] = 'HCAirTemperatureOp1';


#avgdata = averageValues(rawdata[, c('SolarRadiation', 'Precipitation','HCAirTemperatureOp1')], 
 #                       list(date=rawdata$date), rawdata);

#fc = calf(measurements_table, avgdata)
#fm = PredYield(ydat, tdat);
#fm = PredYieldDummy(ydat, tdat);
m = calAllEq(measurements_table, avgdata)


