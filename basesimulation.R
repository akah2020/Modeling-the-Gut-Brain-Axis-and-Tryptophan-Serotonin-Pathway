library(deSolve)   # load the ode solver library

## Define the right-hand side of the ode
## Note that the function's arguments must be in this order, although
## you can call them other names if you wish
food = function(t){
   #just to try something out: on for first two time units, then off 
  #if (t<=2){ food = 1}
  #else { food = 0}
  #food = (sin(0.8*t)+1)^2 # 8 hour period; k = 2*pi/8 = 0.75
  #what if we ate at 7 am, 12 pm, and 6 pm over three days with increasing meal sizes as the day goes on (2, 4, and 6)
  B = 2
  L = 4
  D = 6
  tt = t %% 24
food = ((tt>=7)&(tt<7.5))*B + ((tt>=12)&(tt<12.5))*L + ((tt>=18)&(tt<19))*D
 
}

BBBode = function(t,x,parms){
  d1 = parms[1]
  d2 = parms[2]
  d3 = parms[3]
  d4 = parms[4]
  d5 = parms[5]
  d6 = parms[6]
  K1 = parms[7]
  K2 = parms[8]
  K3 = parms[9]
  K4 = parms[10]
  K5 = parms[11]
  Ki = parms[12]
  p1 = parms[13]
  p2 = parms[14]
  V1 = parms[15]
  V2 = parms[16]
  V3 = parms[17]
  V4 = parms[18]
  TPHB = parms[19]
  TPHG = parms[20]
  TB = x[1]
  TG = x[2]
  PB = x[3]
  PG = x[4]
  HB = x[5]
  HG = x[6]
  s = food(t) # describes food going in as a function of time
  
  dTBdt = V1*TG/(K1+TG) - (d2*TB) - (((0.01*V3*TB*TPHB)/((K3+TB+(TB^2/Ki))*(K4 + TPHB)))*(1.5-(HB^2/(0.0007682 + HB^2))))
  dTGdt = s - (V1*TG)/(K1 + TG) - (d1*TG) - (((0.02*V2*TG*TPHG)/((K2+TG+(TG^2/Ki))*(K4 + TPHG)))*(1.5-(HG^2/(0.0007682 + HG^2))))
  dPBdt = (((0.01*V3*TB*TPHB)/((K3+TB+(TB^2/Ki))*(K4 + TPHB)))*(1.5-(HB^2/(0.0007682 + HB^2)))) - ((V4*PB)/(K5+PB))
  dPGdt = (((0.02*V2*TG*TPHG)/((K2+TG+(TG^2/Ki))*(K4 + TPHG)))*(1.5-(HG^2/(0.0007682 + HG^2)))) - ((V4*PG)/(K5+PG))
  dHBdt = ((V4*PB)/(K5+PB)) - d5*HB
  dHGdt = ((V4*PG)/(K5+PG)) - d5*HG
  
  # make 0.02 param -- Overall, about 2 percent of ingested tryptophan is used for the synthesis of serotonin (Best 2010)
  # 0.01 of trp is used for serotonin, other trp used for other reactions (Richard 2009)
  # make all numbers into parameters 
  
  xdot = c(dTBdt,dTGdt, dPBdt, dPGdt, dHBdt, dHGdt)
  list(xdot)   
}

## Set the parameters
d1 = 0.277 # we got this number using half life? the half life of tryptophan is 2 hrs?
d2 = 1
d3 = 1
d4 = 1
d5 = 1
d6 = 1
K1 = 64
K2 = 40
K3 = 40
K4 = 20
K5 = 160
Ki = 1000
p1 = 1
p2 = 1
V1 = 400*0.01 # trying to get serotonin in gut 9 times serotonin in brain since 90% of serotonin synthesis is done in gut
V2 = 400
V3 = 400
V4 = 400

TPHB = p2/d3
TPHG = p1/d4
# USE STEADY STATE FOR TPHB and TPHG

# put them in a vector
parms = c(d1, d2, d3, d4, d5, d6, K1, K2, K3, K4, K5, Ki, p1, p2, V1, V2, V3, V4, TPHB, TPHG)

## Set the initial condition
TB0 = 0   # change this? start at where it is at end of the day
TG0 = 0
PB0 = 0
PG0 = 0
HB0 = 0
HG0 = 0
## Set the time interval for the integration
t0 = 0
tf = 24 #multiple days 
dt =.001

## Solve the ode
sol = ode(func = BBBode, y = c(TB0,TG0,PB0,PG0,HB0,HG0), times = seq(t0,tf,by = dt),parms = parms)
## The first column of the solution is the time values, the next columns are the state variables

# run days
L = length(sol[,1])
Xinit=sol[L,2:7]
tf = 72
sol = ode(func = BBBode, y = Xinit, times = seq(t0,tf,by = dt),parms = parms)
# use this as initial values 


## Plot the first solution
plot(sol[,1],sol[,2],type = 'l', lwd = 6, col = 'blue',
     xlab = 'Time',ylab = 'x(t)',ylim=c(0,7),
     main = 'Tryptophan in the Gut and Brain (Base)')
lines(sol[,1],sol[,3], lwd = 6, col='red')
legend('topright', legend = c('Gut','Brain'),col = c('red','blue'), lwd = 6, lty = c(1,1))

## Plot the second solution
plot(sol[,1],sol[,4],type = 'l', lwd = 6, col = 'blue',
     xlab = 'Time',ylab = 'x(t)',ylim=c(0,0.02),
     main = '5-HTP in the Gut and Brain (Base)')
lines(sol[,1],sol[,5], lwd = 6, col='red')
legend('topright', legend = c('Gut','Brain'),col = c('red','blue'), lwd = 6, lty = c(1,1))

## Plot the third solution
plot(sol[,1],sol[,6],type = 'l', lwd = 6, col = 'blue',
     xlab = 'Time',ylab = 'x(t)',ylim=c(0,0.05),
     main = 'Serotonin in the Gut and Brain (Base)')
lines(sol[,1],sol[,7], lwd = 6, col='red')
legend('topright', legend = c('Gut','Brain'),col = c('red','blue'), lwd = 6, lty = c(1,1))
# observe that in our model the brain is 10 times smaller
# is the gut sufficiently ahead of the brain, by how much, is this time difference reasonable 

# how much trp in the brain is used to make serotonin as opposed to other things, even relative to gut
# justify changes (ie 0.01) in front of V3
# determining d and p params
# how much trp gets from gut to brain, how much is used up

# can we mimic some kinds of pathologies that they describe
# different diets, conditions
## what happens if we change some of the parameters we guessed 
# look at other things we modeled, understand them better
