
library(OpenOptions)

#### Set Parameters ########

S = 200;
K = 196;
r = .1;
vol = .3
time = 2
q = .05

#### Test the functions ####

Call = BSMCall(S, K, r, vol, time, q)
Put = BSMPut(S, K, r, vol, time, q)
IV_Call = BSMIVCall(Call, S, K, r, time, q)
IV_Put = BSMIVPut(Put, S, K, r, time, q)

Delta_Call = BSMDeltaCall(S, K, r, vol, time, q)
Delta_Put = BSMDeltaPut(S, K, r, vol, time, q)

Epsilon_Call = BSMEpsilonCall(S, K, r, vol, time, q)
Epsilon_Put = BSMEpsilonPut(S, K, r, vol, time, q)

Theta_Call = BSMThetaCall(S, K, r, vol, time, q)
Theta_Put = BSMThetaPut(S, K, r, vol, time, q)

Rho_Call = BSMRhoCall(S, K, r, vol, time, q)
Rho_Put = BSMRhoPut(S, K, r, vol, time, q)

Charm_Call = BSMCharmCall(S, K, r, vol, time, q)
Charm_Put = BSMCharmPut(S, K, r, vol, time, q)

Vega = BSMVega(S, K, r, vol, time, q)
Gamma = BSMGamma(S, K, r, vol, time, q)
Vanna = BSMVanna(S, K, r, vol, time, q)

Lambda_Call = BSMLambdaCall(S, K, r, vol, time, q)
Lambda_Put = BSMLambdaPut(S, K, r, vol, time, q)

Lambda_Call - S/Call*Delta_Call
Lambda_Put - S/Put*Delta_Put

cat(paste0(
  "Call: ", round(Call,6) ,"\n",
  "Put: ", round(Put,6) ,"\n",
  "IV_Call: ", round(IV_Call,6) ,"\n",
  "IV_Put: ", round(IV_Put,6) ,"\n",
  "Delta_Call: ", round(Delta_Call,6) ,"\n",
  "Delta_Put: ", round(Delta_Put,6) ,"\n",
  "Rho_Call: ", round(Rho_Call,6) ,"\n",
  "Rho_Put: ", round(Rho_Put,6) ,"\n",
  "Theta_Call: ", round(Theta_Call,6) ,"\n",
  "Theta_Put: ", round(Theta_Put,6) ,"\n",
  "Epsilon_Call: ", round(Epsilon_Call,6) ,"\n",
  "Epsilon_Put: ", round(Epsilon_Put,6) ,"\n",
  "Lambda_Call: ", round(Lambda_Call,6) ,"\n",
  "Lambda_Put: ", round(Lambda_Put,6) ,"\n",
  "Charm_Call: ", round(Charm_Call,6) ,"\n",
  "Charm_Put: ", round(Charm_Put,6) ,"\n",
  "Vega: ", round(Vega,6), "\n",
  "Vanna: ", round(Vanna,6), "\n",
  "Gamma: ", round(Gamma,6), "\n"
))


#### Call Greeks #####

cat(paste0(
  "Call: ", round(Call,6) ,"\n",
  "Delta_Call: ", round(Delta_Call,6) ,"\n",
  "Gamma: ", round(Gamma,6), "\n",
  "Vega: ", round(Vega,6), "\n",
  "Vanna: ", round(Vanna,6), "\n",
  "Rho_Call: ", round(Rho_Call,6) ,"\n",
  "Theta_Call: ", round(Theta_Call,6) ,"\n",
  "Epsilon_Call: ", round(Epsilon_Call,6) ,"\n",
  "Lambda_Call: ", round(Lambda_Call,6) ,"\n",
  "Charm_Call: ", round(Charm_Call,6) ,"\n"
))

#### Put Greeks #####

cat(paste0(
  "Put: ", round(Put,6) ,"\n",
  "Delta_Put: ", round(Delta_Put,6) ,"\n",
  "Rho_Put: ", round(Rho_Put,6) ,"\n",
  "Theta_Put: ", round(Theta_Put,6) ,"\n",
  "Epsilon_Put: ", round(Epsilon_Put,6) ,"\n",
  "Lambda_Put: ", round(Lambda_Put,6) ,"\n",
  "Charm_Put: ", round(Charm_Put,6) ,"\n",
  "Vega: ", round(Vega,6), "\n",
  "Vanna: ", round(Vanna,6), "\n",
  "Gamma: ", round(Gamma,6), "\n"
))

#### Check Greek Letters for Call #########

dt = 1e-07

Diff = (BSMCall(S+dt, K, r, vol, time, q) - BSMCall(S, K, r, vol, time, q))/dt - Delta_Call
cat(paste0("Difference between Delta and -dV/dS: ", round(Diff,7), "\n"))

Diff = (BSMCall(S, K, r, vol + dt, time , q) - BSMCall(S, K, r, vol, time, q))/dt - Vega
cat(paste0("Difference between Gamma and -dV/dvol: ", round(Diff,7), "\n"))

Diff = (BSMDeltaCall(S + dt, K, r, vol, time , q) - BSMDeltaCall(S, K, r, vol, time, q))/dt - Gamma
cat(paste0("Difference between Gamma and -d^2V/dS^2: ", round(Diff,7), "\n"))

Diff = (BSMCall(S, K, r, vol, time - dt, q) - BSMCall(S, K, r, vol, time, q))/dt - Theta_Call
cat(paste0("Difference between Theta and -dV/dt: ", round(Diff,7), "\n"))

Diff = (BSMCall(S, K, r  + dt, vol, time , q) - BSMCall(S, K, r, vol, time, q))/dt - Rho_Call
cat(paste0("Difference between Rho and -dV/dr: ", round(Diff,7), "\n"))

Diff = (BSMCall(S  + dt, K, r, vol, time , q) - BSMCall(S, K, r, vol, time, q))/dt * S/BSMCall(S, K, r, vol, time, q) - Lambda_Call
cat(paste0("Difference between Lambda and -dlogV/dlogS: ", round(Diff,7), "\n"))

Diff = (BSMDeltaCall(S, K, r, vol + dt, time , q) - BSMDeltaCall(S, K, r, vol, time, q))/dt - Vanna
cat(paste0("Difference between Vanna and dDelta/dvol: ", round(Diff,7), "\n"))

Diff = (BSMDeltaCall(S, K, r, vol, time - dt , q) - BSMDeltaCall(S, K, r, vol, time, q))/dt - Charm_Call
cat(paste0("Difference between Charm and -dDelta/dt: ", round(Diff,7), "\n"))

Diff = (BSMCall(S, K, r, vol, time , q  + dt) - BSMCall(S, K, r, vol, time, q))/dt- Epsilon_Call
cat(paste0("Difference between Epsilon and -dV/dq: ", round(Diff,7), "\n"))



#### Check Greek Letters for Put #########

dt = 1e-07

Diff = (BSMPut(S+dt, K, r, vol, time, q) - BSMPut(S, K, r, vol, time, q))/dt - Delta_Put
cat(paste0("Difference between Delta and -dV/dS: ", round(Diff,7), "\n"))

Diff = (BSMPut(S, K, r, vol + dt, time , q) - BSMPut(S, K, r, vol, time, q))/dt - Vega
cat(paste0("Difference between Gamma and -dV/dvol: ", round(Diff,7), "\n"))

Diff = (BSMDeltaPut(S + dt, K, r, vol, time , q) - BSMDeltaPut(S, K, r, vol, time, q))/dt - Gamma
cat(paste0("Difference between Gamma and -d^2V/dS^2: ", round(Diff,7), "\n"))

Diff = (BSMPut(S, K, r, vol, time - dt, q) - BSMPut(S, K, r, vol, time, q))/dt - Theta_Put
cat(paste0("Difference between Theta and -dV/dt: ", round(Diff,7), "\n"))

Diff = (BSMPut(S, K, r  + dt, vol, time , q) - BSMPut(S, K, r, vol, time, q))/dt - Rho_Put
cat(paste0("Difference between Rho and -dV/dr: ", round(Diff,7), "\n"))

Diff = (BSMPut(S  + dt, K, r, vol, time , q) - BSMPut(S, K, r, vol, time, q))/dt * S/BSMPut(S, K, r, vol, time, q) - Lambda_Put
cat(paste0("Difference between Lambda and -dlogV/dlogS: ", round(Diff,7), "\n"))

Diff = (BSMDeltaPut(S, K, r, vol + dt, time , q) - BSMDeltaPut(S, K, r, vol, time, q))/dt - Vanna
cat(paste0("Difference between Vanna and dDelta/dvol: ", round(Diff,7), "\n"))

Diff = (BSMDeltaPut(S, K, r, vol, time - dt , q) - BSMDeltaPut(S, K, r, vol, time, q))/dt - Charm_Put
cat(paste0("Difference between Charm and -dDelta/dt: ", round(Diff,7), "\n"))

Diff = (BSMPut(S, K, r, vol, time , q  + dt) - BSMPut(S, K, r, vol, time, q))/dt- Epsilon_Put
cat(paste0("Difference between Epsilon and -dV/dq: ", round(Diff,7), "\n"))




#### Set Parameters ########

S = 200;
K = 196 + 1:8;
r = .1;
vol = .3
time = 2
q = .05

BSMGreeks(S,K,r,vol,time, q, option_type = "call")
BSMGreeks(S,K,r,vol,time, q, option_type = "put")


#### End ####
rm(list = ls())
