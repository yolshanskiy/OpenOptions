library(OpenOptions)

#### Set Parameters and Check Consistency with BSM ########

S = 810;
K = 800;
r = .05;
vol = .2
time = 0.5
q = .1
N = 5  ## start with a small tree example

#### Calculate All State Prices and Build Stock Price/Dividends tree #######
#### Note that depending on application we might have different assumptions
## about the stock-price process:
## The first way is to assume that the value of the stock has volatility sigma,
## and that it pays dividends after out of the value:
## set nonstandard_dividend_payment = TRUE for this approach.
## the more common approach is to assume that stock price has an increment with volatility sigma and then
## pays approximately S*qdt (more specifically S*(exp(qdt) - 1))
## set nonstandard_dividend_payment = FALSE for this approach.

StatePrices = binomialStatePrices(S,K,r,vol,time,Div = q, N,verbose = FALSE, nonstandard_dividend_payment = FALSE)
StockPaths  = binomialStockPaths(S,K,r,vol,time,Div = q, N, nonstandard_dividend_payment = FALSE)

StatePrices2 = binomialStatePrices(S,K,r,vol,time,Div = q, N, verbose = FALSE, nonstandard_dividend_payment = TRUE)
StockPaths2  = binomialStockPaths(S,K,r,vol,time,Div = q, N, nonstandard_dividend_payment = TRUE)
#### Check that state prices add up to 1 after discounting  ####
#### Check State Prices ####

colSums(StatePrices, na.rm = T) * exp( (0:N/N)*time*r)
colSums(StatePrices2, na.rm = T) * exp( (0:N/N)*time*r)

#### Check that 1-period tree is priced correctly: ####
col = 2
ValueStates = (StockPaths$Price + StockPaths$Dividends)[,col] * StatePrices[,col]
sum(ValueStates, na.rm = T)

ValueStates = (StockPaths2$Price + StockPaths2$Dividends)[,col] * StatePrices2[,col]
sum(ValueStates, na.rm = T)


#### Check that 2-period tree is priced correctly ####

check_correct_pricing <- function(StockPaths, StatePrices, period = 2){
  col = period + 1
  ValueTerminalStates = (StockPaths$Price + StockPaths$Dividends)[,col] * StatePrices[,col]
  PriceOfDividends =   sum(StockPaths$Dividends[,2:(col-1)] * StatePrices[,2:(col-1)], na.rm = T)
  sum(ValueTerminalStates, na.rm = T) + PriceOfDividends
}



#### Check that N-period tree ####
col = N + 1
check_correct_pricing (StockPaths, StatePrices, period = 5)
check_correct_pricing (StockPaths2, StatePrices2, period = 5)

## price by binomial-tree approach
binomialOptionPricer_detailed(S,K,r,vol,time,Div = q, N,verbose = FALSE, nonstandard_dividend_payment = F)
binomialOptionPricer_detailed(S,K,r,vol,time,Div = q, N,verbose = FALSE, nonstandard_dividend_payment = T)

## price by binomial-tree approach
binomialOptionPricer(S,K,r,vol,time,Div = q, N,verbose = FALSE, nonstandard_dividend_payment = F)
binomialOptionPricer(S,K,r,vol,time,Div = q, N,verbose = FALSE, nonstandard_dividend_payment = T)

## price by binomial-tree approach
binomialOptionPricer_fast(S,K,r,vol,time,Div = q, N,verbose = FALSE, nonstandard_dividend_payment = F)
binomialOptionPricer_fast(S,K,r,vol,time,Div = q, N,verbose = FALSE, nonstandard_dividend_payment = T)


## price by state prices assuming dividends paid out of the stock
Price = StockPaths$Price; SP = StatePrices
sum(SP[,col] * pmax(Price[,col] - K, 0), na.rm = T)

## price by state prices assuming dividends do not affect price dynamic
Price = StockPaths2$Price; SP = StatePrices2
sum(SP[,col] * pmax(Price[,col] - K, 0), na.rm = T)

#### Price European Call Option #####

N = 252*2*4; ## update position 4 times per day *2*4
col = N + 1
## standard dividends
StatePrices = binomialStatePrices(S,K,r,vol,time,Div = q, N,verbose = FALSE, nonstandard_dividend_payment = FALSE)
StockPaths  = binomialStockPaths(S,K,r,vol,time,Div = q, N, nonstandard_dividend_payment = FALSE)
## not standard dividends
StatePrices2 = binomialStatePrices(S,K,r,vol,time,Div = q, N, verbose = FALSE, nonstandard_dividend_payment = TRUE)
StockPaths2  = binomialStockPaths(S,K,r,vol,time,Div = q, N, nonstandard_dividend_payment = TRUE)

## price by state prices assuming dividends paid out of the stock
sum(StatePrices[,col] * pmax(StockPaths$Price[,col] - K, 0), na.rm = T)
sum(StatePrices[,col] * pmax(K - StockPaths$Price[,col], 0), na.rm = T)
## price by state prices assuming dividends do not affect price dynamic
sum(StatePrices2[,col] * pmax(StockPaths2$Price[,col] - K, 0), na.rm = T)
sum(StatePrices2[,col] * pmax(K - StockPaths2$Price[,col], 0), na.rm = T)

## price by binomial-tree approach
binomialOptionPricer(S,K,r,vol,time,Div = q, N,verbose = FALSE, nonstandard_dividend_payment = F)
binomialOptionPricer(S,K,r,vol,time,Div = q, N,verbose = FALSE, nonstandard_dividend_payment = T)
## price by binomial-tree approach
binomialOptionPricer(S,K,r,vol,time,Div = q, N,verbose = FALSE, nonstandard_dividend_payment = F, type = "put")
binomialOptionPricer(S,K,r,vol,time,Div = q, N,verbose = FALSE, nonstandard_dividend_payment = T, type = "put")

## price by binomial-tree approach
binomialOptionPricer_fast(S,K,r,vol,time,Div = q, N,verbose = FALSE, nonstandard_dividend_payment = F, robust_warning = TRUE)
binomialOptionPricer_fast(S,K,r,vol,time,Div = q, N,verbose = FALSE, nonstandard_dividend_payment = T, robust_warning = TRUE)
## price by binomial-tree approach
binomialOptionPricer_fast(S,K,r,vol,time,Div = q, N,verbose = FALSE, nonstandard_dividend_payment = F, robust_warning = TRUE, type = "put")
binomialOptionPricer_fast(S,K,r,vol,time,Div = q, N,verbose = FALSE, nonstandard_dividend_payment = T, robust_warning = TRUE, type = "put")


## price by BSM
BSMCall(S,K,r,vol,time,q)
BSMPut(S,K,r,vol,time,q)

#### Visualize Convergence ########

library(data.table)
library(ggplot2)
library(dplyr)
N_MAX = 252
dt =
  data.table(
    N = 1:N_MAX,
    `Non-Standard` = sapply(1:N_MAX, function(N) binomialOptionPricer_fast(S,K,r,vol,time,Div = q, N,verbose = FALSE, nonstandard_dividend_payment = F, robust_warning = F)),
    Standard = sapply(1:N_MAX, function(N) binomialOptionPricer_fast(S,K,r,vol,time,Div = q, N,verbose = FALSE, nonstandard_dividend_payment = T, robust_warning = F))
  )

dt %>%
  melt(id = "N") %>%
  ggplot(aes(N, value, colour = variable)) +
  geom_point() +
  geom_hline(yintercept = BSMCall(S,K,r,vol,time,q), size = 1, linetype = 10) +
  theme_bw() + labs(x = "Number Periods", y = "Option Value", colour = "Type Dividend")

#### American Options ######

## Europeans:
binomialOptionPricer(S,K,r,vol,time,Div = q, N,verbose = FALSE, nonstandard_dividend_payment = F)
binomialOptionPricer(S,K,r,vol,time,Div = q, N,verbose = FALSE, nonstandard_dividend_payment = T)
## depending on dividends rule the answer can differ:
binomialAmericanOptionPricer(S,K,r,vol,time,Div = q, N,verbose = FALSE, nonstandard_dividend_payment = F)
binomialAmericanOptionPricer(S,K,r,vol,time,Div = q, N,verbose = FALSE, nonstandard_dividend_payment = T)


## Europeans:
binomialOptionPricer(S,K,r,vol,time,Div = q, N,verbose = FALSE, nonstandard_dividend_payment = F, type = "put")
binomialOptionPricer(S,K,r,vol,time,Div = q, N,verbose = FALSE, nonstandard_dividend_payment = T, type = "put")
## depending on dividends rule the answer can differ:
binomialAmericanOptionPricer(S,K,r,vol,time,Div = q, N,verbose = FALSE, nonstandard_dividend_payment = F, type = "put")
binomialAmericanOptionPricer(S,K,r,vol,time,Div = q, N,verbose = FALSE, nonstandard_dividend_payment = T, type = "put")


#### Check deep In the money options #####
K = 0.1;  N =252;
## deep in the money option value must be close to:
S*exp(-q*time)-K*exp(-r*time)
##
BSMCall(S,K,r,vol,time,q)
## price by binomial-tree approach
binomialOptionPricer_fast(S,K,r,vol,time,Div = q, N,verbose = FALSE, nonstandard_dividend_payment = F)
binomialOptionPricer_fast(S,K,r,vol,time,Div = q, N,verbose = FALSE, nonstandard_dividend_payment = T)

## Put ##
K = 5000;
## deep in the money option value must be close to:
K*exp(-r*time) - S*exp(-q*time)
##
BSMPut(S,K,r,vol,time,q)
## price by binomial-tree approach
binomialOptionPricer_fast(S,K,r,vol,time,Div = q, N,verbose = FALSE, nonstandard_dividend_payment = F, type = "put")
binomialOptionPricer_fast(S,K,r,vol,time,Div = q, N,verbose = FALSE, nonstandard_dividend_payment = T, type = "put")

#### End ####
rm(list = ls())
