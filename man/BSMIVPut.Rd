\name{BSMIVPut}
\alias{BSMIVPut}
\title{Implied volatility of a European Put option according to Black-Scholes-Merton model}
\usage{BSMIVPut(OptionPrice, S, K, r, time, q, precision, vol_initial)}
\arguments{
  \item{OptionPrice}{Option Price}
  \item{S}{Current stock price}
  \item{K}{Strike price of the option}
  \item{r}{Risk-free interest rate (continuous)}
  \item{time}{Time to expiration of the option, in years}
  \item{q}{Continuous dividend rate}
  \item{precision}{Targeted precision (tolerance) of the answer}
  \item{vol_initial}{The initial IV point}
}
\value{
  Returns implied volatility of a European Put option according to Black-Scholes-Merton model.
}

