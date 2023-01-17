\name{BSMGreeks}
\alias{BSMGreeks}
\title{Greeks of European options according to Black-Scholes-Merton formula}
\usage{BSMGreeks(S, K, r, vol, time, q = 0, option_type)}
\arguments{
  \item{S}{Current stock price}
  \item{K}{Vector of strike prices}
  \item{r}{Risk-free interest rate (continuous)}
  \item{vol}{Volatility of the underlying stock (annulized)}
  \item{time}{Time to expiration of the option, in years}
  \item{q}{Continuous dividend rate}
  \item{option_type}{Option type: either "call" or "put".}
}
\value{
  Returns greeks for a list of European options different by strikes only according to Black-Scholes-Merton.
}
