\name{BSMRhoPut}
\alias{BSMRhoPut}
\title{Black-Scholes-Merton Rho of a European put option}
\usage{BSMRhoPut(S, K, r, vol, time, q)}
\arguments{
  \item{S}{Current stock price}
  \item{K}{Strike price of the option}
  \item{r}{Risk-free interest rate (continuous)}
  \item{vol}{Volatility of the underlying stock (annulized)}
  \item{time}{Time to expiration of the option, in years}
  \item{q}{Continuous dividend rate}
}
\value{
  Returns the Black-Scholes-Merton Rho for put option price.
}

