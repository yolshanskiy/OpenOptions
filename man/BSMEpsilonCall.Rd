\name{BSMEpsilonCall}
\alias{BSMEpsilonCall}
\title{Black-Scholes-Merton Epsilon of a European call option}
\usage{BSMEpsilonCall(S, K, r, vol, time, q)}
\arguments{
  \item{S}{Current stock price}
  \item{K}{Strike price of the option}
  \item{r}{Risk-free interest rate (continuous)}
  \item{vol}{Volatility of the underlying stock (annulized)}
  \item{time}{Time to expiration of the option, in years}
  \item{q}{Continuous dividend rate}
}
\value{
  Returns the Black-Scholes-Merton Epsilon for call option price.
}

