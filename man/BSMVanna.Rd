\name{BSMVanna}
\alias{BSMVanna}
\title{Black-Scholes-Merton vanna of a European option}
\usage{BSMVanna(S, K, r, vol, time, q)}
\arguments{
  \item{S}{Current stock price}
  \item{K}{Strike price of the option}
  \item{r}{Risk-free interest rate (continuous)}
  \item{vol}{Volatility of the underlying stock (annulized)}
  \item{time}{Time to expiration of the option, in years}
  \item{q}{Continuous dividend rate}
}
\value{
  Returns the Black-Scholes-Merton vanna of a European option. Vanna is the same for put and call options.
}

