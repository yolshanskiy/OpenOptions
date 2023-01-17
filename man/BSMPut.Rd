\name{BSMPut}
\alias{BSMPut}
\title{Black-Scholes-Merton Price of a European Put Option }
\usage{BSMPut(S, K, r, vol, time, q)}
\arguments{
  \item{S}{Current stock price}
  \item{K}{Strike price of the option}
  \item{r}{Risk-free interest rate (continuous)}
  \item{vol}{Volatility of the underlying stock (annulized)}
  \item{time}{Time to expiration of the option, in years}
  \item{q}{Continuous dividend rate}
}
\value{
  Returns the Black-Scholes-Merton put option price
}
\description{
  The BSMPut function calculates the Black-Scholes-Merton put option price using the Black-Scholes-Merton formula. The formula takes into account the current stock price, strike price, risk-free interest rate, volatility, time to expiration, and continuous dividend rate.
}
\examples{
  BSMPut(100, 110, 0.05, 0.2, 0.5, 0)
  # Returns the put option price for a stock with a current price of 100, strike price of 110, risk-free interest rate of 0.05, volatility of 0.2, time to expiration of 0.5 years, and continuous dividend rate of 0
}
