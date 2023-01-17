#include <Rcpp.h>
using namespace Rcpp;

// Black-Scholes-Merton call option function
// [[Rcpp::export]]
double BSMCall(double S, double K, double r, double vol, double time, double q = 0 ) {
  double d1, d2, ln;
  ln = log(S / K);
  d1 = (ln + (r - q + (vol * vol) / 2) * time) / (vol * sqrt(time));
  d2 = d1 - vol * sqrt(time);
  return exp(-q * time) * S * R::pnorm(d1, 0.0, 1.0, true, false) - K * exp(-r * time) * R::pnorm(d2, 0.0, 1.0, true, false);
}

// Black-Scholes-Merton put option function
// [[Rcpp::export]]
double BSMPut(double S, double K, double r, double vol, double time, double q = 0) {
  double d1, d2, ln;
  ln = log(S / K);
  d1 = (ln + (r - q + (vol * vol) / 2) * time) / (vol * sqrt(time));
  d2 = d1 - vol * sqrt(time);
  return -exp(-q * time) * S * R::pnorm(-d1, 0.0, 1.0, true, false) + K * exp(-r * time) * R::pnorm(-d2, 0.0, 1.0, true, false);
}

// Implied volatility call option function
// [[Rcpp::export]]
double BSMIVCall(double OptionPrice, double S, double K, double r, double time, double q = 0, double precision = 1E-10) {
  double sig = 0.2;
  double sig_old, ln, d1, d2, tmp_p;
  int iter = 0, MAXITER = 100000;
  do {
    sig_old = sig;
    ln = log(S / K);
    d1 = (ln + (r - q + (sig * sig) / 2) * time) / (sig * sqrt(time));
    d2 = d1 - sig * sqrt(time);
    tmp_p = exp(-q * time) * S * R::pnorm(d1, 0.0, 1.0, true, false) - K * exp(-r * time) * R::pnorm(d2, 0.0, 1.0, true, false);
    sig = sig + (OptionPrice - tmp_p) / (S * (1 / (sqrt(2 * M_PI))) * exp(-0.5 * (d1 * d1)) * sqrt(time));
  } while (fabs(sig - sig_old) > precision and iter < MAXITER);

  return sig;
}

// Implied volatility put option function
// [[Rcpp::export]]
double BSMIVPut(double OptionPrice, double S, double K, double r, double time, double q = 0, double precision = 1E-10) {
  double sig = 0.2;
  double sig_old, ln, d1, d2, tmp_p;

  do {
    sig_old = sig;
    ln = log(S / K);
    d1 = (ln + (r - q + (sig * sig) / 2) * time) / (sig * sqrt(time));
    d2 = d1 - sig * sqrt(time);
    tmp_p = -exp(-q * time) * S * R::pnorm(-d1, 0.0, 1.0, true, false) + K * exp(-r * time) * R::pnorm(-d2, 0.0, 1.0, true, false);
    sig = sig + (OptionPrice - tmp_p) / (S * (1 / (sqrt(2 * M_PI))) * exp(-0.5 * (d1 * d1)) * sqrt(time));
  } while (fabs(sig - sig_old) > precision);

  return sig;
}

// Black-Scholes-Merton call option delta function
// [[Rcpp::export]]
double BSMDeltaCall(double S, double K, double r, double vol, double time, double q = 0) {
  double d1, ln;
  ln = log(S / K);
  d1 = (ln + (r - q + (vol * vol) / 2) * time) / (vol * sqrt(time));
  return exp(-q * time) * R::pnorm(d1, 0.0, 1.0, true, false);
}

// Black-Scholes-Merton put option delta function
// [[Rcpp::export]]
double BSMDeltaPut(double S, double K, double r, double vol, double time, double q = 0) {
  double d1, ln;
  ln = log(S / K);
  d1 = (ln + (r - q + (vol * vol) / 2) * time) / (vol * sqrt(time));
  return -exp(-q * time) * R::pnorm(-d1, 0.0, 1.0, true, false);
}

// Black-Scholes-Merton call option epsilon function
// [[Rcpp::export]]
double BSMEpsilonCall(double S, double K, double r, double vol, double time, double q = 0) {
  double d1, ln;
  ln = log(S / K);
  d1 = (ln + (r - q + (vol * vol) / 2) * time) / (vol * sqrt(time));
  return -exp(-q * time) * time * S * R::pnorm(d1, 0.0, 1.0, true, false);
}

// Black-Scholes-Merton put option epsilon function
// [[Rcpp::export]]
double BSMEpsilonPut(double S, double K, double r, double vol, double time, double q = 0) {
  double d1, ln;
  ln = log(S / K);
  d1 = (ln + (r - q + (vol * vol) / 2) * time) / (vol * sqrt(time));
  return exp(-q * time) * time * S * R::pnorm(-d1, 0.0, 1.0, true, false);
}

// Black-Scholes-Merton vega function
// [[Rcpp::export]]
double BSMVega(double S, double K, double r, double vol, double time, double q = 0) {
  double d1, ln;
  ln = log(S / K);
  d1 = (ln + (r - q + (vol * vol) / 2) * time) / (vol * sqrt(time));
  return exp(-q * time) * S * (1 / (sqrt(2 * M_PI))) * exp(-0.5 * (d1 * d1)) * sqrt(time);
}

// Black-Scholes-Merton gamma function
// [[Rcpp::export]]
double BSMGamma(double S, double K, double r, double vol, double time, double q = 0) {
  double d1, ln;
  ln = log(S / K);
  d1 = (ln + (r - q + (vol * vol) / 2) * time) / (vol * sqrt(time));
  return exp(-q * time) * (1 / (sqrt(2 * M_PI))) * exp(-0.5 * (d1 * d1)) / S / sqrt(time) / vol;
}

// Black-Scholes-Mertoncall option theta function
// [[Rcpp::export]]
double BSMThetaCall(double S, double K, double r, double vol, double time, double q = 0) {
  double d1, d2, ln;
  ln = log(S / K);
  d1 = (ln + (r - q + (vol * vol) / 2) * time) / (vol * sqrt(time));
  d2 = d1 - vol * sqrt(time);
  return -S * exp(-q * time) * (1 / (sqrt(2 * M_PI))) * exp(-0.5 * (d1 * d1)) * 0.5 * vol / sqrt(time) - r * K * exp(-r * time) * R::pnorm(d2, 0.0, 1.0, true, false) + q * exp(-q * time) * S * R::pnorm(d1, 0.0, 1.0, true, false);
}

// Black-Scholes-Merton put option theta function
// [[Rcpp::export]]
double BSMThetaPut(double S, double K, double r, double vol, double time, double q = 0) {
  double d1, d2, ln;
  ln = log(S / K);
  d1 = (ln + (r - q + (vol * vol) / 2) * time) / (vol * sqrt(time));
  d2 = d1 - vol * sqrt(time);
  return -S * exp(-q * time) * (1 / (sqrt(2 * M_PI))) * exp(-0.5 * (d1 * d1)) * 0.5 * vol / sqrt(time) + r * K * exp(-r * time) * R::pnorm(-d2, 0.0, 1.0, true, false) - q * exp(-q * time) * S * R::pnorm(-d1, 0.0, 1.0, true, false);
}

// Black-Scholes-Merton call option rho function
// [[Rcpp::export]]
double BSMRhoCall(double S, double K, double r, double vol, double time, double q = 0) {
  double d1, d2, ln;
  ln = log(S / K);
  d1 = (ln + (r - q + (vol * vol) / 2) * time) / (vol * sqrt(time));
  d2 = d1 - vol * sqrt(time);
  return K * time * exp(-r * time) * R::pnorm(d2, 0.0, 1.0, true, false);
}

// Black-Scholes-Merton put option rho function
// [[Rcpp::export]]
double BSMRhoPut(double S, double K, double r, double vol, double time, double q = 0) {
  double d1, d2, ln;
  ln = log(S / K);
  d1 = (ln + (r - q + (vol * vol) / 2) * time) / (vol * sqrt(time));
  d2 = d1 - vol * sqrt(time);
  return -K * time * exp(-r * time) * R::pnorm(-d2, 0.0, 1.0, true, false);
}

// Black-Scholes-Merton lambda function
// [[Rcpp::export]]
double BSMLambdaCall(double S, double K, double r, double vol, double time, double q = 0) {
  double d1, d2, ln, p1;
  ln = log(S / K);
  d1 = (ln + (r - q + (vol * vol) / 2) * time) / (vol * sqrt(time));
  d2 = d1 - vol * sqrt(time);
  p1 = R::pnorm(d1, 0.0, 1.0, true, false);
  return exp(-q * time) * p1 * S /
    (exp(-q * time) * S * p1 - K * exp(-r * time) * R::pnorm(d2, 0.0, 1.0, true, false));
}

// [[Rcpp::export]]
double BSMLambdaPut(double S, double K, double r, double vol, double time, double q = 0) {
  double d1, d2, ln, p1;
  ln = log(S / K);
  d1 = (ln + (r - q + (vol * vol) / 2) * time) / (vol * sqrt(time));
  d2 = d1 - vol * sqrt(time);
  p1 = R::pnorm(-d1, 0.0, 1.0, true, false);
  return -exp(-q * time) * p1 * S /
    ( -exp(-q * time) * S * p1 + K * exp(-r * time) * R::pnorm(-d2, 0.0, 1.0, true, false) );
}

// Black-Scholes-Merton vanna function
// [[Rcpp::export]]
double BSMVanna(double S, double K, double r, double vol, double time, double q = 0) {
  double d1, d2, ln;
  ln = log(S / K);
  d1 = (ln + (r - q + (vol * vol) / 2) * time) / (vol * sqrt(time));
  d2 = d1 - vol * sqrt(time);
  return -exp(-q * time) * (1 / sqrt(2 * M_PI)) * exp(-0.5 * (d1 * d1)) * d2 / vol;
}

// Black-Scholes-Merton charm function
// [[Rcpp::export]]
double BSMCharmCall(double S, double K, double r, double vol, double time, double q = 0) {
  double d1, d2, ln, p1;
  ln = log(S / K);
  d1 = (ln + (r - q + (vol * vol) / 2) * time) / (vol * sqrt(time));
  d2 = d1 - vol * sqrt(time);
  p1 = R::pnorm(d1, 0.0, 1.0, true, false);
  return q * exp(-q * time) * p1  - exp(-q * time) * (1 / sqrt(2 * M_PI)) * exp(-0.5 * (d1 * d1)) * (r - q - 0.5 * vol * d2 / sqrt(time) ) / vol / sqrt(time);
}

// [[Rcpp::export]]
double BSMCharmPut(double S, double K, double r, double vol, double time, double q = 0) {
  double d1, d2, ln, p1;
  ln = log(S / K);
  d1 = (ln + (r - q + (vol * vol) / 2) * time) / (vol * sqrt(time));
  d2 = d1 - vol * sqrt(time);
  p1 = R::pnorm(-d1, 0.0, 1.0, true, false);
  return -q * exp(-q * time) * p1  - exp(-q * time) * (1 / sqrt(2 * M_PI)) * exp(-0.5 * (d1 * d1)) * (r - q - 0.5 * vol * d2 / sqrt(time) ) / vol / sqrt(time);
}
