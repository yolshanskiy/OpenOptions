#include <Rcpp.h>
using namespace Rcpp;

// Functions to set binomial tree
// [[Rcpp::export]]
Rcpp::List binomialStockPaths(double S, double K, double r, double vol, double time, double Div = 0, int N = 1, bool nonstandard_dividend_payment = false) {

  if (S <= 0)
    Rcpp::stop("S must be positive");
  if (K <= 0)
    Rcpp::stop("K must be positive");
  if (vol <= 0)
    Rcpp::stop("vol must be positive");
  if (time <= 0)
    Rcpp::stop("time must be positive");

  if(nonstandard_dividend_payment){
    double dt = time / N;
    double u = exp(vol * sqrt(dt));
    double d = 1 / u;
    double div_rate = (1 - exp(-dt*Div));

    NumericMatrix StockPaths(N + 1, N + 1), Dividends(N + 1, N + 1);
    StockPaths.fill(NA_REAL);
    Dividends.fill(NA_REAL);
    StockPaths(0,0) = S; Dividends(0,0) = 0;

    // Backwards induction to calculate option values at each time step
    for (int j = 1; j <= N; j++) {
      for (int i = 0; i <= j-1; i++) {
        double value = StockPaths(i,j-1) * u;
        Dividends(i,j) = value * div_rate;
        StockPaths(i,j) = value - Dividends(i,j);
      }
      double value =  StockPaths(j-1,j-1) * d;
      Dividends(j,j) = value * div_rate;
      StockPaths(j,j) = value - Dividends(j,j);
    }
    return Rcpp::List::create(Rcpp::Named("Price") = StockPaths, Rcpp::Named("Dividends") = Dividends);
  }else{
    double dt = time / N;
    double u = exp(vol * sqrt(dt));
    double d = 1 / u;
    double div_rate = (exp(dt*Div) - 1);

    NumericMatrix StockPaths(N + 1, N + 1), Dividends(N + 1, N + 1);
    StockPaths.fill(NA_REAL);
    Dividends.fill(NA_REAL);
    StockPaths(0,0) = S; Dividends(0,0) = 0;

    // Backwards induction to calculate option values at each time step
    for (int j = 1; j <= N; j++) {
      for (int i = 0; i <= j-1; i++) {
        double value = StockPaths(i,j-1) * u;
        Dividends(i,j) = value * div_rate;
        StockPaths(i,j) = value;
      }
      double value =  StockPaths(j-1,j-1) * d;
      Dividends(j,j) = value * div_rate;
      StockPaths(j,j) = value;
    }
    return Rcpp::List::create(Rcpp::Named("Price") = StockPaths, Rcpp::Named("Dividends") = Dividends);
  }
}

// [[Rcpp::export]]
NumericMatrix binomialStatePrices(double S, double K, double r, double vol, double time, double Div = 0, int N = 1, bool verbose = false, bool nonstandard_dividend_payment = false) {
  // Check parameters
  if (S <= 0)
    Rcpp::stop("S must be positive");
  if (K <= 0)
    Rcpp::stop("K must be positive");
  if (vol <= 0)
    Rcpp::stop("vol must be positive");
  if (time <= 0)
    Rcpp::stop("time must be positive");

  if(nonstandard_dividend_payment){
    // if end_period_payment true then state prices are not affected by Dividend

    // Introduce standard intermediate variables
    double dt = time / N;
    double u = exp(vol * sqrt(dt));
    double d = 1 / u;
    double p = (exp( (r) * dt) - d) / (u - d);
    double q = 1 - p;

    double discount = exp(-r * dt);
    double psi_up = discount * p, psi_dn = discount * q;

    if(verbose){
      // Print values of u, d, p, and q
      Rcpp::Rcout << "dt: " << dt << std::endl;
      Rcpp::Rcout << "u: " << u << std::endl;
      Rcpp::Rcout << "d: " << d << std::endl;
      Rcpp::Rcout << "p: " << p << std::endl;
      Rcpp::Rcout << "q: " << q << std::endl;
    }
    if (p <= 0 or q <= 0) {
      Rcpp::stop("No arbitrage fails");
    }
    NumericMatrix StatePrices(N + 1, N + 1);
    StatePrices.fill(NA_REAL); StatePrices(0,0) = 1;
    for (int j = 1; j <= N; j++) {
      StatePrices(0,j) = StatePrices(0,j-1) * psi_up;
      StatePrices(j,j) = StatePrices(j-1,j-1) * psi_dn;
      for (int i = 1; i <= j-1; i++) {
        StatePrices(i,j) = StatePrices(i,j-1) * psi_up + StatePrices(i-1,j-1) * psi_dn;
      }
    }
    return StatePrices;
  }else{
    // if end_period_payment false then dividend is reflected in stock

    // Introduce standard intermediate variables
    double dt = time / N;
    double u = exp(vol * sqrt(dt));
    double d = 1 / u;
    double p = (exp( (r - Div) * dt) - d) / (u - d);
    double q = 1 - p;

    double discount = exp(-r * dt);
    double psi_up = discount * p, psi_dn = discount * q;

    if(verbose){
      // Print values of u, d, p, and q
      Rcpp::Rcout << "dt: " << dt << std::endl;
      Rcpp::Rcout << "u: " << u << std::endl;
      Rcpp::Rcout << "d: " << d << std::endl;
      Rcpp::Rcout << "p: " << p << std::endl;
      Rcpp::Rcout << "q: " << q << std::endl;
    }
    if (p <= 0 or q <= 0) {
      Rcpp::stop("No arbitrage fails");
    }
    NumericMatrix StatePrices(N + 1, N + 1);
    StatePrices.fill(NA_REAL); StatePrices(0,0) = 1;
    for (int j = 1; j <= N; j++) {
      StatePrices(0,j) = StatePrices(0,j-1) * psi_up;
      StatePrices(j,j) = StatePrices(j-1,j-1) * psi_dn;
      for (int i = 1; i <= j-1; i++) {
        StatePrices(i,j) = StatePrices(i,j-1) * psi_up + StatePrices(i-1,j-1) * psi_dn;
      }
    }
    return StatePrices;
  }
}

// Functions to calculate option price using binomial tree model

// [[Rcpp::export]]
double binomialOptionPricer(double S, double K, double r, double vol, double time, double Div = 0,
                            int N = 1, std::string type = "call", bool verbose = false, bool nonstandard_dividend_payment = false) {

  if (S <= 0)
    Rcpp::stop("S must be positive");
  if (K <= 0)
    Rcpp::stop("K must be positive");
  if (vol <= 0)
    Rcpp::stop("vol must be positive");
  if (time <= 0)
    Rcpp::stop("time must be positive");
  if (type != "call" && type != "put")
    Rcpp::stop("Type must be 'call' or 'put'");

  double dt = time / N;
  double u = exp(vol * sqrt(dt));
  double d = 1 / u;
  double p, div_discount;
  if(nonstandard_dividend_payment){
    p = (exp( r * dt) - d) / (u - d);
    div_discount = exp( - Div * dt * N);
  }else{
    p = (exp( (r - Div) * dt) - d) / (u - d);
  }

  double q = 1 - p;
  double discount = exp(-r * dt);

  if(verbose){
    // Print values of u, d, p, and q
    Rcpp::Rcout << "dt: " << dt << std::endl;
    Rcpp::Rcout << "u: " << u << std::endl;
    Rcpp::Rcout << "d: " << d << std::endl;
    Rcpp::Rcout << "p: " << p << std::endl;
    Rcpp::Rcout << "q: " << q << std::endl;
    Rcpp::Rcout << "discount: " << discount << std::endl;
  }


  if (p <= 0 or q <= 0)
    Rcpp::stop("No arbitrage fails");

  NumericMatrix optionValues(N + 1, N + 1);
  optionValues.fill(NA_REAL);

  // Initialize option values at maturity
  for (int i = 0; i <= N; i++) {
    double ST;
    if(nonstandard_dividend_payment){
      ST = S * pow(u, N - i) * pow(d, i) * div_discount;
    }else{
      ST = S * pow(u, N - i) * pow(d, i);
    }

    if (type == "call") {
      optionValues(i, N) = std::max(ST - K, 0.0);
    } else {
      optionValues(i, N) = std::max(K - ST, 0.0);
    }
  }

  // Backwards induction to calculate option values at each time step
  for (int j = N - 1; j >= 0; j--) {
    for (int i = 0; i <= j; i++) {
      if (type == "call") {
        optionValues(i, j) = discount * (p * optionValues(i, j + 1) + q * optionValues(i + 1, j + 1));
      } else {
        optionValues(i, j) = discount * (p * optionValues(i, j + 1) + q * optionValues(i + 1, j + 1));
      }
    }
  }

  return optionValues(0,0);
}

// [[Rcpp::export]]
NumericMatrix binomialOptionPricer_detailed(double S, double K, double r, double vol, double time, double Div = 0,
                            int N = 1, std::string type = "call", bool verbose = false, bool nonstandard_dividend_payment = false) {

  if (S <= 0)
    Rcpp::stop("S must be positive");
  if (K <= 0)
    Rcpp::stop("K must be positive");
  if (vol <= 0)
    Rcpp::stop("vol must be positive");
  if (time <= 0)
    Rcpp::stop("time must be positive");
  if (type != "call" && type != "put")
    Rcpp::stop("Type must be 'call' or 'put'");

  double dt = time / N;
  double u = exp(vol * sqrt(dt));
  double d = 1 / u;
  double p, div_discount;
  if(nonstandard_dividend_payment){
    p = (exp( r * dt) - d) / (u - d);
    div_discount = exp( - Div * dt * N);
  }else{
    p = (exp( (r - Div) * dt) - d) / (u - d);
  }

  double q = 1 - p;
  double discount = exp(-r * dt);

  if(verbose){
    // Print values of u, d, p, and q
    Rcpp::Rcout << "dt: " << dt << std::endl;
    Rcpp::Rcout << "u: " << u << std::endl;
    Rcpp::Rcout << "d: " << d << std::endl;
    Rcpp::Rcout << "p: " << p << std::endl;
    Rcpp::Rcout << "q: " << q << std::endl;
    Rcpp::Rcout << "discount: " << discount << std::endl;
  }


  if (p <= 0 or q <= 0)
    Rcpp::stop("No arbitrage fails");

  NumericMatrix optionValues(N + 1, N + 1);
  optionValues.fill(NA_REAL);

  // Initialize option values at maturity
  for (int i = 0; i <= N; i++) {
    double ST;
    if(nonstandard_dividend_payment){
      ST = S * pow(u, N - i) * pow(d, i) * div_discount;
    }else{
      ST = S * pow(u, N - i) * pow(d, i);
    }

    if (type == "call") {
      optionValues(i, N) = std::max(ST - K, 0.0);
    } else {
      optionValues(i, N) = std::max(K - ST, 0.0);
    }
  }

  // Backwards induction to calculate option values at each time step
  for (int j = N - 1; j >= 0; j--) {
    for (int i = 0; i <= j; i++) {
      if (type == "call") {
        optionValues(i, j) = discount * (p * optionValues(i, j + 1) + q * optionValues(i + 1, j + 1));
      } else {
        optionValues(i, j) = discount * (p * optionValues(i, j + 1) + q * optionValues(i + 1, j + 1));
      }
    }
  }

  return optionValues;
}


// [[Rcpp::export]]
double binomialOptionPricer_fast(double S, double K, double r, double vol, double time, double Div = 0,
                                 int N = 1, std::string type = "call", bool verbose = false, bool nonstandard_dividend_payment = false,
                                 bool robust_warning = true) {
  if(robust_warning) Rcpp::Rcout << "Warning, binomialOptionPricer_fast is not robust for large N. Verify answer with binomialOptionPricer." << std::endl;
  if (S <= 0)
    Rcpp::stop("S must be positive");
  if (K <= 0)
    Rcpp::stop("K must be positive");
  if (vol <= 0)
    Rcpp::stop("vol must be positive");
  if (time <= 0)
    Rcpp::stop("time must be positive");
  if (type != "call" && type != "put")
    Rcpp::stop("Type must be 'call' or 'put'");

  double dt = time / N;
  double u = exp(vol * sqrt(dt));
  double d = 1 / u;
  double p, div_discount = 1;
  if(nonstandard_dividend_payment){
    p = (exp( r * dt) - d) / (u - d);
    div_discount = exp( - Div * dt * N);
  }else{
    p = (exp( (r - Div) * dt) - d) / (u - d);
  }

  double q = 1 - p;
  double discount = exp(-r * time);

  if(verbose){
    // Print values of u, d, p, and q
    Rcpp::Rcout << "dt: " << dt << std::endl;
    Rcpp::Rcout << "u: " << u << std::endl;
    Rcpp::Rcout << "d: " << d << std::endl;
    Rcpp::Rcout << "p: " << p << std::endl;
    Rcpp::Rcout << "q: " << q << std::endl;
    Rcpp::Rcout << "discount: " << discount << std::endl;
  }

  if (p <= 0 or q <= 0) Rcpp::stop("No arbitrage fails");

  double optionValue = 0;
  double ST = S * div_discount * pow(u, N), adj_ST = d / u;
  double rn = pow(p, N), adj_rn = q / p; // risk-neutral probability of the state
  double bin = 1;
  // Rcout << bin << "-" << rn << " ";
  for (int i = 0; i <= N; i++) {
    if (type == "call") {
      if(ST - K > 0)
        optionValue += (ST - K) * rn;
    } else {
      if(K - ST > 0)
        optionValue += (K - ST) * rn;
    }
    ST = ST * adj_ST;
    rn = rn * adj_rn * (N - i) / (i + 1);
    // bin = bin * (N - i) / (i + 1);
    // Rcout << bin << "-" << rn << " ";
  }
  // Rcout << "\n";
  return optionValue * discount;
}

// [[Rcpp::export]]
double binomialAmericanOptionPricer(double S, double K, double r, double vol, double time, double Div = 0,
                                    int N = 1, std::string type = "call", bool verbose = false, bool nonstandard_dividend_payment = false) {

  if (S <= 0) Rcpp::stop("S must be positive");
  if (K <= 0) Rcpp::stop("K must be positive");
  if (vol <= 0) Rcpp::stop("vol must be positive");
  if (time <= 0) Rcpp::stop("time must be positive");
  if (type != "call" && type != "put") Rcpp::stop("Type must be 'call' or 'put'");

  double dt = time / N;
  double u = exp(vol * sqrt(dt));
  double d = 1 / u;
  double p, div_discount;
  if(nonstandard_dividend_payment){
    p = (exp( r * dt) - d) / (u - d);
    div_discount = exp( - Div * dt * N);
  }else{
    p = (exp( (r - Div) * dt) - d) / (u - d);
  }

  double q = 1 - p;
  double discount = exp(-r * dt);

  if(verbose){
    // Print values of u, d, p, and q
    Rcpp::Rcout << "dt: " << dt << std::endl;
    Rcpp::Rcout << "u: " << u << std::endl;
    Rcpp::Rcout << "d: " << d << std::endl;
    Rcpp::Rcout << "p: " << p << std::endl;
    Rcpp::Rcout << "q: " << q << std::endl;
    Rcpp::Rcout << "discount: " << discount << std::endl;
  }


  if (p <= 0 or q <= 0) Rcpp::stop("No arbitrage fails");

  NumericMatrix optionValues(N + 1, N + 1);
  optionValues.fill(NA_REAL);

  // Initialize option values at maturity
  for (int i = 0; i <= N; i++) {
    double ST;
    if(nonstandard_dividend_payment){
      ST = S * pow(u, N - i) * pow(d, i) * div_discount;
    }else{
      ST = S * pow(u, N - i) * pow(d, i);
    }

    if (type == "call") {
      optionValues(i, N) = std::max(ST - K, 0.0);
    } else {
      optionValues(i, N) = std::max(K - ST, 0.0);
    }
  }

  // Backwards induction to calculate option values at each time step
  double mult = d / u;
  for (int j = N - 1; j >= 0; j--) {
    double ST;
    if(nonstandard_dividend_payment){
      ST = S * pow(u, j) * div_discount;
    }else{
      ST = S * pow(u, j);
    }
    for (int i = 0; i <= j; i++) {
      if (type == "call") {
        optionValues(i, j) = std::max(discount * (p * optionValues(i, j + 1) + q * optionValues(i + 1, j + 1)), ST - K);
      } else {
        optionValues(i, j) = std::max(discount * (p * optionValues(i, j + 1) + q * optionValues(i + 1, j + 1)), K - ST);
      }
      ST = ST * mult;
    }
  }

  return optionValues(0,0);
}


// [[Rcpp::export]]
Rcpp::List binomialAmericanOptionPricer_detailed(double S, double K, double r, double vol, double time, double Div = 0,
                                    int N = 1, std::string type = "call", bool verbose = false, bool nonstandard_dividend_payment = false) {

  if (S <= 0)  Rcpp::stop("S must be positive");
  if (K <= 0)  Rcpp::stop("K must be positive");
  if (vol <= 0)  Rcpp::stop("vol must be positive");
  if (time <= 0)  Rcpp::stop("time must be positive");
  if (type != "call" && type != "put")  Rcpp::stop("Type must be 'call' or 'put'");

  double dt = time / N;
  double u = exp(vol * sqrt(dt));
  double d = 1 / u;
  double p, div_discount;
  if(nonstandard_dividend_payment){
    p = (exp( r * dt) - d) / (u - d);
    div_discount = exp( - Div * dt * N);
  }else{
    p = (exp( (r - Div) * dt) - d) / (u - d);
  }

  double q = 1 - p;
  double discount = exp(-r * dt);

  if(verbose){
    // Print values of u, d, p, and q
    Rcpp::Rcout << "dt: " << dt << std::endl;
    Rcpp::Rcout << "u: " << u << std::endl;
    Rcpp::Rcout << "d: " << d << std::endl;
    Rcpp::Rcout << "p: " << p << std::endl;
    Rcpp::Rcout << "q: " << q << std::endl;
    Rcpp::Rcout << "discount: " << discount << std::endl;
  }


  if (p <= 0 or q <= 0) Rcpp::stop("No arbitrage fails");

  NumericMatrix optionValues(N + 1, N + 1), HoldingValue(N + 1, N + 1);
  optionValues.fill(NA_REAL); HoldingValue.fill(NA_REAL);

  // Initialize option values at maturity
  for (int i = 0; i <= N; i++) {
    double ST;
    if(nonstandard_dividend_payment){
      ST = S * pow(u, N - i) * pow(d, i) * div_discount;
    }else{
      ST = S * pow(u, N - i) * pow(d, i);
    }

    if (type == "call") {
      optionValues(i, N) = std::max(ST - K, 0.0);
    } else {
      optionValues(i, N) = std::max(K - ST, 0.0);
    }
  }

  // Backwards induction to calculate option values at each time step
  double mult = d / u;
  for (int j = N - 1; j >= 0; j--) {
    double ST;
    if(nonstandard_dividend_payment){
      ST = S * pow(u, j) * div_discount;
    }else{
      ST = S * pow(u, j);
    }
    for (int i = 0; i <= j; i++) {

      double tmp_contin_value = discount * (p * optionValues(i, j + 1) + q * optionValues(i + 1, j + 1));
      if (type == "call") {
        if(ST - K > tmp_contin_value){
          optionValues(i, j) = ST - K;
          HoldingValue(i, j) = tmp_contin_value;
        }else{
          optionValues(i, j) = tmp_contin_value;
        }
      }else{
        if(K - ST > tmp_contin_value){
          optionValues(i, j) = ST - K;
          HoldingValue(i, j) = tmp_contin_value;
        }else{
          optionValues(i, j) = tmp_contin_value;
        }
      }
      ST = ST * mult;
    }
  }

  return Rcpp::List::create(Rcpp::Named("V") = optionValues, Rcpp::Named("HoldingValue") = HoldingValue);
}
