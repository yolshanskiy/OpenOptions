#' @export
BSMGreeks <- function(S, K, r, vol, time, q = 0, option_type){
  stopifnot_custom <- function(A, message){
    if(any(is.na(A)) || !all(A)) stop(message)
  }

  stopifnot_custom(( tolower(option_type) == "call" |
                tolower(option_type) == "put"),
            "The option_type is inappropriate.")
  stopifnot_custom(S > 0 & is.double(S), "S must be positive.")
  stopifnot_custom(min(K) > 0 & is.double(K), "K must be positive.")
  stopifnot_custom(vol > 0 & is.double(vol), "vol must be positive.")
  stopifnot_custom(time > 0 & is.double(time), "time must be positive.")
  stopifnot_custom(is.double(q), "q must be positive number.")

  out = data.frame(K = K)

  if(tolower(option_type) == "call"){
    out[["Value"]] = sapply(K, function(K) BSMCall(S, K, r, vol, time, q = 0) )
    out[["Delta"]] = sapply(K, function(K) BSMDeltaCall(S, K, r, vol, time, q = 0) )
    out[["Gamma"]] = sapply(K, function(K) BSMGamma(S, K, r, vol, time, q = 0) )
    out[["Vega"]] = sapply(K, function(K) BSMVega(S, K, r, vol, time, q = 0) )
    out[["Vanna"]] = sapply(K, function(K) BSMVanna(S, K, r, vol, time, q = 0) )
    out[["Rho"]] = sapply(K, function(K) BSMRhoCall(S, K, r, vol, time, q = 0) )
    out[["Theta"]] = sapply(K, function(K) BSMThetaCall(S, K, r, vol, time, q = 0) )
    out[["Epsilon"]] = sapply(K, function(K) BSMEpsilonCall(S, K, r, vol, time, q = 0) )
    out[["Lambda"]] = sapply(K, function(K) BSMLambdaCall(S, K, r, vol, time, q = 0) )
    out[["Charm"]] = sapply(K, function(K) BSMCharmCall(S, K, r, vol, time, q = 0) )
  }
  if(tolower(option_type) == "put"){
    out[["Value"]] = sapply(K, function(K) BSMPut(S, K, r, vol, time, q = 0) )
    out[["Delta"]] = sapply(K, function(K) BSMDeltaPut(S, K, r, vol, time, q = 0) )
    out[["Gamma"]] = sapply(K, function(K) BSMGamma(S, K, r, vol, time, q = 0) )
    out[["Vega"]] = sapply(K, function(K) BSMVega(S, K, r, vol, time, q = 0) )
    out[["Vanna"]] = sapply(K, function(K) BSMVanna(S, K, r, vol, time, q = 0) )
    out[["Rho"]] = sapply(K, function(K) BSMRhoPut(S, K, r, vol, time, q = 0) )
    out[["Theta"]] = sapply(K, function(K) BSMThetaPut(S, K, r, vol, time, q = 0) )
    out[["Epsilon"]] = sapply(K, function(K) BSMEpsilonPut(S, K, r, vol, time, q = 0) )
    out[["Lambda"]] = sapply(K, function(K) BSMLambdaPut(S, K, r, vol, time, q = 0) )
    out[["Charm"]] = sapply(K, function(K) BSMCharmPut(S, K, r, vol, time, q = 0) )
  }

  return(out)
}
