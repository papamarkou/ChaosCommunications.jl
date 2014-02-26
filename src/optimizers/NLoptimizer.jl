abstract Optimizer

# NLoptimizer holds the parameters passed to Opt. This is a workaround for the lack of serialization of Opt.
immutable NLoptimizer <: Optimizer
  npars::Int # Number of optimization parameters
  alg::Symbol # Optimization algorithm to use
  lb::VF64OrNothing # Lower bound constraint
  ub::VF64OrNothing # Upper bound constraint
  absftol::Float64 # Stopping criterion: absolute tolerance on function value
  relftol::Float64 # Stopping criterion: relative tolerance on function value
  absxtol::Float64 # Stopping criterion: absolute tolerance on optimization parameters
  relxtol::Float64 # Stopping criterion: relative tolerance on optimization parameters
  nfeval::Int64 # Stopping criterion: maximum number of function evaluations
  maxtime::Float64 # Stopping criterion: maximum optimization time

  function NLoptimizer(npars::Int, alg::Symbol, lb::VF64OrNothing, ub::VF64OrNothing, absftol::Float64,
    relftol::Float64, absxtol::Float64, relxtol::Float64, nfeval::Int64, maxtime::Float64)
    if isnan(absftol) && isnan(relftol) && isnan(absxtol) && isnan(relxtol) && nfeval < 0 && maxtime < 0
      error("A stopping criterion must be specified.")
    end

    new(npars, alg, lb, ub, absftol, relftol, absxtol, relxtol, nfeval, maxtime)
  end
end

NLoptimizer(npars::Int;
  alg::Symbol=:LN_COBYLA,
  lb::VF64OrNothing=nothing,
  ub::VF64OrNothing=nothing,
  absftol::Float64=1e-32,
  relftol::Float64=NaN,
  absxtol::Float64=1e-32,
  relxtol::Float64=NaN,
  nfeval::Int64=1_000,
  maxtime::Float64=-1.0) =
  NLoptimizer(npars, alg, lb, ub, absftol, relftol, absxtol, relxtol, nfeval, maxtime)

function NLoptimizer(;
  alg::Symbol=:LN_COBYLA,
  lb::F64OrNothing=nothing,
  ub::F64OrNothing=nothing,
  absftol::Float64=1e-32,
  relftol::Float64=NaN,
  absxtol::Float64=1e-32,
  relxtol::Float64=NaN,
  nfeval::Int64=1_000,
  maxtime::Float64=-1.0)
  lbvec = lb == nothing ? nothing : [lb]
  ubvec = ub == nothing ? nothing : [ub]  
  NLoptimizer(1, alg, lbvec, ubvec, absftol, relftol, absxtol, relxtol, nfeval, maxtime)
end

function convert(::Type{Opt}, o::NLoptimizer)
  opt = Opt(o.alg, o.npars)

  if o.lb != nothing; lower_bounds!(opt, o.lb); end    
  if o.ub != nothing; upper_bounds!(opt, o.ub); end
  if !isnan(o.absftol) ftol_abs!(opt, o.absftol); end
  if !isnan(o.relftol) ftol_rel!(opt, o.relftol); end
  if !isnan(o.absxtol) xtol_abs!(opt, o.absxtol); end
  if !isnan(o.relxtol) xtol_rel!(opt, o.relxtol); end
  if !(o.nfeval < 0) maxeval!(opt, o.nfeval); end
  if !(o.maxtime < 0) maxtime!(opt, o.maxtime); end

  return opt
end
