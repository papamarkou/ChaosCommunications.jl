opt_success = (:NLOPT_SUCCESS, :NLOPT_STOPVAL_REACHED, :NLOPT_FTOL_REACHED, :NLOPT_XTOL_REACHED, :NLOPT_MAXEVAL_REACHED,
  :NLOPT_MAXTIME_REACHED)

opt_failure = (:NLOPT_FAILURE, :NLOPT_INVALID_ARGS, :NLOPT_OUT_OF_MEMORY, :NLOPT_ROUNDOFF_LIMITED, :NLOPT_FORCED_STOP)

# MCML decoder for non-coherent UTURUniChannelCSK
immutable UTURUniChannelCSKMCMLDecoder <: MCMLDecoder
  carrier::Carrier
  nmc::Int64
  opt::Vector{NLoptimizer} # NLopt Opt types for the two functions to be mazimized
  init::Vector{Float64}

  function UTURUniChannelCSKMCMLDecoder(c::Carrier, nmc::Int64, opt::Vector{NLoptimizer}, init::Vector{Float64})
    @assert length(opt) == 2 "2 optimizers are needed for the MCML decoder of UTURUniChannelCSK."
    @assert length(init) == 2 "The init field of the MCML decoder of UTURUniChannelCSK must be a vector of length 2."
    new(c, nmc, opt, init)
  end
end

UTURUniChannelCSKMCMLDecoder(c::Carrier, nmc::Int64, opt::NLoptimizer, init::Vector{Float64}) =
  UTURUniChannelCSKMCMLDecoder(c, nmc, [opt, opt], init)
UTURUniChannelCSKMCMLDecoder(c::Carrier, nmc::Int64, opt::Vector{NLoptimizer}, init::Float64) =
  UTURUniChannelCSKMCMLDecoder(c, nmc, opt, [init, init])
UTURUniChannelCSKMCMLDecoder(c::Carrier, nmc::Int64, opt::NLoptimizer, init::Float64) =
  UTURUniChannelCSKMCMLDecoder(c, nmc, [opt, opt], [init, init])

# Default optimizer for MCML decoder of non-coherent UTURUniChannelCSK
uturunichannelcsk_default_opt = NLoptimizer(alg=:LD_SLSQP, lb=0.0, ub=100.0, absftol=NaN, relftol=NaN, absxtol=1e-32,
  relxtol=NaN, nfeval=-1, maxtime=-1.0)

UTURUniChannelCSKMCMLDecoder(; carrier::Carrier=LogisticCarrier(5), nmc::Int64=50,
  opt::Union(NLoptimizer, Vector{NLoptimizer})=uturunichannelcsk_default_opt,
  init::Union(Float64, Vector{Float64})=1.0) =
  UTURUniChannelCSKMCMLDecoder(carrier, nmc, opt, init)

# Wrapper for MCML decoder construction of any systems that supports MCML decoding
mcml_types = (:UTURUniChannelCSK,)

function mcml_decoder(; system::Symbol=:UTURUniChannelCSK, args...)
  @assert in(system, mcml_types) "MCML decoder not defined for $system."

  if system == :UTURUniChannelCSK
    UTURUniChannelCSKMCMLDecoder(; args...)
  end
end
