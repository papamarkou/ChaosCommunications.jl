abstract Decoder
abstract DeterministicDecoder <: Decoder
abstract QuasiRandomDecoder <: Decoder
abstract RandomDecoder <: Decoder
# abstract EMDecoder <: QuasiRandomDecoder
# abstract GibbsDecoder <: RandomDecoder
# abstract ParticleFilterDecoder <: RandomDecoder

type CorDecoder <: DeterministicDecoder
end

abstract MCMLDecoder <: QuasiRandomDecoder

opt_success = (:NLOPT_SUCCESS, :NLOPT_STOPVAL_REACHED, :NLOPT_FTOL_REACHED, :NLOPT_XTOL_REACHED, :NLOPT_MAXEVAL_REACHED,
  :NLOPT_MAXTIME_REACHED)

opt_failure = (:NLOPT_FAILURE, :NLOPT_INVALID_ARGS, :NLOPT_OUT_OF_MEMORY, :NLOPT_ROUNDOFF_LIMITED, :NLOPT_FORCED_STOP)

# MCML decoder for non-coherent UTURUniChannelCSK
immutable UTURUniChannelCSKMCMLDecoder <: MCMLDecoder
  carrier::Carrier
  nmc::Int64
  opt::Opt # NLopt Opt type
  init::Vector{Float64}

  function UTURUniChannelCSKMCMLDecoder(c::Carrier, nmc::Int64, opt::Opt, init::Vector{Float64})
    @assert length(init) == 2 "The init field of the MCML decoder of UTURUniChannelCSK must be a vector of length 2."
    new(c, nmc, opt, init)
  end
end

UTURUniChannelCSKMCMLDecoder(c::Carrier, nmc::Int64, opt::Opt, init::Float64) =
  UTURUniChannelCSKMCMLDecoder(c, nmc, opt, [init, init])

# Default optimizer for MCML decoder of non-coherent UTURUniChannelCSK
uturunichannelcsk_default_opt = Opt(:LD_MMA, 1)
lower_bounds!(uturunichannelcsk_default_opt, [0.])
upper_bounds!(uturunichannelcsk_default_opt, [10.])
xtol_rel!(uturunichannelcsk_default_opt, 1e-32)

UTURUniChannelCSKMCMLDecoder(; carrier::Carrier=LogisticCarrier(5), nmc::Int64=50,
  opt::Opt=uturunichannelcsk_default_opt, init=1.) = UTURUniChannelCSKMCMLDecoder(carrier, nmc, opt, init)

# Wrapper for MCML decoder construction of any systems that supports MCML decoding
mcml_types = (:UTURUniChannelCSK,)

function mcml_decoder(; system::Symbol=:UTURUniChannelCSK, args...)
  @assert in(system, mcml_types) "MCML decoder not defined for $system"

  if system == :UTURUniChannelCSK
    UTURUniChannelCSKMCMLDecoder(; args...)
  end
end

function decode{T<:Real}(d::CorDecoder, x::Vector{T}, y::Vector{T})
  dec = dot(x, y)
  if dec > 0
    return 1
  elseif dec < 0
    return -1
  elseif dec == 0
    return -2 
  end
end
