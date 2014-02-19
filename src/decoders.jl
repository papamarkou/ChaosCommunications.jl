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

immutable UTURUniChannelCSKMCMLDecoder <: MCMLDecoder
  carrier::Carrier
  nmc::Int64
end

mcmltypes = (:UTURUniChannelCSK,)

function mcml_decoder(; system::Symbol=:UTURUniChannelCSK, args...)
  @assert in(system, mcmltypes) "MCML decoder not defined for $system"

  if system == :UTURUniChannelCSK
    UTURUniChannelCSKMCMLDecoder(snr; args...)
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
