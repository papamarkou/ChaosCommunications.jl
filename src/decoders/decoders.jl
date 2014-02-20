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
