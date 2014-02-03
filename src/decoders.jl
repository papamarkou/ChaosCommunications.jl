abstract Decoder
abstract DeterministicDecoder <: Decoder
abstract QuasiRandomDecoder <: Decoder
abstract RandomDecoder <: Decoder
# abstract EMDecoder <: QuasiRandomDecoder
# abstract GibbsDecoder <: RandomDecoder
# abstract ParticleFilterDecoder <: RandomDecoder

type CorDecoder <: DeterministicDecoder
end

type MCMLDecoder <: QuasiRandomDecoder
  likelihood::Function
  maximize::Function
  optimizer::Dict
  nmc::Int
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
