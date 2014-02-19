immutable UTURUniChannelCSK <: DigitalUTURUniChannelSystem
  coherent::Bool
  carrier::Carrier
  noise::Distribution
  decoder::Decoder

  function UTURUniChannelCSK(coh::Bool, c::Carrier, n::Distribution, d::Decoder)
    if isa(n, Normal)
      @assert n.μ == 0. "Normal noise with μ!=0 not supported in UTURUniChannelCSK."

      if coh
        @assert isa(d, CorDecoder) "$d decoder not supported in coherent UTURUniChannelCSK."
      else
        @assert isa(d, CorDecoder) || isa(d, MCMLDecoder) "$d decoder not supported in non-coherent UTURUniChannelCSK."
      end
    else
      throw(ArgumentError("$n channel noise not supported in UTURUniChannelCSK."))
    end

    new(coh, c, n, d)
  end
end

UTURUniChannelCSK(coh::Bool, c::Carrier; noise::Distribution=Normal(), decoder::Decoder=CorDecoder()) =
  UTURUniChannelCSK(coh, c, noise, decoder)

UTURUniChannelCSK(; coherent::Bool=false, carrier::Carrier=LogisticCarrier(5), noise::Distribution=Normal(),
  decoder::Decoder=CorDecoder()) = UTURUniChannelCSK(coherent, carrier, noise, decoder)

function uturunichannelcsk_snr2var(snr::Float64; carrier::Carrier=LogisticCarrier(5))
  var(carrier.pdf)/(10.0^(snr/10.0))
end

function uturunichannelcsk_var2snr(v::Float64; carrier::Carrier=LogisticCarrier(5))
  10.0*log10(var(carrier.pdf)/v)
end

channel_var_estimator(d::UTURUniChannelCSKMCMLDecoder, tsum::Float64) = 0.5*tsum/d.carrier.len
channel_var_estimator(d::UTURUniChannelCSKMCMLDecoder, t::Vector{Float64}) = channel_var_estimator(d, sum(t))

function logmcml(d::UTURUniChannelCSKMCMLDecoder, bit::Int, v::Float64, r1::Vector{Float64}, r2::Vector{Float64},
  tsum::Float64)
  -0.5*tsum/v-d.carrier.len(log(v)+log2π)-log(d.nmc)
end
logmcml(d::UTURUniChannelCSKMCMLDecoder, bit::Int, v::Float64, r1::Vector{Float64}, r2::Vector{Float64},
  t::Vector{Float64}) = logmcml(d, bit, v, r1, r2, sum(t))

function decode(d::UTURUniChannelCSKMCMLDecoder, r1::Vector{Float64}, r2::Vector{Float64})
  tp1, tm1 = Array(Float64, d.nmc), Array(Float64, d.nmc)
  x = Array(Float64, d.carrier.len)

  for i = 1:d.nmc
    generate!(d.carrier, x)
    x = x-mean(x)
    tp1[i] = dor(r1-x, r1-x)+dot(r2-x, r2-x)
    tm1[i] = dor(r1+x, r1+x)+dot(r2-x, r2-x)    
  end

  tp1sum = sum(tp1)
  tm1sum = sum(tm1)

  varp1 = channel_var_estimator(d, tp1sum)
  varm1 = channel_var_estimator(d, tm1sum)

  dec = logmcml(d, 1, varp1, r1, r2, tp1sum)-logmcml(d, -1, varm1, r1, r2, tm1sum)
  if dec > 0
    return 1
  elseif dec < 0
    return -1
  elseif dec == 0
    return -2 
  end
end

function sim_sys(s::UTURUniChannelCSK, bit::Int)
  @assert in(bit, [-1, 1]) "Bit must be either -1 or 1."

  x = Array(Float64, s.carrier.len)
  generate!(s.carrier, x)
  x = x-mean(x)
  r1 = bit*x+rand(s.noise, s.carrier.len)
  r2 = s.coherent ? x : x+rand(s.noise, s.carrier.len)
  decode(s.decoder, r1, r2)
end

function sim_ber(s::UTURUniChannelCSK, bit::Int, n::Int64)
  estimation_errors::Int64 = 0
  decoding_failures::Int64 = 0

  for i = 1:n
    bit_estimate = sim_sys(s, bit)
    if in(bit_estimate, [-1, 1])
      if bit_estimate != bit
        estimation_errors += 1
      end
    else
      decoding_failures += 1
    end
  end

  estimation_errors/(n-decoding_failures)
end

function psim_ber(s::UTURUniChannelCSK, bit::Int, n::Int64)
  decoding_failures, estimation_errors =
  (@parallel (+) for i = 1:n
    bit_estimate = sim_sys(s, bit)
    in(bit_estimate, [-1, 1]) ? (bit_estimate == bit ? [0, 0]: [0, 1]) : [1, 0]
  end)

  estimation_errors/(n-decoding_failures)
end
