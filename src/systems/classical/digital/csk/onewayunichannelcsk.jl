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

function uturunichannelcsk_snr2var(snr::Float64; carrier::Carrier=LogisticCarrier(5))
  var(carrier.pdf)/(10.0^(snr/10.0))
end

function uturunichannelcsk_var2snr(v::Float64; carrier::Carrier=LogisticCarrier(5))
  10.0*log10(var(carrier.pdf)/v)
end
