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

function logmcml(d::UTURUniChannelCSKMCMLDecoder, v::Float64, t::Vector{Float64})
  log(sum(exp(-0.5*t/v)))-d.carrier.len*(log(v)+log2π)-log(d.nmc)
end

function gradlogmcml(d::UTURUniChannelCSKMCMLDecoder, v::Float64, t::Vector{Float64})
  expterms = exp(-0.5*t/v)
  (0.5*dot(t, expterms)/(sum(expterms)*v)-d.carrier.len)/v
end

function decode(d::UTURUniChannelCSKMCMLDecoder, r1::Vector{Float64}, r2::Vector{Float64})
  tp1, tm1 = Array(Float64, d.nmc), Array(Float64, d.nmc)
  x = Array(Float64, d.carrier.len)

  for i = 1:d.nmc
    generate!(d.carrier, x)
    x = x-mean(d.carrier)
    dotr2 = dot(r2-x, r2-x)
    tp1[i] = dot(r1-x, r1-x)+dotr2
    tm1[i] = dot(r1+x, r1+x)+dotr2 
  end

  function objective_fp1(v::Vector{Float64}, grad::Vector{Float64})
    if length(grad) > 0
      grad[1] = gradlogmcml(d, v[1], tp1)
    end
    logmcml(d, v[1], tp1)
  end
  max_objective!(d.opt[1], objective_fp1)
  (maxfp1, maxxp1, retp1) = optimize(d.opt[1], [d.init[1]])

  if in(retp1, opt_failure)
    return -3
  end

  function objective_fm1(v::Vector{Float64}, grad::Vector{Float64})
    if length(grad) > 0
      grad[1] = gradlogmcml(d, v[1], tm1)
    end
    logmcml(d, v[1], tm1)
  end
  max_objective!(d.opt[2], objective_fm1)
  (maxfm1, maxxm1, retm1) = optimize(d.opt[2], [d.init[2]])

  if in(retm1, opt_failure)
    return -3
  end

  dec = maxfp1-maxfm1

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
