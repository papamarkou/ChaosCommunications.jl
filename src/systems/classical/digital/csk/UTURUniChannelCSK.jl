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

uturunichannelcsk_snrdb2var(s::Float64; carrier::Carrier=LogisticCarrier(5)) = var(carrier)/(10.0^(s/10.0))

uturunichannelcsk_var2snrdb(v::Float64; carrier::Carrier=LogisticCarrier(5)) = 10.0*log10(var(carrier)/v)

uturunichannelcsk_ebn0db2var(eb::Float64; carrier::Carrier=LogisticCarrier(5)) =
  carrier.len*var(carrier)/(10.0^(eb/10.0))

uturunichannelcsk_var2ebn0db(v::Float64; carrier::Carrier=LogisticCarrier(5)) = 10.0*log10(carrier.len*var(carrier)/v)

uturunichannelcsk_ebn0db2snrdb(eb::Float64; carrier::Carrier=LogisticCarrier(5)) = eb-10.0*log10(carrier.len)

uturunichannelcsk_snrdb2ebn0db(s::Float64; carrier::Carrier=LogisticCarrier(5)) = s+10.0*log10(carrier.len)

function logmcml(d::UTURUniChannelCSKMCMLDecoder, v::Float64, t::Vector{Float64})
  logsumexp(-0.5*t/v)-d.carrier.len*(log(v)+log2π)-log(d.nmc)
end

function gradlogmcml(d::UTURUniChannelCSKMCMLDecoder, v::Float64, t::Vector{Float64})
  expterms = exp(-0.5*t/v)
  (0.5*dot(t, expterms)/(sum(expterms)*v)-d.carrier.len)/v
end

function decode(d::UTURUniChannelCSKMCMLDecoder, r1::Vector{Float64}, r2::Vector{Float64})
  tp1, tm1 = Array(Float64, d.nmc), Array(Float64, d.nmc)
  x = Array(Float64, d.carrier.len)
  local maxfp1, maxxp1, retp1, maxfm1, maxxm1, retm1
  opt = Array(Opt, 2)

  for i = 1:d.nmc
    generate!(d.carrier, x)
    subtract!(x, mean(d.carrier))
    dotr2 = sumsqdiff(r2, x)
    tp1[i] = sumsqdiff(r1, x)+dotr2
    tm1[i] = sumsq(r1+x)+dotr2
  end

  function objective_fp1(v::Vector{Float64}, grad::Vector{Float64})
    if length(grad) > 0
      grad[1] = gradlogmcml(d, v[1], tp1)
    end
    logmcml(d, v[1], tp1)
  end

  opt[1] = convert(Opt, d.opt[1])
  max_objective!(opt[1], objective_fp1)

  try
    (maxfp1, maxxp1, retp1) = optimize(opt[1], [d.init[1]])
  catch
    return -3
  end

  if in(retp1, opt_failure)
    return -4
  end

  function objective_fm1(v::Vector{Float64}, grad::Vector{Float64})
    if length(grad) > 0
      grad[1] = gradlogmcml(d, v[1], tm1)
    end
    logmcml(d, v[1], tm1)
  end

  opt[2] = convert(Opt, d.opt[2])
  max_objective!(opt[2], objective_fm1)

  try
    (maxfm1, maxxm1, retm1) = optimize(opt[2], [d.init[2]])
  catch
    return -5
  end

  if in(retm1, opt_failure)
    return -6
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
  subtract!(x, mean(s.carrier))
  r1 = bit*x+rand(s.noise, s.carrier.len)
  r2 = s.coherent ? x : x+rand(s.noise, s.carrier.len)
  decode(s.decoder, r1, r2)
end

function uturunichannelcsk_gen_sys(ebn0db::Ranges{Float64}, sprlen::Ranges{Int}; coherent::Bool=true,
  carrier::Carrier=LogisticCarrier(), decoder::Decoder=CorDecoder())
  systems = Array(UTURUniChannelCSK, ebn0db.len, sprlen.len)
  local d

  for i = 1:ebn0db.len
    for j = 1:sprlen.len
      ckwargs = Dict()
      for k in typeof(carrier).names
        ckwargs[k] = k != :len ? carrier.(k) : sprlen[j] 
      end
      c = typeof(carrier)(; ckwargs...)

      v = ebn0db2var(ebn0db[i]; system=:UTURUniChannelCSK, carrier=c)

      if isa(decoder, CorDecoder)
        d = decoder
      elseif isa(decoder, MCMLDecoder)
        dkwargs = Dict()
        for k in typeof(decoder).names
          dkwargs[k] = k != :carrier ? decoder.(k) : c 
        end
        d = mcml_decoder(system=:UTURUniChannelCSK; dkwargs...)
      else
        error("$decoder is not a supported type of decoder for the uturunichannelcsk_gen_sys function.")
      end

      systems[i, j] = UTURUniChannelCSK(coherent, c; noise=Normal(0., sqrt(v)), decoder=d)
    end
  end

  return systems
end

function sim_ber(s::Vector{UTURUniChannelCSK}, bit::Int, n::Int64)
  slen = length(s)
  nbiterrors, ndecfails, bers = Array(Int64, slen), Array(Int64, slen), Array(Float64, slen)

  for i = 1:slen
    try
      nbiterrors[i], ndecfails[i], bers[i] = sim_ber(s[i], bit, n)
    catch
      nbiterrors[i] = ndecfails[i] = bers[i] = NaN
    end
  end

  return nbiterrors, ndecfails, bers
end

function sim_ber(s::Matrix{UTURUniChannelCSK}, bit::Int, n::Int64)
  nrows, ncols = size(s)
  nbiterrors, ndecfails, bers = Array(Int64, nrows, ncols), Array(Int64, nrows, ncols), Array(Float64, nrows, ncols)

  for i = 1:nrows
    for j = 1:ncols
      try
        nbiterrors[i, j], ndecfails[i, j], bers[i, j] = sim_ber(s[i, j], bit, n)
      catch
        nbiterrors[i, j] = ndecfails[i, j] = bers[i, j] = NaN
      end
    end
  end

  return nbiterrors, ndecfails, bers
end

function psim_ber(s::Vector{UTURUniChannelCSK}, bit::Int, n::Int64; ptype::Symbol=:pmap)
  slen = length(s)
  nbiterrors, ndecfails, bers = Array(Int64, slen), Array(Int64, slen), Array(Float64, slen)

  for i = 1:slen
    try
      nbiterrors[i], ndecfails[i], bers[i] = psim_ber(s[i], bit, n; ptype=ptype)
    catch
      nbiterrors[i] = ndecfails[i] = bers[i] = NaN
    end
  end

  return nbiterrors, ndecfails, bers
end

function psim_ber(s::Matrix{UTURUniChannelCSK}, bit::Int, n::Int64; ptype::Symbol=:pmap)
  nrows, ncols = size(s)
  nbiterrors, ndecfails, bers = Array(Int64, nrows, ncols), Array(Int64, nrows, ncols), Array(Float64, nrows, ncols)

  for i = 1:nrows
    for j = 1:ncols
      try
        nbiterrors[i, j], ndecfails[i, j], bers[i, j] = psim_ber(s[i, j], bit, n; ptype=ptype)
      catch
        nbiterrors[i, j] = ndecfails[i, j] = bers[i, j] = NaN
      end
    end
  end

  return nbiterrors, ndecfails, bers
end

function ber_lb(s::UTURUniChannelCSK; btype::Symbol=:jensen)
  if btype == :jensen
    if isa(s.noise, Normal) && s.noise.μ == 0.0
      cdf(Normal(), -sqrt(s.carrier.len*var(s.carrier))/s.noise.σ)
    else
      error("Jensen's BER lower bound defined for UTURUniChannelCSK only for Normal noise with zero mean.")
    end
  else
    error("$btype type of BER lower bound not known for UTURUniChannelCSK.")
  end
end
