function gen_sys(ebn0db::Ranges{Float64}, sprlen::Ranges{Int}; system::Symbol=:UTURUniChannelCSK, args...)
  if system == :UTURUniChannelCSK
    uturunichannelcsk_gen_sys(ebn0db, sprlen; args...)
  else
    error("sim_ber was called with invalid arguments.")
  end
end

gen_sys(ebn0db::Float64, sprlen::Ranges{Int}; system::Symbol=:UTURUniChannelCSK, args...) =
  gen_sys(ebn0db:ebn0db, sprlen; system=system, args...)

gen_sys(ebn0db::Ranges{Float64}, sprlen::Int; system::Symbol=:UTURUniChannelCSK, args...) =
  gen_sys(ebn0db, sprlen:sprlen; system=system, args...)

gen_sys(ebn0db::Float64, sprlen::Int; system::Symbol=:UTURUniChannelCSK, args...) =
  gen_sys(ebn0db:ebn0db, sprlen:sprlen; system=system, args...)

function sim_ber(s::DigitalSystem, bit::Int, n::Int64)
  nbiterrors::Int64 = 0
  ndecfails::Int64 = 0

  for i = 1:n
    bit_estimate = sim_sys(s, bit)
    if in(bit_estimate, [-1, 1])
      if bit_estimate != bit
        nbiterrors += 1
      end
    else
      ndecfails += 1
    end
  end

  nbiterrors, ndecfails, nbiterrors/(n-ndecfails)
end

function sim_ber{T<:DigitalSystem}(s::Vector{T}, bit::Int, n::Int64)
  slen = length(s)
  output = cell(slen)

  for i = 1:slen
    output[i] = sim_ber(s[i], bit, n)
  end

  return output
end

function sim_ber{T<:DigitalSystem}(s::Matrix{T}, bit::Int, n::Int64)
  nrows, ncols = size(s)
  output = cell(nrows, ncols)

  for i = 1:nrows
    for j = 1:ncols
      output[i, j] = sim_ber(s[i, j], bit, n)
    end
  end

  return output
end

function pfor_sim_ber(s::DigitalSystem, bit::Int, n::Int64)
  ndecfails, nbiterrors =
  (@parallel (+) for i = 1:n
    bit_estimate = sim_sys(s, bit)
    in(bit_estimate, [-1, 1]) ? (bit_estimate == bit ? [0, 0]: [0, 1]) : [1, 0]
  end)

  nbiterrors, ndecfails, nbiterrors/(n-ndecfails)
end

function pmap_sim_ber(s::DigitalSystem, bit::Int, n::Int64)
  nbiterrors::Int64 = 0
  ndecfails::Int64 = 0
  sim_function(bit::Int) = sim_sys(s, bit)

  bit_estimates = pmap(sim_function, fill(bit, n))

  for i = 1:n
    if in(bit_estimates[i], [-1, 1])
      if bit_estimates[i] != bit
        nbiterrors += 1
      end
    else
      ndecfails += 1
    end
  end

  nbiterrors, ndecfails, nbiterrors/(n-ndecfails)
end

function psim_ber(s::DigitalSystem, bit::Int, n::Int64; ptype::Symbol=:pmap)
  if ptype == :pmap
    pmap_sim_ber(s, bit, n)
  elseif ptype == :pfor
    pfor_sim_ber(s, bit, n)
  else
    error("$ptype parallel implementation not known.")
  end
end

function psim_ber{T<:DigitalSystem}(s::Vector{T}, bit::Int, n::Int64; ptype::Symbol=:pmap)
  slen = length(s)
  output = cell(slen)

  for i = 1:slen
    output[i] = psim_ber(s[i], bit, n; ptype=ptype)
  end

  return output
end

function psim_ber{T<:DigitalSystem}(s::Matrix{T}, bit::Int, n::Int64; ptype::Symbol=:pmap)
  nrows, ncols = size(s)
  output = cell(nrows, ncols)

  for i = 1:nrows
    for j = 1:ncols
      output[i, j] = psim_ber(s[i, j], bit, n; ptype=ptype)
    end
  end

  return output
end

function ber_lb{T<:DigitalSystem}(s::Vector{T}; args...)
  slen = length(s)
  output = Array(Float64, slen)

  for i = 1:slen
    output[i] = ber_lb(s[i]; args...)
  end

  return output
end

function ber_lb{T<:DigitalSystem}(s::Matrix{T}; args...)
  nrows, ncols = size(s)
  output = Array(Float64, nrows, ncols)

  for i = 1:nrows
    for j = 1:ncols
      output[i, j] = ber_lb(s[i, j]; args...)
    end
  end

  return output
end
