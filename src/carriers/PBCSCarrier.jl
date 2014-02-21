immutable PBCSCarrier <: QuasiRandomCarrier
  c::Float64 # center of underlying circle
  len::Int
end

PBCSCarrier(len::Int) = PBCSCarrier(0.0, len)

function invcdf_pbcs(c::PBCSCarrier, x::Float64)
  if ((0<=x) && (x<0.5))
    -sqrt(1-2*x)+c.c
  elseif ((0.5<=x) && (x<=1))
    sqrt(2*x-1)+c.c
  else
    error("Input out of domain of inverse cdf of PBCS")
  end
end

function generate!(c::PBCSCarrier, x::Vector{Float64})
  local a::Int
  dunif = Uniform(0, 1)
  dbern = Bernoulli()

  for i = 1:2:c.len
    x[i] = rand(dunif)
  end

  for i = 2:2:c.len
    a = rand(dbern) == 0 ? -1 : 1
    x[i] = c.c+a*sqrt(1-(x[i-1]-c.c)^2)
  end

  return x
end

generate(c::PBCSCarrier) = generate!(c, Array(Float64, c.len))

mean(c::PBCSCarrier) = c.c

var(c::PBCSCarrier) = 0.5
