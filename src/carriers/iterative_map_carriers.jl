initialize(c::IterativeMapCarrier) = rand(c.pdf)

function generate!(c::IterativeMapCarrier, x::RealVector)
  x[1] = initialize(c)

  for i = 2:c.len
    x[i] = iterate(c, x[i-1])
  end

  return x
end

generate(c::IterativeMapCarrier) = generate!(c, Array(Real, c.len))

mean(c::IterativeMapCarrier) = mean(c.pdf)

immutable LogisticCarrier <: IterativeMapCarrier
  pdf::Distribution
  len::Int
end

LogisticCarrier(len::Int) = LogisticCarrier(Beta(0.5, 0.5), len) 

function iterate(c::LogisticCarrier, x::Float64)
  4*x*(1-x)
end

immutable BernoulliCarrier <: IterativeMapCarrier
  l::Float64 # minimum of domain
  u::Float64 # maximum of domain
  nc::Float64 # non-centrality parameter
  pdf::Distribution
  len::Int

  function BernoulliCarrier(l::Float64, u::Float64, nc::Float64, pdf::Distribution, len::Int)
    @assert l<nc && nc<u "Non-centrality parameter must be in the interior of domain of Bernoulli map."
    new(l, u, nc, pdf, len)
  end
end

BernoulliCarrier(l::Float64, u::Float64, nc::Float64, len::Int) = BernoulliCarrier(l, u, nc, Uniform(l, u), len)
BernoulliCarrier(l::Float64, u::Float64, len::Int) = BernoulliCarrier(l, u, 0.5*(l-u), len)
BernoulliCarrier(len::Int) = BernoulliCarrier(0., 1., 0.5, Uniform(l, u), len)

function iterate(c::BernoulliCarrier, x::Float64)
  if c.l<=x && x<c.nc
    return ((c.u-c.l)*x+(c.nc-c.u)*c.l)/(c.nc-c.l)
  elseif c.nc<=x && x<=c.u
    return ((c.u-c.l)*x+(c.l-c.nc)*c.u)/(c.u-c.nc)
  else
    error("Input out of domain of Bernoulli map")
  end
end
