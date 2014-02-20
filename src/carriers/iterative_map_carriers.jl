initialize(c::IterativeMapCarrier) = rand(c.density)

function generate!(c::IterativeMapCarrier, x::RealVector)
  x[1] = initialize(c)

  for i = 2:c.len
    x[i] = c.iterate(x[i-1])
  end

  return x
end

generate(c::IterativeMapCarrier) = generate!(c, Array(Real, c.len))

mean(c::IterativeMapCarrier) = mean(c.density)

immutable LogisticCarrier <: IterativeMapCarrier
  density::Distribution
  len::Int
  iterate::Function
end

LogisticCarrier(density::Distribution, len::Int) = LogisticCarrier(density, len, logistic)
LogisticCarrier(len::Int) = LogisticCarrier(Beta(0.5, 0.5), len, logistic)

immutable FSPWL2B <: FSPWL # 2-branch fully-stretching piece-wise linear maps
  l::Float64 # minimum of domain
  u::Float64 # maximum of domain
  nc::Float64 # non-centrality parameter
  density::Distribution
  len::Int
  iterate::Function

  function FSPWL2B(l::Float64, u::Float64, nc::Float64, density::Distribution, len::Int, iterate::Function)
    @assert l<nc && nc<u "Non-centrality parameter must be in the interior of domain of Bernoulli map."
    new(l, u, nc, density, len, mtype, iterate)
  end
end

fspwl2b_types = (:bernoulli, :nbernoulli, :tent, :valley)

function FSPWL2B(l::Float64, u::Float64, nc::Float64, len::Int, mtype::Symbol)
  @assert in(mtype, fspwl2b_types) "$mtype is not a supported type of FSPWL2B map."

  if mtype == :bernoulli
    return FSPWL2B(l, u, nc, Uniform(l, u), len, (x::Float64)->bernoulli(x, l, u, nc))
  elseif mtype == :nbernoulli
    return FSPWL2B(l, u, nc, Uniform(l, u), len, (x::Float64)->nbernoulli(x, l, u, nc))
  elseif mtype == :tent
    return FSPWL2B(l, u, nc, Uniform(l, u), len, (x::Float64)->tent(x, l, u, nc))
  elseif mtype == :valley
    return FSPWL2B(l, u, nc, Uniform(l, u), len, (x::Float64)->valley(x, l, u, nc))    
  end
end

BernoulliCarrier(l::Float64, u::Float64, nc::Float64, len::Int) = FSPWL2B(l, u, nc, len, :bernoulli)
BernoulliCarrier(l::Float64, u::Float64, len::Int) = BernoulliCarrier(l, u, 0.5*(l-u), len)
BernoulliCarrier(len::Int) = BernoulliCarrier(0., 1., 0.5, len)

NBernoulliCarrier(l::Float64, u::Float64, nc::Float64, len::Int) = FSPWL2B(l, u, nc, len, :nbernoulli)
NBernoulliCarrier(l::Float64, u::Float64, len::Int) = NBernoulliCarrier(l, u, 0.5*(l-u), len)
NBernoulliCarrier(len::Int) = NBernoulliCarrier(0., 1., 0.5, len)

TentCarrier(l::Float64, u::Float64, nc::Float64, len::Int) = FSPWL2B(l, u, nc, len, :tent)
TentCarrier(l::Float64, u::Float64, len::Int) = TentCarrier(l, u, 0.5*(l-u), len)
TentCarrier(len::Int) = TentCarrier(0., 1., 0.5, len)

ValleyCarrier(l::Float64, u::Float64, nc::Float64, len::Int) = FSPWL2B(l, u, nc, len, :valley)
ValleyCarrier(l::Float64, u::Float64, len::Int) = ValleyCarrier(l, u, 0.5*(l-u), len)
ValleyCarrier(len::Int) = ValleyCarrier(0., 1., 0.5, len)
