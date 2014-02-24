initialize(c::IterativeMapCarrier) = rand(c.density)

function generate!(c::IterativeMapCarrier, x::Vector{Float64})
  x[1] = initialize(c)

  for i = 2:c.len
    x[i] = c.iterate(x[i-1])
  end

  return x
end

generate(c::IterativeMapCarrier) = generate!(c, Array(Float64, c.len))

mean(c::IterativeMapCarrier) = mean(c.density)

var(c::IterativeMapCarrier) = var(c.density)

immutable LogisticCarrier <: IterativeMapCarrier
  density::Distribution
  len::Int
  iterate::Function
end

LogisticCarrier(len::Int) = LogisticCarrier(Beta(0.5, 0.5), len, logistic)
LogisticCarrier(; density::Distribution=Beta(0.5, 0.5), len::Int=5, iterate::Function=logistic) =
  LogisticCarrier(density, len, iterate)

immutable FSPWL2BCarrier <: FSPWLCarrier # 2-branch fully-stretching piece-wise linear maps
  l::Float64 # minimum of domain
  u::Float64 # maximum of domain
  nc::Float64 # non-centrality parameter
  density::Distribution
  len::Int
  iterate::Function

  function FSPWL2BCarrier(l::Float64, u::Float64, nc::Float64, density::Distribution, len::Int, iterate::Function)
    @assert l<nc && nc<u "Non-centrality parameter must be in the interior of domain of Bernoulli map."
    new(l, u, nc, density, len, iterate)
  end
end

fspwl2b_types = (:bernoulli, :nbernoulli, :tent, :valley)

function FSPWL2BCarrier(l::Float64, u::Float64, nc::Float64, len::Int, mtype::Symbol)
  @assert in(mtype, fspwl2b_types) "$mtype is not a supported type of FSPWL2BCarrier map."

  if mtype == :bernoulli
    return FSPWL2BCarrier(l, u, nc, Uniform(l, u), len, (x::Float64)->bernoulli(x, l, u, nc))
  elseif mtype == :nbernoulli
    return FSPWL2BCarrier(l, u, nc, Uniform(l, u), len, (x::Float64)->nbernoulli(x, l, u, nc))
  elseif mtype == :tent
    return FSPWL2BCarrier(l, u, nc, Uniform(l, u), len, (x::Float64)->tent(x, l, u, nc))
  elseif mtype == :valley
    return FSPWL2BCarrier(l, u, nc, Uniform(l, u), len, (x::Float64)->valley(x, l, u, nc))    
  end
end

FSPWL2BCarrier(; l::Float64=0.0, u::Float64=1.0, nc::Float64=0.5, density::Distribution=Uniform(0.0, 1.0), len::Int=5,
  iterate::Function=(x::Float64)->bernoulli(x, l, u, nc)) = FSPWL2BCarrier(l, u, nc, density, len, iterate)

BernoulliCarrier(l::Float64, u::Float64, nc::Float64, len::Int) = FSPWL2BCarrier(l, u, nc, len, :bernoulli)
BernoulliCarrier(l::Float64, u::Float64, len::Int) = BernoulliCarrier(l, u, 0.5*(l-u), len)
BernoulliCarrier(len::Int) = BernoulliCarrier(0.0, 1.0, 0.5, len)

NBernoulliCarrier(l::Float64, u::Float64, nc::Float64, len::Int) = FSPWL2BCarrier(l, u, nc, len, :nbernoulli)
NBernoulliCarrier(l::Float64, u::Float64, len::Int) = NBernoulliCarrier(l, u, 0.5*(l-u), len)
NBernoulliCarrier(len::Int) = NBernoulliCarrier(0.0, 1.0, 0.5, len)

TentCarrier(l::Float64, u::Float64, nc::Float64, len::Int) = FSPWL2BCarrier(l, u, nc, len, :tent)
TentCarrier(l::Float64, u::Float64, len::Int) = TentCarrier(l, u, 0.5*(l-u), len)
TentCarrier(len::Int) = TentCarrier(0.0, 1.0, 0.5, len)

ValleyCarrier(l::Float64, u::Float64, nc::Float64, len::Int) = FSPWL2BCarrier(l, u, nc, len, :valley)
ValleyCarrier(l::Float64, u::Float64, len::Int) = ValleyCarrier(l, u, 0.5*(l-u), len)
ValleyCarrier(len::Int) = ValleyCarrier(0.0, 1.0, 0.5, len)

immutable CircularCarrier <: IterativeMapCarrier
  nc::Float64 # non-centrality parameter
  density::Distribution
  len::Int
  iterate::Function

  function CircularCarrier(nc::Float64, density::Distribution, len::Int, iterate::Function)
    @assert -1.0<nc && nc<1.0 "Non-centrality parameter must be in the interior of domain of circular map."
    new(nc, density, len, iterate)
  end
end

CircularCarrier(nc::Float64, len::Int) = CircularCarrier(nc, VDist(nc), len, (x::Float64)->circular(x, nc))
CircularCarrier(len::Int) = CircularCarrier(0.42, VDist(0.42), len, (x::Float64)->circular(x, 0.42))
CircularCarrier(; nc::Float64=0.42, density::Distribution=VDist(0.42), len::Int=5,
  iterate::Function=(x::Float64)->circular(x, 0.42)) = CircularCarrier(nc, density, len, iterate)
