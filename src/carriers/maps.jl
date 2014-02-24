logistic(x::Float64) = 4*x*(1-x)

function bernoulli(x::Float64, l::Float64, u::Float64, nc::Float64)
  if l<=x && x<nc
    return ((u-l)*x+(nc-u)*l)/(nc-l)
  elseif nc<=x && x<=u
    return ((u-l)*x+(l-nc)*u)/(u-nc)
  else
    error("Input out of domain of Bernoulli map.")
  end
end

function nbernoulli(x::Float64, l::Float64, u::Float64, nc::Float64)
  if l<=x && x<nc
    return ((l-u)*x+(nc*u-l^2))/(nc-l)
  elseif nc<=x && x<=u
    return ((l-u)*x+(u^2-l*nc))/(u-nc)
  else
    error("Input out of domain of negative Bernoulli map.")
  end
end

function tent(x::Float64, l::Float64, u::Float64, nc::Float64)
  if l<=x && x<nc
    ((u-l)*x+(nc-u)*l)/(nc-l)
  elseif nc<=x && x<=u
    ((l-u)*x+(u^2-l*nc))/(u-nc)
  else
    error("Input out of domain of tent map.")
  end
end

function valley(x::Float64, l::Float64, u::Float64, nc::Float64)
  if l<=x && x<nc
    ((l-u)*x+(nc*u-l^2))/(nc-l)
  elseif nc<=x && x<=u
    ((u-l)*x+(l-nc)*u)/(u-nc)
  else
    error("Input out of domain of valley map.")
  end
end

function circular(x::Float64, nc::Float64)
  if -1.<=x && x<-sqrt(nc)
    y = -sqrt((1.-x^2)/(1.-nc))
  elseif -sqrt(nc)<=x && x<sqrt(nc)
    y = sqrt(1.-x^2/nc)
  elseif sqrt(nc)<=x && x<=1.
    y = -sqrt((1.-x^2)/(1.-nc))
  else
    error("Input out of domain of circular map.")
  end
end

immutable VDist <: ContinuousUnivariateDistribution
    nc::Float64

    function VDist(nc::Real)
      @assert -1.0<nc && nc<1.0 "Non-centrality parameter must be in (-1, 1) in the case circular map."
      new(float64(nc))
    end

    VDist() = VDist(0.0)
end

@continuous_distr_support VDist -1.0 1.0

function pdf(d::VDist, x::Real)
  if -1.0<=x && x<=0.0
    -2.0*(1.0-d.nc)*x 
  elseif 0.0<x && x<=1.0
    2*d.nc*x
  else
    0.0
  end
end

function rand(d::VDist)
  u = rand()

  if 0.0<=u && u<1.0-d.nc
    -sqrt((u+d.nc-1.0)/(d.nc-1.0))
  elseif 1.0-d.nc<=u && u<=1.0
    sqrt((u+d.nc-1.0)/d.nc)
  end
end

mean(d::VDist) = (4.0*d.nc-2.0)/3.0

var(d::VDist) = 0.5-mean(d)^2
