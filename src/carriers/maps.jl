logistic(x::Float64) = 4*x*(1-x)

function bernoulli(x::Float64, l::Float64, u::Float64, nc::Float64)
  if l<=x && x<nc
    return ((u-l)*x+(nc-u)*l)/(nc-l)
  elseif nc<=x && x<=u
    return ((u-l)*x+(l-nc)*u)/(u-nc)
  else
    error("Input out of domain of Bernoulli map")
  end
end

function nbernoulli(x::Float64, l::Float64, u::Float64, nc::Float64)
  if l<=x && x<nc
    return ((l-u)*x+(nc*u-pow(l, 2)))/(nc-l)
  elseif nc<=x && x<=u
    return ((l-u)*x+(pow(u, 2)-l*nc))/(u-nc)
  else
    error("Input out of domain of negative Bernoulli map")
  end
end

function tent(x::Float64, l::Float64, u::Float64, nc::Float64)
  if l<=x && x<nc
    ((u-l)*x+(nc-u)*l)/(nc-l)
  elseif nc<=x && x<=u
    ((l-u)*x+(pow(u, 2)-l*nc))/(u-nc)
  else
    error("Input out of domain of tent map")
  end
end

function valley(x::Float64, l::Float64, u::Float64, nc::Float64)
  if l<=x && x<nc
    ((l-u)*x+(nc*u-pow(l, 2)))/(nc-l)
  elseif nc<=x && x<=u
    ((u-l)*x+(l-nc)*u)/(u-nc)
  else
    error("Input out of domain of valley map")
  end
end
