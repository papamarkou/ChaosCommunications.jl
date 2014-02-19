type IterativeMapCarrier <: QuasiRandomCarrier
  map::Function
  pdf::Distribution
  len::Int
end

initialize(c::IterativeMapCarrier) = rand(c.pdf)

function generate!(c::IterativeMapCarrier, x::RealVector)
  x[1] = initialize(c)

  for i = 2:c.len
    x[i] = c.map(x[i-1])
  end

  return x
end

generate(c::IterativeMapCarrier) = generate!(c, Array(Real, c.len))

mean(c::IterativeMapCarrier) = mean(c.pdf)

LogisticCarrier(len::Int) = IterativeMapCarrier(logistic, Beta(0.5, 0.5), len) 

function logistic(x::Float64)
  4*x*(1-x)
end
