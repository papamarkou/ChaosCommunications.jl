function sim_ber(s::DigitalSystem, bit::Int, n::Int64)
  estimation_errors::Int64 = 0
  decoding_failures::Int64 = 0

  for i = 1:n
    bit_estimate = sim_sys(s, bit)
    if in(bit_estimate, [-1, 1])
      if bit_estimate != bit
        estimation_errors += 1
      end
    else
      decoding_failures += 1
    end
  end

  estimation_errors, decoding_failures, estimation_errors/(n-decoding_failures)
end

function sim_ber(s::Vector{DigitalSystem}, bit::Int, n::Int64)
  slen = length(s)
  output = cell(slen)

  for i = 1:slen
    output[i] = sim_ber(s[i], bit, n)
  end
end

function psim_ber(s::DigitalSystem, bit::Int, n::Int64)
  decoding_failures, estimation_errors =
  (@parallel (+) for i = 1:n
    bit_estimate = sim_sys(s, bit)
    in(bit_estimate, [-1, 1]) ? (bit_estimate == bit ? [0, 0]: [0, 1]) : [1, 0]
  end)

  estimation_errors, decoding_failures, estimation_errors/(n-decoding_failures)
end

function psim_ber(s::Vector{DigitalSystem}, bit::Int, n::Int64)
  slen = length(s)
  output = cell(slen)

  for i = 1:slen
    output[i] = psim_ber(s[i], bit, n)
  end
end
