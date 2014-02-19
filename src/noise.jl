snr_types = (:UTURUniChannelCSK,)

function snr2var(snr::Float64; system::Symbol=:UTURUniChannelCSK, args...)
  @assert in(system, snr_types) "SNR not defined for $system"

  if system == :UTURUniChannelCSK
    uturunichannelcsk_snr2var(snr; args...)
  end
end

function var2snr(v::Float64; system::Symbol=:UTURUniChannelCSK, args...)
  @assert in(system, snr_types) "SNR not defined for $system"

  if system == :UTURUniChannelCSK
    uturunichannelcsk_var2snr(v; args...)
  end
end
