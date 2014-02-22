snrdb_types = (:UTURUniChannelCSK,)

function snrdb2var(s::Float64; system::Symbol=:UTURUniChannelCSK, args...)
  @assert in(system, snrdb_types) "SNR not defined for $system"

  if system == :UTURUniChannelCSK
    uturunichannelcsk_snrdb2var(s; args...)
  end
end

function var2snrdb(v::Float64; system::Symbol=:UTURUniChannelCSK, args...)
  @assert in(system, snrdb_types) "SNR not defined for $system"

  if system == :UTURUniChannelCSK
    uturunichannelcsk_var2snrdb(v; args...)
  end
end

function ebn0db2var(eb::Float64; system::Symbol=:UTURUniChannelCSK, args...)
  @assert in(system, snrdb_types) "SNR not defined for $system"

  if system == :UTURUniChannelCSK
    uturunichannelcsk_ebn0db2var(eb; args...)
  end
end

function var2ebn0db(v::Float64; system::Symbol=:UTURUniChannelCSK, args...)
  @assert in(system, snrdb_types) "SNR not defined for $system"

  if system == :UTURUniChannelCSK
    uturunichannelcsk_var2ebn0db(v; args...)
  end
end

function ebn0db2snrdb(eb::Float64; system::Symbol=:UTURUniChannelCSK, args...)
  @assert in(system, snrdb_types) "SNR not defined for $system"

  if system == :UTURUniChannelCSK
    uturunichannelcsk_ebn0db2snrdb(eb; args...)
  end
end

function snrdb2ebn0db(s::Float64; system::Symbol=:UTURUniChannelCSK, args...)
  @assert in(system, snrdb_types) "SNR not defined for $system"

  if system == :UTURUniChannelCSK
    uturunichannelcsk_snrdb2ebn0db(s; args...)
  end
end
