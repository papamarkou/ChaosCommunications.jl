using Distributions, ChaosCommunications

carrier = LogisticCarrier(5)

noise_var = snr2var(1., system=:UTURUniChannelCSK, carrier=carrier)

system = UTURUniChannelCSK(true, carrier, noise=Normal(0., sqrt(noise_var)))

@time sim_ber(system, 1, 1000000)
# @time psim_ber(system, 1, 1000000)
