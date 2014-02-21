using Distributions, ChaosCommunications

carrier = LogisticCarrier(4)

noise_var = snr2var(5., system=:UTURUniChannelCSK, carrier=carrier)

system = UTURUniChannelCSK(true, carrier, noise=Normal(0., sqrt(noise_var)))

@time sim_ber(system, 1, 100000)
# @time psim_ber(system, 1, 100000)
