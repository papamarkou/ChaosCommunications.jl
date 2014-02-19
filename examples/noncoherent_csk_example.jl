using Distributions, ChaosCommunications

carrier = LogisticCarrier(4)

noise_var = snr2var(1., system=:UTURUniChannelCSK, carrier=carrier)

decoder = mcml_decoder(system=:UTURUniChannelCSK, carrier=carrier)

system = UTURUniChannelCSK(false, carrier, noise=Normal(0., sqrt(noise_var)), decoder=decoder)

@time sim_ber(system, 1, 100000)
# @time psim_ber(system, 1, 1000000)
