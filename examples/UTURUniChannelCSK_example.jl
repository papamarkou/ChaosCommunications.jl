using Distributions, NLopt, ChaosCommunications

carrier = LogisticCarrier(4)

noise_var = ebn0db2var(5., system=:UTURUniChannelCSK, carrier=carrier)

println("    Testing correlation decoder of coherent single-user CSK...")

system = UTURUniChannelCSK(true, carrier, noise=Normal(0., sqrt(noise_var)))

sim_ber(system, 1, 100000)
# psim_ber(system, 1, 100000)

println("    Testing correlation decoder of non-coherent single-user CSK...")

system = UTURUniChannelCSK(false, carrier, noise=Normal(0., sqrt(noise_var)))

sim_ber(system, 1, 100000)

println("    Testing MCML decoder of non-coherent single user-CSK...")

decoder = mcml_decoder(system=:UTURUniChannelCSK, carrier=carrier)

system = UTURUniChannelCSK(false, carrier, noise=Normal(0., sqrt(noise_var)), decoder=decoder)

sim_ber(system, 1, 1000)
