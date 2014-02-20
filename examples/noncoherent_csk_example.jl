using Distributions, NLopt, ChaosCommunications

println("    Testing correlation decoder of non-coherent single user CSK...")

carrier = LogisticCarrier(4)

noise_var = snr2var(10., system=:UTURUniChannelCSK, carrier=carrier)

system = UTURUniChannelCSK(false, carrier, noise=Normal(0., sqrt(noise_var)))

@time sim_ber(system, 1, 100000)

println("    Testing MCML decoder of non-coherent single user CSK...")

opt = Opt(:LD_SLSQP, 1)
lower_bounds!(opt, [0.])
upper_bounds!(opt, [10.])
ftol_abs!(opt, 1e-32)

decoder = mcml_decoder(system=:UTURUniChannelCSK, nmc=50, opt=opt, carrier=carrier)

system = UTURUniChannelCSK(false, carrier, noise=Normal(0., sqrt(noise_var)), decoder=decoder)

@time sim_ber(system, 1, 1000)
