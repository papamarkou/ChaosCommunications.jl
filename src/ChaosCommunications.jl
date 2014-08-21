module ChaosCommunications

using NumericExtensions
using Distributions
using StatsBase
using NLopt
using Dates
using ProgressMeter

import Base.convert, Base.mean, Base.var, Base.rand
import Distributions.log2π, Distributions.@continuous_distr_support

export
  # Types
  Optimizer,
  NLoptimizer,
  InfoForm,
  Classical,
  Quantum,
  Analog,
  Digital,
  InteractionForm,
  TransmitterForm,
  ReceiverForm,
  UniTransmitter,
  MultiTransmitter,
  UniReceiver,
  MultiReceiver,
  OneWay,
  TwoWay,
  ChannelForm,
  UniChannel,
  MultiChannel,
  System,
  ClassicalSystem,
  QuantumSystem,
  AnalogSystem,
  DigitalSystem,
  DigitalOneWayUniChannelSystem,
  DigitalUTURUniChannelSystem,
  UTURUniChannelCSK,
  VDist,
  Carrier,
  RandomProcessCarrier,
  QuasiRandomCarrier,
  IterativeMapCarrier,
  FSPWLCarrier,
  FSPWL2BCarrier,
  LogisticCarrier,
  BernoulliCarrier, # EMulated type, it is in fact a function
  NBernoulliCarrier, # EMulated type, it is in fact a function
  TentCarrier, # EMulated type, it is in fact a function
  ValleyCarrier, # EMulated type, it is in fact a function
  CircularCarrier,
  PBCSCarrier,
  Decoder,
  DeterministicDecoder,
  QuasiRandomDecoder,
  RandomDecoder,
  CorDecoder,
  MCMLDecoder,

  # Functions
  logistic,
  bernoulli,
  nbernoulli,
  tent,
  valley,
  circular,
  initialize,
  generate!,
  generate,
  snrdb2var,
  var2snrdb,
  ebn0db2var,
  var2ebn0db,
  ebn0db2snrdb,
  snrdb2ebn0db,
  mcml_decoder,
  decode,
  gen_sys,
  sim_sys,
  sim_ber,
  psim_ber,
  ber_lb

abstract InfoForm
abstract Classical <: InfoForm
abstract Quantum <: InfoForm
immutable Analog <: Classical end
immutable Digital <: Classical end

abstract InteractionForm
abstract TransmitterForm
abstract ReceiverForm
immutable UniTransmitter <: TransmitterForm end
immutable MultiTransmitter <: TransmitterForm end
immutable UniReceiver <: ReceiverForm end
immutable MultiReceiver <: ReceiverForm end
abstract OneWay{T<:TransmitterForm, R<:ReceiverForm} <: InteractionForm
abstract TwoWay <: InteractionForm

abstract ChannelForm
immutable UniChannel <: ChannelForm end
immutable MultiChannel <:ChannelForm end

abstract System{Info<:InfoForm, I<:InteractionForm, C<:ChannelForm}

typealias ClassicalSystem{I<:InteractionForm, C<:ChannelForm} System{Classical, I, C}
typealias QuantumSystem{I<:InteractionForm, C<:ChannelForm} System{Quantum, I, C}
typealias AnalogSystem{I<:InteractionForm, C<:ChannelForm} System{Analog, I, C}
typealias DigitalSystem{I<:InteractionForm, C<:ChannelForm} System{Digital, I, C}

typealias DigitalOneWayUniChannelSystem{T<:TransmitterForm, R<:ReceiverForm} DigitalSystem{OneWay{T, R}, UniChannel}
# UTUR := UniTransmitter, UniReceiver
typealias DigitalUTURUniChannelSystem DigitalSystem{OneWay{UniTransmitter, UniReceiver}, UniChannel}

typealias F64OrNothing Union(Float64, Nothing)
typealias VF64OrNothing Union(Vector{Float64}, Nothing)
typealias FunctionOrNothing Union(Function, Nothing)

# type MinMaxError <: Exception
#     msg::String
# end

include(joinpath("carriers", "maps.jl"))
include(joinpath("optimizers", "NLoptimizer.jl"))

include(joinpath("carriers", "carriers.jl"))
include(joinpath("carriers", "iterative_map_carriers.jl"))
include(joinpath("carriers", "PBCSCarrier.jl"))
include(joinpath("decoders", "decoders.jl"))
include(joinpath("decoders", "mcml_decoders.jl"))
include("noise.jl")
include(joinpath("systems", "classical", "digital", "digital_systems.jl"))
include(joinpath("systems", "classical", "digital", "csk", "UTURUniChannelCSK.jl"))
end # module
