module ChaosCommunications

using NumericExtensions
using Distributions
using StatsBase
using NLopt

import Base.mean, Base.var, Base.rand
import Distributions.log2π, Distributions.@continuous_distr_support

export
  # Types
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
  snr2var,
  var2snr,
  mcml_decoder,
  decode,
  sim_sys,
  sim_ber,
  psim_ber

abstract InfoForm
abstract Classical <: InfoForm
abstract Quantum <: InfoForm
type Analog <: Classical end
type Digital <: Classical end

abstract InteractionForm
abstract TransmitterForm
abstract ReceiverForm
type UniTransmitter <: TransmitterForm end
type MultiTransmitter <: TransmitterForm end
type UniReceiver <: ReceiverForm end
type MultiReceiver <: ReceiverForm end
abstract OneWay{T<:TransmitterForm, R<:ReceiverForm} <: InteractionForm
abstract TwoWay <: InteractionForm

abstract ChannelForm
type UniChannel <: ChannelForm end
type MultiChannel <:ChannelForm end

abstract System{Info<:InfoForm, I<:InteractionForm, C<:ChannelForm}

typealias ClassicalSystem{I<:InteractionForm, C<:ChannelForm} System{Classical, I, C}
typealias QuantumSystem{I<:InteractionForm, C<:ChannelForm} System{Quantum, I, C}
typealias AnalogSystem{I<:InteractionForm, C<:ChannelForm} System{Analog, I, C}
typealias DigitalSystem{I<:InteractionForm, C<:ChannelForm} System{Digital, I, C}

typealias DigitalOneWayUniChannelSystem{T<:TransmitterForm, R<:ReceiverForm} DigitalSystem{OneWay{T, R}, UniChannel}
# UTUR := UniTransmitter, UniReceiver
typealias DigitalUTURUniChannelSystem DigitalSystem{OneWay{UniTransmitter, UniReceiver}, UniChannel}

typealias FunctionOrNothing Union(Function, Nothing)

# type MinMaxError <: Exception
#     msg::String
# end

include(joinpath("carriers", "maps.jl"))
include(joinpath("carriers", "carriers.jl"))
include(joinpath("carriers", "iterative_map_carriers.jl"))
include(joinpath("carriers", "PBCSCarrier.jl"))
include(joinpath("decoders", "decoders.jl"))
include(joinpath("decoders", "mcml_decoders.jl"))
include("noise.jl")
include(joinpath("systems", "classical", "digital", "digital_systems.jl"))
include(joinpath("systems", "classical", "digital", "csk", "UTURUniChannelCSK.jl"))
end # module
