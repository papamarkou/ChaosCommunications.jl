abstract RandomProcessCarrier
abstract QuasiRandomCarrier
abstract IterativeMapCarrier <: QuasiRandomCarrier

typealias Carrier Union(Distribution, RandomProcessCarrier, QuasiRandomCarrier)
