abstract RandomProcessCarrier
abstract QuasiRandomCarrier
abstract IterativeMapCarrier <: QuasiRandomCarrier
abstract FSPWL <: IterativeMapCarrier # Fully-stretching piece-wise linear maps

typealias Carrier Union(Distribution, RandomProcessCarrier, QuasiRandomCarrier)
