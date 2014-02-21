abstract RandomProcessCarrier
abstract QuasiRandomCarrier
abstract IterativeMapCarrier <: QuasiRandomCarrier
abstract FSPWLCarrier <: IterativeMapCarrier # Fully-stretching piece-wise linear maps

typealias Carrier Union(Distribution, RandomProcessCarrier, QuasiRandomCarrier)
