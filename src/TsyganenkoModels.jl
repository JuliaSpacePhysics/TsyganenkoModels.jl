module TsyganenkoModels

using Dates: AbstractTime
using GeoCotrans
using LinearAlgebra: dot
export t89, t96, t01, ts04, dipole_tilt
export T89, T96, T01, TS04

include("models.jl")
include("dipole.jl")
include("Birkeland_current.jl")
include("ring_current.jl")
include("tail.jl")
include("t89.jl")
include("t96.jl")
include("t01.jl")
include("ts04.jl")

for f in (:t89, :t96, :t01, :ts04)
    @eval $f(x, y, z, t::AbstractTime, args...; kw...) = $f(x, y, z, dipole_tilt(t), args...; kw...)
    @eval @inline function $f(r, args...; kw...)
        @assert length(r) == 3
        return GSM($f(r[1], r[2], r[3], args...; kw...))
    end
end

for T in (:T89, :T96, :T01, :TS04)
    # support namedtuple input
    @eval $T(nt::NamedTuple) = $T(; nt...)
    # support functor interface
    @eval (m::$T)(x, y, z, t::AbstractTime) = m(x, y, z, dipole_tilt(t))
    @eval (m::$T)(r, args...; kw...) = begin
        @assert length(r) == 3
        m(r[1], r[2], r[3], args...; kw...)
    end
end

end
