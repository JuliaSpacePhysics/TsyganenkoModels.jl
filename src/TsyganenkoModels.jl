module TsyganenkoModels

using Dates: AbstractTime
using GeoCotrans
using LinearAlgebra: dot
export t89, t96, t01, ts04, dipole_tilt

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

end
