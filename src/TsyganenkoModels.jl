module TsyganenkoModels

using Dates: AbstractTime
using GeoCotrans
export t89, t96, t01, dipole_tilt
using LinearAlgebra: dot

include("dipole.jl")
include("t89.jl")
include("t96.jl")
include("t01.jl")

end
