module TsyganenkoModels

using Dates: AbstractTime
using GeoCotrans
export t89, t96, dipole_tilt
using LinearAlgebra: dot

include("t89.jl")
include("t96.jl")
include("dipole.jl")

end
