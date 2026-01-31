module Geopack

using PythonCall
using PythonCall: pynew, pyimport
using StaticArrays: SVector
using Dates: AbstractTime, datetime2unix

const geopack = pynew()

function recalc(ut)
    return _float((@pyconst pyimport("geopack.geopack").recalc)(ut))
end

recalc(time::AbstractTime) = recalc(datetime2unix(time))

_tuple(x) = pyconvert(Tuple, x)
_float(x) = pyconvert(Float64, x)
_sv3(x) = SVector{3}(x)
_sv3(py::Py) = SVector{3}(_float(py[0]), _float(py[1]), _float(py[2]))

function dip(xgsm, ygsm, zgsm)
    py = (@pyconst pyimport("geopack.geopack").dip)(xgsm, ygsm, zgsm)
    return _tuple(py)
end

function igrf_gsm(xgsm, ygsm, zgsm)
    py = (@pyconst pyimport("geopack.geopack").igrf_gsm)(xgsm, ygsm, zgsm)
    return _sv3(py)
end


for f in (:t89, :t96, :t01, :t04)
    mod = "geopack.$f"
    @eval function $f(args...)
        py = (@pyconst pyimport($mod).$f)(args...)
        return _sv3(py)
    end
end

function load_igrf(ut)
    py = (@pyconst pyimport("geopack.geopack").load_igrf)(ut)
    pyType = PyArray{Float64, 1, true, true, Float64}
    return (pyType(py[0]; copy = false), pyType(py[1]; copy = false))
end

function __init__()
    PythonCall.pycopy!(geopack, pyimport("geopack"))

    certifi = pyimport("certifi")
    ENV["SSL_CERT_FILE"] = pyconvert(String, certifi.where())
    return
end

function reload()
    (@pyconst pyimport("importlib").reload)(geopack)
end

end