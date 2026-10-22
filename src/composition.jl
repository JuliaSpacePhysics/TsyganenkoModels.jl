struct TsyIGRF{T, N} <: CompositeFieldModel
    models::T
    name::N
end

getcsys(::TsyIGRF) = (GSM(), Cartesian3())

"""
    TsyIGRF(model = T89(iopt=3); name=nothing)

Create a composite Earth magnetic field model combining internal IGRF and external Tsyganenko field `model`.

# Example
```julia
model = TsyIGRF(T89(iopt=3))

t = Date(2020, 6, 21)  # Summer solstice 2020
B = model(6.6, 0.0, 0.0, t)  # GSM (x, y, z) [Re] -> (Bx, By, Bz) [nT]
```
"""
function TsyIGRF(model = T89(iopt = 3); name = nothing)
    return TsyIGRF((IGRF(), model), name)
end

function evaluate_model(m::TsyIGRF, pos, t, in, csys, out; kw...)
    return mapreduce(+, m.models) do model
        model_csys = getcsys(model)
        evaluate_model(model, pos, t, in, model_csys, out; kw...)
    end
end
