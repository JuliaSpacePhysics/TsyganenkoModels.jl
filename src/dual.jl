# Minimal forward-mode dual number for the exact Jacobians the Fortran originals take by
# finite differences. Covers only the operations the differentiated kernels use.
struct Dual{T, N} <: Real
    v::T
    d::NTuple{N, T}
end

Dual{T, N}(x::Real) where {T, N} = Dual{T, N}(x, ntuple(_ -> zero(T), Val(N)))
# Values are at least Float64 (see `_duals`), so the kernels' Float64 literals never widen T.
Base.promote_rule(::Type{Dual{T, N}}, ::Type{<:Union{AbstractFloat, Integer}}) where {T, N} = Dual{T, N}

function _duals(xs::Vararg{Any, N}) where {N}
    vs = Base.front(promote(xs..., 1.0))
    return ntuple(i -> Dual(vs[i], ntuple(j -> j == i ? one(vs[i]) : zero(vs[i]), Val(N))), Val(N))
end

@inline _chain(x::Dual, f, df) = Dual(f, df .* x.d)
# Where `df` can be infinite, a zero partial stays zero, e.g. `sqrt(max(c, 0.0))` once
# clamped.
@inline _chain0(x::Dual, f, df) = Dual(f, map(p -> ifelse(iszero(p), p, df * p), x.d))

Base.:+(x::Dual, y::Dual) = Dual(x.v + y.v, x.d .+ y.d)
Base.:-(x::Dual, y::Dual) = Dual(x.v - y.v, x.d .- y.d)
Base.:-(x::Dual) = Dual(-x.v, map(-, x.d))
Base.:*(x::Dual, y::Dual) = Dual(x.v * y.v, y.v .* x.d .+ x.v .* y.d)
function Base.:/(x::Dual, y::Dual)
    q = x.v / y.v
    return Dual(q, (x.d .- q .* y.d) ./ y.v)
end
# Without direct Dual–constant methods, `@fastmath` promotes constants to zero-partial
# duals, wasting most of the arithmetic.
for (op, fop) in ((:+, :add_fast), (:-, :sub_fast), (:*, :mul_fast), (:/, :div_fast))
    @eval Base.FastMath.$fop(x::Dual, y::Dual) = $op(x, y)
    @eval Base.FastMath.$fop(x::Dual, c::AbstractFloat) = $op(x, c)
    @eval Base.FastMath.$fop(c::AbstractFloat, x::Dual) = $op(c, x)
end
Base.:+(x::Dual, c::AbstractFloat) = Dual(x.v + c, x.d)
Base.:+(c::AbstractFloat, x::Dual) = x + c
Base.:-(x::Dual, c::AbstractFloat) = Dual(x.v - c, x.d)
Base.:-(c::AbstractFloat, x::Dual) = Dual(c - x.v, map(-, x.d))
Base.:*(x::Dual, c::AbstractFloat) = Dual(x.v * c, x.d .* c)
Base.:*(c::AbstractFloat, x::Dual) = x * c
Base.:/(x::Dual, c::AbstractFloat) = Dual(x.v / c, x.d ./ c)
Base.:/(c::AbstractFloat, x::Dual) = (q = c / x.v; _chain(x, q, -q / x.v))

@inline function _pow(x::Dual, p::Real)
    y = x.v^p
    return _chain0(x, y, iszero(x.v) ? p * x.v^(p - 1) : p * y / x.v)
end
@inline Base.literal_pow(::typeof(^), x::Dual, ::Val{p}) where {p} =
    _chain(x, Base.literal_pow(^, x.v, Val(p)), p * Base.literal_pow(^, x.v, Val(p - 1)))
Base.:^(x::Dual, p::Integer) = _pow(x, p)
Base.:^(x::Dual, p::AbstractFloat) = _pow(x, p)
Base.FastMath.pow_fast(x::Dual, p::AbstractFloat) = _pow(x, p)
Base.sqrt(x::Dual) = (s = sqrt(x.v); _chain0(x, s, inv(2s)))
Base.exp(x::Dual) = (e = exp(x.v); _chain(x, e, e))
Base.log(x::Dual) = _chain(x, log(x.v), inv(x.v))
function Base.sincos(x::Dual)
    s, c = sincos(x.v)
    return _chain(x, s, c), _chain(x, c, -s)
end
Base.sin(x::Dual) = sincos(x)[1]
Base.cos(x::Dual) = sincos(x)[2]

Base.:<(x::Dual, y::Dual) = x.v < y.v
Base.isless(x::Dual, y::Dual) = isless(x.v, y.v)
