@inline function _switch(f, sigma, s0, dsig, ps, x, y, z, oimf; q0)
    return if sigma < (s0 + dsig)
        b = f()
        if sigma < (s0 - dsig)
            b
        else # Interpolation region near magnetopause
            fint = 0.5 * (1 - (sigma - s0) / dsig)
            fext = 1.0 - fint
            q = dipole(ps, x, y, z; q0)
            @. (b + q) * fint + oimf * fext - q
        end
    else  # Outside magnetosphere
        oimf .- dipole(ps, x, y, z; q0)
    end
end
