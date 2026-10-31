# Development notes

Data Flow

User call: model(x, y, z, ps)
    └─ evalmodel(m::T01, ...) → t01(x,y,z,ps,params)
        └─ t01_extall(...)
            ├─ dipole_shield()          # Chapman-Ferraro (magnetopause currents)
            ├─ deformed()           # Tail current sheet
            ├─ birk_tot()           # Birkeland (Region 1 & 2)
            ├─ full_rc()            # Ring current ← ring_current.jl
            └─ linear combination   # weighted sum via model params

Key Design Patterns

- Functor interface — (m::T01)(x, y, z, ps) dispatches to evalmodel, then the raw function
- _switch — smooth magnetopause boundary: inside → full model, outside → just IMF penetration
- Component functions shared across T01/TS04: full_rc, birk_tot, deformed, dipole_shield all live in the parent module and are imported

## TS04 vs T01 Comparison

TS04 supersedes T01 with improved storm-time dynamics.
Many sub-functions are identical and shared: birk_tot, full_rc, dipole,
deformed, dipole_shield, warped, unwarped, taildisk, tail_shield.
Key Differences:


| Parameter        | T01                          | TS04                           |
|------------------|------------------------------|--------------------------------|
| Input indices    | g1, g2                       | w1-w6 (storm time integrals)   |
| Coefficients     | 43 elements                  | 69 elements                    |
| dsig (MP layer)  | 0.003                        | 0.005                          |
| rh0 (hinging)    | a[40] ≈ 9.0                  | Fixed 7.5                      |
| g (warping)      | a[41] ≈ 28.2                 | Fixed 35.0                     |
| IMF factor       | theta-dependent (a[24]+a[25]*θ) | Fixed a[20]                 |
| Tail d           | Fixed a[28]                  | a[36]*exp(-w1/a[37])+a[69]     |
| deltady          | a[29]                        | Fixed 4.7                      |
| dxshift          | g2-dependent                 | dst-dependent (|dst|^0.37)     |
| xkappa (Birk)    | g2-dependent                 | dst-dependent (znam/20)^a[k]   |
| phi (RC)         | 0.5π*tanh(|dst|/a[34])       | Fixed a[38]                    |
| Birk terms       | r11, r12, r21, r22           | Only r11, r21                  |

Amplitude Formulas:
T01:  tamp1 = a[2] + a[3]*dlp1 + a[4]*g1 + a[5]*dst
TS04: tamp1 = a[2] + a[3]*dlp1 + a[4]*a[39]*w1/√(w1²+a[39]²) + a[5]*dst
The w-index terms use saturation functions: a*w/√(w²+a²) → 1 as w→∞
This captures the nonlinear magnetospheric response to prolonged driving.
