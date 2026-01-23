# Development notes

## TS04 vs T01 Comparison

TS04 is an evolution of T01 with improved storm-time dynamics.
Many sub-functions are identical and shared: birk_tot, full_rc, dipole,
deformed, shlcar3x3, warped, unwarped, taildisk, shlcar5x5.
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
#