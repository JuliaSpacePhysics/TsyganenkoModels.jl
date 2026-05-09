# Development notes

Repo Structure

TsyganenkoModels.jl/
├── src/
│   ├── TsyganenkoModels.jl     # Module root; exports T89–TS04, TsyIGRF, IGRF, dipole_tilt
│   ├── models.jl               # Abstract type TsyganenkoModel; struct defs T89/T96/T01/TS04
│   ├── composition.jl          # TsyIGRF = IGRF + Tsyganenko, evaluate_model dispatch
│   │
│   ├── t89.jl / t96.jl / t01.jl / ts04.jl    # Model entry points; functor (m::T01)(x,y,z,ps)
│   ├── t01_funcs.jl / ts04_funcs.jl           # extall() — assemble all field components
│   ├── t01_consts.jl / t96_consts.jl / ts04_consts.jl  # Coefficient arrays (T01_A, TS04_A, …)
│   │
│   ├── ring_current.jl         # full_rc, src_prc, ap_rc/apprc, rc_shield — C_SY, C_PR
│   ├── Birkeland_current.jl    # birk_tot, one_cone, birk_shl — SH coefficient matrices
│   ├── tail.jl                 # deformed, warped, taildisk — TAIL_A1/A2, TAIL_F/B/C
│   ├── dipole.jl               # dipole_tilt, dipole, shlcar3x3 — SHLCAR_A
│   └── helpers.jl              # _switch (magnetopause boundary logic)
│
├── test/
│   ├── runtests.jl             # Aqua, dipole_tilt, all 4 models, field-line tracing
│   └── test_geopack.jl         # Cross-validation vs Python geopack
│
├── lib/Geopack.jl/             # Julia wrapper for Python geopack (reference)
├── python/                     # PyPI wrapper package (tsyganenkomodels-jl)
└── docs/                       # Docusaurus site

---
Data Flow

User call: model(x, y, z, ps)
    └─ evalmodel(m::T01, ...) → t01(x,y,z,ps,params)
           └─ T01Impl.extall(...)
                  ├─ shlcar3x3()          # Chapman-Ferraro (magnetopause currents)
                  ├─ deformed()           # Tail current sheet
                  ├─ birk_tot()           # Birkeland (Region 1 & 2)
                  ├─ full_rc()            # Ring current ← ring_current.jl
                  └─ linear combination   # weighted sum via model params

TS04 has the same structure but adds 6 activity weights (w1–w6) and slightly different fixed parameters (rh0=7.5, g=35 vs T01's parameter-fitted values).

---
Key Design Patterns

- Each model is a submodule (module T01Impl) — isolates constants, avoids namespace pollution
- Functor interface — (m::T01)(x, y, z, ps) dispatches to evalmodel, then the raw function
- _switch — smooth magnetopause boundary: inside → full model, outside → just IMF penetration
- Component functions shared across T01/TS04: full_rc, birk_tot, deformed, shlcar3x3 all live in the parent module and are imported
- ring_current.jl is the most complex component — the only one with both a symmetric solver (_rc_symm via gradient of Aφ) and a shielding field (rc_shield with 86-element coefficient vectors)

## Findings by Category

### Design

Four implementations of the same harmonic shielding pattern: birk_shl, rc_shield, shlcar3x3, shlcar5x5 all implement "exponential decay × trig in y × trig in z, tilted, coefficient-weighted". They differ in grid size (3×3, 5×5, 2×3×3×2×2), basis type, and whether sps multiplies the second sum. Conceptually the same family — auditing one tells you nothing about the others.

_rc_symm computes a vector-potential gradient via finite differences. ForwardDiff through the elliptic/vector-potential path was slower in sampled benchmarks, so analytic derivatives would be needed for a useful speedup.

### Performance

fialcos loop is always 1 or 2 iterations (n = mode ∈ {1, 2}) — the loop body is expensive (branches, trig, tangent calls). Unrolling as mode == 1 ? ... : ... would eliminate branch overhead and expose more to the compiler.

sincos underuse: fialcos calls tan(tetanp*0.5) and tan(tetanm*0.5) without sincos.

Bumper.jl used only in t96_helpers.jl (3 spots: birk1tot_02 and r2inner). The T01/TS04 hot paths don't use it — they happen to avoid allocations by returning scalar tuples, which is good but accidental.

### Simplicity / Readability

shlcar3x3 has a silent special case at k=3 (lines 83–88 of dipole.jl) using a different basis formula with no comment ferent type of harmonic basis function — worth one line.

birk_1n2's phi convention uses atan(-zsc, xsc) (note the minus) — solar magnetic coordinates have z→-z relative to GSM

t96_helpers.jl is 1306 lines — a T96-specific monolith with its own Birkeland, ring current, and dipole shielding that components. These aren't wrong (T96 uses different parameterizations), but the file has no internal section structure.

apprc:107-108 — two lines of chained arithmetic longer than the screen. Each computes alpha_s and gamma_s via a sum ofeach term onto its own line would make the fitting structure visible.

### Summary Table

┌─────────────┬───────────────────────────────────────────┬──────────────────────────┐
│    Area     │                   Issue                   │          Impact          │
├─────────────┼───────────────────────────────────────────┼──────────────────────────┤
│ Performance │ sincos half-used                          │ minor                    │
├─────────────┼───────────────────────────────────────────┼──────────────────────────┤
│ Simplicity  │ apprc two-liner expressions               │ readability              │
└─────────────┴───────────────────────────────────────────┴──────────────────────────┘

## TS04 vs T01 Comparison

TS04 supersedes T01 with improved storm-time dynamics.
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
