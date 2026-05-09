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
- ring_current.jl is the most complex component — the only one with both a symmetric solver (_rc_symm via finite-diff gradient of Aφ) and a shielding field (rc_shield with 86-element coefficient vectors)

## Findings by Category

### Design

birk_shl and rc_shield are the same function — identical 5-level nested loop (m:2, i:3, k:3, n:2, nn:2 = 72 terms), same coefficient indexing (a[72+i], a[75+k], a[85], a[86]), same l counter. The only difference is rc_shield premultiplies basis vectors by fac_sc = (x_sc + 1)^3. One function with a scale parameter eliminates ~50 lines of duplication.

Four implementations of the same harmonic shielding pattern: birk_shl, rc_shield, shlcar3x3, shlcar5x5 all implement "exponential decay × trig in y × trig in z, tilted, coefficient-weighted". They differ in grid size (3×3, 5×5, 2×3×3×2×2), basis type, and whether sps multiplies the second sum. Conceptually the same family — auditing one tells you nothing about the others.

_rc_symm, prc_quad, and one_cone all compute Jacobians via finite differences — this cross-cutting pattern appears 3 times in 3 different files with slightly different step sizes (1e-4, 1e-4, 1e-6). No shared helper.

Cardano cubic solver duplicated in ap_rc:68-75 and apprc:109-116 — identical 7 lines.

prox proximity correction pattern identical in both ap_rc and apprc.

### Performance

one_cone (called 4× per field point via birk_tot) fires 4 FD calls to r_s and theta_s (lines 92–95), costing 8 function evaluations for what is a 2×2 Jacobian. Replacing with ForwardDiff.jl dual numbers would cost ~2 function evaluations. Straightforward win since r_s/theta_s are pure scalar functions.

prc_quad calls br_prc_q and bt_prc_q 6× total (3 FD pairs). Each has 18/17 basis terms. ForwardDiff or analytic derivatives would halve this.

fialcos loop is always 1 or 2 iterations (n = mode ∈ {1, 2}) — the loop body is expensive (branches, trig, tangent calls). Unrolling as mode == 1 ? ... : ... would eliminate branch overhead and expose more to the compiler.

Missing @inbounds on all the inner accumulation loops (taildisk, shlcar5x5, rc_shield, birk_shl, shlcar3x3). Arrays are fixed-size constants so bounds checks are pure waste.

Missing @fastmath on pure numerical kernels like ap_rc, apprc, br_prc_q, bt_prc_q, r_s, theta_s. No special values, no IEEE requirements — would unlock FMA and reassociation.

sincos underuse: prc_quad:138 does sin(tp), cos(tp), sin(tm), cos(tm) separately (inconsistent with line 130 which uses sincos). one_cone:92-95 does the same. fialcos calls tan(tetanp*0.5) and tan(tetanm*0.5) without sincos.

Bumper.jl used only in t96_helpers.jl (3 spots: birk1tot_02 and r2inner). The T01/TS04 hot paths don't use it — they happen to avoid allocations by returning scalar tuples, which is good but accidental.

### Simplicity / Readability

l counter anti-pattern in all four shielding functions:
l = 0
for ...
    l += 1; bx += a[l] * ...
If loop order ever changes, all indices silently shift. Also blocks @inbounds. Direct flat-index arithmetic (e.g. a[(mes the layout auditable against the original Fortran and is @inbounds-safe.

shlcar3x3 has a silent special case at k=3 (lines 83–88 of dipole.jl) using a different basis formula with no comment ferent type of harmonic basis function — worth one line.

birk_1n2's phi convention uses atan(-zsc, xsc) (note the minus) — solar magnetic coordinates have z→-z relative to GSM

t96_helpers.jl is 1306 lines — a T96-specific monolith with its own Birkeland, ring current, and dipole shielding that components. These aren't wrong (T96 uses different parameterizations), but the file has no internal section structure.

apprc:107-108 — two lines of chained arithmetic longer than the screen. Each computes alpha_s and gamma_s via a sum ofeach term onto its own line would make the fitting structure visible.

### Summary Table

┌─────────────┬───────────────────────────────────────────┬──────────────────────────┐
│    Area     │                   Issue                   │          Impact          │
├─────────────┼───────────────────────────────────────────┼──────────────────────────┤
│ Design      │ birk_shl ≡ rc_shield (scale factor apart) │ -50 LoC                  │
├─────────────┼───────────────────────────────────────────┼──────────────────────────┤
│ Design      │ Cardano solver duplicated                 │ -7 LoC                   │
├─────────────┼───────────────────────────────────────────┼──────────────────────────┤
│ Performance │ one_cone 4 FD evals → ForwardDiff         │ ~2–4× for Birkeland      │
├─────────────┼───────────────────────────────────────────┼──────────────────────────┤
│ Performance │ prc_quad 6 FD evals → ForwardDiff         │ ~2× for PRC quadrupole   │
├─────────────┼───────────────────────────────────────────┼──────────────────────────┤
│ Performance │ Missing @inbounds on inner loops          │ ~5–15%                   │
├─────────────┼───────────────────────────────────────────┼──────────────────────────┤
│ Performance │ Missing @fastmath on kernels              │ ~10–20%                  │
├─────────────┼───────────────────────────────────────────┼──────────────────────────┤
│ Performance │ sincos half-used                          │ minor                    │
├─────────────┼───────────────────────────────────────────┼──────────────────────────┤
│ Simplicity  │ l counter pattern × 4                     │ fragile, blocks inbounds │
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
