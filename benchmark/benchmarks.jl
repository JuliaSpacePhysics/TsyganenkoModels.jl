using BenchmarkTools
using TsyganenkoModels
using Dates

const SUITE = BenchmarkGroup()

# Representative inputs
const r_gsm = (-5.1, 0.3, 2.8)
const ps = -0.533585131
const param = (; pdyn = 2.0, dst = -87.0, byimf = 2.0, bzimf = -5.0)
const x, y, z = r_gsm
const time = DateTime(2001, 1, 1, 2, 3, 4)

# ── Model-level ────────────────────────────────────────────────────────────────
models = SUITE["models"] = BenchmarkGroup()
models["T89"] = @benchmarkable T89(2)($r_gsm, $ps)
models["T96"] = @benchmarkable T96($param)($r_gsm, $ps)
models["T01"] = @benchmarkable T01($param)($r_gsm, $ps)
models["TS04"] = @benchmarkable TS04($param)($r_gsm, $ps)

# ── Shared field components ────────────────────────────────────────────────────
components = SUITE["components"] = BenchmarkGroup()

# Ring current: symmetric + partial + shielding
components["full_rc"] = @benchmarkable TsyganenkoModels.full_rc($ps, $x, $y, $z, 1.0, 1.05, 0.95)

# Ring current sub-pieces
components["rc_symm"] = @benchmarkable TsyganenkoModels.rc_symm($x, $y, $z)
components["prc_quad"] = @benchmarkable TsyganenkoModels.prc_quad($x, $y, $z)

# Birkeland: 4 region×mode evaluations
components["birk_tot"] = @benchmarkable TsyganenkoModels.birk_tot($ps, $x, $y, $z, 1.0, 1.0)

# Tail: deformed current sheet (two modes)
const tail_state = (; dxshift1 = 4.0, dxshift2 = 0.0, d = 3.0, deltady = 2.0, g = 28.2)
components["deformed"] = @benchmarkable TsyganenkoModels.deformed($ps, $x, $y, $z, 9.0, $tail_state)
components["shlcar5x5"] = @benchmarkable TsyganenkoModels.shlcar5x5(TsyganenkoModels.TAIL_A1, $x, $y, $z, 4.0)

# Dipole shielding
components["shlcar3x3"] = @benchmarkable TsyganenkoModels.shlcar3x3($x, $y, $z, $ps)

# ── Numerical kernels ──────────────────────────────────────────────────────────
kernels = SUITE["kernels"] = BenchmarkGroup()

# Vector potential kernels (called from _rc_symm via finite differences, 4x each)
const r_test, sint_test, cost_test = 4.0, 0.6, 0.8
kernels["ap_rc"] = @benchmarkable TsyganenkoModels.ap_rc($r_test, $sint_test, $cost_test)
kernels["apprc"] = @benchmarkable TsyganenkoModels.apprc($r_test, $sint_test, $cost_test)

# PRC quadrupole kernels (called from prc_quad via finite differences, 6x each)
kernels["br_prc_q"] = @benchmarkable TsyganenkoModels.br_prc_q($r_test, $sint_test, $cost_test)
kernels["bt_prc_q"] = @benchmarkable TsyganenkoModels.bt_prc_q($r_test, $sint_test, $cost_test)

# Birkeland coordinate mapping (called from one_cone via finite differences, 8x each)
const a11 = TsyganenkoModels.A11
# Single Birkeland cone (FD-heavy)
kernels["one_cone"] = @benchmarkable TsyganenkoModels.one_cone($a11, $x, $y, $z, 1, 0.06)

# RC shielding (72-term harmonic sum)
kernels["rc_shield_sy"] = @benchmarkable TsyganenkoModels.rc_shield(
    TsyganenkoModels.C_SY, $ps, 0.05, $x, $y, $z
)
