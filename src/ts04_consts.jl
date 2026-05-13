# TS04 fitted coefficients.
# Reference: Tsyganenko & Sitnov (2005), JGR, doi:10.1029/2004JA010798.
# Each entry mirrors a term in the published formula; Fortran A[i] indices noted for audit.
# Saturation form: sat(w, β) = β·w / √(w² + β²), with per-driver scales w_sat[1..6].
const TS04_C = (
    amp_cf    = 1.0,                                                       # a[1]
    tamp1     = (c =  5.44118, dlp = 0.891995, w =  9.09684),              # a[2:4]; a[5]=0
    tamp2     = (c = -7.18972, dlp = 12.2700,  w = -4.89408),              # a[6:8]; a[9]=0
    src       = (c =  0.870536, w = 1.36081),                              # a[10:11]; a[12]=0
    prc       = (c =  0.688650, w = 0.602330),                             # a[13:14]; a[15]=0
    r11       = (c =  0.316346,  w =  1.22728),                            # a[16:17]
    r21       = (c = -0.0363620, w = -0.405821),                           # a[18:19]
    factimf   = 0.452536,                                                  # a[20]
    dlp1_exp  = 0.755831,                                                  # a[21]
    dlp2_exp  = 0.215662,                                                  # a[22]
    xappa_exp = 0.152759,                                                  # a[23]
    dxshift1  = (c =  5.96235, znam = 23.2036),                            # a[24:25]; c − znam/|dst*|^0.37
    dxshift2  = (c = 11.2994,  znam = 69.9596),                            # a[26:27]
    sc_sy     = (c = 0.989596, exp = -0.0132131),                          # a[28:29]
    sc_pr     = (c = 0.985681, exp =  0.0344212),                          # a[30:31]
    xkappa1   = (c = 1.02389, exp = 0.207867),                             # a[32:33]
    xkappa2   = (c = 1.51220, exp = 0.0682715),                            # a[34:35]
    d         = (a = 1.84714, b = 1.76977, off = 0.619441),                # a[36], a[37], a[69]
    phi       = 1.37690,                                                   # a[38]
    w_sat     = (0.696350, 0.343280, 3.28846, 111.293, 5.82287, 4.39664),  # a[39:44]
    rh0       = 7.5,
    g         = 35.0,
    deltady   = 4.7,
)
