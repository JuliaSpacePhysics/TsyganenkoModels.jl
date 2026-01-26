# Tsyganenko 89 Magnetospheric Magnetic Field Model (T89c/d)

# Fixed constants
const A02 = 25.0
const XLW2 = 170.0
const YN = 30.0
const RT = 30.0
const XD = 0.0
const XLD2 = 40.0
const SXC = 4.0
const XLWC2 = 50.0
const DXL = 20.0

# Parameter array: 30 parameters × 7 Kp levels
# Kp levels: 1=(0,0+), 2=(1-,1,1+), 3=(2-,2,2+), 4=(3-,3,3+), 5=(4-,4,4+), 6=(5-,5,5+), 7=(≥6-)
const PARAM = [
    # Kp=1
    -116.53 -10719.0 42.375 59.753 -11363.0 1.7844 30.268 -0.035372 -0.066832 0.016456 -1.3024 0.0016529 0.0020293 20.289 -0.025203 224.91 -9234.8 22.788 7.8813 1.8362 -0.27228 8.8184 2.8714 14.468 32.177 0.01 0.0 7.0459 4.0 20.0;
    # Kp=2
    -55.553 -13198.0 60.647 61.072 -16064.0 2.2534 34.407 -0.038887 -0.094571 0.027154 -1.3901 0.001346 0.0013238 23.005 -0.030565 55.047 -3875.7 20.178 7.9693 1.4575 0.89471 9.4039 3.5215 14.474 36.555 0.01 0.0 7.0787 4.0 20.0;
    # Kp=3
    -101.34 -13480.0 111.35 12.386 -24699.0 2.6459 38.948 -0.03408 -0.12404 0.029702 -1.4052 0.0012103 0.0016381 24.49 -0.037705 -298.32 4400.9 18.692 7.9064 1.3047 2.4541 9.7012 7.1624 14.288 33.822 0.01 0.0 6.7442 4.0 20.0;
    # Kp=4
    -181.69 -12320.0 173.79 -96.664 -39051.0 3.2633 44.968 -0.046377 -0.16686 0.048298 -1.5473 0.0010277 0.0031632 27.341 -0.050655 -514.1 12482.0 16.257 8.5834 1.0194 3.6148 8.6042 5.5057 13.778 32.373 0.01 0.0 7.3195 4.0 20.0;
    # Kp=5
    -436.54 -9001.0 323.66 -410.08 -50340.0 3.9932 58.524 -0.038519 -0.26822 0.074528 -1.4268 -0.0010985 0.0096613 27.557 -0.056522 -867.03 20652.0 14.101 8.3501 0.72996 3.8149 9.2908 6.4674 13.729 28.353 0.01 0.0 7.4237 4.0 20.0;
    # Kp=6
    -707.77 -4471.9 432.81 -435.51 -60400.0 4.6229 68.178 -0.088245 -0.21002 0.11846 -2.6711 0.0022305 0.01091 27.547 -0.05408 -424.23 1100.2 13.954 7.5337 0.89714 3.7813 8.2945 5.174 14.213 25.237 0.01 0.0 7.0037 4.0 20.0;
    # Kp=7
    -1190.4 2749.9 742.56 -1110.3 -77193.0 7.6727 102.05 -0.096015 -0.74507 0.11214 -1.3614 0.0015157 0.022283 23.164 -0.074146 -2219.1 48253.0 12.714 7.6777 0.57138 2.9633 9.3909 9.7263 11.123 21.558 0.01 0.0 4.4518 4.0 20.0
]'  # Transpose to get 30×7 matrix

"""
    t89(𝐫, ps, iopt; cache=nothing) -> GSM(Bx, By, Bz)
    t89(x, y, z, ps, iopt)
    t89(x, y, z, t::AbstractTime, iopt)

Compute GSM components of the magnetic field [nT] produced by extraterrestrial current systems
in the geomagnetosphere using the Tsyganenko 89 model, given the position `𝐫` or `x`, `y`, `z` in GSM coordinates [Earth radii], the geodipole tilt angle `ps` [radians] / time `t`.
"""
function t89(x, y, z, ps, iopt; cache = nothing)
    @assert 1 ≤ iopt ≤ 7 "iopt must be between 1 and 7"

    # Get parameters for this Kp level
    if isnothing(cache) || cache.iopt != iopt
        cache = t89_init(iopt)
    end
    return GSM(_t89_compute(x, y, z, ps, cache))
end

(m::T89)(x, y, z, ps) = t89(x, y, z, ps, m.iopt; cache = m.cache)

function t89_init(iopt::Int)
    # Extract parameters for this Kp level
    A = @view PARAM[:, iopt]

    # Parse parameters
    AK1, AK2, AK3, AK4, AK5 = A[1], A[2], A[3], A[4], A[5]
    AK6, AK7, AK8, AK9, AK10 = A[6], A[7], A[8], A[9], A[10]
    AK11, AK12, AK13, AK14, AK15 = A[11], A[12], A[13], A[14], A[15]
    AK16, AK17 = A[16], A[17]
    DX = A[18]
    ADR, D0, DD, RC = A[19], A[20], A[21], A[22]
    G, AT = A[23], A[24]
    P = A[25]
    DEL = A[26]
    Q, SX, GAM = A[27], A[28], A[29]
    DYC = A[30]
    DT = D0

    # Precompute constants
    DYC2 = DYC^2
    HA02 = 0.5 * A02
    RDX2M = -1.0 / DX^2
    RDX2 = -RDX2M
    RDYC2 = 1.0 / DYC2
    HLWC2M = -0.5 * XLWC2
    DRDYC2 = -2.0 * RDYC2
    DRDYC3 = 2.0 * RDYC2 * sqrt(RDYC2)
    HXLW2M = -0.5 * XLW2
    HXLD2M = -0.5 * XLD2
    W1 = -0.5 / DX
    W2 = W1 * 2.0
    W4 = -1.0 / 3.0
    W3 = W4 / DX
    W5 = -0.5
    W6 = -3.0
    RDXL = 1.0 / DXL
    HRDXL = 0.5 * RDXL
    A6H = AK6 * 0.5
    A9T = AK9 / 3.0
    YNP = 1 / π / YN * 0.5
    YND = 2.0 * YN

    AK610 = AK6 * W1 + AK10 * W5
    AK711 = AK7 * W2 - AK11
    AK812 = AK8 * W2 + AK12 * W6
    AK913 = AK9 * W3 + AK13 * W4

    return (;
        iopt,
        AK1, AK2, AK3, AK4, AK5, AK6, AK7, AK8, AK9, AK10,
        AK11, AK12, AK13, AK14, AK15, AK16, AK17,
        DX, ADR, D0, DD, RC, G, AT, DT, DEL, P, Q, SX, GAM, DYC,
        HA02, RDX2M, RDX2, RDYC2, HLWC2M, DRDYC2, DRDYC3, HXLW2M,
        HXLD2M, W1, W2, W3, W4, W5, W6, RDXL, HRDXL,
        A6H, A9T, YNP, YND, AK610, AK711, AK812, AK913, DYC2,
    )
end

function _t89_compute(x, y, z, ps, c)
    # Trigonometric functions of tilt angle
    SPS, CPS = sincos(ps)
    # Coordinate transformations
    x2 = x * x
    y2 = y * y
    z2 = z * z
    TPS = SPS / CPS
    HTP = TPS * 0.5
    XSM = x * CPS - z * SPS
    ZSM = x * SPS + z * CPS

    # Calculate tail current sheet shape function ZS
    XRC = XSM + c.RC
    XRC16 = XRC^2 + 16.0
    SXRC = sqrt(XRC16)
    y4 = y2 * y2
    Y410 = y4 + 1.0e4
    SY4 = SPS / Y410
    GSY4 = c.G * SY4
    ZS1 = HTP * (XRC - SXRC)
    DZSX = -ZS1 / SXRC
    ZS = ZS1 - GSY4 * y4
    D2ZSGY = -SY4 / Y410 * 4.0e4 * y2 * y
    DZSY = c.G * D2ZSGY

    # Ring current contribution
    XSM2 = XSM^2
    DSQT = sqrt(XSM2 + A02)
    FA0 = 0.5 * (1.0 + XSM / DSQT)
    DDR = c.D0 + c.DD * FA0
    DFA0 = c.HA02 / DSQT^3
    ZR = ZSM - ZS
    TR = sqrt(ZR^2 + DDR^2)
    RTR = 1.0 / TR
    RO2 = XSM2 + y2
    ADRT = c.ADR + TR
    ADRT2 = ADRT^2
    FK = 1.0 / (ADRT2 + RO2)
    DSFC = sqrt(FK)
    FC = FK^2 * DSFC
    FACXY = 3.0 * ADRT * FC * RTR
    XZR = XSM * ZR
    YZR = y * ZR
    DBXDP = FACXY * XZR
    DER25 = FACXY * YZR
    XZYZ = XSM * DZSX + y * DZSY
    FAQ = ZR * XZYZ - DDR * c.DD * DFA0 * XSM
    DBZDP = FC * (2.0 * ADRT2 - RO2) + FACXY * FAQ
    DER15 = DBXDP * CPS + DBZDP * SPS
    DER35 = DBZDP * CPS - DBXDP * SPS

    # Tail current sheet contribution
    DELY2 = c.DEL * y2
    D = c.DT + DELY2

    # Tail thickness modulation by dipole tilt
    ADSL = 0.0
    XGHS = 0.0
    if abs(c.GAM) ≥ 1.0e-6
        XXD = XSM - XD
        RQD = 1.0 / (XXD^2 + XLD2)
        RQDS = sqrt(RQD)
        H = 0.5 * (1.0 + XXD * RQDS)
        HS = -c.HXLD2M * RQD * RQDS
        GAMH = c.GAM * H
        D = D + GAMH
        XGHS = XSM * c.GAM * HS
        ADSL = -D * XGHS
    end

    D2 = D^2
    T = sqrt(ZR^2 + D2)
    XSMX = XSM - c.SX
    RDSQ2 = 1.0 / (XSMX^2 + XLW2)
    RDSQ = sqrt(RDSQ2)
    V = 0.5 * (1.0 - XSMX * RDSQ)
    DVX = c.HXLW2M * RDSQ * RDSQ2
    OM = sqrt(sqrt(XSM2 + 16.0) - XSM)
    OMS = -OM / (OM * OM + XSM) * 0.5
    RDY = 1.0 / (c.P + c.Q * OM)
    OMSV = OMS * V
    RDY2 = RDY^2
    FY = 1.0 / (1.0 + y2 * RDY2)
    W = V * FY
    YFY1 = 2.0 * FY * y2 * RDY2
    FYPR = YFY1 * RDY
    FYDY = FYPR * FY
    DWX = DVX * FY + FYDY * c.Q * OMSV
    YDWY = -V * YFY1 * FY
    DDY = 2.0 * c.DEL * y
    ATT = c.AT + T
    S1 = sqrt(ATT^2 + RO2)
    F5 = 1.0 / S1
    F7 = 1.0 / (S1 + ATT)
    F1 = F5 * F7
    F3 = F5^3
    F9 = ATT * F3
    FS = ZR * XZYZ - D * y * DDY + ADSL
    XDWX = XSM * DWX + YDWY
    RTT = 1.0 / T
    WT = W * RTT
    BRRZ1 = WT * F1
    BRRZ2 = WT * F3
    DBXC1 = BRRZ1 * XZR
    DBXC2 = BRRZ2 * XZR

    TLT2 = ps^2

    DER21 = BRRZ1 * YZR
    DER22 = BRRZ2 * YZR
    WTFS = WT * FS
    DBZC1 = W * F5 + XDWX * F7 + WTFS * F1
    DBZC2 = W * F9 + XDWX * F1 + WTFS * F3
    DER11 = DBXC1 * CPS + DBZC1 * SPS
    DER12 = DBXC2 * CPS + DBZC2 * SPS
    DER31 = DBZC1 * CPS - DBXC1 * SPS
    DER32 = DBZC2 * CPS - DBXC2 * SPS
    DER16 = (DER11, DER21, DER31) .* TLT2
    DER17 = (DER12, DER22, DER32) .* TLT2

    # Closure currents contribution
    ZPL = z + RT
    ZMN = z - RT
    ROGSM2 = x2 + y2
    SPL = sqrt(ZPL^2 + ROGSM2)
    SMN = sqrt(ZMN^2 + ROGSM2)
    XSXC = x - SXC
    RQC2 = 1.0 / (XSXC^2 + XLWC2)
    RQC = sqrt(RQC2)
    FYC = 1.0 / (1.0 + y2 * c.RDYC2)
    WC = 0.5 * (1.0 - XSXC * RQC) * FYC
    DWCX = c.HLWC2M * RQC2 * RQC * FYC
    DWCY = c.DRDYC2 * WC * FYC * y
    SZRP = 1.0 / (SPL + ZPL)
    SZRM = 1.0 / (SMN - ZMN)
    XYWC = x * DWCX + y * DWCY
    WCSP = WC / SPL
    WCSM = WC / SMN
    FXYP = WCSP * SZRP
    FXYM = WCSM * SZRM
    FPL = (x * FXYP, y * FXYP, WCSP + XYWC * SZRP)
    FMN = (-x * FXYM, -y * FXYM, WCSM + XYWC * SZRM)

    # Chapman-Ferraro sources
    EX = exp(x / c.DX)
    EC = EX * CPS
    ES = EX * SPS
    ECZ = EC * z
    ESZ = ES * z
    ESZY2 = ESZ * y2
    ESZZ2 = ESZ * z2
    ECZ2 = ECZ * z
    ESY = ES * y

    # C.-F. field components
    SX1 = c.AK6 * ECZ + c.AK7 * ES + c.AK8 * ESY * y + c.AK9 * ESZ * z
    SY1 = c.AK10 * ECZ * y + c.AK11 * ESY + c.AK12 * ESY * y2 + c.AK13 * ESY * z2
    SZ1 = c.AK14 * EC + c.AK15 * EC * y2 + c.AK610 * ECZ2 + c.AK711 * ESZ +
        c.AK812 * ESZY2 + c.AK913 * ESZZ2

    # Combine all contributions
    BCL = @. c.AK3 * (FPL + FMN) + c.AK4 * SPS * (FPL - FMN)
    BT = @. c.AK1 * (DER11, DER21, DER31) + c.AK2 * (DER12, DER22, DER32)
    return @. BT + BCL + c.AK5 * (DER15, DER25, DER35) + c.AK16 * DER16 + c.AK17 * DER17 + (SX1, SY1, SZ1)
end
