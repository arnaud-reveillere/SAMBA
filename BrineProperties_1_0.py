#Copyright (C) 2010  Karl Bandilla
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see http://www.gnu.org/licenses/.

import math as math

def GetBrineDensityViscosity(T, P, XS):
    '''compute brine density (kg/m3) and viscosity (Pa s) for a given temperature in C (T),
       pressure in Pa (P) and salinity as mass fraction (XS)'''
    return GetBrineDensity(T, P, XS), GetBrineViscosity(T, P, XS)
      
def GetBrineViscosity(T, P, XS):
    '''compute brine viscosity (Pa s) for a given temperature in C (T),
       pressure in Pa (P) and salinity as mass fraction (XS)'''
    #compute brine density based on the approach described in
    #S. L. Phillips, A. Igbene, J. A. Fair, H. Ozbek and M .Tavana (1981)
    #A Technical Databook for Geothermal Energy Utilization. Lawrence
    #Berkeley Laboratory Report LBL-12810
    
    #compute saturation pressure of pure water
    Psatw = SAT(T)
    #molecular weight of NaCl
    AMSalt = 58.448
    xncgl = 0.0
    #convert to molal concentration
    smol = XS / AMSalt / (1.0 - XS - xncgl) * 1000.0
    #compute ratio of brine viscosity to pure water viscosity using
    #equation (1) from Phillips et al. (1981)
    ratio = 1.0 + 0.0816 * smol + 0.0122 * smol**2 + 0.000128 * smol**3 \
            + 0.000629 * T * (1.0 - math.exp(-0.7 * smol))
    dp = P - Psatw
    if dp < 1000.0:
        dp = 1000.0
    #get pure water density for viscosity calculation
    dw0 = GetPureWaterDensity(T, Psatw + dp)
    return ratio * GetPureWaterViscosity(T, dw0)

def GetBrineDensity(T, P, XS):
    '''compute brine density (kg/m3) for a given temperature in C (T),
       pressure in Pa (P) and salinity as mass fraction (XS)'''
    # compute brine density based on the approach described in
    #A. Battistelli, C. Calore and K. Pruess (1997) The simulator
    #THOUGH2/EWASG for modelling geothermal reservoirs with brines  
    #and non-condensible gas. Geothermics 26(4), pp. 437-464 and using
    #J.L. Haas (1976) Physical properties of the coexisting phases and
    #Thermochemical properties of the H2) component in boiling NaCl
    #solutions. USGS Bulletin 1421-A, Washington, DC, 73 pp.
    
    #atomic weight of NaCl
    AMSalt = 58.448
    #atomic weight of H2)
    AMWater = 18.016
    #compute saturation pressure of pure water
    Psatw = SAT(T)
    #compute density of pure water
    dws = GetPureWaterDensity(T, Psatw)
    xncgl = 0.0
    #convert to molal concentration
    smol = XS / AMSalt / (1.0 - XS - xncgl) * 1000.0
    #convert to specific volume plus unit conversion kg/m3 to g/cm3
    v0 = 1000.0 / dws
    #equation (10) from Haas (1976) with vc = 3.1975e0
    #there is a slight discrepancy from the paper for the last term
    #according to the paper the last term should be
    #(v0 / (3.1975e0 - v0))**2
    sk = (-13.644e0 + 13.97e0 * v0) * (3.1975e0 / (3.1975e0 - v0))**2
    #equation (8) from Haas (1976)
    fi = -167.219e0 + 448.55e0 * v0 - 261.07e0 * v0**2 + sk * math.sqrt(smol)
    #equation (10) from Haas (1976) plus unit conversion
    dbs = (1000.0 + smol * AMSalt) / (1000.0 * v0 + smol * fi) * 1000.0
    #compute critical temperature
    tc = GetCriticalTemperature(XS)
    #equation (26) from Battistelli et al. (1997)
    tau = 1.0 - (T + 273.15) / (tc + 273.15)
    #convert to mole fraction
    xmol = smol / (1000.0 / AMWater + smol)
    #equation (25) from Battistelli et al. (1997)
    c = -1.6534e-10 / (tau**1.25 - 5.6 * xmol**1.5 + 0.005)
    #get saturated vapor pressure of brine
    Psat = SATB(T, XS)
    #equation (24) from Battistelli et al. (1997)
    return dbs / (1.0 + c * (P - Psat))

    
def GetPureWaterViscosity(T, DX):
    #dynamic viscosity for pure water (based on the Eight International
    #Conference on the Properties of Steam, Release on Dynamic viscosity
    #of Water Substance, as published in the ASME Steam Tables, Appendix 6

    #constants from appendix C equation (4)
    a0 = 0.0181583
    a1 = 0.0177624
    a2 = 0.0105287
    a3 = -0.0036744
    #constants from appendix C table a
    b00 = 0.501938
    b01 = 0.235622
    b02 = -0.274637
    b03 = 0.145831
    b04 = -0.0270448
    b10 = 0.162888
    b11 = 0.789393
    b12 = -0.743539
    b13 = 0.263129
    b14 = -0.0253093
    b20 = -0.130356
    b21 = 0.673665
    b22 = -0.959456
    b23 = 0.347247
    b24 = -0.0267758
    b30 = 0.907919
    b31 = 1.207552
    b32 = -0.687343
    b33 = 0.213486
    b34 = -0.0822904
    b40 = -0.551119
    b41 = 0.0670665
    b42 = -0.497089
    b43 = 0.100754
    b44 = 0.0602253
    b50 = 0.146543
    b51 = -0.0843370
    b52 = 0.195286
    b53 = -0.032932
    b54 = -0.0202595

    #convert from C to K
    tk = T + 273.15
    #make dimensionless T*/T
    tr = 647.27 / tk
    #T/T*
    tri = 1.0 / tr
    #T*/T - 1
    tr1 = tr - 1.0
    #pure water density
    dx1 = DX
    if dx1 < 0.0:
        dx1 = DW0
    #dimensionless density rho/rho*
    dr = dx1 / 317.763
    #rho/rho* - 1
    dr1 = dr - 1.0
    #summation in square brackets of equation (2)
    amu0 = a0 + a1 * tr + a2 * tr**2 + a3 * tr**3
    #equation (2)
    amu0 = math.sqrt(tri) / amu0 * 1.e-6
    #equation (1) summation for i = 0
    sun = b00 + b01 * dr1 + b02 * dr1**2 + b03 * dr1**3 + b04 * dr1**4
    #equation (1) summation for i = 1
    sun += (b10 + b11 * dr1 + b12 * dr1**2 + b13 * dr1**3 + b14 * dr1**4) * tr1
    #equation (1) summation for i = 2
    sun += (b20 + b21 * dr1 + b22 * dr1**2 + b23 * dr1**3 + b24 * dr1**4) * tr1**2
    #equation (1) summation for i = 3
    sun += (b30 + b31 * dr1 + b32 * dr1**2 + b33 * dr1**3 + b34 * dr1**4) * tr1**3
    #equation (1) summation for i = 4
    sun += (b40 + b41 * dr1 + b42 * dr1**2 + b43 * dr1**3 + b44 * dr1**4) * tr1**4
    #equation (1) summation for i = 5
    sun += (b50 + b51 * dr1 + b52 * dr1**2 + b53 * dr1**3 + b54 * dr1**4) * tr1**5
    #equation (1)
    return amu0 * math.exp(dr * sun)

def SAT(T):
    #saturation pressure for pure water (based on the 1967 IFC Formulation
    #For Industrial Use, A Formulation of the Thermodynamic Properties
    #of Ordinary Water, as published in the ASME Steam Tables, Appendix 1

    #constants from section 7.1.5 Saturation Line
    k1 = -7.691234564
    k2 = -2.608023696e1
    k3 = -1.681706546e2
    k4 = 6.423285504e1
    k5 = -1.189646225e2
    k6 = 4.167117320
    k7 = 2.097506760e1
    k8 = 1.e9
    k9 = 6

    #convert to dimensionless temperature from 2.1 a)
    th = (T + 273.15) / 647.3

    #compute dimensionless beta K from section 5
    #short form of 1-theta 
    x = 1.0 - th
    #compute sumation term
    sc = k5 * x + k4
    sc = sc * x + k3
    sc = sc * x + k2
    sc = sc * x + k1
    sc = sc * x
    betaK = math.exp(sc / (th * (1.0 + k6 * x + k7 * x**2)) \
                  - x / (k8 * x**2 + k9))
    #convert beta K to pressure
    return betaK * 2.212e7

def SATB(T, XS):
    #vapor pressure of brine (from J.L. Haas: Physical properties of the
    #coexisting phases and Thermochemical properties of the H2) component
    #in boiling NaCl solutions. USGS Bulletin 1421-A, Washington, DC, 73 pp.
    #1976) 
    AMSalt = 58.448
    xg = 0.0
    #convert to molal concentration
    smol = XS / AMSalt / (1.0 - XS- xg) * 1000.0
    #equation (4)
    a = 1.0 + 5.93582e-6 * smol - 5.19386e-5 * smol**2 \
        + 1.23156e-5 * smol**3
    #equation (5)
    b = 1.15420e-6 * smol + 1.41254e-7 * smol**2 - 1.92476e-8 * smol**3 \
        - 1.70717e-9 * smol**4 + 1.05390e-10 * smol**5
    #convert temperature from C to K
    tk = T + 273.15
    if tk <= 0.0:
        print 'broken in SATB'
    #equation (3) with m = (a+bTx)^-1 and c=0, and convert back to C
    t0 = math.exp(math.log(tk) / (a + b * tk)) - 273.15
    #use the equivalent temperature (t0) to find pressure using pure water
    return SAT(t0)

def GetCriticalTemperature(XS):
    #find critical temperature for brine. Numerical solution using
    #Newton iterations of equation (27) in A. Battistelli, C. Calore
    #and K. Pruess (1997) The simulator THOUGH2/EWASG for modelling
    #geothermal reservoirs with brines and non-condensible gas.
    #Geothermics 26(4), pp. 437-464
    
    #intial guess
    xnew = 374.1
    #not sure why
    xpc = XS * 100.0
    #Newton iterations
    for i in range(1, 31):
        #equation (27) * 100, not sure where 0.1852385e-6 comes from
        fx = -92.682482 + 0.1852385e-6 + 0.43077335 * xnew \
              - 6.2561155e-4 * xnew**2 + 3.6441625e-7 * xnew**3- xpc
        #derivative of fx
        fdx = 0.43077335 - 6.2561155e-4 * 2.0 * xnew \
             + 3.6441625e-7 * 3.0 * xnew**2
        xold = xnew
        #compute change
        dx = -1.0 * fx / fdx
        #update critical temperature
        xnew = xold + dx
        #convergence criterium
        if abs(dx) < 0.1:
            break
    tc = xnew
    if tc < 374.15:
        tc = 374.15
    return tc

def GetPureWaterDensity(T, Psat):
    #pure water density (based on the 1967 IFC Formulation For Industrial 
    #Use, A Formulation of the Thermodynamic Properties of Ordinary 
    #Water, as published in the ASME Steam Tables, Appendix 1
    
    #constants from section 7.1.1: sa are lower case a
    sa1 = 8.438375405e-1
    sa2 = 5.362162162e-4
    sa3 = 1.720000000e0
    sa4 = 7.342278489e-2
    sa5 = 4.975858870e-2
    sa6 = 6.537154300e-1
    sa7 = 1.150e-6
    sa8 = 1.51080e-5
    sa9 = 1.41880e-1
    sa10 = 7.002753165e0
    sa11 = 2.995284926e-4
    sa12 = 2.040E-1

    #constants from section 7.1.1: a are upper case a
    a11 = 7.982692717e0
    a12 = -2.616571843E-2
    a13 = 1.522411790E-3
    a14 = 2.284279054E-2
    a15 = 2.421647003E2
    a16 = 1.269716088E-10
    a17 = 2.074838328E-7
    a18 = 2.174020350E-8
    a19 = 1.105710498E-9
    a20 = 1.293441934E1
    a21 = 1.308119072E-5
    a22 = 6.047626338E-14

    #convert to dimensionless temperature from 2.1 a)
    tkr = (T + 273.15) / 647.3
    #convert to dimensionless temperature from 2.1 a)
    pnmr = Psat / 2.212e7
    #from section 9.1: Sub-region 1, Reduced volume
    #preliminary terms Y and Z
    y = 1.0 - sa1 * tkr**2 - sa2 / tkr**6
    zp = sa3 * y**2 - 2.0 *sa4 * tkr + 2.0 * sa5 * pnmr
    z = y + math.sqrt(zp)
    #first term of Chi 
    par1 = a11 * sa5 * z**(-5.0 / 17.0)
    #second term of Chi
    par2 = a12 + a13 * tkr + a14 * tkr**2 + a15 * (sa6 - tkr)**10 \
           + a16 * (sa7 + tkr**19)**-1
    #third term of Chi
    par3 = (a17 + 2.0 * a18 * pnmr + 3.0 * a19 * pnmr**2) \
           * (sa8 + tkr**11)**-1
    #fourth term of Chi
    par4 = a20 * tkr**18 * (sa9 + tkr**2) \
           * (-3.0 * (sa10 + pnmr)**(-4) + sa11)
    #fifth term of Chi
    par5 = 3.0 * a21 * (sa12 - tkr) * pnmr**2
    #sixth term of Chi
    par6 = 4.0 * a22 * tkr**(-20) * pnmr**3
    #add up terms to get Chi
    Chi = par1 + par2 - par3 - par4 + par5 + par6
    #convert from from dimensionless specific volume to regular
    #specific volume from 2.1 a) 
    v = Chi * 3.17e-3
    #convert from specific volume to density
    return 1.0 / v

def GetBatzleWangDensity(T, P, XS):
    #equation (27a)
    rhow = 1.0 + 1e-6 * (-80.0 * T - 3.3 * T**2 + 0.00175 * T**3 \
                         + 489.0 * P - 2.0 * T * P + 0.016 * T**2 * P \
                         - 1.3e-5 * T**3 * P - 0.333 * P**2 \
                         - 0.002 * T * P**2)
    #equation (27b)
    rhob = rhow + XS * (0.668 + 0.44 * XS + 1e-6 * (300.0 * P \
                        - 2400.0 *P * XS + T * (80.0 + 3.0 * T \
                        - 3300 * XS - 13.0 * P + 47.0 * P * XS)))
    return rhob
    
