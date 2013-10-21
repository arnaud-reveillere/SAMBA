# SAMBA is a computer program that solves the problem of leakage between connected aquifers using a semi-analytic method
# It is running on Windows 7 32 bits with python-2.7.1 + numpy-1.6.0 + scipy-0.9.0
# It has been developed under BRGM research project CO2 storage risks management
# 
# Copyright (C) 2012 Arnaud REVEILLERE
# Contact: a.reveillere at brgm.fr / arnaud.reveillere @ centraliens.net
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see http://www.gnu.org/licenses/.


from numpy import array, arange, append, pi, log, sign, interp
from scipy.special import exp1
from BrinePropertiesLinearized_1_0 import *


# List of notations used in this document

# prefix     meaning
# D          Difference (e.g. DP = Delta P = Pressure difference)

# suffix     meaning
# D          Dimensionless
# b          property of the bottom aquifer, either uniform in the aquifer (e.g. salinity Xs_b) or evaluated at the middle of the layer (e.g. pressure P_b)
# t          property of the bottom aquifer, either uniform in the aquifer (e.g. salinity Xs_t) or evaluated at the middle of the layer (e.g. pressure P_t)
# l          property of the leak
# pc         property of the porous column
# wb         property of the wellbore
# lb         property evaluated at the bottom of the leak (e.g. T_lb)
# Xs         constant salinty
# TXs        constant salinity and temperature
# lifted     property of the lifted brine
# lifting    property of the lifting brine




def LeakageSolution(h_b, w_b, k_b, Cr_b,                        # Bottom aquifer properties
                    h_t, w_t, k_t, Cr_t, MixingCoef,            # Top aquifer properties 
                    T_lb, gradT, P_lb,                          # T,P initial conditions
                    r_l, h_l, h_pc, w_pc, k_pc,                 # leak properties. Set h_pc = 0, w_pc = 1 and any value to k_pc (e.g. k_pc = 1) if there is no porous column
                    Xs_b, Xs_t, Xs0_lifted, gradXs_lifted,      # Brines salinities as salt mass fraction
                    d, Q0, t_inj, t_sim,                        # Leak to injection distance, injection flow rate and duration, simulation time
                    IsThermalEq = True, IsTubingPLossIncluded = False, ComputePressures = False, ComputeDrivingAndCounteringForces = False, ExportInTextFile = False, discretization = 500):
    
    '''Compute the time (s), the leakage flow rate (m3/s) and the position of the lifting/lifted brines interface (m) based on the following parameters defining the problem:
    - bottom aquifer height h_b (m), porosity w_b, permeability k_b (m2) and pore compressibility Cr_b (Pa-1);
    - top aquifer    height h_t (m), porosity w_t, permeability k_t (m2) and pore compressibility Cr_t (Pa-1); Part of the lifting brin mixing in the top aquifer MixingCoef;
    - initial temperature T_lb (Celsius) and pressure P_lb (Pa) at the bottom of the leak; temperature gradient gradT (Celsius/m);
    - leak radius r_l (m), height h_l (m); porous column height h_pc (m), porosity w_pc and permeability k_pc (m2). Set h_pc=0 and w_pc = 1 if there is no pc.
    - bottom and top aquifer salinities Xs_b and Xs_t; and salinity at the bottom of the leak Xs0_lifted and salinity gradient of the lifted brine gradXs_lifted (m-1)
    - leak to injection distance d (m) injection flow rate Q0 (m3/s), injection duration t_inj (s) and simulation duration t_sim (s)
    if ComputePressures = True, it will also compute the pressure increase (Pa) at the bottom of the leak, at the top of the porous column and at the top of the leak
    if ComputeDrivingAndCounteringForces = True, it will compute the pressure terms (Pa) injection-induced over-pressurization, density difference, porous column and aquifers flow resistance and tubing pressure loss
    ExportInTextFile should be False or the recording text file name in string format (e.g. 'MySimulation.txt')'''

    # Compute dimensionless groupings
    #####################################################################
    
    [beta, kappa, phi, tau, delta, Dxi, Drho0, DP1bk_pc, DP1bk_wb, mu_pc_lifting, mu_pc_lifted, W0, W_lb, W_lt, LambdaTable] = GetDimensionlessGroupings(h_b, w_b, k_b, Cr_b, h_t, w_t, k_t, Cr_t, T_lb, gradT, P_lb, r_l, h_l, h_pc, w_pc, k_pc, Xs_b, Xs_t, Xs0_lifted, gradXs_lifted, d, Q0, IsThermalEq, IsTubingPLossIncluded)
    # Note that if there is no porous column, we have delta =  DP1bk_pc = mu_pc_lifting = mu_pc_lifted = 0. ; kappa = w_pc = 1. These parameters have no influence.
    
    # Semi analytic resolution
    #####################################################################
    
    # If pressure and/or driving&countering forces are demanded, local variables used during this computation are required
    if ComputePressures or ComputeDrivingAndCounteringForces: 
        [t, Ql, zi, DP1_pc, DP1_wb, TubingPLoss, mu_pc_av, MonotonicQl, ExtremaTimes, gamma] = SemiAnalyticalResolution(beta, kappa, phi, tau, delta, w_pc, Dxi, Drho0, DP1bk_pc, DP1bk_wb, mu_pc_lifting, mu_pc_lifted, MixingCoef, W0, W_lb, W_lt, t_inj, t_sim, LambdaTable, True, discretization)
        
        # Compute pressures if demanded
        if ComputePressures: 
            [DP_lbD, DP_cD, DP_ltD] = GetPressures(t, DP1_wb, TubingPLoss, MonotonicQl, ExtremaTimes, gamma, W0, W_lb, W_lt, t_inj)
        
        # Compute driving and countering forces if demanded
        if ComputeDrivingAndCounteringForces: 
            [InjectionForce, DensityDifference, PorousColumnResistance, AquifersResistance, TubingPressureLoss] = GetDrivingAndCounteringForces(t, DP1_pc, DP1_wb, TubingPLoss, mu_pc_av, MonotonicQl, ExtremaTimes, gamma, kappa, W0, W_lb, W_lt, t_inj)
    
    # Base case: local variables used during this computation are not interesting and not returned.
    else:
        [t, Ql, zi] = SemiAnalyticalResolution(beta, kappa, phi, tau, delta, w_pc, Dxi, Drho0, DP1bk_pc, DP1bk_wb, mu_pc_lifting, mu_pc_lifted, MixingCoef, W0, W_lb, W_lt, t_inj, t_sim, LambdaTable, False, discretization)
        
        
    # Compute dimensional results
    ########################################################################
    
    Ql = GetDimensionalFlow(Ql, Q0)           # leakage rate in m3/s
    zi = GetDimensionalZ(zi, h_l)             # lifted / lifting brines interface position in m    
    
    if ComputePressures:
        DP_lb = GetDimensionalPressure(DP_lbD, Q0, h_b, k_b, P_lb, Xs_b, T_lb, gradT)
        DP_c  = GetDimensionalPressure(DP_cD , Q0, h_b, k_b, P_lb, Xs_b, T_lb, gradT)
        DP_lt = GetDimensionalPressure(DP_ltD, Q0, h_b, k_b, P_lb, Xs_b, T_lb, gradT)

    if ComputeDrivingAndCounteringForces: 
        InjectionForce         = GetDimensionalPressure(InjectionForce        , Q0, h_b, k_b, P_lb, Xs_b, T_lb, gradT)
        DensityDifference      = GetDimensionalPressure(DensityDifference     , Q0, h_b, k_b, P_lb, Xs_b, T_lb, gradT)
        PorousColumnResistance = GetDimensionalPressure(PorousColumnResistance, Q0, h_b, k_b, P_lb, Xs_b, T_lb, gradT)
        AquifersResistance     = GetDimensionalPressure(AquifersResistance    , Q0, h_b, k_b, P_lb, Xs_b, T_lb, gradT)
        TubingPressureLoss     = GetDimensionalPressure(TubingPressureLoss    , Q0, h_b, k_b, P_lb, Xs_b, T_lb, gradT)
        
        
        
    # Export results in a text file
    ###########################################################################
    
    # if ExportInTextFile = False, the output text file is not demanded
    # otherwise, it is demanded and ExportInTextFile must be a the text file name, in string format.
    if ExportInTextFile != False:
        f = open(ExportInTextFile, 'w')
        f.write('Time (s)\tUpwards flow rate in the passive well (m3/s)\tLifting-lifted brines interface position in the passive well (m)')
        
        if ComputePressures:
            f.write('\tPressure under the leak P_lb(t) - P_lb(t=0) (Pa)\tPressure at the top of the porous column P_c(t) - P_c(t=0) (Pa)\tPressure at the top of the leak P_lt(t) - P_lt(t=0) (Pa)')
        
        if ComputeDrivingAndCounteringForces:
            f.write('\tInjection pressurization (Pa)\tLifting-lifted brines density difference (Pa)\tPorous column flow resistance (Pa)\tAquifers flow resistance (Pa)\tPressure Loss in tubings (Pa)')
        f.write('\n')
        
        for n in range(len(t)):
            f.write('%s\t%s\t%s' % (t[n], Ql[n], zi[n]))
            
            if ComputePressures:
                f.write('\t%s\t%s\t%s' % (DP_lb[n], DP_c[n], DP_lt[n]))
            
            if ComputeDrivingAndCounteringForces:
                f.write('\t%s\t%s\t%s\t%s\t%s' % (InjectionForce[n], DensityDifference[n], PorousColumnResistance[n], AquifersResistance[n], TubingPressureLoss[n]))
            f.write('\n')
        
        f.close()
        
    # Return results
    ###########################################################################
    
    # Return time, flow rate, interface position; the pressures; the Heaviside approximation parameters
    if ComputePressures and not ComputeDrivingAndCounteringForces:
        return t, Ql, zi, DP_lb, DP_c, DP_lt, MonotonicQl, gamma
    
    # Return time, flow rate, interface position; the driving and countering pressures; the Heaviside approximation parameters        
    elif ComputeDrivingAndCounteringForces and not ComputePressures:
        return t, Ql, zi, InjectionForce, DensityDifference, PorousColumnResistance, AquifersResistance, TubingPressureLoss, MonotonicQl, gamma
    
    # Return time, flow rate, interface position; the pressures; the driving and countering forces; the Heaviside approximation parameters   
    elif ComputeDrivingAndCounteringForces and ComputePressures:
        return t, Ql, zi, DP_lb, DP_c, DP_lt, InjectionForce, DensityDifference, PorousColumnResistance, AquifersResistance, TubingPressureLoss, MonotonicQl, gamma
    
    # Base case: return time, flow rate and interface position.
    else:
        return t, Ql, zi
    

    
    
def GetDimensionlessGroupings(h_b, w_b, k_b, Cr_b,                  # Bottom aquifer properties
                              h_t, w_t, k_t, Cr_t,                  # Top aquifer properties
                              T_lb, gradT, P_lb,                    # T,P initial conditions
                              r_l, h_l, h_pc, w_pc, k_pc,           # leak properties
                              Xs_b, Xs_t, Xs0_lifted, gradXs_lifted,# Brines salinities as salt mass fraction
                              d, Q0,                                # Leak to injection distance, injection flow rate
                              IsThermalEq = True, IsTubingPLossIncluded = False):
    '''Get all dimensionless groupings and functions defined in Reveillere (2013) based on the following parameters defining the problem:
    - bottom aquifer height h_b (m), porosity w_b, permeability k_b (m2) and pore compressibility Cr_b (Pa-1);
    - top aquifer    height h_t (m), porosity w_t, permeability k_t (m2) and pore compressibility Cr_t (Pa-1);
    - initial temperature T_lb (Celsius) and pressure P_lb (Pa) at the bottom of the leak; temperature gradient gradT (Celsius/m);
    - leak radius r_l (m), height h_l (m); porous column height h_pc (m), porosity w_pc and permeability k_pc (m2). Set h_pc=0 and w_pc = 1 if there is no pc.
    - bottom and top aquifer salinities Xs_b and Xs_t; and salinity at the bottom of the leak Xs0_lifted and salinity gradient of the lifted brine gradXs_lifted (m-1)
    - leak to injection distance d (m) injection flow rate Q0 (m3/s)'''

    # Bottom aquifer brine properties
    T_b = T_lb + gradT*h_b/2
    P_b = P_lb + GetBrineDensity(T_b, P_lb, Xs_b)*9.81*h_b/2
    rho_b = GetBrineDensity(T_b, P_b, Xs_b)
    mu_b = GetBrineViscosity(T_b , P_lb, Xs_b)
    Cb_b = GetBrineCompressibility(T_b, P_b, Xs_b)

    # Properties of the lifted brine in the leak (salinity defined by Xs0_lifted and gradXs_lifted) 
    [rho0_lifted, xi] = GetBrineDensityLinearFit(h_l, T_lb, gradT, P_lb, Xs0_lifted, gradXs_lifted) 
    
    # Properties of the lifting brine in the leak, of constant salinity Xs_b
    # isothermal case
    if IsThermalEq == True:
        [rho0_lifting, xi_Xs_or_TXs] = GetBrineDensityLinearFit(h_l, T_lb, gradT, P_lb, Xs_b)
    # adiabatic case
    else:
        [rho0_lifting, xi_Xs_or_TXs] = GetBrineDensityLinearFit(h_l, T_b , 0.   , P_lb, Xs_b)
    
    # wellbore height
    h_wb = h_l - h_pc
    
    # Static pressure differences between top/bottom of wb when filled with leak native brine
    DP0_pc = (rho0_lifted + xi*h_pc/2)*9.81*h_pc
    DP0_wb = (rho0_lifted + xi*(h_l+h_pc)/2)*9.81*h_wb
    
    # dimensionless pressure increase due to lifting / lifted brines density difference in the wb
    DP1bk_pc = 4*pi*h_b*k_b/(Q0*mu_b)   * ((rho0_lifting - rho0_lifted) + (xi_Xs_or_TXs - xi)*h_pc/2 )*9.81*h_pc
    DP1bk_wb = 4*pi*h_b*k_b/(Q0*mu_b)   * ((rho0_lifting - rho0_lifted) + (xi_Xs_or_TXs - xi)*(h_l+h_pc)/2 )*9.81*h_wb
    
    # viscosity and permeability in the porous column for the Darcy flow equation
    if h_pc > 0:
        mu_pc_lifting = GetAverageBrineViscosity(h_pc, T_lb, gradT, P_lb, Xs_b)/mu_b                        # flow of the bottom reservoir brine
        mu_pc_lifted  = GetAverageBrineViscosity(h_pc, T_lb, gradT, P_lb, Xs0_lifted, gradXs_lifted)/mu_b   # flow of the brine initially filling the leak
        kappa = r_l**2*k_pc/(4*h_b*k_b*h_pc)
    else:
        mu_pc_lifting = mu_pc_lifted = 0. ; kappa = 1.
        
    # Top aquifer brine properties (Temperature, Pressure, density, viscosity and compressibility in SI units)
    T_t = T_lb - gradT*(h_l + h_t/2)
    P_lt = P_lb - DP0_wb # top of the leak
    P_t = P_lt - GetBrineDensity(T_t, P_lt, Xs_t)*9.81*h_t/2
    rho_t = GetBrineDensity(T_t, P_t, Xs_t)
    mu_t = GetBrineViscosity(T_t, P_t, Xs_t)
    Cb_t = GetBrineCompressibility(T_t, P_t, Xs_t) 

    # dimensionless groupings
    beta = 4*pi*h_b*k_b*9.81*rho0_lifting*h_l/(Q0*mu_b)
    phi = xi_Xs_or_TXs*h_l/(2*rho0_lifting)
    Dxi = beta*(xi_Xs_or_TXs - xi)*h_l/rho0_lifting
    Drho0 = beta*(rho0_lifting-rho0_lifted)/rho0_lifting
    delta = h_pc/h_l

    # characteristic time
    tau = pi*r_l**2*h_l/Q0
    
    # Compute the lambda coefficient of the tubing pressure loss due to frictional forces for the bottom brine (salt mass fraction Xs_b) flowing in a very highly corroded tubing (rugosity 0.005 m) in adiabatic conditions (constant temperature T_b)
    # Build the a table of dimensionless flow -> dimensionless pressure table [[0, ..., 1], [Pressure loss(0 m3/s), ... , Pressure loss(Q0 m3/s)]]
    if IsTubingPLossIncluded == False:
        LambdaTable = [[0., 1.], [0.,0.]] # interpolation will return 0 to any flow rate
    
    else:
        # Dimensional table. Set the coefficient to 1 for Q = 0.
        LambdaTable = [[0]+[2**(-i) for i in range(10,-1,-1)], []]
        for Qv in LambdaTable[0]:
            LambdaTable[1].append(GetPressureLossCoefficient(r_l, 0.005, Q0*Qv, rho_b, mu_b) * k_b*h_b*rho_b*h_wb*Q0/(pi*mu_b*r_l**5))

        
    # bottom aquifer dimensionless well function
    def W0(t):
        '''compute the dimensionless pressure increase in the bottom aquifer created by the injection, at the distance d'''
        return exp1(d**2/(4*t*GetTSratio(w_b, k_b, Cr_b, mu_b, Cb_b)))
    
    def W_lb(t):
        '''compute the dimensionless well pressure decrease in the bottom aquifer created by the leakage, at the leak location'''
        return exp1(r_l**2/(4*t*GetTSratio(w_b, k_b, Cr_b, mu_b, Cb_b)))
    
    # top aquifer dimensionless well function
    def W_lt(t):
        '''compute the dimensionless well pressure increase in the top aquifer created by the leakage, at the leak location'''
        return rho_t*mu_t*h_b*k_b/(rho_b*mu_b*h_t*k_t) * exp1(r_l**2/(4*t*GetTSratio(w_t, k_t, Cr_t, mu_t, Cb_t)))
    
    # return dimensionless groupings and the three well functions
    return beta, kappa, phi, tau, delta, Dxi, Drho0, DP1bk_pc, DP1bk_wb, mu_pc_lifting, mu_pc_lifted, W0, W_lb, W_lt, LambdaTable


        
        

def SemiAnalyticalResolution(beta, kappa, phi, tau, delta, w_pc, Dxi, Drho0, DP1bk_pc, DP1bk_wb, mu_pc_lifting, mu_pc_lifted, MixingCoef, W0, W_lb, W_lt, t_inj, t_sim, LambdaTable, ReturnLocalVar = False, discretization = 500):
    '''Compute the time (s), and dimensionless leakage flow rate and position of the lifting/lifted brines interface as a function of:
    - the dimensionless groupings beta, kappa, phi, tau, delta, w_pc, Dxi, Drho0, DP1bk_pc, DP1bk_wb,  MixingCoef, (cf. Reveillere, 2013 for definitions)
    - the dimensionless well functions W0, W_lb, W_lt; average dimensionless viscosities of the brines flowing in the pc: mu_pc_lifting, mu_pc_lifted (cf. Reveillere, 2013 for definitions)
    - the injection duration t_inj (s) and the simulation duration t_sim (s); the table of the dimensionless \Lambda(QlD) parameter, the number of steps in the time discretization
    if ReturnLocalVar = True, it also returns the following dimensionless variables: pressure increase due to lifting-lifted brines density difference in the pc and in the wb; 
    tubings pressure loss, average viscosity in the porous column; leakage monotonic components, their starting time and their gamma parameters'''

    ###################################################################################################################
    # Initialization
    ###################################################################################################################
    
    # Negligible quantity
    negligible = 1.e-5
    
    # nb of time steps over which the presence of an extrema is tested. 
    # 10 was tested successively for avoiding finding too numerous "numeric" extrema (i.e. oscillations around the "real" solution) due to the explicit resolution. 
    Nb_dt_ExtremaTest = 4
    
    # Time discretization
    dt = t_sim/discretization
    t = arange(dt, t_sim*(1+negligible), dt)

    # Initialization of the main variables
    zi = array([0.])        # dimensionless Interface position
    Ql = array([0.])        # total leakage volumetric flow rate at z = z_lb
    
    # Initialization of internal variables
    DP1_wb = array([0.])               # dimensionless pressure change in the wellbore when lifting brine replaces lifted one due to their density difference
    DP1_pc = array([0.])               # dimensionless pressure change in the porous column when lifting brine replaces lifted one due to their density difference
    TubingPLoss = array([0.])          # dimensionless pressure loss due to the dynamic leakage in the tubing (not a wellbore anymore)
    mu_pc_av = array([mu_pc_lifted])   # Average viscosity in the porous column used in the Darcy equation. It depends on which brine is flowing; initially lifting brine is flowing
    V_mixed = 0                        # Volume that has leaked into the top aquifer, mixed in it and can not go back in the leak if the flow reverses
    
    # Initialization of the variables used for the Heaviside step function approximation (cf. Reveillere, 2013).
    # A new monotonic leakage components will be constructed at every new extremum; the sum of these monotonic leakage components equals the total leakage rate
    nb_extrema = 0                     # nb of local extrema of the leakage flow rate. Global extrema are not included. nb of leakage monotonic components = nb_extrema + 1
    gamma = [array([negligible])]      # list of gamma parameters of every component. This gamma corresponds to 1 minus the one introduced in Nordbotten et al. 2004 for the approximation of the convolution integral
    ExtremaTimes = [0]                 # list of the times of the extrema (i.e. beginning time steps of monotonic flow rates; the first one starts at 0)
    MonotonicQl = [array([0.])]        # list of dimensionless leakage monotonic components at z = z_lb

        
    ###################################################################################################################
    # Resolution
    ###################################################################################################################
    
    for n  in range(discretization-1):
      
        # 1st step. Compute the total leakage rate Ql at time step n+1. Leakage flow rates are at time n+1, all other elements being at time step n. 
        # Note that the slowly varying gamma parameter, which is part of the leakage function, is an exception to that rule and is taken at time n.
        # Consequently compute the tubing pressure drop, which is a function of Ql, at time n+1.
        #################################################################################################################
        
        # Compute the injection driving over-pressure
        if t[n] <= t_inj:
            InjectionTerm = W0(t[n])
        else:
            InjectionTerm = W0(t[n]) - W0(t[n] - t_inj)

        # Compute constant leakage components, if any: Qlk, k = 0..(nb_extrema-1)
        ConstantLeakageTerm = 0
        Ql = append(Ql, 0.)

        
        for k in range(nb_extrema):
            # append the constant monotonic leakage components
            MonotonicQl[k] = append(MonotonicQl[k], MonotonicQl[k][n])
            # Evaluate the constant leakage term. It is evaluated at time n+1 in order to have all leakage components evaluated the same way.
            ConstantLeakageTerm += MonotonicQl[k][n+1] * ( mu_pc_av[n]/kappa + W_lt(gamma[k][n]*(t[n+1]-ExtremaTimes[k])) + W_lb(gamma[k][n]*(t[n+1]-ExtremaTimes[k])))
            
            # add it to the total leakage rate
            Ql[n+1] += MonotonicQl[k][n+1]

        # Compute the leakage rate of the only fluctuating monotonic components (i.e. the last added, k=nb_extrema. Others are constant)
        MonotonicQl[nb_extrema] = append(   MonotonicQl[nb_extrema], (InjectionTerm - ConstantLeakageTerm - DP1_pc[n] - DP1_wb[n] - TubingPLoss[n])/ \
        ( mu_pc_av[n]/kappa + W_lb(gamma[nb_extrema][n]*(t[n+1]-ExtremaTimes[nb_extrema])) + W_lt(gamma[nb_extrema][n]*(t[n+1]-ExtremaTimes[nb_extrema])))   )
   
        Ql[n+1] += MonotonicQl[nb_extrema][n+1]
        
        # Test if there is a local extremum
        # The condition that the leakage has been monotonic over the last Nb_dt_ExtremaTest time steps avoids finding too many local extremum due to poor time discretization
        if n > Nb_dt_ExtremaTest + 1 and t[n] > ExtremaTimes[-1] + Nb_dt_ExtremaTest*dt and sign((Ql[n+1] - Ql[n])*(Ql[n] - Ql[n-Nb_dt_ExtremaTest])) == -1:
            IsLocExtremum = True
        else:
            IsLocExtremum = False
        
        # Introduce a new monotonic leakage function in case of new local extremum or of major flow regime change (e.g. stopping the injection)
        if IsLocExtremum or  t_inj < t[n] <= t_inj+dt:
            ExtremaTimes.append(t[n])                                 # The extrema occurred at t[n] (s). Record it.
            MonotonicQl[nb_extrema][n+1] = MonotonicQl[nb_extrema][n] # Correct the previously computed time-step n+1 of the last monotonic rate. Set it constant
            NewMonotonicQl = [0]*(n+1) + [Ql[n+1] - Ql[n]]            # Introduce a new monotonic component with a null flow rate record before t[n], and non-null at t[n+1]
            MonotonicQl.append(array(NewMonotonicQl))                 # add it to the list of monotonic components
            gamma.append(array([1-negligible]*(n+1)))                 # add its new gamma parameter, initially = negligible
            nb_extrema += 1
            
            
        
        
        # Compute the pressure loss in the tubing
        if Ql[n+1] >1:
            print("Leakage flow supperior than injected flow rate (not possible) -> refine time discretization or modify numerical scheme. Tubing pressure losses are explicit.")
            break
        TubingPLoss = append(TubingPLoss, sign(Ql[n+1])*interp(abs(Ql[n+1]), LambdaTable[0], LambdaTable[1])*Ql[n+1]**2)

            
            
        # 2nd step : compute lifting - lifted brine interface position zi at n+1 based on the leakage history record from 0 to n+1 (included).
        # Consequently compute functions of zi at n+1: brine average viscosity in the porous column and pressure increases due to brine density difference in the leak.
        #################################################################################################################    
        
        # In adiabatic case, phi<0 and some cases can result in errors without that condition
        if .25 + phi*(sum(Ql)*dt/tau-V_mixed)/w_pc <=0:
            h_pc = 820.
            zi = append(zi, 1.)
            mu_pc_av = append(mu_pc_av, mu_pc_lifting)
            DP1_wb = append(DP1_wb, DP1bk_wb)
            DP1_pc = append(DP1_pc, DP1bk_pc)
            
        else: 
            # compute zi supposing it hasn't reached the top of the porous column
            zi_tmp = (-.5 + (.25 + phi*(sum(Ql)*dt/tau-V_mixed)/w_pc)**.5)/phi

            # 0 < zi <= h_pc: interface in the porous column
            if 0 < zi_tmp <= delta:
                zi = append(zi, zi_tmp)
                # Update viscosity and pressures in the porous column
                mu_pc_av =  append(mu_pc_av, zi[n+1]/delta * mu_pc_lifting + (delta-zi[n+1])/delta * mu_pc_lifted)
                DP1_pc = append(DP1_pc, Drho0*zi[n+1] + Dxi*zi[n+1]*(delta-zi[n+1]/2) )
                DP1_wb = append(DP1_wb, 0.)
                
            # zi > h_cp: interface has reached the top of the porous column
            elif zi_tmp > delta:
                # update viscosity and pressure in the porous column
                mu_pc_av = append(mu_pc_av, mu_pc_lifting)
                DP1_pc = append(DP1_pc, DP1bk_pc )
                
                # compute zi supposing it hasn't reached the top of the leak
                c = delta*(w_pc-1)*(1+phi*delta)-(sum(Ql)*dt/tau-V_mixed)
                zi_tmp2 = (-(1-phi*delta) + ((1-phi*delta)**2 -4*phi*c)**.5)/(2*phi)

                # h_cp <= zi < h_l: interface in the wb
                if  zi_tmp2 < 1:
                    zi = append(zi, zi_tmp2)
                    # update pressure in the wellbore
                    DP1_wb = append(DP1_wb, Drho0*(zi[n+1] - delta) + Dxi*(zi[n+1] - delta)*(1-(zi[n+1] - delta)/2) )
                    
                # zi >= h_l: interface has reached the top of the leak
                else:
                    zi = append(zi, 1.)
                    # update pressure in the wellbore
                    DP1_wb = append(DP1_wb, DP1bk_wb )
                    
                    # calculate the volume of lifting brine that has leaked into the top aquifer and can go back in the leak if the flow reverses
                    # It is estimated as MixingCoef x the total leakage
                    if Ql[n+1]>0:
                        V_mixed += MixingCoef*Ql[n+1]*dt/tau
                        
            # zi <=0, the interface has reached back the bottom of the leak
            else: 
                zi = append(zi, 0.)
                mu_pc_av = append(mu_pc_av, mu_pc_lifted)
                DP1_pc = append(DP1_pc, 0.)
                DP1_wb = append(DP1_wb, 0.)
            
            
        # 3rd step. Update the gamma parameter at n+1 based on history record from 0 to n+1 (included) for each monotonic leakage
        ###############################################################################################################
        
        for k in range(nb_extrema+1):
            gamma[k] = append(   gamma[k], max(negligible, dt*sum(MonotonicQl[k])/((t[n+1]-ExtremaTimes[k])*MonotonicQl[k][n+1]))   )
    
    
    ###################################################################################################################
    # Return results
    ###################################################################################################################
    
    # Everything is dimensionless except the time (s)
    if nb_extrema > 20:
        print("Warning: the leakage rate as been split up as a sum of %s monotonic leakage rates, which is abnormly high. Increasing time discretization is recomanded" % nb_extrema)
    if ReturnLocalVar == True:
        return t, Ql, zi, DP1_pc, DP1_wb, TubingPLoss, mu_pc_av, MonotonicQl, ExtremaTimes, gamma
    else:
        return t, Ql, zi

        
        
        
        
def GetPressures(t, DP1_wb, TubingPLoss, MonotonicQl, ExtremaTimes, gamma, W0, W_lb, W_lt, t_inj):
    '''Compute the dimensionless pressure increase at the bottom of the leak, at the top of the porous column and at the top of the leak based on:
    - the array of time steps time t (s), the duration of the injection t_inj
    - the following outputs of the semi-analytic resolution: pressure increase in the wellbore DP1_wb; pressure loss in the tubing;
      leakage components MonotonicQl, their starting time ExtremaTimes and their gamma parameters
    - the well functions W0, Wlb and Wlt'''

    
    # Pressure increase at the top of the leak      DP_lt(t) = P_lt(t) - P_lt(t=0)
    # initialization 
    DP_lt = array([0.]*len(t))
    
    # Pressure increase under the leak              DP_lb(t) = P_lb(t) - P_lb(t=0)
    #initialization 
    DP_lb = array([0.]*len(t))
    
    # Add the over-pressurization due to the injection
    for n in range(len(t)-1):
        if t[n] <= t_inj:
            DP_lb[n] += W0(t[n])
        else:
            DP_lb[n] += W0(t[n]) - W0((t[n]-t_inj))
            
            
    # Add pressure increase/decrease due to the monotonic leakage components
    for k in range(len(ExtremaTimes)):
        for n in range(len(t)-1):
            # Only compute the influence of the leakage component k after its start at tk.
            if t[n] > ExtremaTimes[k]:
                DP_lt[n] += MonotonicQl[k][n]*W_lt(gamma[k][n]*(t[n]-ExtremaTimes[k]))
                DP_lb[n] += -MonotonicQl[k][n]*W_lb(gamma[k][n]*(t[n]-ExtremaTimes[k]))
    
    # Pressure increase in between the pc and the wb   DP_c(t) = P_c(t) - P_c(t=0)
    DP_c  = DP_lt + DP1_wb + TubingPLoss
    

        
    return DP_lb, DP_c, DP_lt
    
    
def GetDrivingAndCounteringForces(t, DP1_pc, DP1_wb, TubingPLoss, mu_pc_av, MonotonicQl, ExtremaTimes, gamma, kappa, W0, W_lb, W_lt, t_inj):
    '''Compute the pressure terms (injection-induced over-pressurization, density difference, porous column and aquifers flow resistance, tubing pressure loss) over time based on:
    - the list of time steps time t (s), the duration of the injection t_inj
    - the following outputs of the semi-analytic resolution: pressure increase in the porous column DP1_pc and the wellbore DP1_wb; average viscosity in the pc over time; 
      leakage components MonotonicQl, their starting time ExtremaTimes and their gamma parameters
    - the dimensionless parameter kappa and well functions W0, Wlb and Wlt'''

    
    ## First compute the non Ql-proportional forces.  Due to the explicit numerical resolution method, they are taken from time 0 to n-1
    # Overpressure created by the injection.
    InjectionForce = W0(t[:-1])
    for n in range(len(t)-1):
        if t[n] > t_inj:
            InjectionForce[n] +=  - W0((t[n]-t_inj))
    
    # Brines density difference mechanism
    DensityDifference = DP1_wb[:-1] + DP1_pc[:-1]
    
    # Tubing Pressure loss
    TubingPressureLoss = TubingPLoss[:-1]

    # Set the value of these forces at n equal to the value at n-1
    InjectionForce = append(InjectionForce, InjectionForce[-1])
    DensityDifference = append(DensityDifference, DensityDifference[-1])
    TubingPressureLoss = append(TubingPressureLoss, TubingPressureLoss[-1])

    ## Then compute the Ql-proportional countering forces.  Due to the explicit numerical resolution method, they are taken from time 1 to n
    
    # aquifers and porous column resistance
    # initialisation
    AquifersResistance = array([0.]*len(t))
    PorousColumnResistance = array([0.]*len(t))
    
    # successively add every leakage monotonic component
    for k in range(len(ExtremaTimes)):
        for n in range(1, len(t)):
            # Only compute the influence of the leakage component k after its start at tk.
            if t[n] > ExtremaTimes[k]:
                AquifersResistance[n-1] += MonotonicQl[k][n]*(W_lt(gamma[k][n]*(t[n]-ExtremaTimes[k]))+W_lb(gamma[k][n]*(t[n]-ExtremaTimes[k])))
                PorousColumnResistance[n-1] += MonotonicQl[k][n]*mu_pc_av[n]/kappa
    
    # Set the value of these forces at time 0 to the one at time 1
    AquifersResistance[-1] = AquifersResistance[-2]
    PorousColumnResistance[-1] = PorousColumnResistance[-2]
    

    
    return InjectionForce, -DensityDifference, -PorousColumnResistance, -AquifersResistance, -TubingPressureLoss




def GetDimensionalPressure(PD, Q0, h_b, k_b, P_lb, Xs_b, T_lb, gradT):
    '''Compute the pressure (Pa) base on the dimensionless pressure PD, injection flow rate Q0 (m3/s) and:
    - bottom aquifer height h_b (m), permeability k_b (m2) and salt mass fraction Xs_b
    - the temperature T_lb (Celsius) and pressure P_lb (Pa) at the bottom of the leak; temperature gradient gradT (Celsius/m);'''
    T_b = T_lb + gradT*h_b/2
    mu_b = GetBrineViscosity(T_b , P_lb, Xs_b)
    
    return Q0*mu_b*PD/(4*pi*h_b*k_b)


def GetDimensionalFlow(QD, Q0):
    '''Compute the flow rate(m3/s) based on the dimensionless flow Q and on the injection flow Q0 (m3/s).'''
    return QD*Q0

    
def GetDimensionalZ(zD, h_l, z_lb=0):
    '''Compute the vertical unit z (m) based on the dimensionless one zD, on the bottom of the leak z_lb (m) and its height h_l (m).'''
    return zD*h_l+z_lb
    
    
def Nordbotten_et_al_2004(kappa, ViscosityD, W0, W_lb, W_lt, gamma = 0.92):
    '''Compute the dimensionless leakage function QlD(t), with t in s, according to Nordbotten et al., 2004 solution based on the dimensionless 
        parameter kappa (Darcy equation in the porous column), brine viscosity and well functions W0, W_lb, W_lt (cf. definitions in Reveillere, 2013).'''
    def f(t):
        return W0(t) / (ViscosityD/kappa + W_lb((1-gamma)*t) + W_lt((1-gamma)*t))
    return f
 
    
    
def Modified_Nordbotten_et_al_2004(kappa, ViscosityD, W0, W_lb, W_lt, tinj, discretization = 500):
    '''Compute the dimensionless leakage over time (s) according to Nordbotten et al., 2004 solution modified by introducing a varying gamma parameter.
    It is based on the dimensionless parameter kappa (Darcy equation in the porous column), brine viscosity and  well functions W0, W_lb, W_lt (cf. definitions in Reveillere, 2013).'''
    
    # time discretization
    dt = tinj/discretization
    t = arange(dt, tinj*1.00000001, dt)

    Ql = array([0.]) # leakage rate
    g = array([0.92]) # gamma in Nordbotten et al. 2004 for the brine leakage

    for n  in range(0,len(t)-1):
      
        #1st step
        Ql = append(Ql, W0(t[n])/(ViscosityD/kappa + W_lb((1-g[n])*t[n]) + W_lt((1-g[n])*t[n])) )

        #2nd step    
        gg = dt*sum(Ql)/(t[n+1]*Ql[n+1])
        if gg >= 1. or gg <0.1:
            g = append(g, g[-1])
        else:
            g = append(g, gg)

    return t, Ql

def GetPressureLossInTubing(length, radius, rugosity, FlowRate, density, viscosity):
    '''Compute the pressure loss of brine flowing in a tubing at a flow rate FlowRate (m3/s).
    Based on tubing length (m), diameter (m) and rugosity (m); and brine density (km.m-3) and viscosity (Pa.s)
    Based on the Darcy Weisbach and the Colebrook equations: en.wikipedia.org/wiki/Darcy-Weissbach_equation'''
    
    if FlowRate == 0:
        return 0
        
    else:
        
        # Compute the pressure loss coefficient
        l =  GetPressureLossCoefficient(radius, rugosity, FlowRate, density, viscosity)
        return l*density/4 * length/radius* (FlowRate/(pi*radius**2))**2
        
        
def GetPressureLossCoefficient(radius, rugosity, FlowRate, density, viscosity, negligible = 1.e-8, nmax = 1000):
    '''Compute the pressure loss coefficient for brine flowing in a tubing at a flow rate FlowRate (m3/s).
    Based on tubing length (m), diameter (m) and rugosity (m); and brine density (km.m-3) and viscosity (Pa.s)
    Based on the Darcy Weisbach and the Colebrook equations: en.wikipedia.org/wiki/Darcy-Weissbach_equation'''
    
    if FlowRate == 0:
        return 0
        
    else:
        
        # Reynolds number
        Re = 2*FlowRate*density/(pi*radius*viscosity)
        
        if Re < 30:
            l = 1.
            
        # Case 1: laminar flow.
        elif Re < 46*radius/rugosity:
            l = 64/Re

        # Case 2: any significant flow. Use the Colebrook-White expression, and apply a basic algorithm
        else:
            
            # initialization
            x1 = 1.
            x2 = 2.
            n = 0
            
            # resolution
            while abs(1-x1/x2) > negligible and n < nmax:
                x1 = x2
                x2 = -.87*log(2.51/Re * x1 + rugosity/(7.42*radius))
                n+=1

            # Colebrook-White parameter
            l = x2**-2
            
        return l

"""
References

Nordbotten, J.M., Celia, M.A., Bachu, S., 2004: Analytical solutions for leakage rates through abandoned wells.
Water Resour. Res. 40, W04204

Reveillere, A., 2013: Semi-analytical Solution for Brine Leakage Through Passive AbandonedWells Taking Account of Brine Density Differences.
Transp Porous Med. DOI 10.1007/s11242-013-0221-3
"""
