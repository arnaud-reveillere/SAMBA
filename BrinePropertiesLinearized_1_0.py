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

from numpy import arange, array, append, average
from scipy import polyfit
from BrineProperties_1_0 import *


def GetBrineDensityLinearFit(h, T0, gradT, P0, Xs0, gradXs=0, ReturnAll = False, discretization = 100):
    
    '''Compute (rho0, xi), the brine density linear fit rho(z) = rho0 + xi*(z-z0), in kg.m-3, for a given
    temperature T(z) = T0 + gradT+(z-z0) in Celcius, Pressure P(z0) = P0 in Pa and salinity Xs(z) = Xs0 + gradXs*(z-z0) as mass fraction.
    T0 = T(z0), P0 = P(z0), Xs0 = Xs(z0)'''
    
    # variables definitions
    dz = h/discretization
    z = arange(0., h*1.00001, dz)
    rho = array([]) # density

    # initialization at z = 0
    T = T0; P = P0; Xs = Xs0;
    
    # incrementation. We move towards the surface.
    for zz in z:
        rho = append(rho, GetBrineDensity(T, P, Xs))
        T  += -gradT*dz
        Xs += -gradXs*dz
        P  += -rho[-1]*9.81*dz

    # Linear fit
    [xi, rho0] = polyfit(z, rho, 1)

    if ReturnAll == True:
        return rho0, xi, rho, z
    else:
        return rho0, xi


def GetAverageBrineViscosity(h, T0, gradT, P0, Xs0, gradXs=0, ReturnAll = False, discretization = 100):
    
    '''returns the average brine viscosity on the height of interest for a given temperature T(z) = T0 + gradT+(z-z0) in Celcius,
    Pressure P(z0) = P0 in Pa and salinity  Xs(z) = Xs0 + gradXs*(z-z0) as mass fraction, Xs0 = Xs(0)'''    

    # variables definitions
    dz = h/discretization
    z = arange(0., h*1.00001, dz)
    rho = array([]) # density
    mu = array([]) # viscosity
    
    # initialization
    T = T0; P = P0; Xs = Xs0;
    
    # induction
    for zz in z:
        rho = append(rho, GetBrineDensity(T, P, Xs))
        mu = append( mu, GetBrineViscosity(T, P, Xs))
        T  += -gradT*dz
        Xs += -gradXs*dz
        P  += -rho[-1]*9.81*dz

    if ReturnAll == True:
        return average(mu), mu, z
    else:
        return average(mu)

def GetBrineCompressibility(T, P, Xs, dP = 1.e5):
    
    '''Compute the brine compressibility for a given temperature T in Celsius,
    Pressure P in Pa, salinity  Xs as mass fraction and a typical pressure increase dP in Pa'''
     
    rho0 = GetBrineDensity(T, P, Xs)
    rho1 = GetBrineDensity(T, P+dP, Xs)
    return (rho1 - rho0)/(rho0*dP)
    
    
def GetTransmissivity(h, k, rho, mu, g = 9.81):
    '''Compute the transmissivity (m2/s) for a given aquifer of height h (m),
    permeability k (m2), filled with brine of density rho (kg/m3) and viscosity mu (Pa.s)'''
    return k*rho*g*h/mu
    
    
def GetStorability(h, w, Cr, rho, Cb, g = 9.81):
    '''Compute the storability for a given aquifer of height h (m), porosity w (m2),
    pore compressibility Cr (Pa-1), filled with brine of density rho (kg/m3) and compressibility Cb (Pa-1)'''
    return rho*g*w*h*(Cr + Cb)
    
def GetTSratio(w, k, Cr, mu, Cb):
    '''Compute the Transmissivity/Storability ratio (m2/s) for a given aquifer of porosity w (m2), permeability k (m2),
    pore compressibility Cr (Pa-1), filled with brine of viscosity mu (Pa.s) and compressibility Cb (Pa-1)'''
    return k/(mu*w*(Cr + Cb))

    
def GetMolality(Xs):
    '''Compute the molality (g/kgH2O) for a given salt mass fraction Xs(no unit)'''
    # molality = (solute mass)/(solvent mass)
    # salt mass fraction Xs = (solute mass)/(solution mass)
    return 1000*Xs/(1-Xs)
    

def GetSaltMassFraction(Molality):
    '''Compute the salt mass fraction (no unit) for a given Molality (g/kgH2O)'''
    # molality = (solute mass)/(solvent mass)
    # salt mass fraction Xs = (solute mass)/(solution mass)
    return Molality/(1000+Molality)




def GetMassInjectionRates(Qv, T, P, Xs):
    '''Compute the H2O and NaCl mass injection rate (kg/s) for a volumetric injection flow rate Qv (m3/s)
    at a temperature T (Celsius), pressure P (Pa) and salt mass fraction Xs'''
    
    rho = GetBrineDensity(T, P, Xs)
    
    # injection mass flows
    Qm_total = Qv*rho
    Qm_NaCl = Qm_total * Xs
    Qm_H2O = Qm_total - Qm_NaCl
    
    return Qm_H2O, Qm_NaCl


    
