# Progamm is running on Windows 7 32 bits with python-2.7.1 + numpy-1.6.0 + scipy-0.9.0 + matplotlib-1.0.1
# contact: Arnaud Reveillere a.reveillere@brgm.fr / arnaud.reveillere@gmail.com

# Import some matplotlib function (only required for the plots)
from pylab import figure, xlabel, ylabel, plot, legend, subplot, savefig, axis

# Add the repository where the SAMBA model is located (absolute or relative path)
from sys import path
path.append('../')  
from SAMBA import *

###########################################################
# Parameters describing the problem
###########################################################

# Bottom aquifer
h_b = 30.                   # m    height
w_b = 0.12                  #      porosity
k_b = 6.67e-13              # m2   permeability
Cr_b = 4.5e-10              # Pa-1 pore compressibility

# Top aquifer
h_t = 80.                   # m    height
w_t = 0.2                   #      porosity
k_t = 1.e-12                # m2   permeability
Cr_t = 4.5e-10              # Pa-1 pore compressibility
MixingCoef = 0.             #      part of the lifting brine mixing in the top aquifer

# T,P initial conditions
T_lb = 65.                  # Celcius   temperature at the bottom of the leak
gradT = 0.03                # Celcius/m geothermal gradient 
P_lb = 146.e5               # Pa        pressure at the bottom of the leak

# Leak properties
r_l =.1102                  # m    leak radius
h_l = 820.                  # m    total height
h_pc = 820.                 # m    porous column height; 0 < h_pc <= h_l
w_pc = 0.15                 #      porous column porosity
k_pc = 1.e-11               # m2   porous column permeability

# Injection parameters
d = 3025.                   # m    injection - leak distance
Q0 = 200./3600              # m3/s injection volumetric flow rate
t_inj = 5*3.1557e7          # s    injection duration
t_sim = 20*3.1557e7         # s    simulation duration

# Brine salinities (Salt mass fractions, no unit)
Xs_b = 0.07                 #      bottom aquifer
Xs_t = 0.001                #      top aquifer
Xs0_lifted = 0.001          #      bottom of the leak
gradXs_lifted = 0.          # m-1  salinity gradient


###########################################################
# Semi-analytical resolution
###########################################################

# Semi-analytical resolution
[t, Ql, zi, DP_lb, DP_c, DP_lt, 
InjectionForce, DensityDifference, PorousColumnResistance, AquifersResistance, TubingPressureLoss,
MonotonicQl, gamma
]=LeakageSolution(h_b, w_b, k_b, Cr_b,                    # Bottom aquifer properties
                  h_t, w_t, k_t, Cr_t, MixingCoef,        # Top aquifer properties
                  T_lb, gradT, P_lb,                      # T,P initial conditions
                  r_l, h_l, h_pc, w_pc, k_pc,             # leak properties
                  Xs_b, Xs_t, Xs0_lifted, gradXs_lifted,  # Brines salinities as salt mass fraction
                  d, Q0, t_inj, t_sim,                    # Leak to injection distance, injection flow rate and duration, simulation time
                  True, False, True, True,                # Isothermal leakage, tubing P losses not included, Compute driving P
                  "Example_1_SAMBA_User_Guide.txt")       # Export results in text file

###########################################################
# Plots
###########################################################

# Modify some units.
Ql = 3.1557e7*Ql             # m3 per year
t = t/3.1557e7               # years

# Plot the flow rate, the lifting - lifted brines interface position and the the over-pressures (in MPa)
figure(1)

subplot(311)
plot(t, Ql)
ylabel('$Q_L\ (m^3\ per\ year)$')

subplot(312)
plot(t, zi)
ylabel('$z_I\ (m)$')

subplot(313)
plot(t, 1.e-6*DP_lt, label = r'$P_{L+}$')
plot(t, 1.e-6*DP_lb, label = r'$P_{L-}$')
ylabel(r'$P-P(t=0)\ (MPa)$')
legend(loc = 'best')
xlabel('$\mathrm{Time}\  (years)$')
savefig('Example_1_SAMBA_User_Guide_FlowRate&Pressures.png', dpi = 200)


# Plot driving and countering pressures over time
figure(2)

plot(t, 1.e-6*InjectionForce, label = "$\mathrm{Injection\ pressurization}$")
plot(t, 1.e-6*DensityDifference, label = "$\mathrm{Brines\ density\ difference}$")
plot(t, 1.e-6*PorousColumnResistance, label = "$\mathrm{Porous\ column\ resistance}$")
plot(t, 1.e-6*AquifersResistance, label = "$\mathrm{Aquifers\ flow\ resistance}$")
ylabel('$\mathrm{Driving\ and\ countering\ pressures}\ (MPa)$')
legend(loc = 'best')
xlabel('$\mathrm{Time}\  (years)$')
savefig("Example_1_SAMBA_User_Guide_Drinving&CounteringPressures.png", format = 'png', dpi = 200)


# Plot the leakage flow rate decomposition in monotic components and the gamma paramter
figure(3)

subplot(211)
k=0
for Qlm in MonotonicQl:
    plot(t, 3.1557e7*GetDimensionalFlow(Qlm, Q0), label = '$Q_{L,%s}$' % k)
    k += 1
ylabel('$Q_{L,k}\ (m^3\ per\ year)$')
legend(loc = 'center right')

subplot(212)
k=0
for g in gamma:
    plot(t, g, label = '$\gamma_{%s}$' % k)
    k += 1
ylabel('$\gamma_k$')
legend(loc = 'center right')
xlabel('$\mathrm{Time}\ (year)$')
savefig("Example_1_SAMBA_User_Guide_ConvolutionIntegralsApproximation.png", format = 'png', dpi = 200)
