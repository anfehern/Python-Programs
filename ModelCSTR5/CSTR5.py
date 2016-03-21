# -*- coding: utf-8 -*-
"""
Created on Sun Mar 13 21:45:10 2016

@author: Hexen

Info: Nonlinear model of two CSTRs in series.
      Henson, M.A. and Seborg, D.E., Feedback Linearizing Control
      Chap. 4 of Nonlinear Process Control
      Edited by Hensen, M.A. and Seborg, D.E., Prentice Hall (1997)
"""

# Importing Libraries
import numpy as np
import pandas as pd
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Defining Global Varables for Control Purposes
global u

# Model ODE System of two CSTRs in series
def CSTR5(y, t):

    # Input (1) - Coolant Flowrate (L/min)
    qc = u;
    
    # States (4):
    # Concentration of A in Reactor #1 (mol/L)
    Ca1 = y[0];
    # Temperature of Reactor #1 (K)
    T1 = y[1]
    # Concentration of A in Reactor #2 (mol/L)
    Ca2 = y[2]
    # Temperature of Reactor #2 (K)
    T2 = y[3]
    
    # Parameters
    q = 100                 # Flowrate (L/min) 
    Caf = 1.0               # Feed Concentration of A (mol/L)
    Tf = 350.0              # Feed Temperature (K)
    Tcf = 350.0             # Coolant Temperature (K)
    V1 = 100                # Volume of Reactor #1 (L)
    V2 = 100                # Volume of Reactor #2 (L)
    # UA1 or UA2 - Overall Heat Transfer Coefficient (J/min-K)    
    UA1 = 1.67e5            
    UA2 = 1.67e5
    # Pre-exponential Factor for A->B Arrhenius Equation
    k0 = 7.2e10             
    # EoverR - E/R (K) - Activation Energy (J/mol) / Gas Constant (J/mol-K)    
    EoverR = 1e4            
    # Heat of Reaction - Actually (-dH) for an exothermic reaction (J/mol)
    dH = 4.78e4
    rho = 1000              # Density of Fluid (g/L)
    rhoc = 1000             # Density of Coolant Fluid (g/L)
    Cp = 0.239              # Heat Capacity of Fluid (J/g-K)
    Cpc = 0.239             # Heat Capacity of Coolant Fluid (J/g-K)
    
    # Dynamic Balances    
    # dCa1/dt
    dCa1dt = q*(Caf-Ca1)/V1 - k0*Ca1*np.exp(-EoverR/T1)
    # dT1/dt
    dT1dt = q*(Tf-T1)/V1 + (dH*k0/(rho*Cp))*Ca1*np.exp(-EoverR/T1) + \
       (rhoc*Cpc/(rho*Cp*V1)) * qc * (1-np.exp(-UA1/(qc*rhoc*Cpc))) * (Tcf-T1)
    # dCa2/dt
    dCa2dt = q*(Ca1-Ca2)/V2 - k0*Ca2*np.exp(-EoverR/T2)
    # dT2/dt
    dT2dt  = q*(T1-T2)/V2 + (dH*k0/(rho*Cp))*Ca2*np.exp(-EoverR/T2) + \
       (rhoc*Cpc/(rho*Cp*V2))  * qc * (1-np.exp(-UA2/(qc*rhoc*Cpc))) * \
       (T1 - T2 + np.exp(-UA1/(qc*rhoc*Cpc))*(Tcf-T1))
    
    # Results
    return np.array([dCa1dt, dT1dt, dCa2dt, dT2dt])
    

# Steady State - Initial Condition for the Control
u_ss = 100.0

# Steady State - Initial Conditions for the States
CA1_ss = 0.088227891617
T1_ss = 441.219326816202
CA2_ss = 0.005292690885
T2_ss = 449.474581253729
x_ss = [CA1_ss, T1_ss, CA2_ss, T2_ss]

# Final Time (sec)
tf = 16

# Open Loop - Step Change
u = u_ss * 1.075

# Plotting System of ODEs
tspan=np.linspace(0,16,1000)
ty=odeint(CSTR5,x_ss,tspan)

f1 = plt.figure()
plt.plot(tspan,ty[:,0])
plt.plot(tspan,ty[:,2])
plt.legend(['$C_{A,1}$','$C_{A,2}$'])
plt.xlabel('z (m)')
plt.ylabel('Dependent variable')

f2 = plt.figure()
plt.plot(tspan,ty[:,1])
plt.plot(tspan,ty[:,3])
plt.legend(['$T_1$','$T_2$'])
plt.xlabel('z (m)')
plt.ylabel('Dependent variable')


