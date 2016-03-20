# -*- coding: utf-8 -*-
"""
Created on Sun Dec 27 15:50:15 2015

@author: The Mystery Machine

This is my first attempt to learn to write Python code for applications in
Chemical Engineering. Code copied from IPython Notebooks. 
"""

# Importing libraries
import numpy as np
import pandas as pd
from scipy.optimize import fsolve
from scipy.optimize import curve_fit
from scipy.integrate import odeint

from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression

import matplotlib.pyplot as plt

# Define function to evaluate
def f(V):
    return R*T/(V-b)-alph*a/(V*(V+b))-P

# Define variables for evaluate this function
R=0.08206       # Gas constant (L-atm/gmol-K)
P_c=72.9
T_c=304.2
w_bar=0.225
T=300.
P=100.
a=0.42747*(R**2.0*T_c**2.0/P_c)
b=0.08664*(R*T_c/P_c)
m=0.48508+1.55171*w_bar-0.1561*w_bar**2
alph=(1+m*(1-np.sqrt(T/T_c)))**2

# Solving function
print('initial f(V0)= %f' %f(0.07))
print('a = %f' %a)
print('b = %f' %b)
print('m = %f' %m)
print('alpha = %f' %(alph))
print('V= %f' %(fsolve(f,0.07))) # Print out solution result





