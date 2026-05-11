import numpy as np

# User Input
# Units: mm, MPa, tonnes/mm^3

# model geometry
L = 100.0
w_f = 12.7

# Material properties
#FRP
E1_FRP = 222000.0
E2_FRP = 7010.0
nu12_FRP = 0.314
G12_FRP = 4661.0
G13_FRP = G12_FRP
G23_FRP = 2540.0
rho_FRP = 1.2e-9

"""E1_FRP = 200000.0
E2_FRP = 200000.0
nu12_FRP = 0.3
G12_FRP = E1_FRP/(2.0*(1.0 + nu12_FRP))
G13_FRP = G12_FRP
G23_FRP = G12_FRP
rho_FRP = 7.85e-9"""

# TPU
rho_m = 1.0e-9
C10_m = 30.972
D1_m = 0.0

# tape
rho_t = 1.2e-9
E_t = 1.0
nu_t = 0.3

# FRP layup
t_t = 0.13
t_FRP = 0.04
layup = [0, 90, 0]

# mesh controls
meshSize = w_f/5.0

# prestress
prestress = 6.0

# out of plane pulling displacement, selected to achieve same rotation as a simulation that is known to work
s_alpha = np.sqrt(2)*8.0/12.7
uz_pull = w_f/np.sqrt(2)*s_alpha

# job setup
cpus = 1