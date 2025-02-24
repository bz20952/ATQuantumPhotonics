import numpy as np

rho = 0.006  # Resistivity of n-doped silicon at dopant denisty of 10E16
L = 25E-6  # Actuation length
b = 0.16E-6  # Waveguide height
h = 0.36E-6  # Waveguide width
A_w = b*h  # Waveguide cross-sectional area
epsilon = 8.854E-12  # Permittivity of free space
d = 100E-9  # Unperturbed slot width
A_p = b*L  # Plate area normal to electrostatic force

C_air = epsilon * A_p / d
R_si = rho * L / A_w
tau = C_air * R_si

print('Resistance:', R_si/1E6, 'MOhm')
print('Capacitance:', C_air, 'F')
print('Time constant:', tau, 's')

f_c = 1/(2 * np.pi * tau)

print('Cutoff frequency:', f_c/1E6, 'MHz')

