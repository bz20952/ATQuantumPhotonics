import numpy as np
from vals import b, density, E, zeta


h = 0.34E-6  # Waveguide width
A_w = b*h  # Waveguide cross-sectional area

L = 40E-6
V = A_w*L
I = (b*h**3)/12


def flex_freq(beta: float):

    return ((beta**2)/(2*np.pi))*np.sqrt((E*I)/(density*A_w*L**4))

def damp(f: float):

    return f*np.sqrt(1-zeta**2)    

# Assume fixed-fixed
print('Fixed-fixed:')
f_flex = flex_freq(4.730)
print('Undamped:', f_flex/1E6, 'MHz')
f_flex_damped = damp(f_flex)
print('Damped:', f_flex_damped/1E6, 'MHz')

print('\n')

# Assume free-free
print('Free-free:')
f_flex = flex_freq(3.516)
print('Undamped:', f_flex/1E6, 'MHz')
f_flex_damped = damp(f_flex)
print('Damped:', f_flex_damped/1E6, 'MHz')

print('\n')

# Assume simply supported
print('Simply supported:')
f_flex = flex_freq(np.pi)
print('Undamped:', f_flex/1E6, 'MHz')
f_flex_damped = damp(f_flex)
print('Damped:', f_flex_damped/1E6, 'MHz')