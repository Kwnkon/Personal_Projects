import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.special import jn

#%% Task 1 Initial Value Problem

'''
d2J0/dz2 = - 1/z dJ0/dz - J0
         = rhs 


f1 = J0
f2 = dJ0/dz



f1h = f1 + dt/2 * f2

f2h = f2 + dt/2 * rhs

initial conditions state:
J0(0) = 1 => f1(0) = 1
dJ0/dz|z=0 = 0 => f2(0) = 0

'''

# Func to calculate the 2nd differential
def rhs(f1, f2, z):
    return -(1/z)*f2- f1

z = 0.001 # this is close to zero but not exactly as that would cause problems in the rhs func
dz = 0.01
f1 = 1 # f1(0) = 1 as in inital conditions
f2 = 0 # f2(0) = 0 as in initial conditions
    
z_vals = []
f1_vals = []

# Loop to find and plot estimates of the bessel function
while z < 30:
    
    z_vals.append(z)
    f1_vals.append(f1)
    
    f1h = f1 + dz/2 * f2
    f2h = f2 + dz/2 * rhs(f1,f2,z)
    
    f1 += dz * f2h
    f2 += dz * rhs(f1h,f2h, z + dz/2)
    
    z += dz
    #plt.plot(z,f1,'r.')
    #plt.pause(0.01)

# Inbuilt ODE solver func to compare with explicit RK2 method
def solver(z, y):
    f1, f2 = y # for some reason the solve_ivp function prefers these inputs being packed into 1 variable ¯\_(ツ)_/¯
    return [f2, - (1/z) * f2 - f1]

z_range = np.linspace(0.001, 30, 1000) # make it be as large a range as before and will also be used in a sec
sol = solve_ivp(solver, [0.001,30], [1,0], t_eval=z_range)

z_results = sol.t
f1_results = sol.y[0]

# Lastly the exact Bessel func generator to compare with above
bessel = jn(0, z_range)

# Plot results
plt.plot(z_vals, f1_vals, 'b--',label='Runga-Kutta method')
plt.plot(z_results, f1_results, 'r--', label='ODE Solver')
plt.plot(z_range, bessel, 'g-', label='Exact Bessel Function')
plt.xlabel('z')
plt.ylabel(r'$J_0(z)$')
plt.title('Approximation of Bessel Function')
plt.legend()

'''
The differnce betweeen each approximation method is
not very clear on first inspection but if you zoom
in you will be able to see the difference in
estimation methods for the Bessel function
'''