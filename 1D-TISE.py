import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from matplotlib.widgets import Slider
'''
let's call hbar2 / 2 m a2 = 1

so V = inf      for x<= 0 
     = - V0   for 0<x<=a
     = 0        for x>a

-  a2 φ'' + V φ = E φ

V0 = 0

for x <= 0:
    φ = 0

for 0 < x <= a:
    - a2 φ'' - 10 φ = E φ
    φ΄'' = ((E - 10) / a2) φ

    so λ = (E - 10) / a2
    
for x > a:
    - a2 φ'' = Ε φ
    φ'' = - Ε/a2 φ


'''
# Set up values that will be used later
dx = 0.01
a = 1
V0 = 10
x_range = np.arange(0, 5, 0.1)

# The Potential function as described in the quesiton
def V(x):
    if x <= 0:
        return np.inf
    elif 0 < x <= a:
        return - V0
    else:
        return 0

# This is creating an array for the potential values to be stored that we need to plot it
potential = []

for i in x_range:
    potential.append(V(i))




# This is the function that does the main work of solving the schroedinger equation
def solve(E):
    
    x_vals = []
    phi_vals = []
    
    # Defining this here cause it will be used later and it is easier to just call a function
    def rhs(phi, x, E):
        return (V(x) - E) * phi
    
    # Boundary conditions for this problem
    x = 0.001 # this is small but not zero to avoid the infinite potential
    phi = 0 # φ(0) = 0 due to potential well wall
    dphi = 1 # φ'(0) = 1 it doesn't actually matter what value this takes as long as it is non-zero

    # This loop uses Runge-Kutta method to numerically solve the equation
    while x <= 5:
        # these two just create arrays to plot later
        x_vals.append(x)
        phi_vals.append(phi)
        
        # half step method for RK2
        phi_h = phi + dx/2 * dphi
        dphi_h = dphi + dx/2 * rhs(phi, x, E)
        
        phi += dx * dphi_h
        dphi += dx * rhs(phi,x, E)
        
        x += dx
    
    # We need φ(x) = 0 as x -> oo but we can't do that, obviously, so we will 
    # check whether the function is near zero at x = 3a
    # This has the problem that unbounded states at high energies that happen
    # to have phi ~ 0 at x = 3a will be considered 'bound' with this method
    # specifically this will occur for 2.05 < E < 2.42
    
    i = np.searchsorted(x_vals, 3) # This finds the index, i, for φ(3a)
    phi3a = phi_vals[i]
    
    # This checks whether the value of φ(3a) is near zero and set the colour
    # to green if yes or red if not
    if np.abs(phi3a) < 0.05:
        colour = 'g'
    else:
        colour = 'r'        
        
        
    return x_vals, phi_vals, colour
    

# Creates the subplots for the potential and wavefunction graphs. The comma a the end of the plots are necessary for the slider
fig, ax = plt.subplots()
V_plot, = ax.plot(x_range, potential, label = 'V(x)')
phi_plot, = ax.plot([], [], label = 'φ(x)')
ax.set_xlabel('x')
ax.set_ylabel('V(x) and φ(x)')


# This sets the axes limits to be reasonable. The choice of -10.5 is for clarity of the lower limit of the potential
ax.set_xlim(0,3)
ax.set_ylim(-10.5,10)
plt.legend()

# These text boxes inform the user what the colour of the wavefunctions mean (although it is probably already intuitive)
# and explain why they might see an unphysical bound state for high energy
plt.text(0.3, 9, 'Correct Guess in Green', color='g')
plt.text(1.4,9, 'Incorrect Guess in Red', color='r')
plt.text(0.7, 5, "Note: You may notice 'bound states' \nat higher energy values than expected. \nThey are present due to only checking\n the wavefunction until x = 3a", color='blue', fontsize=9)

# This creates the slider we want to use
ax_slider = plt.axes([0.2, 0.005, 0.6, 0.03])
slider = Slider(ax_slider, label = 'E', valmin = -10, valmax = 3, valinit = -4, valstep=0.01)

# This function allows us to update the state of the plots with the slider input change
def update(val):
    E = slider.val
    xs, phis, colour = solve(E) # this solves the TISE for the desired E value
    phi_plot.set_data(xs, phis) # this updates the data
    phi_plot.set_color(colour)  # this updates the colour
    fig.canvas.draw_idle()      # and this refreshes the plot so we can see the new graph
    
# This links the slider to the function to make it all work
slider.on_changed(update) 

# Sets up an initial plot of the wavefunction so the user sees an initial prediction and can adjust based on that
update(-4)

# This actually shows the plots
plt.show()
