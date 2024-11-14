import numpy as np
import matplotlib.pyplot as plt

# Set up the figure with a wider layout and style adjustments
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))
plt.subplots_adjust(wspace=0.4)

# Define x-axis values over a wider range for better visualization of potential shapes
x = np.linspace(-5, 5, 400)

# Define new Hamiltonian potentials with a well-like shape
H1 = 3 * np.tanh(x)**2 - 2  # A double-well type potential
H2 = 2 * np.tanh(x)**2 - 1   # A similar well with adjusted depth
H3 = np.tanh(x)**2           # A softer single-well shape

def V(a):
    return a**2 - (a*(a+1)) / ((np.cosh(x))**2)
# Define wavefunctions that correspond to each Hamiltonian
# Adjusting the forms to be consistent with typical bound-state solutions in potential wells
psi01 = (1 / np.cosh(x)) ** 2
psi11 = np.tanh(x) / (np.cosh(x)) ** 2
psi21 = (3 * np.tanh(x)**2 - 1) / ((np.cosh(x)) ** 2)

psi02 = 1 / (np.cosh(x)) ** 2
psi12 = np.tanh(x) / np.cosh(x)

psi03 = 1 / np.cosh(x)



# Customize plot styles for a cleaner, more professional look with labels and grid lines
ax1.plot(x, V(3), label=r'$V_1$', color='blue', linestyle='--')
ax1.plot(x, psi01, label=r'$\psi_0^{(1)}$', color='purple')
ax1.plot(x, psi11 + 5, label=r'$\psi_1^{(1)}$', color='green')
ax1.plot(x, psi21 + 8, label=r'$\psi_2^{(1)}$', color='red')
ax1.set_title(r"$H_1$ System")
ax1.set_xlabel("x")
ax1.set_ylabel("Amplitude")
ax1.legend(loc='upper right')
ax1.grid(True)

ax2.plot(x, V(2) + 5, label=r'$V_2$', color='blue', linestyle='--')
ax2.plot(x, psi02 + 5, label=r'$\psi_0^{(2)}$', color='purple')
ax2.plot(x, psi12 + 8, label=r'$\psi_1^{(2)}$', color='green')
ax2.set_title(r"$H_2$ System")
ax2.set_xlabel("x")
ax2.legend(loc='upper right')
ax2.grid(True)

ax3.plot(x, V(1) + 8, label=r'$V_3$', color='blue', linestyle='--')
ax3.plot(x, psi03 + 8, label=r'$\psi_0^{(3)}$', color='purple')
ax3.set_title(r"$H_3$ System")
ax3.set_xlabel("x")
ax3.legend(loc='upper right')
ax3.grid(True)

# Set the overall title
plt.suptitle("Partner Hamiltonian Modelling for Rosen-Morse Potential", fontsize=16)
plt.show()
