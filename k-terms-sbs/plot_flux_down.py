from turtle import color
import numpy as np
import matplotlib.pyplot as plt

height, fdown1, fdown2 = np.loadtxt(fname='_Fdown_', skiprows=1, unpack=True)

fig, ax = plt.subplots()

ax.plot(fdown1, height, color='r', label='line-by-line')
ax.plot(fdown2, height, color='b', label='fast calculations')
ax.set_xlabel(r'downward flux, $\mathregular{W/m^{2}}$')
ax.set_ylabel(r'height, km')
ax.set_ylim(bottom=0)
ax.grid(which='major', axis='both', color='gray', alpha=0.5)
ax.legend()

plt.title(r'K1 thermal downward fluxes considering SO$_2$ for 400 - 600 $\mathregular{cm^{-1}}$ spectral range (Voigt line shape)')
plt.tight_layout()
plt.show()