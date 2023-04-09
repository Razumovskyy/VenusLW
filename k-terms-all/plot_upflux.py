import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

height, fup1, fup2 = np.loadtxt(fname='CO2_Fup_3_', skiprows=1, unpack=True)

fig, ax = plt.subplots()

ax.plot(fup1, height, color='red', label='line-by-line')
ax.plot(fup2, height, color='blue', label='fast calculations')
ax.set_xlabel(r'upward flux, $\mathregular{W/m^{2}}$', fontsize=20)
ax.set_ylabel(r'height, km', fontsize=20)
ax.set_ylim(bottom=0, top=100)
ax.grid(which='major', axis='both', color='gray', alpha=0.5)
ax.legend(loc='best', fontsize=14)

#plt.title(r'Aggregate thermal upward fluxes in the 400 - 600 $\mathregular{cm^{-1}}$ spectral range (CO$_2$, Tonkov)', fontsize=13)

# add the inset
ax_inset = inset_axes(ax, width="100%", height="100%", bbox_to_anchor=(0.55, 0.4, 0.4, 0.4), bbox_transform=ax.transAxes)
discrepancy = 100 * (abs(fup1 - fup2)) / fup1
ax_inset.plot(discrepancy, height, color='black')
ax_inset.set_xlabel(r'discrepancy, $\%$', fontweight='bold')
ax_inset.set_xlim(left=0)
ax_inset.set_ylim(bottom=0, top=120)

ticks = ax_inset.get_yticks()
ticks = ticks[ticks != 0]
ax_inset.set_yticks(ticks)

ax_inset.grid(which='major', axis='both', color='gray', alpha=0.5)

plt.tight_layout()
plt.show()