import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

height, cool1, cool2 = np.loadtxt(fname='_COOLING', skiprows=1, unpack=True)
fig, ax = plt.subplots()

ax.plot(1013.25*cool1, height, color='r', label='line-by-line')
ax.plot(1013.25*cool2, height, color='b', label='fast calculations')
ax.set_xlabel(r'cooling rate, $\mathregular{K/day}$', fontsize=16)
ax.set_ylabel(r'height, km', fontsize=16)
ax.set_ylim(bottom=0, top=100)
ax.grid(which='major', axis='both', color='gray', alpha=0.5)
ax.legend(loc='best')

plt.title(r'K1 cooling rates considering SO$_{2}$ for 400 - 600 $\mathregular{cm^{-1}}$ spectral range (Voigt line shape)', fontsize=16)

ax_inset = inset_axes(ax, width="100%", height="100%", bbox_to_anchor=(0.8, 0.15, 0.2, 0.2), bbox_transform=ax.transAxes)
#discrepancy = np.where((1013.25*abs(cool1) > 0.1) & (1013.25*abs(cool2) > 0.1), 1013.25*abs(cool1 - cool2), 0.0)
discrepancy= 1013.25*abs(cool1-cool2)
ax_inset.plot(discrepancy, height, color='black')
ax_inset.set_xlabel(r'discrepancy, $\bf{\mathregular{K/day}}$', fontweight='bold')
ax_inset.set_xlim(left=0)
ax_inset.set_ylim(bottom=0, top=100)

ticks = ax_inset.get_yticks()
ticks = ticks[ticks != 0]
ax_inset.set_yticks(ticks)

ax_inset.grid(which='major', axis='both', color='gray', alpha=0.5)


# Inset plot
# fig, ax = plt.subplots()

# discrepancy = np.where((1013.25*cool1 > 0.1) & (1013.25*cool2 > 0.1), 1013.25*abs(cool1 - cool2), 0.0)

# ax.plot(discrepancy, height, color='black', linewidth=3)
# ax.set_xlabel(r'discrepancy, $\bf{\mathregular{K/day}}$', fontweight='bold')
# ax.set_xlim(left=0)
# ax.set_ylim(bottom=0, top=120)
# ticks = ax.get_yticks()
# ticks = ticks[ticks != 0]
# ax.set_yticks(ticks)
# ax.grid(which='major', axis='both', color='gray', alpha=0.5)


plt.tight_layout()
plt.show()