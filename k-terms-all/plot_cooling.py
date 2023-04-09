import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

height, cool1, cool2 = np.loadtxt(fname='CO2_COOLING_3_', skiprows=1, unpack=True)

fig, ax = plt.subplots()

ax.plot(cool1, height, color='r', label='line-by-line')
ax.plot(cool2, height, color='b', label='fast calculations')
ax.set_xlabel(r'cooling rate, $\mathregular{K/day}$', fontsize=20)
#ax.set_ylabel(r'height, km', fontsize=20)
ax.set_ylim(bottom=0, top=100)
ax.set_yticklabels([])
ax.grid(which='major', axis='both', color='gray', alpha=0.5)
ax.legend(loc='best', fontsize=14)

#plt.title(r'Aggregate cooling rates in the 400 - 600 $\mathregular{cm^{-1}}$ spectral range (CO2, Tonkov)', fontsize=16)


# Discrepancy inset plot
ax_inset = inset_axes(ax, width="100%", height="100%", bbox_to_anchor=(0.2, 0.2, 0.5, 0.5), bbox_transform=ax.transAxes)
discrepancy = np.where((abs(cool1) > 0.1) & (abs(cool2) > 0.1), abs(cool1 - cool2), 0.0)
ax_inset.plot(discrepancy, height, color='black')
ax_inset.set_xlabel(r'discrepancy, $\bf{\mathregular{K/day}}$', fontweight='bold')
ax_inset.set_xlim(right=1)
ax_inset.set_ylim(bottom=0, top=100)

ticks = ax_inset.get_yticks()
ticks = ticks[ticks != 0]
ax_inset.set_yticks(ticks)

ax_inset.grid(which='major', axis='both', color='gray', alpha=0.5)



# Inset as a second plot
# fig, ax = plt.subplots()

# discrepancy = np.where((abs(cool1) > 0.1) & (abs(cool2) > 0.1), abs(cool1 - cool2), 0.0)

# ax.plot(discrepancy, height, color='black', linewidth=3)
# ax.set_xlabel(r'discrepancy, $\bf{\mathregular{K/day}}$', fontweight='bold')
# ax.set_xlim(right=2)
# ax.set_ylim(bottom=0, top=100)
# ticks = ax.get_yticks()
# ticks = ticks[ticks != 0]
# ax.set_yticks(ticks)
# ax.grid(which='major', axis='both', color='gray', alpha=0.5)


plt.tight_layout()
plt.show()