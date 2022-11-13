import numpy as np
import matplotlib.pyplot as plt

height, cool1, cool2 = np.loadtxt(fname='CO2_COOLING_3_', skiprows=1, unpack=True)

fig, ax = plt.subplots()

ax.plot(cool1, height, color='r')
ax.plot(cool2, height, color='b')
plt.grid(visible=True, which='major', color='g')
plt.xlim(-1, 3)
plt.ylim(50, 120)
plt.xlabel('Cooling rate [K/day]')
plt.ylabel('Height [km]')
plt.show()